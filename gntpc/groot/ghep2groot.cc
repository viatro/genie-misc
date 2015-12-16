#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using std::cerr;
using std::cout;
using std::endl;

using namespace genie;

void Convert(const TString&, const TString&);

int main(int argc, char ** argv) {
    if ( (argc < 2) || (argc > 3) || (TString(argv[1]) == TString("-h")) || (TString(argv[1]) == TString("--help")) ) {
        cout << "USAGE:\t" << argv[0] << " INPUT_FILE.root [OUTPUT_FILE.root]" << endl;
    } else if (argc == 3) {
        Convert(argv[1], argv[2]);
    } else {
        TString s = ((TObjString*)(TString(argv[1]).Tokenize("/")->Last()))->GetName();
        TString old_ext = "ghep.root";
        TString new_ext = "groot.root";
        s.Replace(s.Length() - old_ext.Length(), old_ext.Length(), new_ext);
        Convert(argv[1], s);
    }
    return 0;
}

void Convert(const TString& InpFileName, const TString& OutFileName) {
    int    brIev           =  0;   // Event number
    int    brNParticles    =  0;   // Nu. of primary particles
    TString brScattering   = "";
    TString brInteraction  = "";
    double brXSec          =  0;   // Cross section for selected event (femtobarns: 1E-39 cm2)
    double brDXSec         =  0;   // Cross section for selected event kinematics (femtobarns: 1E-39 cm2 /{K^n})
    double brProb          =  0;   // Probability for that event (given cross section, path lengths, etc)
    double brWeight        =  0;   // Event weight
        
    TClonesArray * parr = new TClonesArray("TParticle", 25);
    TClonesArray &ar = *parr;
    TLorentzVector *vertex = new TLorentzVector(0,0,0,0);
    
    TFile fout(OutFileName.Data(),"recreate");
    TTree * T = new TTree("T", "GENIE to be G4 tracked");
    
    T->Branch("iev",    &brIev,             "iev/I"     );
    T->Branch("n",      &brNParticles,      "n/I"       );
    T->Branch("sctr",   &brScattering                   );
    T->Branch("intr",   &brInteraction                  );
    T->Branch("xsec",   &brXSec,            "xsec/D"    );
    T->Branch("dxsec",  &brDXSec,           "dxsec/D"   );
    T->Branch("prob",   &brProb,            "prob/D"    );
    T->Branch("weight", &brWeight,          "weight/D"  );
    T->Branch("vtx",    &vertex                         );
    T->Branch("p",      &parr                           );
    
    // Open the ROOT file and get the TTree & its header
    TFile fin(InpFileName.Data(),"READ");
    TTree * ghep_tree = 0;
    ghep_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
    if (!ghep_tree) {
        cerr << "Null input GHEP event tree" << endl;
        return;
    }
    
    // Get the MC record
    NtpMCEventRecord *mcrec = 0;
    ghep_tree->SetBranchAddress("gmcrec", &mcrec);
    if (!mcrec) {
        cerr << "Null MC record" << endl;
        return;
    }
    
    TVector3 polz(0,0,0);
    
    // Event loop
    for(Long64_t iev = 0; iev < ghep_tree->GetEntries(); iev++) {
        ghep_tree->GetEntry(iev);
        brIev = iev;
        
        //NtpMCRecHeader rec_header = mcrec->hdr;
        EventRecord &  event      = *(mcrec->event);
        
        if (event.IsUnphysical()) {
            cerr << "Event " << iev << "is unphysical. Skipping" << endl; 
            mcrec->Clear();
            continue;
        }
        
        const ProcessInfo &  proc_info  = event.Summary()->ProcInfo();
        brScattering = proc_info.ScatteringTypeAsString();
        brInteraction = proc_info.InteractionTypeAsString();
        
        
        // Get event weight, probability, cross-sections
        brXSec   = (1E+38/units::cm2) * event.XSec();
        brDXSec  = (1E+38/units::cm2) * event.DiffXSec();
        brProb   = event.Probability();
        brWeight = event.Weight();
        
        vertex = event.Vertex();
        
        //Particles loop
        int ip = 0;
        GHepParticle *p = 0;
        TIter event_iter(&event);
        while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
            if (!p) continue;
            if (pdg::IsPseudoParticle(p->Pdg()) && p->Status() != kIStFinalStateNuclearRemnant )    continue;
            //p[GeV],x[fm]
            //p->Px(), p->Py(), p->Pz(), p->E(), p->Vx(), p->Vy(), p->Vz(), p->Vt() );
            TParticle* particle = new(ar[ip]) TParticle(p->Pdg(), p->Status(), p->FirstMother(),p->LastMother(),p->FirstDaughter(),p->LastDaughter(), *p->P4(), *p->X4() );  
            if (p->PolzIsSet()) {
                p->GetPolarization(polz);
                particle->SetPolarisation(polz);
            }
            ip++;            
        }
        brNParticles = ip;
        //parr->Compress();
        T->Fill();
        mcrec->Clear();
        parr->Delete();//Clear();
    } // event loop
    
    fin.Close();
    
    fout.Write();
    fout.Close();
}