#include <iostream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TFolder.h>
#include <TBits.h>
#include <TString.h>
#include <TObjString.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpWriter.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/SystemUtils.h"

using std::cerr;
using std::cout;
using std::endl;

using namespace genie;
using namespace genie::constants;

const int kNPmax = 250;

void Convert(const TString&, const TString&);
void GetCommandLineArgs (int argc, char ** argv);

int main(int argc, char ** argv) {
    
    GetCommandLineArgs(argc, argv);
    
    utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
    //utils::app_init::RandGen(gOptRanSeed);
    
    GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
    
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
    int    nParticles    =  0;   // Nu. of primary particles
    double Weight        =  0;   // Event weight
    double Prob          =  0;   // Probability for that event (given cross section, path lengths, etc)
    double XSec          =  0;   // Cross section for selected event (femtobarns: 1E-39 cm2)
    double DXSec         =  0;   // Cross section for selected event kinematics (femtobarns: 1E-39 cm2 /{K^n})
    TClonesArray *parr = new TClonesArray("TParticle",100);
    TClonesArray &ar = *parr;
    TLorentzVector *vertex = new TLorentzVector(0,0,0,0);
    
    TFile fout(OutFileName.Data(),"recreate");
    TTree * T = new TTree("T","GENIE to be G4 tracked");
    T->Branch("nParticles",&nParticles,"nParticles/I");
    T->Branch("weight",&Weight,"weight/D");
    T->Branch("prob",&Prob,"prob/D");
    T->Branch("xsec",&XSec,"xsec/D");
    T->Branch("dxsec",&DXSec,"dxsec/D");
    T->Branch("vertex",&vertex,64);
    T->Branch("particles",&parr/*,1024000*/);
    
    // Open the ROOT file and get the TTree & its header
    TFile fin(InpFileName.Data(),"READ");
    TTree *           ghep_tree = 0;
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
        
        //NtpMCRecHeader rec_header = mcrec->hdr;
        EventRecord &  event      = *(mcrec->event);
        
        if (event.IsUnphysical()) {
            cerr << "Event " << iev << "is unphysical. Skipping" << endl; 
            mcrec->Clear();
            continue;
        }
        
        // Get event weight-probability-xsections
        Weight = event.Weight();
        Prob = event.Probability();
        XSec = (1E+39/units::cm2) * event.XSec();
        DXSec = (1E+39/units::cm2) * event.DiffXSec();
        
        vertex = event.Vertex();
        
        //Particles loop
        int ip = 0;
        GHepParticle *p = 0;
        TIter event_iter(&event);
        while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
            if (!p) continue; //assert(p);
            //if( (p->Status() != kIStStableFinalState) || (!pdg::IsParticle( p->Pdg() )) ) continue;
            //if( (p->Status() != kIStStableFinalState) && (p->Status() != kIStFinalStateNuclearRemnant) ) continue;
            //if ( pdg::IsPseudoParticle( p->Pdg() ) && (p->Status() != kIStFinalStateNuclearRemnant) ) continue;
            TParticle* particle = new(ar[ip]) TParticle(p->Pdg(), p->Status(), p->FirstMother(),p->LastMother(),p->FirstDaughter(),p->LastDaughter(), *p->P4(), *p->X4() );  //p[GeV],x[fm]
            //p->Px(), p->Py(), p->Pz(), p->E(), p->Vx(), p->Vy(), p->Vz(), p->Vt() );
            if (p->PolzIsSet()) {
                p->GetPolarization(polz);
                particle->SetPolarisation(polz);
            }
            ip++;            
        }
        nParticles = ip;
        //parr->Compress();
        T->Fill();
        mcrec->Clear();
        parr->Delete();//Clear();
    } // event loop
    
    fin.Close();
    
    fout.Write();
    fout.Close();
}

void GetCommandLineArgs(int argc, char ** argv)
{
  // Common run options. 
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app
  CmdLnArgParser parser(argc,argv);
}
//