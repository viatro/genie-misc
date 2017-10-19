//* made from gNtpConv.cxx */

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TFolder.h>
#include <TBits.h>
#include <TObjString.h>
#include <TMath.h>
#include <TDatabasePDG.h>

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
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
#include "Utils/T2KEvGenMetaData.h"

using std::string;
using std::ostringstream;
using std::ofstream;
using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::vector;

using namespace genie;
using namespace genie::constants;

//func prototypes
void   ConvertToGEST            (void);
void   GetCommandLineArgs        (int argc, char ** argv);
void   PrintSyntax               (void);
string DefaultOutputFile         (void);
bool   CheckRootFilename         (string filename);
//
int    BaryonNumber(int pdg);

//input options (from command line arguments):
string     gOptInpFileName;         ///< input file name
string     gOptOutFileName;         ///< output file name
int        gOptVersion;             ///< output file format version
Long64_t   gOptN;                   ///< number of events to process
bool       gOptCopyJobMeta = false; ///< copy MC job metadata (gconfig, genv TFolders)
long int   gOptRanSeed;             ///< random number seed

//genie version used to generate the input event file 
int gFileMajorVrs = -1;
int gFileMinorVrs = -1;
int gFileRevisVrs = -1;

//consts
const int kNPmax = 250;

//
TDatabasePDG *pdg_db = new TDatabasePDG();

//____________________________________________________________________________________
int main(int argc, char ** argv) {
    GetCommandLineArgs(argc, argv);
    
    utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
    utils::app_init::RandGen(gOptRanSeed);
    
    GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
    
    // Call the conversion function
    ConvertToGEST();
    
    return 0;
}
//____________________________________________________________________________________
void ConvertToGEST(void) {
    // Define branch variables
    //
    int     brIev             = 0;      // Event number 
    int     brNeutrino        = 0;      // Neutrino pdg code
    int     brFSPrimLept      = 0;      // Final state primary lepton pdg code
    int     brTarget          = 0;      // Nuclear target pdg code (10LZZZAAAI)
    int     brTargetZ         = 0;      // Nuclear target Z (extracted from pdg code above)
    int     brTargetA         = 0;      // Nuclear target A (extracted from pdg code above)
    int     brHitNuc          = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
    int     brHitQrk          = 0;      // Hit quark pdg code        (set for DIS events only)
    bool    brFromSea         = false;  // Hit quark is from sea     (set for DIS events only)
    int     brResId           = 0;      // Produced baryon resonance (set for resonance events only)
    TString brProcInfo        = "";
    TString brScatteringType  = "";
    TString brInteractionType = "";
    double  brXSec            = 0;      // Cross section for selected event (1E-38 cm2)
    double  brDXSec           = 0;      // Cross section for selected event kinematics (1E-38 cm2 /{K^n})
    double  brEvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
    double  brEv              = 0;      // Neutrino energy @ LAB
    double  brPxv             = 0;      // Neutrino px @ LAB
    double  brPyv             = 0;      // Neutrino py @ LAB
    double  brPzv             = 0;      // Neutrino pz @ LAB
    double  brEn              = 0;      // Initial state hit nucleon energy @ LAB
    double  brPxn             = 0;      // Initial state hit nucleon px @ LAB
    double  brPyn             = 0;      // Initial state hit nucleon py @ LAB
    double  brPzn             = 0;      // Initial state hit nucleon pz @ LAB
    int     brN               = 0;      // Nu. of final state particles in hadronic system
    int     brPdg        [kNPmax];      // Pdg code of k^th final state particle in hadronic system
    double  brE          [kNPmax];      // Energy     of k^th final state particle in hadronic system @ LAB
    double  brPx         [kNPmax];      // Px         of k^th final state particle in hadronic system @ LAB
    double  brPy         [kNPmax];      // Py         of k^th final state particle in hadronic system @ LAB
    double  brPz         [kNPmax];      // Pz         of k^th final state particle in hadronic system @ LAB
    double  brP          [kNPmax];      // P          of k^th final state particle in hadronic system @ LAB
    double  brKinE       [kNPmax];      // Kin.energy of k^th final state particle in hadronic system @ LAB
    double  brCosth      [kNPmax];      // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
    int     brPdgb            = 0;      // Nuclear remnant (hadronic blob) PDG code
    double  brEb              = 0;      // Nuclear remnant (hadronic blob) energy @ LAB
    double  brPxb             = 0;      // Nuclear remnant (hadronic blob) px @ LAB
    double  brPyb             = 0;      // Nuclear remnant (hadronic blob) py @ LAB
    double  brPzb             = 0;      // Nuclear remnant (hadronic blob) pz @ LAB
    double  brPb              = 0;      // Nuclear remnant (hadronic blob) p  @ LAB
    double  brKinEb           = 0;      // Nuclear remnant (hadronic blob) kin.energy @ LAB
    double  brVtxX               ;      // Vertex x in detector coord system (SI)
    double  brVtxY               ;      // Vertex y in detector coord system (SI)
    double  brVtxZ               ;      // Vertex z in detector coord system (SI)
    double  brVtxT               ;      // Vertex t in detector coord system (SI)
    double  brSumKE              ;      // Sum of kinetic energies of all final state particles
    
    // Open output file & create output summary tree & create the tree branches
    //
    LOG("gntpc", pNOTICE) 
            << "*** Saving summary tree to: " << gOptOutFileName;
    TFile fout(gOptOutFileName.c_str(),"recreate");
    
    TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");
    
    // Create tree branches
    //
    s_tree->Branch("iev",          &brIev,           "iev/I"         );
    s_tree->Branch("neu",          &brNeutrino,      "neu/I"         );
    s_tree->Branch("fspl",         &brFSPrimLept,    "fspl/I"        );
    s_tree->Branch("tgt",          &brTarget,        "tgt/I"         );
    s_tree->Branch("Z",            &brTargetZ,       "Z/I"           );
    s_tree->Branch("A",            &brTargetA,       "A/I"           );
    s_tree->Branch("hitnuc",       &brHitNuc,        "hitnuc/I"      );
    s_tree->Branch("hitqrk",       &brHitQrk,        "hitqrk/I"      );
    s_tree->Branch("resid",        &brResId,         "resid/I"       );
    s_tree->Branch("proc_info",    &brProcInfo                       );
    s_tree->Branch("scattering",   &brScatteringType                 );
    s_tree->Branch("interaction",  &brInteractionType                );
    s_tree->Branch("sea",          &brFromSea,       "sea/O"         );
    s_tree->Branch("xsec",         &brXSec,          "xsec/D"        );
    s_tree->Branch("dxsec",        &brDXSec,         "dxsec/D"       );
    s_tree->Branch("EvRF",         &brEvRF,          "EvRF/D"        );
    s_tree->Branch("Ev",           &brEv,            "Ev/D"          );
    s_tree->Branch("pxv",          &brPxv,           "pxv/D"         );
    s_tree->Branch("pyv",          &brPyv,           "pyv/D"         );
    s_tree->Branch("pzv",          &brPzv,           "pzv/D"         );
    s_tree->Branch("En",           &brEn,            "En/D"          );
    s_tree->Branch("pxn",          &brPxn,           "pxn/D"         );
    s_tree->Branch("pyn",          &brPyn,           "pyn/D"         );
    s_tree->Branch("pzn",          &brPzn,           "pzn/D"         );
    s_tree->Branch("n",            &brN,             "n/I"           );
    s_tree->Branch("pdg",           brPdg,           "pdg[n]/I "     );
    s_tree->Branch("E",             brE,             "E[n]/D"        );
    s_tree->Branch("px",            brPx,            "px[n]/D"       );
    s_tree->Branch("py",            brPy,            "py[n]/D"       );
    s_tree->Branch("pz",            brPz,            "pz[n]/D"       );
    s_tree->Branch("p",             brP,             "p[n]/D"        );
    s_tree->Branch("KE",            brKinE,          "KE[n]/D"       );
    s_tree->Branch("costh",         brCosth,         "costh[n]/D"    );
    s_tree->Branch("pdgb",         &brPdgb,          "pdgb/I"        );
    s_tree->Branch("Eb",           &brEb,            "Eb/D"          );
    s_tree->Branch("pxb",          &brPxb,           "pxb/D"         );
    s_tree->Branch("pyb",          &brPyb,           "pyb/D"         );
    s_tree->Branch("pzb",          &brPzb,           "pzb/D"         );
    s_tree->Branch("pb",           &brPb,            "pb/D"          );
    s_tree->Branch("KEb",          &brKinEb,         "KEb/D"         );
    s_tree->Branch("vtxx",         &brVtxX,          "vtxx/D"        );
    s_tree->Branch("vtxy",         &brVtxY,          "vtxy/D"        );
    s_tree->Branch("vtxz",         &brVtxZ,          "vtxz/D"        );
    s_tree->Branch("vtxt",         &brVtxT,          "vtxt/D"        );
    s_tree->Branch("sumKE",        &brSumKE,         "sumKE/D"       );
    
    // Open the ROOT file and get the TTree & its header
    TFile fin(gOptInpFileName.c_str(),"READ");
    TTree *           er_tree = 0;
    NtpMCTreeHeader * thdr    = 0;
    er_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
    thdr    = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );
    if (!er_tree) {
        LOG("gntpc", pERROR) << "Null input GHEP event tree";
        return;
    }
    LOG("gntpc", pINFO) << "Input tree header: " << *thdr;
    
    // Get the mc record
    NtpMCEventRecord * mcrec = 0;
    er_tree->SetBranchAddress("gmcrec", &mcrec);
    if (!mcrec) {
        LOG("gntpc", pERROR) << "Null MC record";
        return;
    }
    
    // Figure out how many events to analyze
    Long64_t nmax = (gOptN<0) ? 
            er_tree->GetEntries() : TMath::Min( er_tree->GetEntries(), gOptN );
    if (nmax<0) {
        LOG("gntpc", pERROR) << "Number of events = 0";
        return;
    }
    
    LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";
    
    TLorentzVector pdummy(0,0,0,0);
    
    // Event loop
    for(Long64_t iev = 0; iev < nmax; iev++) {
        er_tree->GetEntry(iev);
    
        NtpMCRecHeader rec_header = mcrec->hdr;
        EventRecord &  event      = *(mcrec->event);
    
        LOG("gntpc", pINFO) << rec_header;
        LOG("gntpc", pINFO) << event;
    
        // Go further only if the event is physical
        bool is_unphysical = event.IsUnphysical();
        if(is_unphysical) {
            LOG("gntpc", pINFO) << "Skipping unphysical event";
            mcrec->Clear();
            continue;
        }
    
        // Clean-up arrays
        //
        for(int j=0; j<kNPmax; j++) {
            brPdg   [j] =  0;
            brE     [j] =  0;
            brPx    [j] =  0;
            brPy    [j] =  0;
            brPz    [j] =  0;
            brP     [j] =  0;
            brKinE  [j] =  0;
            brCosth [j] =  0;
        }
    
        // Computing event characteristics
        //
    
        //input particles
        GHepParticle * neutrino = event.Probe();
        GHepParticle * target = event.TargetNucleus();
        GHepParticle * fspl = event.FinalStatePrimaryLepton();
        GHepParticle * hitnucl = event.HitNucleon();
    
        int tgtZ = 0;
        int tgtA = 0;
        if(pdg::IsIon(target->Pdg())) {
            tgtZ = pdg::IonPdgCodeToZ(target->Pdg());
            tgtA = pdg::IonPdgCodeToA(target->Pdg());
        } 
        if(target->Pdg() == kPdgProton   ) { tgtZ = 1; tgtA = 1; }    
        if(target->Pdg() == kPdgNeutron  ) { tgtZ = 0; tgtA = 1; }    
    
        // Summary info
        const Interaction * interaction = event.Summary();
        const InitialState & init_state = interaction->InitState();
        const ProcessInfo &  proc_info  = interaction->ProcInfo();
        //const Kinematics &   kine       = interaction->Kine();
        const XclsTag &      xcls       = interaction->ExclTag();
        const Target &       tgt        = init_state.Tgt();
    
        // Vertex in detector coord system
        TLorentzVector * vtx = event.Vertex();
    
        // Process id
        //bool is_qel     = proc_info.IsQuasiElastic();
        bool is_res     = proc_info.IsResonant();
        bool is_dis     = proc_info.IsDeepInelastic();
        bool is_coh     = proc_info.IsCoherent();
        //bool is_dfr     = proc_info.IsDiffractive();
        bool is_imd     = proc_info.IsInverseMuDecay();
        bool is_imdanh  = proc_info.IsIMDAnnihilation();
        //bool is_singlek = proc_info.IsSingleKaon();    
        bool is_nuel    = proc_info.IsNuElectronElastic();
        //bool is_em      = proc_info.IsEM();
        //bool is_weakcc  = proc_info.IsWeakCC();
        //bool is_weaknc  = proc_info.IsWeakNC();
        //bool is_mec     = proc_info.IsMEC();
    
        if(!hitnucl) { assert(is_coh || is_imd || is_imdanh || is_nuel); }
    
        // Hit quark - set only for DIS events
        int  qrk  = (is_dis) ? tgt.HitQrkPdg() : 0;     
        bool seaq = (is_dis) ? tgt.HitSeaQrk() : false; 
    
        // Resonance id ($GENIE/src/BaryonResonance/BaryonResonance.h) -
        // set only for resonance neutrinoproduction
        int resid = (is_res) ? EResonance(xcls.Resonance()) : -99;
    
        // (qel or dis) charm production?
        //bool charm = xcls.IsCharmEvent();
    
        double xsec   = (1E+38/units::cm2) * event.XSec();
        double dxsec  = (1E+38/units::cm2) * event.DiffXSec();
    
        const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
        //const TLorentzVector & k2 = *(fspl->P4());                          // l 4-p (k2)
        const TLorentzVector & p1 = (hitnucl) ? *(hitnucl->P4()) : pdummy; // N 4-p (p1)      
    
        //double M  = kNucleonMass; 
        //TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
        //double Q2 = -1 * q.M2();                       // momemtum transfer
        //double v  = (hitnucl) ? q.Energy()       : -1; // v (E transfer to the nucleus)
        //double x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
        //double y  = (hitnucl) ? v/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1
        //double W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
        //double W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
        //double t  = 0;
    
        // Get v 4-p at hit nucleon rest-frame
        TLorentzVector k1_rf = k1;         
        if(hitnucl) {
            k1_rf.Boost(-1.*p1.BoostVector());
        }
        
        //
        TObjArrayIter piter(&event);
        GHepParticle * p = 0;
        int ip=-1;
        
        // Extract the final state system originating from the hadronic vertex 
        // (after the intranuclear rescattering step)
        //
    
        LOG("gntpc", pDEBUG) << "Extracting final state hadronic system";
    
        vector<int> final_had_syst;
        final_had_syst.push_back(event.FinalStatePrimaryLeptonPosition());
        
        while ( (p = (GHepParticle *) piter.Next()) ) {
            ip++;
            // don't count final state lepton as part hadronic system
            if (event.Particle(ip)->FirstMother()==0) continue;
            int pdgc = p->Pdg();
            int ist  = p->Status();
            if (pdg::IsPseudoParticle(pdgc)) continue;
            if (ist==kIStStableFinalState  &&  pdgc < 2000000000)  final_had_syst.push_back(ip);
        }//particle-loop
    
        //if (count(final_had_syst.begin(), final_had_syst.end(), -1) > 0) {
        //    mcrec->Clear();
        //    continue;
        //}
    
        //
        // All information has been assembled -- Start filling up the tree branches
        //
        brIev        = (int) iev;      
        brNeutrino   = neutrino->Pdg();      
        brFSPrimLept = fspl->Pdg();
        brTarget     = target->Pdg(); 
        brTargetZ    = tgtZ;
        brTargetA    = tgtA;   
        brHitNuc     = (hitnucl) ? hitnucl->Pdg() : 0;      
        brHitQrk     = qrk;     
        brFromSea    = seaq;  
        brResId      = resid;
        brProcInfo   = proc_info.AsString();
        brScatteringType = proc_info.ScatteringTypeAsString();
        brInteractionType = proc_info.InteractionTypeAsString();
        brXSec       = xsec;
        brDXSec      = dxsec;
        brEvRF       = k1_rf.Energy();      
        brEv         = k1.Energy();      
        brPxv        = k1.Px();  
        brPyv        = k1.Py();  
        brPzv        = k1.Pz();  
        brEn         = (hitnucl) ? p1.Energy() : 0;
        brPxn        = (hitnucl) ? p1.Px()     : 0;
        brPyn        = (hitnucl) ? p1.Py()     : 0;
        brPzn        = (hitnucl) ? p1.Pz()     : 0;
        
        brSumKE = 0;
        brN = final_had_syst.size();
        for (int j = 0; j < brN; j++) {
            p = event.Particle(final_had_syst[j]);
    
            int    fpdg = p->Pdg();     
            double fE   = p->Energy();     
            double fKE  = p->KinE();     
            double fpx  = p->Px();     
            double fpy  = p->Py();     
            double fpz  = p->Pz();     
            double fp   = TMath::Sqrt(fpx*fpx + fpy*fpy + fpz*fpz);
            double fcth = TMath::Cos( p->P4()->Vect().Angle(k1.Vect()) );
    
            brPdg  [j] = fpdg;
            brE    [j] = fE;
            brPx   [j] = fpx;
            brPy   [j] = fpy;
            brPz   [j] = fpz;
            brP    [j] = fp;
            brKinE [j] = fKE;
            brCosth[j] = fcth;
    
            brSumKE += fKE;
        }
        
        brPdgb    = 0;
        brEb      = 0;
        brKinEb   = 0;
        brPxb     = 0;
        brPyb     = 0;
        brPzb     = 0;
        brPb      = 0;
        
        piter.Reset();
        GHepParticle * hadr_blob = 0;
        double charge = 0;
        int baryon_number = 0;
        
        while ( (p = (GHepParticle *) piter.Next())) {
            int pdgc = p->Pdg();
            int ist  = p->Status();
            if (pdg::IsPseudoParticle(pdgc)) {
                if (ist == kIStFinalStateNuclearRemnant) hadr_blob = p;
                continue;
            }
            if (ist == kIStInitialState) {
                if (pdg::IsIon(pdgc)) {
                    charge += pdg::IonPdgCodeToZ(pdgc);
                    baryon_number += pdg::IonPdgCodeToA(pdgc);
                } else {
                   charge += TMath::Nint(p->Charge()/3);
                   baryon_number += BaryonNumber(pdgc);
                }
            } else if (ist == kIStStableFinalState) {
                if (pdg::IsIon(pdgc)) {
                    charge -= pdg::IonPdgCodeToZ(pdgc);
                    baryon_number -= pdg::IonPdgCodeToA(pdgc);
                } else {
                   charge -= TMath::Nint(p->Charge()/3);
                   baryon_number -= BaryonNumber(pdgc);
                }
            }
        }
        
        if (hadr_blob) {
            /*if (int(charge) == 0) {
                hadr_blob->SetPdgCode(2112);
                hadr_blob->SetEnergy(hadr_blob->Energy()/baryon_number);
                hadr_blob->SetPx(hadr_blob->Px()/baryon_number);
                hadr_blob->SetPy(hadr_blob->Py()/baryon_number);
                hadr_blob->SetPz(hadr_blob->Pz()/baryon_number);
                for (int ib = 0; ib < baryon_number; ++ib) {
                    brPdg  [brN]  = 2112;
                    brE    [brN]  = hadr_blob->Energy();
                    brPx   [brN]  = hadr_blob->Px();
                    brPy   [brN]  = hadr_blob->Py();
                    brPz   [brN]  = hadr_blob->Pz();
                    brP    [brN]  = TMath::Sqrt(brPxb*brPxb + brPyb*brPyb + brPzb*brPzb);
                    brKinE [brN]  = hadr_blob->KinE(true);
                    if (brKinE[brN] < 0) brKinE[brN] = 0;
                    brCosth[brN] = 0;
                    
                    brN++;
                }
            } else {*/
                brPdgb   = pdg::IonPdgCode(baryon_number, charge);
                brEb     = hadr_blob->Energy();
                brKinEb  = hadr_blob->KinE();
                brPxb    = hadr_blob->Px();
                brPyb    = hadr_blob->Py();
                brPzb    = hadr_blob->Pz();
                brPb     = TMath::Sqrt(brPxb*brPxb + brPyb*brPyb + brPzb*brPzb);
            //}
        }
    
        brVtxX = vtx->X();   
        brVtxY = vtx->Y();   
        brVtxZ = vtx->Z();   
        brVtxT = vtx->T();
    
        s_tree->Fill();
    
        mcrec->Clear();

    } // event loop


    // Copy MC job metadata (gconfig and genv TFolders)
    if(gOptCopyJobMeta) {
        TFolder * genv    = (TFolder*) fin.Get("genv");
        TFolder * gconfig = (TFolder*) fin.Get("gconfig");
        fout.cd();       
        genv    -> Write("genv");
        gconfig -> Write("gconfig");
    }
    
    fin.Close();
    
    fout.cd();
    s_tree->Write("",TObject::kOverwrite);
    //fout.Write();
    fout.Close();
}

int BaryonNumber(int pdg) {
    if (abs(pdg) < 1000000000) {
        TParticlePDG *pPDG = pdg_db->GetParticle(pdg);
        if (!pPDG) return 0;
        if ( ! TString(pPDG->ParticleClass()).EqualTo("Baryon") ) return 0;
        if (pdg > 0) return  1;
        if (pdg < 0) return -1;
        return 0;
    }
    if (abs(pdg) >= 2000000000) return 0;
    else return (pdg / 10) % 1000;
}
//____________________________________________________________________________________
// FUNCTIONS FOR PARSING CMD-LINE ARGUMENTS 
//____________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // Common run options. 
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app
  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('i') ) {
    LOG("gntpc", pINFO) << "Reading input filename";
    gOptInpFileName = parser.ArgAsString('i');
  } else {
    LOG("gntpc", pFATAL)
         << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // check input GENIE ROOT file
  bool inpok = !(gSystem->AccessPathName(gOptInpFileName.c_str()));
  if (!inpok) {
    LOG("gntpc", pFATAL)
        << "The input ROOT file ["
        << gOptInpFileName << "] is not accessible";
    gAbortingInErr = true;
    exit(2);
  }

  // get output file name 
  if( parser.OptionExists('o') ) {
    LOG("gntpc", pINFO) << "Reading output filename";
    gOptOutFileName = parser.ArgAsString('o');
  } else {
    LOG("gntpc", pINFO)
         << "Unspecified output filename - Using default";
    gOptOutFileName = DefaultOutputFile();
  }

  // get number of events to convert
  if( parser.OptionExists('n') ) {
    LOG("gntpc", pINFO) << "Reading number of events to analyze";
    gOptN = parser.ArgAsInt('n');
  } else {
    LOG("gntpc", pINFO)
         << "Unspecified number of events to analyze - Use all";
    gOptN = -1;
  }

  // check whether to copy MC job metadata (only if output file is in ROOT format)
  gOptCopyJobMeta = parser.OptionExists('c');

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gntpc", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gntpc", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  LOG("gntpc", pNOTICE) << "Input filename  = " << gOptInpFileName;
  LOG("gntpc", pNOTICE) << "Output filename = " << gOptOutFileName;
  LOG("gntpc", pNOTICE) << "Conversion to format = " << gOptRanSeed 
                        << ", vrs = " << gOptVersion;
  LOG("gntpc", pNOTICE) << "Number of events to be converted = " << gOptN;
  LOG("gntpc", pNOTICE) << "Copy metadata? = " << ((gOptCopyJobMeta) ? "Yes" : "No");
  LOG("gntpc", pNOTICE) << "Random number seed = " << gOptRanSeed;

  LOG("gntpc", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________________
string DefaultOutputFile(void)
{
  // filename extension - depending on file format
  string ext = "gest.root";
  string inpname = gOptInpFileName;
  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }

  // remove ghep.
  size_t pos = inpname.find("ghep.");
  if(pos != string::npos) {
    inpname.erase(pos, pos+4);
  }

  ostringstream name;
  name << inpname << ext;

  return gSystem->BaseName(name.str().c_str());
}
//____________________________________________________________________________________
//____________________________________________________________________________________
void PrintSyntax(void)
{
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/stdapp/gNtpConv.cxx");
  string cmd      = "less " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//____________________________________________________________________________________
