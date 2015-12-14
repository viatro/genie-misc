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
void   ConvertToSGST             (void);
void   GetCommandLineArgs        (int argc, char ** argv);
void   PrintSyntax               (void);
string DefaultOutputFile         (void);
bool   CheckRootFilename         (string filename);

//format enum
typedef enum EGNtpcFmt {
  kConvFmt_undef = 0,
  kConvFmt_gst,
  kConvFmt_sgst,
  kConvFmt_gxml,
  kConvFmt_ghep_mock_data,
  kConvFmt_rootracker,
  kConvFmt_rootracker_mock_data,
  kConvFmt_t2k_rootracker,
  kConvFmt_numi_rootracker,
  kConvFmt_t2k_tracker,
  kConvFmt_nuance_tracker,
  kConvFmt_ghad,
  kConvFmt_ginuke
} GNtpcFmt_t;

//input options (from command line arguments):
string     gOptInpFileName;         ///< input file name
string     gOptOutFileName;         ///< output file name
GNtpcFmt_t gOptOutFileFormat;       ///< output file format id
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
//____________________________________________________________________________________
int main(int argc, char ** argv) {
    GetCommandLineArgs(argc, argv);
    
    utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
    utils::app_init::RandGen(gOptRanSeed);
    
    GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
    
    // Call the conversion function
    ConvertToSGST();
    
    return 0;
}
//____________________________________________________________________________________
void ConvertToSGST(void) {
  // Some constants
  const double e_h = 1.3; // typical e/h ratio used for computing mean `calorimetric response'

  // Define branch variables
  //
  int     brIev         = 0;      // Event number 
  int     brNeutrino    = 0;      // Neutrino pdg code
  int     brFSPrimLept  = 0;      // Final state primary lepton pdg code
  int     brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int     brTargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
  int     brTargetA     = 0;      // Nuclear target A (extracted from pdg code above)
  int     brHitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
  int     brHitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
  bool    brFromSea     = false;  // Hit quark is from sea     (set for DIS events only)
  int     brResId       = 0;      // Produced baryon resonance (set for resonance events only)
  TString brProcInfo    = "";
  TString brScatteringType = "";
  TString brInteractionType = "";
  double  brWeight      = 0;      // Event weight
  double  brKineXs      = 0;      // Bjorken x as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double  brKineYs      = 0;      // Inelasticity y as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double  brKineTs      = 0;      // Energy transfer to nucleus at COH events as was generated during kinematical selection
  double  brKineQ2s     = 0;      // Momentum transfer Q^2 as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double  brKineWs      = 0;      // Hadronic invariant mass W as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double  brKineX       = 0;      // Experimental-like Bjorken x; neglects fermi momentum / off-shellness 
  double  brKineY       = 0;      // Experimental-like inelasticity y; neglects fermi momentum / off-shellness 
  double  brKineT       = 0;      // Experimental-like energy transfer to nucleus at COH events 
  double  brKineQ2      = 0;      // Experimental-like momentum transfer Q^2; neglects fermi momentum / off-shellness
  double  brKineW       = 0;      // Experimental-like hadronic invariant mass W; neglects fermi momentum / off-shellness 
  double  brEvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
  double  brEv          = 0;      // Neutrino energy @ LAB
  double  brPxv         = 0;      // Neutrino px @ LAB
  double  brPyv         = 0;      // Neutrino py @ LAB
  double  brPzv         = 0;      // Neutrino pz @ LAB
  double  brEn          = 0;      // Initial state hit nucleon energy @ LAB
  double  brPxn         = 0;      // Initial state hit nucleon px @ LAB
  double  brPyn         = 0;      // Initial state hit nucleon py @ LAB
  double  brPzn         = 0;      // Initial state hit nucleon pz @ LAB
  double  brEl          = 0;      // Final state primary lepton energy @ LAB
  double  brPxl         = 0;      // Final state primary lepton px @ LAB
  double  brPyl         = 0;      // Final state primary lepton py @ LAB
  double  brPzl         = 0;      // Final state primary lepton pz @ LAB
  double  brPl          = 0;      // Final state primary lepton p  @ LAB
  double  brCosthl      = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
  int     brNf          = 0;      // Nu. of final state particles in hadronic system
  int     brPdgf  [kNPmax];       // Pdg code of k^th final state particle in hadronic system
  double  brEf    [kNPmax];       // Energy     of k^th final state particle in hadronic system @ LAB
  double  brPxf   [kNPmax];       // Px         of k^th final state particle in hadronic system @ LAB
  double  brPyf   [kNPmax];       // Py         of k^th final state particle in hadronic system @ LAB
  double  brPzf   [kNPmax];       // Pz         of k^th final state particle in hadronic system @ LAB
  double  brPf    [kNPmax];       // P          of k^th final state particle in hadronic system @ LAB
  double  brCosthf[kNPmax];       // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
  int     brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int     brPdgi[kNPmax];         // Pdg code of k^th particle in 'primary' hadronic system 
  int     brResc[kNPmax];         // FSI code of k^th particle in 'primary' hadronic system 
  double  brEi  [kNPmax];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
  double  brPxi [kNPmax];         // Px       of k^th particle in 'primary' hadronic system @ LAB
  double  brPyi [kNPmax];         // Py       of k^th particle in 'primary' hadronic system @ LAB
  double  brPzi [kNPmax];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
  double  brVtxX;                 // Vertex x in detector coord system (SI)
  double  brVtxY;                 // Vertex y in detector coord system (SI)
  double  brVtxZ;                 // Vertex z in detector coord system (SI)
  double  brVtxT;                 // Vertex t in detector coord system (SI)
  double  brSumKEf;               // Sum of kinetic energies of all final state particles

  // Open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE) 
         << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");

  TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");

  // Create tree branches
  //
  s_tree->Branch("iev",           &brIev,           "iev/I"         );
  s_tree->Branch("neu",           &brNeutrino,      "neu/I"         );
  s_tree->Branch("fspl",          &brFSPrimLept,    "fspl/I"        );
  s_tree->Branch("tgt",           &brTarget,        "tgt/I"         );
  s_tree->Branch("Z",             &brTargetZ,       "Z/I"           );
  s_tree->Branch("A",             &brTargetA,       "A/I"           );
  s_tree->Branch("hitnuc",        &brHitNuc,        "hitnuc/I"      );
  s_tree->Branch("hitqrk",        &brHitQrk,        "hitqrk/I"      );
  s_tree->Branch("resid",         &brResId,         "resid/I"       );
  s_tree->Branch("proc_info",     &brProcInfo                       );
  s_tree->Branch("scattering",    &brScatteringType                 );
  s_tree->Branch("interaction",   &brInteractionType                );
  s_tree->Branch("sea",           &brFromSea,       "sea/O"         );
  s_tree->Branch("wght",          &brWeight,        "wght/D"        );
  s_tree->Branch("xs",            &brKineXs,        "xs/D"          );
  s_tree->Branch("ys",            &brKineYs,        "ys/D"          );
  s_tree->Branch("ts",            &brKineTs,        "ts/D"          );
  s_tree->Branch("Q2s",           &brKineQ2s,       "Q2s/D"         );
  s_tree->Branch("Ws",            &brKineWs,        "Ws/D"          );
  s_tree->Branch("x",             &brKineX,         "x/D"           );
  s_tree->Branch("y",             &brKineY,         "y/D"           );
  s_tree->Branch("t",             &brKineT,         "t/D"           );
  s_tree->Branch("Q2",            &brKineQ2,        "Q2/D"          );
  s_tree->Branch("W",             &brKineW,         "W/D"           );
  s_tree->Branch("EvRF",          &brEvRF,          "EvRF/D"        );
  s_tree->Branch("Ev",            &brEv,            "Ev/D"          );
  s_tree->Branch("pxv",           &brPxv,           "pxv/D"         );
  s_tree->Branch("pyv",           &brPyv,           "pyv/D"         );
  s_tree->Branch("pzv",           &brPzv,           "pzv/D"         );
  s_tree->Branch("En",            &brEn,            "En/D"          );
  s_tree->Branch("pxn",           &brPxn,           "pxn/D"         );
  s_tree->Branch("pyn",           &brPyn,           "pyn/D"         );
  s_tree->Branch("pzn",           &brPzn,           "pzn/D"         );
  s_tree->Branch("El",            &brEl,            "El/D"          );
  s_tree->Branch("pxl",           &brPxl,           "pxl/D"         );
  s_tree->Branch("pyl",           &brPyl,           "pyl/D"         );
  s_tree->Branch("pzl",           &brPzl,           "pzl/D"         );
  s_tree->Branch("pl",            &brPl,            "pl/D"          );
  s_tree->Branch("cthl",          &brCosthl,        "cthl/D"        );
  s_tree->Branch("ni",           &brNi,             "ni/I"          );
  s_tree->Branch("pdgi",          brPdgi,           "pdgi[ni]/I"    );
  s_tree->Branch("resc",          brResc,           "resc[ni]/I"    );
  s_tree->Branch("Ei",            brEi,             "Ei[ni]/D"      );
  s_tree->Branch("pxi",           brPxi,            "pxi[ni]/D"     );
  s_tree->Branch("pyi",           brPyi,            "pyi[ni]/D"     );
  s_tree->Branch("pzi",           brPzi,            "pzi[ni]/D"     );
  s_tree->Branch("nf",           &brNf,             "nf/I"          );
  s_tree->Branch("pdgf",          brPdgf,           "pdgf[nf]/I "   );
  s_tree->Branch("Ef",            brEf,             "Ef[nf]/D"      );
  s_tree->Branch("pxf",           brPxf,            "pxf[nf]/D"     );
  s_tree->Branch("pyf",           brPyf,            "pyf[nf]/D"     );
  s_tree->Branch("pzf",           brPzf,            "pzf[nf]/D"     );
  s_tree->Branch("pf",            brPf,             "pf[nf]/D"      );
  s_tree->Branch("cthf",          brCosthf,         "cthf[nf]/D"    );
  s_tree->Branch("vtxx",         &brVtxX,           "vtxx/D"        );
  s_tree->Branch("vtxy",         &brVtxY,           "vtxy/D"        );
  s_tree->Branch("vtxz",         &brVtxZ,           "vtxz/D"        );
  s_tree->Branch("vtxt",         &brVtxT,           "vtxt/D"        );
  s_tree->Branch("sumKEf",       &brSumKEf,         "sumKEf/D"      );

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
         brPdgi   [j] =  0;     
         brResc   [j] = -1;     
         brEi     [j] =  0;     
         brPxi    [j] =  0;     
         brPyi    [j] =  0;     
         brPzi    [j] =  0;     
         brPdgf   [j] =  0;     
         brEf     [j] =  0;     
         brPxf    [j] =  0;     
         brPyf    [j] =  0;     
         brPzf    [j] =  0;     
         brPf     [j] =  0;     
         brCosthf [j] =  0;     
    }

    // Computing event characteristics
    //

    //input particles
    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    GHepParticle * target = event.Particle(1);
    assert(target);
    GHepParticle * fspl = event.FinalStatePrimaryLepton();
    assert(fspl);
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
    const Kinematics &   kine       = interaction->Kine();
    const XclsTag &      xcls       = interaction->ExclTag();
    const Target &       tgt        = init_state.Tgt();

    // Vertex in detector coord system
    TLorentzVector * vtx = event.Vertex();

    // Process id
    bool is_qel     = proc_info.IsQuasiElastic();
    bool is_res     = proc_info.IsResonant();
    bool is_dis     = proc_info.IsDeepInelastic();
    bool is_coh     = proc_info.IsCoherent();
    bool is_dfr     = proc_info.IsDiffractive();
    bool is_imd     = proc_info.IsInverseMuDecay();
    bool is_imdanh  = proc_info.IsIMDAnnihilation();
    bool is_singlek = proc_info.IsSingleKaon();    
    bool is_nuel    = proc_info.IsNuElectronElastic();
    bool is_em      = proc_info.IsEM();
    bool is_weakcc  = proc_info.IsWeakCC();
    bool is_weaknc  = proc_info.IsWeakNC();
    bool is_mec     = proc_info.IsMEC();

    if(!hitnucl) { assert(is_coh || is_imd || is_imdanh || is_nuel); }
  
    // Hit quark - set only for DIS events
    int  qrk  = (is_dis) ? tgt.HitQrkPdg() : 0;     
    bool seaq = (is_dis) ? tgt.HitSeaQrk() : false; 

    // Resonance id ($GENIE/src/BaryonResonance/BaryonResonance.h) -
    // set only for resonance neutrinoproduction
    int resid = (is_res) ? EResonance(xcls.Resonance()) : -99;

    // (qel or dis) charm production?
    bool charm = xcls.IsCharmEvent();

    // Get event weight
    double weight = event.Weight();

    // Access kinematical params _exactly_ as they were selected internally
    // (at the hit nucleon rest frame; 
    // for bound nucleons: taking into account fermi momentum and off-shell kinematics)
    //
    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = (is_coh || is_dfr) ? kine.t (get_selected) : -1;
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    LOG("gntpc", pDEBUG) 
         << "[Select] Q2 = " << Q2s << ", W = " << Ws 
         << ", x = " << xs << ", y = " << ys << ", t = " << ts;

    // Calculate the same kinematical params but now as an experimentalist would 
    // measure them by neglecting the fermi momentum and off-shellness of bound nucleons
    //

    const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
    const TLorentzVector & k2 = *(fspl->P4());                          // l 4-p (k2)
    const TLorentzVector & p1 = (hitnucl) ? *(hitnucl->P4()) : pdummy; // N 4-p (p1)      

    double M  = kNucleonMass; 
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
    double Q2 = -1 * q.M2();                       // momemtum transfer
    double v  = (hitnucl) ? q.Energy()       : -1; // v (E transfer to the nucleus)
    double x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
    double y  = (hitnucl) ? v/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1
    double W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
    double W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
    double t  = 0;

    // Get v 4-p at hit nucleon rest-frame
    TLorentzVector k1_rf = k1;         
    if(hitnucl) {
         k1_rf.Boost(-1.*p1.BoostVector());
    }

//    if(is_mec){
//      v = q.Energy();
//      x = 0.5*Q2/(M*v);
//      y = v/k1.Energy();
//      W2 = M*M + 2*M*v - Q2;
//      W = TMath::Sqrt(W2);
//    }

    LOG("gntpc", pDEBUG) 
         << "[Calc] Q2 = " << Q2 << ", W = " << W 
         << ", x = " << x << ", y = " << y << ", t = " << t;

    // Extract more info on the hadronic system
    // Only for QEL/RES/DIS/COH/MEC events
    //
    bool study_hadsyst = (is_qel || is_res || is_dis || is_coh || is_mec || is_singlek);    
    
    //
    TObjArrayIter piter(&event);
    GHepParticle * p = 0;
    int ip=-1;

    //
    // Extract the final state system originating from the hadronic vertex 
    // (after the intranuclear rescattering step)
    //

    LOG("gntpc", pDEBUG) << "Extracting final state hadronic system";

    vector<int> final_had_syst;
    while( (p = (GHepParticle *) piter.Next()) && study_hadsyst)
    {
        ip++;
        // don't count final state lepton as part hadronic system 
        //if(!is_coh && event.Particle(ip)->FirstMother()==0) continue;
        if(event.Particle(ip)->FirstMother()==0) continue;
        if(pdg::IsPseudoParticle(p->Pdg())) continue;
        int pdgc = p->Pdg();
        int ist  = p->Status();
        if(ist==kIStStableFinalState) {
            if (pdgc == kPdgGamma || pdgc == kPdgElectron || pdgc == kPdgPositron)  {
               int igmom = p->FirstMother();
               if(igmom!=-1) {
               // only count e+'s e-'s or gammas not from decay of pi0
               if(event.Particle(igmom)->Pdg() != kPdgPi0) { final_had_syst.push_back(ip); }
               }
            } else {
               final_had_syst.push_back(ip);
            }
        }
        // now add pi0's that were decayed as short lived particles
        else if(pdgc == kPdgPi0){
    int ifd = p->FirstDaughter();
    int fd_pdgc = event.Particle(ifd)->Pdg();
    // just require that first daughter is one of gamma, e+ or e-  
    if(fd_pdgc == kPdgGamma || fd_pdgc == kPdgElectron || fd_pdgc == kPdgPositron){
        final_had_syst.push_back(ip);
    }
        }
    }//particle-loop

    if( count(final_had_syst.begin(), final_had_syst.end(), -1) > 0) {
        mcrec->Clear();
    continue;
    }

    //
    // Extract info on the primary hadronic system (before any intranuclear rescattering)
    // looking for particles with status_code == kIStHadronInTheNucleus 
    // An exception is the coherent production and scattering off free nucleon targets 
    // (no intranuclear rescattering) in which case primary hadronic system is set to be 
    // 'identical' with the final  state hadronic system
    //

    LOG("gntpc", pDEBUG) << "Extracting primary hadronic system";
    
    ip = -1;
    TObjArrayIter piter_prim(&event);

    vector<int> prim_had_syst;
    if(study_hadsyst) {
        // if coherent or free nucleon target set primary states equal to final states
        if(!pdg::IsIon(target->Pdg()) || (is_coh)) {
         vector<int>::const_iterator hiter = final_had_syst.begin();
         for( ; hiter != final_had_syst.end(); ++hiter) {
             prim_had_syst.push_back(*hiter);
         }
        } 
        // otherwise loop over all particles and store indices of those which are hadrons
        // created within the nucleus
        else {
    while( (p = (GHepParticle *) piter_prim.Next()) ){
        ip++;      
        int ist_comp  = p->Status();
        if(ist_comp==kIStHadronInTheNucleus) {
        prim_had_syst.push_back(ip); 
        }
    }//particle-loop   
    //
    // also include gammas from nuclear de-excitations (appearing in the daughter list of the 
    // hit nucleus, earlier than the primary hadronic system extracted above)
    for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
        if(i<0) continue;
        if(event.Particle(i)->Status()==kIStStableFinalState) { prim_had_syst.push_back(i); }
    }      
        }//freenuc?
    }//study_hadsystem?
    
    if( count(prim_had_syst.begin(), prim_had_syst.end(), -1) > 0) {
        mcrec->Clear();
    continue;
    }

    //
    // Al information has been assembled -- Start filling up the tree branches
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
    brWeight     = weight;      
    brKineXs     = xs;      
    brKineYs     = ys;      
    brKineTs     = ts;      
    brKineQ2s    = Q2s;            
    brKineWs     = Ws;      
    brKineX      = x;      
    brKineY      = y;      
    brKineT      = t;      
    brKineQ2     = Q2;      
    brKineW      = W;      
    brEvRF       = k1_rf.Energy();      
    brEv         = k1.Energy();      
    brPxv        = k1.Px();  
    brPyv        = k1.Py();  
    brPzv        = k1.Pz();  
    brEn         = (hitnucl) ? p1.Energy() : 0;      
    brPxn        = (hitnucl) ? p1.Px()     : 0;      
    brPyn        = (hitnucl) ? p1.Py()     : 0;      
    brPzn        = (hitnucl) ? p1.Pz()     : 0;            
    brEl         = k2.Energy();      
    brPxl        = k2.Px();      
    brPyl        = k2.Py();      
    brPzl        = k2.Pz();      
    brPl         = k2.P();
    brCosthl     = TMath::Cos( k2.Vect().Angle(k1.Vect()) );

    // Primary hadronic system (from primary neutrino interaction, before FSI)
    brNi = prim_had_syst.size();
    for(int j=0; j<brNi; j++) {
        p = event.Particle(prim_had_syst[j]);
        assert(p);
        brPdgi[j] = p->Pdg();     
        brResc[j] = p->RescatterCode();     
        brEi  [j] = p->Energy();     
        brPxi [j] = p->Px();     
        brPyi [j] = p->Py();     
        brPzi [j] = p->Pz();

        LOG("gntpc", pINFO) 
        << "Counting in primary hadronic system: idx = " << prim_had_syst[j]
        << " -> " << p->Name();
    }

    // Final state (visible) hadronic system
    brSumKEf     = fspl->KinE();

    brNf = final_had_syst.size();
    for(int j=0; j<brNf; j++) {
        p = event.Particle(final_had_syst[j]);
        assert(p);

        int    hpdg = p->Pdg();     
        double hE   = p->Energy();     
        double hKE  = p->KinE();     
        double hpx  = p->Px();     
        double hpy  = p->Py();     
        double hpz  = p->Pz();     
        double hp   = TMath::Sqrt(hpx*hpx + hpy*hpy + hpz*hpz);
        double hm   = p->Mass();     
        double hcth = TMath::Cos( p->P4()->Vect().Angle(k1.Vect()) );

        brPdgf  [j] = hpdg;
        brEf    [j] = hE;
        brPxf   [j] = hpx;
        brPyf   [j] = hpy;
        brPzf   [j] = hpz;
        brPf    [j] = hp;
        brCosthf[j] = hcth;

        brSumKEf += hKE;

        LOG("gntpc", pINFO) 
        << "Counting in f/s system from hadronic vtx: idx = " << final_had_syst[j]
        << " -> " << p->Name();
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

  fout.Write();
  fout.Close();
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
  
  // get output file format
  gOptOutFileFormat = kConvFmt_sgst;

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
  string ext="";
  if      (gOptOutFileFormat == kConvFmt_gst                  ) { ext = "gst.root";         }
  else if (gOptOutFileFormat == kConvFmt_sgst                 ) { ext = "sgst.root";        }
  else if (gOptOutFileFormat == kConvFmt_gxml                 ) { ext = "gxml";             }
  else if (gOptOutFileFormat == kConvFmt_ghep_mock_data       ) { ext = "mockd.ghep.root";  }
  else if (gOptOutFileFormat == kConvFmt_rootracker           ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_rootracker_mock_data ) { ext = "mockd.gtrac.root"; }
  else if (gOptOutFileFormat == kConvFmt_t2k_rootracker       ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_numi_rootracker      ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_t2k_tracker          ) { ext = "gtrac.dat";        }
  else if (gOptOutFileFormat == kConvFmt_nuance_tracker       ) { ext = "gtrac_legacy.dat"; }
  else if (gOptOutFileFormat == kConvFmt_ghad                 ) { ext = "ghad.dat";         }
  else if (gOptOutFileFormat == kConvFmt_ginuke               ) { ext = "ginuke.root";      }

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
