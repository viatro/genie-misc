#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
//#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/CmdLnArgParser.h"

#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/PointGeomAnalyzer.h"

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void ParseFluxHst       (string);
void PrintSyntax        (void);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 1000;    // n-events to generate
NtpMCFormat_t kDefOptNtpFormat = kNFGHEP; // ntuple format
Long_t        kDefOptRunNu     = 0;       // default run number

//User-specified options:
int             gOptNevents;      // n-events to generate
double          gOptNuEnergyMin;  // min neutrino energy in spectrum
double          gOptNuEnergyMax;  // max neutrino energy in spectrum
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
map<int,TH1D*>  gOptFluxHst;      // flux histos (nu pdg  -> spectrum)
Long_t          gOptRunNu;        // run number
bool            gOptWeighted;     // 
long int        gOptRanSeed;      // random number seed
string          gOptInpXSecFile;  // cross-section splines

int main(int argc, char ** argv) {
    GetCommandLineArgs(argc,argv);
    
    // Initialization of random number generators, cross-section table, 
    // messenger thresholds, cache file
    utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
    utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
    utils::app_init::RandGen(gOptRanSeed);
    utils::app_init::XSecTable(gOptInpXSecFile, false);
    
    // Set GHEP print level
    GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
    
    //
    // Generate neutrino events
    
    // Get flux and geom drivers
    GeomAnalyzerI * geom_driver = new geometry::PointGeomAnalyzer(gOptTgtMix);
    
    flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
    
    TVector3 bdir (0,0,1);
    TVector3 bspot(0,0,0);
    
    flux->SetNuDirection      (bdir);
    flux->SetBeamSpot         (bspot);
    flux->SetTransverseRadius (-1);
    
    map<int,TH1D*>::iterator it = gOptFluxHst.begin();
    for (; it != gOptFluxHst.end(); ++it) {
        int    pdg_code = it->first;
        TH1D * spectrum = it->second;
        flux->AddEnergySpectrum(pdg_code, spectrum);
        // once the histogram has been added to the GCylindTH1Flux driver
        // it is owned by the driver and it is up to the the driver
        // to clean up (i.e. delete it).  
        // remove it from this map to avoid double deletion.
        it->second = 0;
    }
    
    GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
    
    // Create the monte carlo job driver
    GMCJDriver * mcj_driver = new GMCJDriver;
    mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
    mcj_driver->SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
    mcj_driver->UseFluxDriver(flux_driver);
    mcj_driver->UseGeomAnalyzer(geom_driver);
    mcj_driver->Configure();
    mcj_driver->UseSplines();
    if (!gOptWeighted) mcj_driver->ForceSingleProbScale();
    
    // Initialize an Ntuple Writer to save GHEP records into a TTree
    NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
    ntpw.Initialize();
    
    // Create an MC Job Monitor
    GMCJMonitor mcjmonitor(gOptRunNu);
    mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());
    
    // Generate events / print the GHEP record / add it to the ntuple
    int ievent = 0;
    while ( ievent < gOptNevents) {
    
        LOG("gevgen_mix", pNOTICE) << " *** Generating event............ " << ievent;
    
        // generate a single event for neutrinos coming from the specified flux
        EventRecord * event = mcj_driver->GenerateEvent();
    
        LOG("gevgen_mix", pNOTICE) << "Generated Event GHEP Record: " << *event;
    
        // add event at the output ntuple, refresh the mc job monitor & clean-up
        ntpw.AddEventRecord(ievent, event);
        mcjmonitor.Update(ievent,event);
        ievent++;
        delete event;
    }
    
    // Save the generated MC events
    ntpw.Save();
    
    //Clean up
    delete flux_driver;
    delete geom_driver;
    delete mcj_driver;
    
    /*map<int,TH1D*>::iterator */it = gOptFluxHst.begin();
    for(; it != gOptFluxHst.end(); ++it) {
        TH1D * spectrum = it->second;
        if(spectrum) delete spectrum;
    }
    gOptFluxHst.clear();
    
    LOG("gevgen_t2k", pNOTICE) << "Done!";
    
    return 0;
}

void GetCommandLineArgs(int argc, char ** argv) {
    LOG("gevgen_mix", pINFO) << "Parsing command line arguments";
    
    // Common run options. Set defaults and read.
    RunOpt::Instance()->EnableBareXSecPreCalc(true);
    RunOpt::Instance()->ReadFromCommandLine(argc,argv);
    
    // Parse run options for this app
    CmdLnArgParser parser(argc,argv);
    
    // help?
    bool help = parser.OptionExists('h');
    if (help) {
        PrintSyntax();
        exit(0);
    }
    
    // number of events
    if( parser.OptionExists('n') ) {
        LOG("gevgen_mix", pINFO) << "Reading number of events to generate";
        gOptNevents = parser.ArgAsInt('n');
    } else {
        LOG("gevgen_mix", pINFO)
        << "Unspecified number of events to generate - Using default";
        gOptNevents = kDefOptNevents;
    }
    
    // run number
    if( parser.OptionExists('r') ) {
        LOG("gevgen_mix", pINFO) << "Reading MC run number";
        gOptRunNu = parser.ArgAsLong('r');
    } else {
        LOG("gevgen_mix", pINFO) << "Unspecified run number - Using default";
        gOptRunNu = kDefOptRunNu;
    }
    
    // neutrino energy
    if( parser.OptionExists('e') ) {
        LOG("gevgen_mix", pINFO) << "Reading neutrino energy";
        string nue = parser.ArgAsString('e');
        
        // split the comma separated list
        vector<string> nurange = utils::str::Split(nue, ",");
        if (nurange.size() != 2) {
            LOG("gevgen_mix", pFATAL) << "Wrong energy range - Exiting";
            PrintSyntax();
            exit(1);
        }
        gOptNuEnergyMin = atof(nurange[0].c_str());
        gOptNuEnergyMax = atof(nurange[1].c_str());
        assert(gOptNuEnergyMax > gOptNuEnergyMin && gOptNuEnergyMin >= 0);
    } else {
        LOG("gevgen_mix", pFATAL) << "Unspecified neutrino energy - Exiting";
        PrintSyntax();
        exit(1);
    }
    
    if ( parser.OptionExists('f') ) {
        LOG("gevgen_mix", pINFO) << "Getting input flux";
        ParseFluxHst(parser.ArgAsString('f'));
    } else {
        LOG("gevgen_mix", pFATAL) << "No flux info was specified - Exiting";
        PrintSyntax();
        exit(1);
    }
    
    // generate weighted events option (only relevant if using a flux)
    gOptWeighted = parser.OptionExists('w');
    
    // target mix (their PDG codes with their corresponding weights)
    if( parser.OptionExists('t') ) {
        LOG("gevgen_mix", pINFO) << "Reading target mix";
        string stgtmix = parser.ArgAsString('t');
        gOptTgtMix.clear();
        vector<string> tgtmix = utils::str::Split(stgtmix,",");
        if (tgtmix.size()==1) {
            int    pdg = atoi(tgtmix[0].c_str());
            double wgt = 1.0;
            gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
        } else {
            vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
            for ( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
                string tgt_with_wgt = *tgtmix_iter;
                string::size_type open_bracket  = tgt_with_wgt.find("[");
                string::size_type close_bracket = tgt_with_wgt.find("]");
                string::size_type ibeg = 0;
                string::size_type iend = open_bracket;
                string::size_type jbeg = open_bracket+1;
                string::size_type jend = close_bracket-1;
                int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend).c_str());
                double wgt = atof(tgt_with_wgt.substr(jbeg,jend).c_str());
                LOG("gevgen_mix", pNOTICE) << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
                gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
            }//tgtmix_iter
        }//>1
    } else {
        LOG("gevgen_mix", pFATAL) << "Unspecified target PDG code - Exiting";
        PrintSyntax();
        exit(1);
    }
    
    // random number seed
    if( parser.OptionExists("seed") ) {
        LOG("gevgen_mix", pINFO) << "Reading random number seed";
        gOptRanSeed = parser.ArgAsLong("seed");
    } else {
        LOG("gevgen_mix", pINFO) << "Unspecified random number seed - Using default";
        gOptRanSeed = -1;
    }
    
    // input cross-section file
    if( parser.OptionExists("cross-sections") ) {
        LOG("gevgen_mix", pINFO) << "Reading cross-section file";
        gOptInpXSecFile = parser.ArgAsString("cross-sections");
    } else {
        LOG("gevgen_mix", pINFO) << "Unspecified cross-section file";
        gOptInpXSecFile = "";
    }
    
    //
    // print-out the command line options
    //
    LOG("gevgen_mix", pNOTICE) << "\n" << utils::print::PrintFramedMesg("gevgen job configuration");
    
    PDGLibrary * pdglib = PDGLibrary::Instance();
    
    LOG("gevgen_mix", pNOTICE) << "Target code (PDG) & weight fraction (in case of multiple targets): ";
    ostringstream tgtinfo;
    map<int,double>::const_iterator titer;
    for (titer = gOptTgtMix.begin(); titer != gOptTgtMix.end(); ++titer) {
        int    tgtpdgc = titer->first;
        double wgt     = titer->second;
        tgtinfo << tgtpdgc << " (weight fraction = " << wgt << ") | ";
    }
    
    ostringstream fluxinfo;
    map<int,TH1D*>::const_iterator hiter;
    for(hiter = gOptFluxHst.begin(); hiter != gOptFluxHst.end(); ++hiter) {
          int    pdg_code = hiter->first;
          TH1D * spectrum = hiter->second;
          TParticlePDG * p = pdglib->Find(pdg_code);
          if(p) {
            string name = p->GetName();
            fluxinfo << "(" << name << ") -> " << spectrum->GetName() << " | ";
          }//p?
    }
    
    LOG("gevgen_mix", pNOTICE) 
        << "\n - Run number: " << gOptRunNu
        << "\n - Number of events requested: " << gOptNevents
        << "\n - Neutrino energy range: [" << gOptNuEnergyMin << ", " << gOptNuEnergyMax << "] GeV"
        << "\n - Generate weighted events? " << gOptWeighted
        << "\n - Random number seed: " << gOptRanSeed
        << "\n - Using cross-section file: " << gOptInpXSecFile
        << "\n - Target   @ " << tgtinfo.str()
        << "\n - Flux     @ " << fluxinfo.str();
    
    LOG("gevgen_mix", pNOTICE) << "\n";
    
    LOG("gevgen_mix", pNOTICE) << *RunOpt::Instance();
    
}

//____________________________________________________________________________

void ParseFluxHst(string flux) {
    // Using flux from histograms
    
    vector<string> fluxv = utils::str::Split(flux,",");
    string fluxfile = fluxv[0];
    bool accessible_flux_file = !(gSystem->AccessPathName(fluxfile.c_str()));
    if (!accessible_flux_file) {
        LOG("gevgen_mix", pFATAL) 
        << "Can not access flux file: " << fluxfile;
        PrintSyntax();
        exit(1);
    }
    // Extract energy spectra for all specified neutrino species
    TFile flux_file(fluxfile.c_str(), "read");
    for(unsigned int inu=1; inu<fluxv.size(); ++inu) {
        string nutype_and_histo = fluxv[inu];
        string::size_type open_bracket  = nutype_and_histo.find("[");
        string::size_type close_bracket = nutype_and_histo.find("]");
        if (open_bracket == string::npos || close_bracket==string::npos) {
            LOG("gevgen_mix", pFATAL) << "You made an error in specifying the flux histograms"; 
            PrintSyntax();
            exit(1);
        }
        string::size_type ibeg = 0;
        string::size_type iend = open_bracket;
        string::size_type jbeg = open_bracket+1;
        string::size_type jend = close_bracket;
        string nutype = nutype_and_histo.substr(ibeg,iend-ibeg);
        string histo  = nutype_and_histo.substr(jbeg,jend-jbeg);
        // access specified histogram from the input root file
        TH1 * ihst = (TH1*) flux_file.Get(histo.c_str()); 
        if (!ihst) {
            LOG("gevgen_mix", pFATAL) << "Can not find histogram: " << histo << " in flux file: " << fluxfile;
            PrintSyntax();
            exit(1);
        }
        
        // create a local copy of the input histogram
        TString origname = ihst->GetName();
        TString tmpname; tmpname.Form("%s_", origname.Data());
        TH1D * spectrum = new TH1D( tmpname.Data(), ihst->GetTitle(), ihst->GetNbinsX(),  
                                    ihst->GetXaxis()->GetXmin(), ihst->GetXaxis()->GetXmax());
        spectrum->SetDirectory(0);
        for (int ibin = 1; ibin <= ihst->GetNbinsX(); ++ibin) {
            if (ihst->GetBinLowEdge(ibin) < gOptNuEnergyMax && ihst->GetBinLowEdge(ibin) + ihst->GetBinWidth(ibin) > gOptNuEnergyMin) {
                spectrum->SetBinContent(ibin, ihst->GetBinContent(ibin));
                LOG("gevgen_mix", pDEBUG) << "adding => " << ibin << ": " << ihst->GetBinContent(ibin);
            }
        }
        // get rid of original
        delete ihst;
        // rename copy
        spectrum->SetName(origname.Data());
        
        // convert neutrino name -> pdg code
        int pdg = atoi(nutype.c_str());
        if (!pdg::IsNeutrino(pdg) && !pdg::IsAntiNeutrino(pdg)) {
            LOG("gevgen_mix", pFATAL) << "Unknown neutrino type: " << nutype; 
            PrintSyntax();
            exit(1);
        }
        // store flux neutrino code / energy spectrum
        LOG("gevgen_mix", pDEBUG) << "Adding energy spectrum for flux neutrino: pdg = " << pdg;
        gOptFluxHst.insert(map<int, TH1D*>::value_type(pdg, spectrum));
    }//inu
    
    if (gOptFluxHst.size() < 1) {
        LOG("gevgen_mix", pFATAL) << "You have not specified any flux histogram!";
        PrintSyntax();
        exit(1);
    }
    
    flux_file.Close();
}

void PrintSyntax(void) {
    LOG("gevgen_mix", pNOTICE)
        << "\n\n" << "Syntax:" << "\n"
        << "\n      gevgen_mix [-h]"
        << "\n              [-r run#]"
        << "\n               -n nev"
        << "\n               -e energy range"
        << "\n               -t target_pdg "
        << "\n               -f flux"
        << "\n              [-w]"
        << "\n              [--seed random_number_seed]"
        << "\n              [--cross-sections xml_file]"
        << "\n              [--event-generator-list list_name]"
        << "\n              [--message-thresholds xml_file]"
        << "\n              [--unphysical-event-mask mask]"
        << "\n              [--event-record-print-level level]"
        << "\n              [--mc-job-status-refresh-rate  rate]"
        << "\n              [--cache-file root_file]"
        << "\n";
}