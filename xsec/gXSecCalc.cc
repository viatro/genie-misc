#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "stream_demultiplexer.h"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GEVGDriver.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGCodeList.h"
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"

using namespace std;
using namespace genie;

// Prototypes:
void          GetCommandLineArgs (int argc, char ** argv);
void          PrintSyntax        (void);
PDGCodeList * GetPdgCodes        (string);

// User-specified options:
PDGCodeList * gOptNuPdgCodeList  = 0;
PDGCodeList * gOptTgtPdgCodeList = 0;
long int      gOptRanSeed        = -1;  // random number seed
string        gOptInpXSecFile    = "";  // input cross-section file
ofstream      gOptOutputFile         ;

stream_demultiplexer sout;

int main(int argc, char ** argv) {
    sout.add_stream(cout);
    // Parse command line arguments
    GetCommandLineArgs(argc,argv);
    
    // Init
    utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
    utils::app_init::RandGen(gOptRanSeed);
    utils::app_init::XSecTable(gOptInpXSecFile, false);
    
    PDGCodeList::const_iterator nuiter;
    PDGCodeList::const_iterator tgtiter;
    int nupdg, tgtpdg;
    map<int, map<int, GEVGDriver *> >  driver_storage;
    
    for (tgtiter = gOptTgtPdgCodeList->begin(); tgtiter != gOptTgtPdgCodeList->end(); ++tgtiter) {
        tgtpdg = *tgtiter;
        for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
            nupdg  = *nuiter;
            
            GEVGDriver * driver = new GEVGDriver();
            driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
            InitialState init_state(tgtpdg, nupdg);
            driver->Configure(init_state);
            if (!gOptInpXSecFile.empty()) driver->UseSplines();
            
            driver_storage[tgtpdg][nupdg] = driver;
        }
    }
    
    double Enu = 0;
    for (tgtiter = gOptTgtPdgCodeList->begin(); tgtiter != gOptTgtPdgCodeList->end(); ++tgtiter) {
        tgtpdg = *tgtiter;
        
        cout << "=====================================================================================================================" << endl;
        sout << "#Target = " << tgtpdg;
        for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
            nupdg  = *nuiter;
            sout << TString::Format("\t%-17s", TString::Format("Neutrino = %d", nupdg).Data());
        }
        sout << endl;
        sout << TString::Format("%-20s", "#Energy(neutrino) [MeV]");
        for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
            sout << TString::Format("\t%-17s", "xsec_total [cm^2]");
        }
        sout << endl;
        
        Enu = 0;
        sout << TString::Format("%20.15f", Enu/units::MeV);
        for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
            nupdg  = *nuiter;
            sout << TString::Format("\t%17.15E", driver_storage[tgtpdg][nupdg]->XSecSum(TLorentzVector(0, 0, Enu, Enu))/units::cm2);
        }
        sout << endl;
        
        for (int p = -4000; p < 0; ++p) {
            Enu = TMath::Power(10, double(p)/1000);
            sout << TString::Format("%20.15f", Enu/units::MeV);
            for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
                nupdg  = *nuiter;
                sout << TString::Format("\t%17.15E", driver_storage[tgtpdg][nupdg]->XSecSum(TLorentzVector(0, 0, Enu, Enu))/units::cm2);
            }
            sout << endl;
        }
        
        Enu = 1;
        sout << TString::Format("%20.15f", Enu/units::MeV);
        for (nuiter = gOptNuPdgCodeList->begin(); nuiter != gOptNuPdgCodeList->end(); ++nuiter) {
            nupdg  = *nuiter;
            sout << TString::Format("\t%17.15E", driver_storage[tgtpdg][nupdg]->XSecSum(TLorentzVector(0, 0, Enu, Enu))/units::cm2);
        }
        sout << "\n" << endl;
        
    }
}

void GetCommandLineArgs(int argc, char ** argv) {
    
    // Common run options. Set defaults and read.
    RunOpt::Instance()->EnableBareXSecPreCalc(false); //default: true
    RunOpt::Instance()->ReadFromCommandLine(argc,argv);
    
    // Parse run options for this app
    CmdLnArgParser parser(argc,argv);
    
    if (parser.OptionExists('o')) {
        gOptOutputFile.open(parser.ArgAsString('o').c_str());
        sout.add_stream(gOptOutputFile);
    }
    /*
    // number of knots
    if (parser.OptionExists('n'))  gOptNKnots = parser.ArgAsInt('n');
    else gOptNKnots = -1; // using default
    
    // max spline energy (if < max of validity range)
    if (parser.OptionExists('e'))  gOptMaxE = parser.ArgAsDouble('e');
    else gOptMaxE = -1; // using default
    */
    
    // comma-separated neutrino PDG code list
    if (parser.OptionExists('p')) {
        gOptNuPdgCodeList = GetPdgCodes(parser.ArgAsString('p'));
        if(!gOptNuPdgCodeList || gOptNuPdgCodeList->size() == 0 ) {
            cerr << "Empty neutrino PDG code list" << endl;
            PrintSyntax();
            exit(2);
        }
    } else {
        cerr << "Unspecified neutrino PDG code list - Exiting" << endl;
        PrintSyntax();
        exit(1);
    }
    
    // comma-separated target PDG code list or input geometry file
    if (parser.OptionExists('t')) {
        gOptTgtPdgCodeList = GetPdgCodes(parser.ArgAsString('t'));
        if(!gOptTgtPdgCodeList || gOptTgtPdgCodeList->size() == 0 ) {
            cerr << "Empty target PDG code list" << endl;
            PrintSyntax();
            exit(3);
        }
    } else {
        cerr << "Unspecified target PDG code list - Exiting" << endl;
        PrintSyntax();
        exit(1);
    }
    
    // random number seed
    if (parser.OptionExists("seed"))  gOptRanSeed = parser.ArgAsLong("seed");
    else gOptRanSeed = -1; // using default
    
    // input cross-section file
    if (parser.OptionExists("cross-sections"))  gOptInpXSecFile = parser.ArgAsString("cross-sections");
    else gOptInpXSecFile = "";
    
    //
    // print the command-line options 
    //
    cout << "\n"
        << utils::print::PrintFramedMesg("gxscalc job configuration")
        << "\n Neutrino PDG codes : " << *gOptNuPdgCodeList
        << "\n Target PDG codes : " << *gOptTgtPdgCodeList
        << "\n Input cross-section file : " << gOptInpXSecFile
        << "\n Random number seed : " << gOptRanSeed
        << endl;
    
    cout << *RunOpt::Instance();
}
//____________________________________________________________________________
PDGCodeList * GetPdgCodes(string spdglist) {
    // split the comma separated list
    vector<string> pdgvec = utils::str::Split(spdglist,  ",");
    // fill in the PDG code list
    PDGCodeList * list = new PDGCodeList;
    vector<string>::const_iterator iter;
    for(iter = pdgvec.begin(); iter != pdgvec.end(); ++iter) {
        list->push_back( atoi(iter->c_str()) );
    }
    return list;
}
//____________________________________________________________________________
void PrintSyntax(void) {
    cerr
        << "\n\n" << "Syntax:" << "\n"
        << "   gxscalc -p nupdg -t tgtpdg "
        << " [-n nknots] [-e max_energy] "
        << " [--seed seed_number]"
        << " [--input-cross-section xml_file]"
        << " [--event-generator-list list_name]"
        << " [--message-thresholds xml_file]\n"
        << endl;
}
//____________________________________________________________________________
