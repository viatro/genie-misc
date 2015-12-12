#include "TROOT.h"
#include "TFile.h"
//#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include <vector>
#include <algorithm>
//#include <set>
#include <map>
#include <iostream>
#include <iomanip>
//#include <cmath>

using namespace std;

TString elementName[] = {"",
    "H",                                                                                                 "He", 
    "Li", "Be",                                                            "B",  "C",  "N",  "O",  "F",  "Ne", 
    "Na", "Mg",                                                            "Al", "Si", "P",  "S",  "Cl", "Ar", 
    "K",  "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo","Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", 
    "Cs", "Ba", 
                "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
                "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
    "Fr", "Ra", 
                "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", 
                "Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"
};

TString PdgName(Int_t pdg) {
    if (pdg < 1000000000 && pdg > -1000000000) {
        if (!TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetPDG()) return TString::Itoa(pdg, 10);
        return TString(TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetName());
    }
    else if (pdg == 2000000002) return "HadrBlob";
    else if (pdg <= -1000000000) return TString("!!! ANTI-") + PdgName(-pdg) +  TString(" !!!");
    else {
        Int_t ionA = (pdg / 10) % 1000;
        Int_t ionZ = (pdg / 10000) % 1000;
        if (ionZ == 0) {
            if (ionA > 1) return TString::Itoa(ionA,10) + TString(".neutron");
            else if (ionA == 1) return "neutron";
            else return "";
        }
        if (ionZ>ionA) return TString(" !!! ") + elementName[ionZ] + TString::Itoa(ionA,10) + TString(" (Z==") + TString::Itoa(ionZ,10) + TString(" > A==") + TString::Itoa(ionA,10) + TString(") !!! ");
        else return elementName[ionZ] + TString::Itoa(ionA,10);//TString::Format("%s%d",elementName[ionZ],ionA);
    }
}

Double_t PdgCharge(Int_t pdg) {
    if (pdg < 1000000000 && pdg > -1000000000) {
        TParticlePDG *pPDG = TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetPDG();
        if (!pPDG) return 0.;
        return pPDG->Charge()/3;
    }
    else return ((pdg / 10000) % 1000);
}

TString PdgClass(Int_t pdg) {
    if (pdg < 1000000000) {
        TParticlePDG *pPDG = TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetPDG();
        if (!pPDG) return "NotInDB";
        return pPDG->ParticleClass();
    }
    if (pdg >= 2000000000) return "GeniePseudoParticle";
    else return "Ion";
}

Int_t PdgBaryonNumber(Int_t pdg) {
    if (pdg < 1000000000) {
        TParticlePDG *pPDG = TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetPDG();
        if (!pPDG) return 0;
        if ( ! TString(pPDG->ParticleClass()).EqualTo("Baryon") ) return 0;
        if (pdg > 0) return  1;
        if (pdg < 0) return -1;
        return 0;
    }
    if (pdg >= 2000000000) return 0;
    else return (pdg / 10) % 1000;
}

Long64_t PdgZAtoIon(Int_t Z, Int_t A) {
    Long64_t pdg = 1000000000L + 10000*abs(Z) + 10*A;
    return Z>=0 ? pdg : -pdg;
}

typedef struct Reaction_t {
    Int_t              neu = 0;
    Int_t              tgt = 0;
    Int_t              fspl = 0;
    map<Int_t, UInt_t> fParticles;
    Int_t              nuclear_remnant = 0;
    
    TString AsStringPdgs() const {
        TString s = "";
        s += TString::Format("%d + %d --> %d + ", neu, tgt, fspl);
        for (auto p: fParticles) {
            if (p.second == 1) s += TString::Format("%d + ", p.first);
            else s += TString::Format("%d*%d + ", p.second, p.first);
        }
        if (nuclear_remnant == 0) s.Remove(s.Length() - 2, 2);
        else s += TString::Itoa(nuclear_remnant, 10);
        return s;
    }
    TString AsStringNames() const {
        TString s = "";
        s += PdgName(neu) + " + " + PdgName(tgt) + " --> " + PdgName(fspl) + " + ";
        for (auto p: fParticles) {
            if (p.second == 1) s += PdgName(p.first) + " + ";
            else s += TString::Itoa(p.second,10) + "." + PdgName(p.first) + " + ";
        }
        if (nuclear_remnant == 0) s.Remove(s.Length() - 2, 2);
        else s += PdgName(nuclear_remnant);
        return s;
    }
    void Reset() {
        neu = 0;
        tgt = 0;
        fspl = 0;
        fParticles.clear();
        nuclear_remnant = 0;
    }
} Reaction_t;

inline bool operator < (const Reaction_t &l, const Reaction_t &r) { // !!! To use struct as key in map
    if (l.neu  < r.neu ) return true;
    if (l.neu  > r.neu ) return false;
    if (l.tgt  < r.tgt ) return true;
    if (l.tgt  > r.tgt ) return false; 
    if (l.fspl < r.fspl) return true;
    if (l.fspl > r.fspl) return false;
    if (l.fParticles.size() < r.fParticles.size()) return true;
    if (l.fParticles.size() > r.fParticles.size()) return false;
    for (auto lit = l.fParticles.begin(), rit = r.fParticles.begin(); lit != l.fParticles.end(); ++lit, ++rit) {
        if (lit->first < rit->first) return true;
        if (lit->first > rit->first) return false;
        if (lit->second < rit->second) return true;
        if (lit->second > rit->second) return false;
    }
    return false; 
}

typedef map<Reaction_t, UInt_t> TallyReactionsMap_t;

TallyReactionsMap_t GetTallyReactionsMap(TChain* gst) {
    
    Long64_t nentries = gst->GetEntries();
        
    TallyReactionsMap_t tally;
       
    Reaction_t reaction;
    
    Int_t           neu;
    Int_t           fspl;
    Int_t           tgt;
    Int_t           Z;
    Int_t           A;
    Int_t           nf;
    Int_t           pdgf[250];
    
    gst->SetBranchStatus("*",0);
    gst->SetBranchStatus("neu", 1);
    gst->SetBranchStatus("fspl",1);
    gst->SetBranchStatus("tgt", 1);
    gst->SetBranchStatus("Z",   1);
    gst->SetBranchStatus("A",   1);
    gst->SetBranchStatus("nf",  1);
    gst->SetBranchStatus("pdgf*",1);
    
    gst->SetBranchAddress("neu",  &neu );
    gst->SetBranchAddress("fspl", &fspl);
    gst->SetBranchAddress("tgt",  &tgt );
    gst->SetBranchAddress("Z",    &Z   );
    gst->SetBranchAddress("A",    &A   );
    gst->SetBranchAddress("nf",   &nf  );
    gst->SetBranchAddress("pdgf", pdgf );
    
    Double_t charge;
    Int_t baryon_number;
    
    for (Long64_t i = 0; i < nentries; ++i) {
        
        gst->GetEntry(i);
        reaction.Reset();
        charge = Z;
        baryon_number = A;
        
        reaction.neu  = neu;
        reaction.tgt  = tgt;
        reaction.fspl = fspl;
        charge -= PdgCharge(fspl);
        
        for (Int_t j = 0; j < nf; ++j) {
            ++(reaction.fParticles[pdgf[j]]);
            charge -= PdgCharge(pdgf[j]);
            baryon_number -= PdgBaryonNumber(pdgf[j]);
        }
        
        //++(reaction.fParticles[PdgZAtoIon(charge,baryon_number)]);
        //if ( charge == 0 && baryon_number == 0 ) reaction.nuclear_remnant = 0;
        /*else */if ( reaction.fParticles.size() == 0 ) reaction.nuclear_remnant = neu;
        else if ( ! (charge == 0 && baryon_number == 0) ) reaction.nuclear_remnant = PdgZAtoIon(charge,baryon_number);
        
        ++(tally[reaction]);
if ((i+1) % 100000 == 0) cout << "\r" << fixed << setprecision(1) << i/1e6 << " M " << flush;//{printf("\r%.1f M",i/1.e6); fflush(stdout);}
    }
cout << "\r";// << endl;
    return tally;
}

void PrintSortedTally(TallyReactionsMap_t& tally) {
    multimap<UInt_t, Reaction_t,greater<UInt_t>> sortedtally;
    UInt_t n = 0;
    for (auto p: tally) sortedtally.insert(pair<UInt_t, Reaction_t>(p.second,p.first));
    for (auto p: sortedtally) n += p.first;
    UInt_t cumul = 0;
    
    cout << setfill(' ') << setw(log10(n)+1) << right << n << " : "
         << setfill(' ') << setw(10) << right << fixed << setprecision(6) << 100. << " % : "
         << setfill(' ') << setw(10) << right << fixed << setprecision(6) << 100. << " %"
         << endl;
    for (auto p: sortedtally) {
        cumul += p.first;
        cout << setfill(' ') << setw(log10(n)+1) << right << p.first << " : "
             << setfill(' ') << setw(10) << right << fixed << setprecision(6) << p.first*100./n << " % : "
             << setfill(' ') << setw(10) << right << fixed << setprecision(6) << cumul*100./n << " % :  \t"
             << p.second.AsStringNames() << endl;
    }
    map<Int_t,UInt_t> fp_counter;
    for (auto s: sortedtally) {
        auto reaction = s.second;
        fp_counter[reaction.fspl] += s.first;
        fp_counter[reaction.nuclear_remnant] += s.first;
        for (auto p : reaction.fParticles) {
            fp_counter[p.first] += p.second * s.first;
        }
    }
    cout << "Counter of final particles:" << endl;
    for (auto p : fp_counter) cout << setfill(' ') << setw(20) << right << PdgName(p.first) << " : " << setfill(' ') << setw(log10(n)+1) << p.second << endl;
}

void ViewHelp() {
    cout << "USAGE: ./gstTally FileName.root [FileName_2.root FileName_3.root ...]" << endl;
}

void gstTally(vector<TString>& filenames) {
    
    //TFile* file = TFile::Open(filename);
    //TTree *gst = (TTree*)file->Get("gst");
    TChain *gst = new TChain("gst");
    for (TString & filename : filenames) gst->Add(filename);
    
    TallyReactionsMap_t tally = GetTallyReactionsMap(gst);
    PrintSortedTally(tally);
    
//cout << "Tally() done."  << endl;
}

#ifndef __CINT__
int main(int argc, char *argv[]) {
    if (argc < 2) {
        ViewHelp();
        return 1;
    }
    vector<TString> filenames;
    TString arg = "";
    for (int i = 1; i < argc; ++i) {
        arg = TString(argv[i]);
        if ( arg.EndsWith(".root") ) filenames.push_back(arg);
        else if (arg == "--help" || arg == "-h" || arg == "-?" || arg == "-help") {
            ViewHelp();
            return 0;
        } else {
            cout << "Unexpected argument : " << arg << endl;
            ViewHelp();
            return 1;
        }
    }
    gstTally(filenames);
    return 0;
}
#endif