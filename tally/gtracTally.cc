#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include <vector>
#include <algorithm>
#include <set>
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
    if (pdg < 1000000000) {
        if (!TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetPDG()) return "unknown";
        return TString(TParticle(pdg,0,0,0,0,0,0,0,0,0,0,0,0,0).GetName());
    }
    else {
        if (pdg == 2000000002) return "HadrBlob";
        Int_t ionA = (pdg / 10) % 1000;
        Int_t ionZ = (pdg / 10000) % 1000;
        if (ionZ == 0) return (ionA == 1) ? "neutron" : TString::Itoa(ionA,10) + TString(".neutron");
        return elementName[ionZ] + TString::Itoa(ionA,10);//TString::Format("%s%d",elementName[ionZ],ionA);
    }
}

Double_t PdgCharge(Int_t pdg) {
    if (pdg < 1000000000) {
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

Int_t PdgZAtoIon(Int_t Z, Int_t A) {
    return 1000000000 + 10000*Z + 10*A;
}

/*Int_t PdgZAtoIon(Double_t Z, Int_t A) {
    return 1000000000 + 10000*(Int_t)round(Z) + 10*A;
}*/

//typedef map<Int_t, UInt_t> fParticles_t;

typedef struct Reaction_t {
    vector<Int_t>      iParticles;
    //multiset<Int_t>    iParticles;
    map<Int_t, UInt_t> fParticles;
    TString AsStringPdgs() const {
        TString s = "";
        for (auto ip: iParticles) s += TString::Format("%d + ", ip);
        s.Remove(s.Length() - 2, 2);
        s += "--> ";
        for (auto p: fParticles) {
            if (p.second == 1) s += TString::Format("%d + ", p.first);
            else s += TString::Format("%d*%d + ", p.second, p.first);
        }
        s.Remove(s.Length() - 2, 2);
        return s;
    }
    TString AsStringNames() const {
        TString s = "";
        for (auto ip: iParticles) s += PdgName(ip) + " + ";
        s.Remove(s.Length() - 2, 2);
        s += "--> ";
        for (auto p: fParticles) {
            if (p.second == 1) s += PdgName(p.first) + " + ";
            else s += TString::Itoa(p.second,10) + "." + PdgName(p.first) + " + ";
        }
        s.Remove(s.Length() - 2, 2);
        return s;
    }
    void Reset() {
        iParticles.clear();
        fParticles.clear();
    }
} Reaction_t;

inline bool operator < (const Reaction_t &l, const Reaction_t &r) { // !!! To use struct as key in map
    if (l.iParticles.size() < r.iParticles.size()) return true;
    if (l.iParticles.size() > r.iParticles.size()) return false;
    else for (size_t i = 0; i < l.iParticles.size(); i++) {
        if (l.iParticles.at(i) < r.iParticles.at(i)) return true;
        if (l.iParticles.at(i) > r.iParticles.at(i)) return false;
    }
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

/*bool operator < (const Reaction_t &l, const Reaction_t &r) { // !!! To use struct as key in map
    if (l.iParticles.size() < r.iParticles.size()) return true;
    if (l.iParticles.size() > r.iParticles.size()) return false;
    else for (int i = 0; i < l.iParticles.size(); i++) {
        if (l.iParticles.at(i) < r.iParticles.at(i)) return true;
        if (l.iParticles.at(i) > r.iParticles.at(i)) return false;
    }
    return (l.fParticles < r.fParticles);
}*/

typedef map<Reaction_t, UInt_t> TallyReactionsMap_t;

TallyReactionsMap_t GetTallyReactionsMap(TTree* gRooTracker) {
    
    Long64_t nentries = gRooTracker->GetEntries();
        
    TallyReactionsMap_t tally;
       
    Reaction_t reaction;
    
    Int_t           StdHepN;
    Int_t           StdHepPdg[250];
    Int_t           StdHepStatus[250];
/*    Int_t           StdHepFd[100];
    Int_t           StdHepLd[100];
    Int_t           StdHepFm[100];
    Int_t           StdHepLm[100];*/
    
    gRooTracker->SetBranchAddress("StdHepN",      &StdHepN    );
    gRooTracker->SetBranchAddress("StdHepPdg",    StdHepPdg   );
    gRooTracker->SetBranchAddress("StdHepStatus", StdHepStatus);
/*    gRooTracker->SetBranchAddress("StdHepFd",     StdHepFd    );
    gRooTracker->SetBranchAddress("StdHepLd",     StdHepLd    );
    gRooTracker->SetBranchAddress("StdHepFm",     StdHepFm    );
    gRooTracker->SetBranchAddress("StdHepLm",     StdHepLm    );*/
    
    Double_t charge;
    Int_t baryon_number;
    
    for (Long64_t i = 0; i < nentries; ++i) {
        
        gRooTracker->GetEntry(i);
        reaction.Reset();
        charge = 0;
        baryon_number = 0;
        
        for (Int_t j = 0; j < StdHepN; ++j) {
            if (StdHepStatus[j] == 0) {
                reaction.iParticles.push_back(StdHepPdg[j]);
                charge += PdgCharge(StdHepPdg[j]);
                baryon_number += PdgBaryonNumber(StdHepPdg[j]);
            }
            else if (StdHepStatus[j] == 1) {
                if (StdHepPdg[j] >= 2000000000) continue;
                charge -= PdgCharge(StdHepPdg[j]);
                baryon_number -= PdgBaryonNumber(StdHepPdg[j]);
                ++(reaction.fParticles[StdHepPdg[j]]);
            }
            else if (StdHepStatus[j] == 15) {
                ++(reaction.fParticles[PdgZAtoIon(charge,baryon_number)]);
            }
            
            /*
            if ((StdHepPdg[j] < 2000000000) || (StdHepPdg[j] == 2000000002)) {
                if (StdHepStatus[j] == 0) reaction.iParticles.push_back(StdHepPdg[j]);
                //if (StdHepStatus[j] == 0) reaction.iParticles.insert(StdHepPdg[j]);
                else if ((StdHepStatus[j] == 1) || (StdHepStatus[j] == 15)) (reaction.fParticles[StdHepPdg[j]])++;
            }
            */
        }
        ++(tally[reaction]);
if ((i+1) % 100000 == 0) cout << "\r" << fixed << setprecision(1) << i/1e6 << " M " << flush;//{printf("\r%.1f M",i/1.e6); fflush(stdout);}
    }
cout << "\r";// << endl;
    return tally;
}

typedef multimap<UInt_t, Reaction_t,greater<UInt_t>> SortedStringsTallyReactionsMap_t;

void PrintSortedTally(TallyReactionsMap_t& tally) {
    SortedStringsTallyReactionsMap_t sortedtally;
    UInt_t n = 0;
    for (auto p: tally) sortedtally.insert(pair<UInt_t, Reaction_t>(p.second,p.first));
    for (auto p: sortedtally) n += p.first;
    cout << setfill(' ') << setw(log10(n)+1) << right << n << " : " << setfill(' ') << setw(10) << right << fixed << setprecision(6) << 100. << " % :" << endl;
    for (auto p: sortedtally) {
        cout << setfill(' ') << setw(log10(n)+1) << right << p.first << " : " << setfill(' ') << setw(10) << right << fixed << setprecision(6) << p.first*100./n << " % :\t" << p.second.AsStringNames() << endl;
    }
}

void ViewHelp() {
    cout << "USAGE: ./gtracTally FileName.root" << endl;
}

void gtracTally(TString filename) {
    
    TFile* file = TFile::Open(filename);
    TTree *gRooTracker = (TTree*)file->Get("gRooTracker");
    TallyReactionsMap_t tally = GetTallyReactionsMap(gRooTracker);
    PrintSortedTally(tally);
    
//cout << "Tally() done."  << endl;
}

#ifndef __CINT__
int main(int argc, char *argv[]) {
    if ( (argc != 2) || !TString(argv[1]).EndsWith(".root") ) 
        ViewHelp();
    else
        gtracTally(TString(argv[1]));
}
#endif