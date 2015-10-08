#ifndef TYPE_H
#define TYPE_H

#include "std.h"

typedef map<unsigned int, TString>     TMapuiTS;
typedef map<int, TString>              TMapiTS;
typedef map<float, TString>            TMapfTS;
typedef map<int, int>                  TMapii;
typedef map<TString, unsigned int>     TMapTSui;
typedef map<TString, int>              TMapTSi;
typedef map<TString, float>            TMapTSf;
typedef map<TString, TString>          TMapTSTS;
typedef map<TString, string>           TMapTSs;
typedef map<TString, Color_t>          TMapTSTC;
typedef map<TString, Int_t>            TMapTSI;
typedef map<TString, vector<int>* >    TMapTSP2vi;
typedef map<TString, vector<float>* >  TMapTSP2vf;
typedef map<TString, vector<string>* > TMapTSP2vs;
typedef map<TString, TH1*>             TMapTSP2TH1;
typedef map<TString, TH2*>             TMapTSP2TH2;
typedef map<TString, TProfile*>        TMapTSP2TProfile;
typedef map<TString, TFile*>           TMapTSP2TFILE;
typedef map<TString, TTree*>           TMapTSP2TTREE;

#endif