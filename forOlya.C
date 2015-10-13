///////////////////////////////////
//// root -b -l -q forOlya.C++ ////
///////////////////////////////////
#include "postBDTcuts.h"

static const float mTrainingMinGlob = 750;
static const float mTrainingMaxGlob = 2500;
static const float mBlindMinGlob = 1690;
static const float mBlindMaxGlob = 1870;
static const float mSideBandLeftLowerMeVGlob  = 1450;
static const float mSideBandRightUpperMeVGlob = 2110;
static const float mSRminMeVGlob = 1713;
static const float mSRmaxMeVGlob = 1841;
static const float optBDTcut = 0.933;
static const float minBDTcut = -0.9;

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}

void forOlya()
{
	TString ntuple = "flatout"; // or "mvaout"
	TFile* f = (ntuple=="flatout") ? new TFile("flatout.periodall+MC.muons.mva.analysis.n0.j0.loose.root","READ") : new TFile("mvaout.muons.root","READ");;
	TTree* tS0 = (TTree*)f->Get((ntuple=="flatout") ? "flatout_Wtaunu_3mu"      : "fltmva_Wtaunu_3mu");
	TTree* tS  = (TTree*)f->Get((ntuple=="flatout") ? "flatout_Wtaunu_200k_3mu" : "fltmva_Wtaunu_200k_3mu");
	TTree* tD  = (TTree*)f->Get((ntuple=="flatout") ? "flatout_Data"            : "fltmva_Data");

	TString smMinTraining  = tstr(mTrainingMinGlob,0);
	TString smMaxTraining  = tstr(mTrainingMaxGlob,0);
	TString smMinSBleft    = tstr(mSideBandLeftLowerMeVGlob,0);
	TString smMaxSBleft    = tstr(mBlindMinGlob,0);
	TString smMinSBright   = tstr(mBlindMaxGlob,0);
	TString smMaxSBright   = tstr(mSideBandRightUpperMeVGlob,0);
	TString sxSRmin        = tstr(mSRminMeVGlob,0);
	TString sxSRmax        = tstr(mSRmaxMeVGlob,0);
	TString soptBDTcut     = tstr(optBDTcut,4);
	TString sminBDTcut     = tstr(minBDTcut,4);
	
	TString cuts_sig = "";
	TString cuts_bkg = "";
	Float_t ninitS   = -1;
	Float_t npassedS = -1;
	Float_t npassedD = -1;
	bool blinded = true; bool unblinded = false;
	bool doBDT   = true; bool noBDT     = false;
	
	//// loose training selection
	cuts_sig = postBDTcut(smMinTraining,smMinSBleft,smMaxSBright,smMaxTraining,"-1","-1","sig", "","", unblinded,noBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinTraining,smMinSBleft,smMaxSBright,smMaxTraining,"-1","-1","bkg", "","", blinded,noBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose+training only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose+training only (SB data): " << npassedD << endl;
	
	//// loose selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", sxSRmin,sxSRmax, unblinded,noBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkg", sxSRmin,sxSRmax, blinded,noBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose only (SB data): " << npassedD << endl;
	
	//// loose+x>x0 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose+x>x0 only (SB data): " << npassedD << endl;
	
	//// tight selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", sxSRmin,sxSRmax, unblinded,noBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkg", sxSRmin,sxSRmax, blinded,noBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight only (SB data): " << npassedD << endl;
	
	//// tight+x>x0 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x0 only (SB data): " << npassedD << endl;
	
	//// tight+x>x1 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x1 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x1 only (SB data): " << npassedD << endl;
}