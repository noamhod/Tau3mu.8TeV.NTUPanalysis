#include "allFitSystAuto.h"

TString postBDTcut(TString xFullMin, TString xBlindMin, TString xBlindMax, TString xFullMax, TString xBDT, TString sigma, TString type, TString xSRmin="", TString xSRmax="")
{
	TString sidebandscuts = "(m3body>"+xFullMin+". && m3body<"+xFullMax+".)";
	TString blindedcuts   = "(m3body<"+xBlindMin+". || m3body>"+xBlindMax+".)";
	TString srcuts        = (xSRmin!="" && xSRmax!="") ? "(m3body>"+xSRmin+". && m3body<"+xSRmax+".)" : "(1)";
	TString bdtcuts       = "score>"+xBDT;
	TString ncandcuts     = "(@m3body->size()==1)";
	TString onPhi         = "(TMath::Abs(mOS1-1020.)<"+sigma+" || TMath::Abs(mOS2-1020.)<"+sigma+")";
	TString onOmega       = "(TMath::Abs(mOS1-782.)<"+sigma+"  || TMath::Abs(mOS2-782.)<"+sigma+")";
	TString onRho         = "(TMath::Abs(mOS1-770.)<"+sigma+"  || TMath::Abs(mOS2-770.)<"+sigma+")";
	TString onRhoOmegaPhi = "("+onPhi+"||"+onOmega+"||"+onRho+")";
	//TString basecuts      = ncandcuts+" && !"+onRhoOmegaPhi+" && "+sidebandscuts;
	TString basecuts      = ncandcuts+" && "+sidebandscuts;
	
	TString cut = "";
	if     (type=="bkg") cut = basecuts+" && "+bdtcuts+" && "+blindedcuts;
	else if(type=="sig") cut = basecuts+" && "+bdtcuts+" && "+srcuts;
	else _FAT("unknown type: "<<type);
	return cut;
}

void testMVAout()
{
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");

	float mMinSBleft  = 1450;
	float mMaxSBleft  = 1690;
	float mMinSBright = 1870;
	float mMaxSBright = 2290;

	TString sxSRmin    = "1713";
	TString sxSRmax    = "1841";

	TString smMinSBleft  = tstr(mMinSBleft,0);
	TString smMaxSBleft  = tstr(mMaxSBleft,0);
	TString smMinSBright = tstr(mMinSBright,0);
	TString smMaxSBright = tstr(mMaxSBright,0);

	TString basecuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","30","bkg");
	TString basecuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","30","sig",sxSRmin,sxSRmax);

	TString cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9","30","bkg");
	TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9","30","sig",sxSRmin,sxSRmax);

	TEventList* elistS;
	TEventList* elistD;
	Float_t npassedS;
	Float_t npassedD;

	tS->SetEventList(0);
        tS->Draw(">>elistS",basecuts_sig);
        elistS = (TEventList*)gDirectory->Get("elistS");
        npassedS = elistS->GetN();  // number of events to pass cuts
        cout << basecuts_sig << ":  npassedS=" << npassedS << endl;
	elistS->Clear();
	tS->SetEventList(0);
	tS->Draw(">>elistS",cuts_sig);
	elistS = (TEventList*)gDirectory->Get("elistS");
	npassedS = elistS->GetN();  // number of events to pass cuts
	cout << cuts_sig << ":  npassedS=" << npassedS << endl;

	tD->SetEventList(0);
        tD->Draw(">>elistD",basecuts_bkg);
        elistD = (TEventList*)gDirectory->Get("elistD");
        npassedD = elistD->GetN();  // number of events to pass cuts
        cout << basecuts_bkg << ":  npassedD=" << npassedD << endl;
	elistD->Clear();
	tD->SetEventList(0);
	tD->Draw(">>elistD",cuts_bkg);
	elistD = (TEventList*)gDirectory->Get("elistD");
	npassedD = elistD->GetN();  // number of events to pass cuts
	cout << cuts_bkg << ":  npassedD=" << npassedD << endl;
}
