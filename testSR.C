#include "postBDTcuts.h"
#include "const.h"

void testSR()
{
	TFile* f = new TFile("mvaout.muons.root","READ");
	TTree* t = (TTree*)f->Get("fltmva_Data");
	

	// Events in the SR:
	//mass=1720.2867, RunNumber=213157, EventNumber=51326767
	//mass=1744.3396, RunNumber=205017, EventNumber=10521463
	t->Scan("score:m3body:mOS1:mOS2:mSS:pval:iso020:iso030:trkspval:lxysig:a0xysig:pt3body:mettrk:metcal:dphihtcal:dphihttrk:mtcal:mttrk:metsdphi:ht:njets","(RunNumber==213157 && EventNumber==51326767) || (RunNumber==205017 && EventNumber==10521463)");
	
	cout << endl;
	
	TString cuts_bkg = postBDTcut("1450","1690","1870","2110","-0.9","0.933","bkg","1713","1841",false,true,"postTraining","mvaout");
	Double_t npassedD = t->GetEntries(cuts_bkg);
	cout << "Tight+x>x1: " << npassedD << endl;
	t->Scan("score:m3body:mOS1:mOS2:mSS:pval:iso020:iso030:trkspval:lxysig:a0xysig:pt3body:mettrk:metcal:dphihtcal:dphihttrk:mtcal:mttrk:metsdphi:ht:njets",cuts_bkg);
}
