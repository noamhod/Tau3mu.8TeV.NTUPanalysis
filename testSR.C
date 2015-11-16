#include "postBDTcuts.h"
#include "const.h"

void testSR()
{
	TFile* f1 = new TFile("mvaout.muons.root","READ");
	TTree* t1 = (TTree*)f1->Get("fltmva_Data");
	
	TFile* f2 = new TFile("flatout.periodall+MC.muons.mva.analysis.n0.j0.loose.root","READ");
	TTree* t2 = (TTree*)f2->Get("flatout_Data");
	
	TString cuts_bkg = "";

	// Events in the SR:
	//mass=1720.2867, RunNumber=213157, EventNumber=51326767
	//mass=1744.3396, RunNumber=205017, EventNumber=10521463
	t1->Scan("score:m3body:mOS1:mOS2:mSS:pval:iso020:iso030:trkspval:lxysig:a0xysig:pt3body:mettrk:metcal:dphihtcal:dphihttrk:mtcal:mttrk:metsdphi:ht:njets","(RunNumber==213157 && EventNumber==51326767) || (RunNumber==205017 && EventNumber==10521463)");
	cuts_bkg = postBDTcut("1450","1690","1870","2110","-0.9","0.933","bkg","1713","1841",false,true,"postTraining","mvaout");
	t1->Scan("RunNumber:EventNumber:lbn:score:m3body:mOS1:mOS2:mSS:pval:iso020:iso030:trkspval:lxysig:a0xysig:pt3body:mettrk:metcal:dphihtcal:dphihttrk:mtcal:mttrk:metsdphi:ht:njets",cuts_bkg);
	
	cuts_bkg = postBDTcut("1450","1690","1870","2110","-0.9","0.933","bkg","1713","1841",false,true,"postTraining","");
	t2->Scan("evt_RunNumber:evt_EventNumber:evt_lbn:mva_score:vtx_mass:vtx_mOS1:vtx_mOS2:vtx_mSS:vtx_pval:vtx_isolation020:vtx_isolation030:trks_fitprob:geo_lxySig:geo_a0xySig:vtx_pt:met_track_et:met_muons_et:ht_dphimet_muons:ht_dphimet_track:met_muons_mT:met_track_mT:mets_dphi:ht_pt:vtx_pvNtrk:vtx_dRmax:vtx_dRmin:jets_n:vtx_code:mu_pt1:mu_pt2:mu_pt3:mu_eta1:mu_eta2:mu_eta3:mu_phi1:mu_phi2:mu_phi3:mu_pbalsig1:mu_pbalsig2:mu_pbalsig3:mu_sctangsig1:mu_sctangsig2:mu_sctangsig3:mu_pvaltrkfit1:mu_pvaltrkfit2:mu_pvaltrkfit3:mu_type1:mu_type2:mu_type3:mu_trkqoverp1:mu_trkqoverp2:mu_trkqoverp3:mu_srcqoverp1:mu_srcqoverp2:mu_srcqoverp3",cuts_bkg);	
}
