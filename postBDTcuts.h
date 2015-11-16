#ifndef POSTBDTCUTS_H
#define POSTBDTCUTS_H

#include "std.h"
#include "type.h"

int visual = 0;

TString postBDTcut(TString xFullMin, TString xBlindMin, TString xBlindMax, TString xFullMax, TString xBDTMin, TString xBDTOpt, TString type, TString xSRmin="", TString xSRmax="", bool isBlinded=1, bool doBDT=1, TString mode="postTraining", TString ntuple="mvaout", TString jetmetword="", TString trkmetword="", bool resVeto=true)
{
	TString pt3body   = (ntuple=="mvaout") ? "pt3body"   : "vtx_pt";
	TString m3body    = (ntuple=="mvaout") ? "m3body"    : "vtx_mass";
	TString mOS1      = (ntuple=="mvaout") ? "mOS1"      : "vtx_mOS1";
	TString mOS2      = (ntuple=="mvaout") ? "mOS2"      : "vtx_mOS2";
	TString mSS       = (ntuple=="mvaout") ? "mSS"       : "vtx_mSS";
	TString score     = (ntuple=="mvaout") ? "score"     : "mva_score";
	TString pval      = (ntuple=="mvaout") ? "pval"      : "vtx_pval";
	TString iso003    = (ntuple=="mvaout") ? "iso003"    : "vtx_isolation003";
	TString iso020    = (ntuple=="mvaout") ? "iso020"    : "vtx_isolation020";
	TString iso030    = (ntuple=="mvaout") ? "iso030"    : "vtx_isolation030";
	TString trkspval  = (ntuple=="mvaout") ? "trkspval"  : "trks_fitprob";
	TString lxysig    = (ntuple=="mvaout") ? "lxysig"    : "geo_lxySig";
	TString a0xysig   = (ntuple=="mvaout") ? "a0xysig"   : "geo_a0xySig";
	TString mettrk    = (ntuple=="mvaout") ? "mettrk"    : "met_track_et"+trkmetword;
	TString metcal    = (ntuple=="mvaout") ? "metcal"    : "met_muons_et"+jetmetword;
	TString mttrk     = (ntuple=="mvaout") ? "mttrk"     : "met_track_mT"+trkmetword;
	TString mtcal     = (ntuple=="mvaout") ? "mtcal"     : "met_muons_mT"+jetmetword;
	TString dphihttrk = (ntuple=="mvaout") ? "dphihttrk" : "ht_dphimet_track"+trkmetword;
	TString dphihtcal = (ntuple=="mvaout") ? "dphihtcal" : "ht_dphimet_muons"+jetmetword;
	TString njets     = (ntuple=="mvaout") ? "njets"     : "jets_n"+jetmetword;;
	///////////////////////////////////////////////////////////////////////////
	
	TString isLoose       = (ntuple=="mvaout") ? "(1)" : "(pass_loose==1)";
	
	TString ncandcuts     = "(@"+m3body+"->size()==1)";
	TString sidebandscuts = "("+m3body+">"+xFullMin+" && "+m3body+"<"+xFullMax+")";
	
	TString bdtmincuts    = score+">"+xBDTMin;
	TString bdtoptcuts    = score+">"+xBDTOpt;
	TString BDTcuts       = (doBDT) ? "("+bdtmincuts+" && "+bdtoptcuts+")" : "(1)";
	
	TString iso           = "("+iso020+"<0.3 && "+iso030+"<1)";
	TString trksprob      = "("+trkspval+">1.e-9)";
	TString distances     = "("+lxysig+">-10 && "+lxysig+"<50 && "+a0xysig+"<25)";
	TString kinematics    = "("+pt3body+">10000 && ("+mettrk+">10000 && "+mettrk+"<250000 && "+mttrk+">20000) && ("+metcal+">10000 && "+metcal+"<250000 && "+mtcal+">20000))";
	TString lowmass1      = "("+mOS1+">220 && "+mOS2+">220 && "+mSS+">220)";
	TString preTraining   = "("+isLoose+" && "+iso+" && "+trksprob+" && "+distances+" && "+kinematics+" && "+lowmass1+")";	
	
	
	
	TString pvallxytrks   = "("+lxysig+">1 && "+trkspval+">8.e-9 && "+pval+">0.2)";
	TString MTHT          = "("+mttrk+">45000 && "+mtcal+">45000 && "+dphihttrk+">2 && "+dphihtcal+">2)";
	// TString passNjets     = "("+njets+"<2)";
	TString lowmass2      = "("+mOS1+">300 && "+mOS2+">300 && "+mSS+">300)";
	TString resonances    = "!((TMath::Abs("+mOS1+"-1020)<50 || TMath::Abs("+mOS2+"-1020)<50 || TMath::Abs("+mOS1+"-782)<50 || TMath::Abs("+mOS2+"-782)<50) && ("+mettrk+"<35000 || "+metcal+"<35000 || "+pt3body+"<35000))";
	TString passDs        = "!((TMath::Abs("+mOS1+"-1020)<50 || TMath::Abs("+mOS2+"-1020)<50) && TMath::Abs("+m3body+"-1968)<100)"; // 100 because of 2*3body_resolution=~65 + margins because of the pi-mu mass differernce
	// TString postTraining  = pvallxytrks+" && "+MTHT+" && "+passNjets+" && "+resonances+" && "+lowmass2+" && "+passDs;
	TString postTraining  = (resVeto) ? (pvallxytrks+" && "+MTHT+" && "+resonances+" && "+lowmass2+" && "+passDs) : (pvallxytrks+" && "+MTHT+" && "+lowmass2);
	
	TString basecuts = ncandcuts+" && "+sidebandscuts+" && "+BDTcuts;
	if     (mode=="preTraining")  basecuts+=" && "+preTraining;
	else if(mode=="postTraining") basecuts+=" && "+preTraining+" && "+postTraining;
	else _FAT("unknown mode="<<mode);

	TString srcuts        = (xSRmin!="" && xSRmax!="") ? "("+m3body+">"+xSRmin+" && "+m3body+"<"+xSRmax+")" : "(1)";
	TString blindedcuts   = (isBlinded) ? "("+m3body+"<"+xBlindMin+" || "+m3body+">"+xBlindMax+")" : "(1)";

	TString cut = basecuts;
	if     (type=="bkg")   cut += " && "+blindedcuts;
	else if(type=="sig")   cut += " && "+srcuts;
	else if(type=="bkgsr") cut += " && "+srcuts;
	else if(type=="full")  cut += "";
	else _FAT("unknown type: "<<type);
	return cut;
}


bool passPostBDTcut(unsigned int vtx, TMapTSP2vf& vf, TMapTSP2vi& vi, float xFullMin, float xBlindMin, float xBlindMax, float xFullMax, float xBDTMin, float xBDTOpt, TString type, float xSRmin=-1, float xSRmax=-1, bool isBlinded=1, bool doBDT=1, TString mode="postTraining", bool resVeto=true)
{
	float pt3body   = vf["pt3body"]->at(vtx);
	float m3body    = vf["m3body"]->at(vtx);    
	float mSS       = vf["mSS"]->at(vtx);       
	float mOS1      = vf["mOS1"]->at(vtx);      
	float mOS2      = vf["mOS2"]->at(vtx);      
	float score     = vf["score"]->at(vtx);     
	float pval      = vf["pval"]->at(vtx);      
	float iso003    = vf["iso003"]->at(vtx);    
	float iso020    = vf["iso020"]->at(vtx);    
	float iso030    = vf["iso030"]->at(vtx);    
	float trkpval   = vf["trkspval"]->at(vtx);  
	float lxysig    = vf["lxysig"]->at(vtx);    
	float a0xysig   = vf["a0xysig"]->at(vtx);   
	float mettrk    = vf["mettrk"]->at(vtx);    
	float metcal    = vf["metcal"]->at(vtx);    
	float mttrk     = vf["mttrk"]->at(vtx);     
	float mtcal     = vf["mtcal"]->at(vtx);     
	float dphihttrk = vf["dphihttrk"]->at(vtx); 
	float dphihtcal = vf["dphihtcal"]->at(vtx); 
	int   njets     = vi["njets"]->at(vtx);    
	int   isLoose   = vi["pass_loose"]->at(vtx);    
	
	bool ncandcuts     = (vf["m3body"]->size()==1);
	bool passloose     = (isLoose==1);
	bool sidebandscuts = (m3body>xFullMin && m3body<xFullMax);

	bool bdtmincuts    = (score>xBDTMin);
	bool bdtoptcuts    = (score>xBDTOpt);
	bool BDTcuts       = (doBDT) ? (bdtmincuts && bdtoptcuts) : 1;
	
	bool iso           = (iso020<0.3 && iso030<1);
	bool trkspval      = (trkpval>1.e-9);
	bool distances     = (lxysig>-10 && lxysig<50 && a0xysig<25);
	bool kinematics    = (pt3body>10000 && (mettrk>10000 && mettrk<250000 && mttrk>20000) && (metcal>10000 && metcal<250000 && mtcal>20000));
	bool lowmass1      = (mOS1>220. && mOS2>220. && mSS>220.);
	bool preTraining   = (iso && trkspval && distances && kinematics && lowmass1);
	
	bool pvallxytrks  = (lxysig>1 && trkpval>8.e-9 && pval>0.2);
	bool MTHT         = (mttrk>45000 && mtcal>45000 && dphihttrk>2 && dphihtcal>2);
	// bool Njets        = (njets<2);
	bool lowmass2     = (mOS1>300 && mOS2>300 && mSS>300);
	bool onResonances = ((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50 || TMath::Abs(mOS1-782)<50 || TMath::Abs(mOS2-782)<50) && (mettrk<35000 || metcal<35000 || pt3body<35000));
	bool onDs         = ((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50) && TMath::Abs(m3body-1968)<100); // 100 because of 2*3body_resolution=~65 + margins because of the pi-mu mass differernce
	// bool postTraining = (pvallxytrks && MTHT && !onResonances && lowmass2 && Njets && !onDs);
	bool postTraining = (resVeto) ? (pvallxytrks && MTHT && !onResonances && lowmass2 && !onDs) : (pvallxytrks && MTHT && lowmass2);
	
	bool basecuts = (ncandcuts && passloose && sidebandscuts && BDTcuts);
	if     (mode=="preTraining")  basecuts = (basecuts && preTraining);
	else if(mode=="postTraining") basecuts = (basecuts && preTraining && postTraining);
	else _FAT("unknown mode="<<mode);

	bool srcuts        = (xSRmin>0 && xSRmax>0) ? (m3body>xSRmin && m3body<xSRmax)       : 1;
	bool blindedcuts   = (isBlinded)            ? (m3body<xBlindMin || m3body>xBlindMax) : 1;
	

	bool pass = basecuts;
	if     (type=="bkg")  pass = (pass && blindedcuts);
	else if(type=="sig")  pass = (pass && srcuts);
	else if(type=="full") pass = (pass);
	else _FAT("unknown type: "<<type);
	return pass;
}

#endif
