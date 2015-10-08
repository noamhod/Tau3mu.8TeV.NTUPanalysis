#ifndef SELECT_H
#define SELECT_H

#include "std.h"
#include "jets.h"
// #include "enums.h"

bool MCP(unsigned int vtx, int itrk, TMapTSP2vi& vi, TMapTSP2vf& vf)
{
	// !expectBLayerHit OR numberOfBLayerHits > 0,  ===> REMOVED !!!
	// number of pixel hits + number of crossed dead pixel sensors > 0,
	// number of SCT hits + number of crossed dead SCT sensors > 4,
	// Number of pixel holes + number of SCT holes < 3,
	// TRT hits > 5 for 0.1<|eta|<1.9, outlier fraction < 0.9.
	// If TRT hits > 5 for |eta|<=0.1 or |eta|>=1.9, outlier fraction < 0.9  ===> REMOVED !!!

	int nPIXhits = -1; int nDeadPIX     = -1; int nPIXholes  = -1;
	int nSCThits = -1; int nDeadSCT     = -1; int nSCTholes  = -1;
	int nTRThits = -1; int nTRToutliers = -1; int nHtTRThits = -1;
	float eta = -999;
	int type = -1;
	if(itrk==1)
	{
		nPIXhits = vi["mu_nPIXhits1"]->at(vtx); nDeadPIX     = vi["mu_nDeadPIX1"]->at(vtx);     nPIXholes  = vi["mu_nPIXholes1"]->at(vtx);
		nSCThits = vi["mu_nSCThits1"]->at(vtx); nDeadSCT     = vi["mu_nDeadSCT1"]->at(vtx);     nSCTholes  = vi["mu_nSCTholes1"]->at(vtx);
		nTRThits = vi["mu_nTRThits1"]->at(vtx); nTRToutliers = vi["mu_nTRToutliers1"]->at(vtx); nHtTRThits = vi["mu_htTRThits1"]->at(vtx);
		eta = vf["mu_eta1"]->at(vtx);
		type = vi["mu_type1"]->at(vtx);
	}
	else if(itrk==2)
	{
		nPIXhits = vi["mu_nPIXhits2"]->at(vtx); nDeadPIX     = vi["mu_nDeadPIX2"]->at(vtx);     nPIXholes  = vi["mu_nPIXholes2"]->at(vtx);
		nSCThits = vi["mu_nSCThits2"]->at(vtx); nDeadSCT     = vi["mu_nDeadSCT2"]->at(vtx);     nSCTholes  = vi["mu_nSCTholes2"]->at(vtx);
		nTRThits = vi["mu_nTRThits2"]->at(vtx); nTRToutliers = vi["mu_nTRToutliers2"]->at(vtx); nHtTRThits = vi["mu_htTRThits2"]->at(vtx);
		eta = vf["mu_eta2"]->at(vtx);
		type = vi["mu_type2"]->at(vtx);
	}
	else if(itrk==3)
	{
		nPIXhits = vi["mu_nPIXhits3"]->at(vtx); nDeadPIX     = vi["mu_nDeadPIX3"]->at(vtx);     nPIXholes  = vi["mu_nPIXholes3"]->at(vtx);
		nSCThits = vi["mu_nSCThits3"]->at(vtx); nDeadSCT     = vi["mu_nDeadSCT3"]->at(vtx);     nSCTholes  = vi["mu_nSCTholes3"]->at(vtx);
		nTRThits = vi["mu_nTRThits3"]->at(vtx); nTRToutliers = vi["mu_nTRToutliers3"]->at(vtx); nHtTRThits = vi["mu_htTRThits3"]->at(vtx);
		eta = vf["mu_eta3"]->at(vtx);
		type = vi["mu_type3"]->at(vtx);
	}
	else _FAT("no index: "<<itrk);
	if(type<0) _FAT("track type could not be resolved: "<<type);
	
	/////////////////////
	//// the cut is below

	bool passedPIX = ((nPIXhits+nDeadPIX) > 0);
	if(!passedPIX) return false;

	bool passedSCT = ((nSCThits+nDeadSCT) > 4.);
	if(!passedSCT) return false;

	bool passedHOL = ((nPIXholes+nSCTholes) < 3.);
	if(!passedHOL) return false;

	float TRTratio = 0.9;
	float nTRTmin  = 5.;
	float nTRT     = nTRThits+nTRToutliers;
	float nTRTol   = nTRToutliers;
	if(fabs(eta)>0.1 && fabs(eta)<1.9)
	{
		// if(nTRT<=nTRTmin)         { _INF(1,"fail n>5 ("<<nTRT<<")"); return false; }
		// if(nTRTol>=TRTratio*nTRT) { _INF(1,"fail OL fraction ("<<nTRTol/nTRT<<")"); return false; }
		if(nTRT>nTRTmin)
		{
			if(nTRTol>=TRTratio*nTRT) return false;
		}
	}

	return true;
}

bool highThresholdTRTRhits(unsigned int vtx, int itrk, TMapTSP2vi& vi)
{
	int nHtTRThits = -1;
	int nTRThits   = -1;
	int type = -1;
	if(itrk==1)
	{
		nTRThits = vi["mu_nTRThits1"]->at(vtx);
		nHtTRThits = vi["mu_htTRThits1"]->at(vtx);
		type = vi["mu_type1"]->at(vtx);
	}
	else if(itrk==2)
	{
		nTRThits = vi["mu_nTRThits2"]->at(vtx);
		nHtTRThits = vi["mu_htTRThits2"]->at(vtx);
		type = vi["mu_type2"]->at(vtx);
	}
	else if(itrk==3)
	{
		nTRThits = vi["mu_nTRThits3"]->at(vtx);
		nHtTRThits = vi["mu_htTRThits3"]->at(vtx);
		type = vi["mu_type3"]->at(vtx);
	}
	else _FAT("no index: "<<itrk);
	if(type<0) _FAT("track type could not be resolved: "<<type);
	
	/////////////////////
	//// the cut is below
	
	// if(type==2               && nHtTRThits>=5)  return false;
	// if((type==0 || type==1)  && nHtTRThits>=10) return false;
	
	float htTRTfraction = (nTRThits>0) ? (float)nHtTRThits/(float)nTRThits : 0.;
	if(htTRTfraction>0.35) return false;
	
	return true;
}

bool TPaEtaPhiMSlayers(unsigned int vtx, int itrk, TMapTSP2vi& vi)
{
	int nEtaPhiLayers = -1;
	int nPrecisionHits = -1;
	int type = -1;
	if(itrk==1)
	{
		nEtaPhiLayers  = vi["mu_nEtaPhiLayers1"]->at(vtx);
		nPrecisionHits = vi["mu_nPrecisionHits1"]->at(vtx);
		type = vi["mu_type1"]->at(vtx);
	}
	else if(itrk==2)
	{
		nEtaPhiLayers  = vi["mu_nEtaPhiLayers2"]->at(vtx);
		nPrecisionHits = vi["mu_nPrecisionHits2"]->at(vtx);
		type = vi["mu_type2"]->at(vtx);
	}
	else if(itrk==3)
	{
		nEtaPhiLayers  = vi["mu_nEtaPhiLayers3"]->at(vtx);
		nPrecisionHits = vi["mu_nPrecisionHits3"]->at(vtx);
		type = vi["mu_type3"]->at(vtx);
	}
	else _FAT("no index: "<<itrk);
	if(type<0) _FAT("track type could not be resolved: "<<type);
	
	/////////////////////
	//// the cut is below
	
	if(type==0 && nEtaPhiLayers<2 && nPrecisionHits<12)  return false;
	if(type==1 && nEtaPhiLayers<3 && nPrecisionHits<13)  return false;
	
	return true;
}

bool TGCRPChits(unsigned int vtx, int itrk, TMapTSP2vi& vi)
{
	int nRPCphi = -1; int nRPCeta = -1;
	int nTGCphi = -1; int nTGCeta = -1;
	int type = -1;
	if(itrk==1)
	{
		nRPCphi = vi["mu_nRPCPhiHits1"]->at(vtx); nTGCphi = vi["mu_nTGCPhiHits1"]->at(vtx);
		nRPCeta = vi["mu_nRPCEtaHits1"]->at(vtx); nTGCeta = vi["mu_nTGCEtaHits1"]->at(vtx);
		type = vi["mu_type1"]->at(vtx);
	}
	else if(itrk==2)
	{
		nRPCphi = vi["mu_nRPCPhiHits2"]->at(vtx); nTGCphi = vi["mu_nTGCPhiHits2"]->at(vtx);
		nRPCeta = vi["mu_nRPCEtaHits2"]->at(vtx); nTGCeta = vi["mu_nTGCEtaHits2"]->at(vtx);
		type = vi["mu_type2"]->at(vtx);
	}
	else if(itrk==3)
	{
		nRPCphi = vi["mu_nRPCPhiHits3"]->at(vtx); nTGCphi = vi["mu_nTGCPhiHits3"]->at(vtx);
		nRPCeta = vi["mu_nRPCEtaHits3"]->at(vtx); nTGCeta = vi["mu_nTGCEtaHits3"]->at(vtx);
		type = vi["mu_type3"]->at(vtx);
	}
	else _FAT("no index: "<<itrk);
	if(type<0) _FAT("track type could not be resolved: "<<type);
	
	/////////////////////
	//// the cut is below
	
	// type==1==TPA
	
	if(type==1 && (nRPCphi>0 && nTGCphi>0))  return false;
	if(type==1 && (nRPCeta>0 && nTGCeta>0))  return false;
	
	return true;
}

bool isMedium(unsigned int vtx, int itrk, TMapTSP2vi& vi)
{
	bool medium;
	int type = -1;
	if(itrk==1)
	{
		type = vi["mu_type1"]->at(vtx);
		medium = (type==0) ? (vints["mu_isMedium1"]->at(vtx)==1) : true;
	}
	else if(itrk==2)
	{
		type = vi["mu_type2"]->at(vtx);
		medium = (type==0) ? (vints["mu_isMedium2"]->at(vtx)==1) : true;
	}
	else if(itrk==3)
	{
		type = vi["mu_type3"]->at(vtx);
		medium = (type==0) ? (vints["mu_isMedium3"]->at(vtx)==1) : true;
	}
	else _FAT("no index: "<<itrk);
	if(type<0) _FAT("track type could not be resolved: "<<type);
	
	// type==0==MUON
	
	/////////////////////
	//// the cut is below
	
	return medium;
}


bool isRhoOmegaPhi(unsigned int vtx, TMapTSP2vi& vi, TMapTSP2vf& vf)
{
	int type1 = vi["mu_type1"]->at(vtx);
	int type2 = vi["mu_type1"]->at(vtx);
	int type3 = vi["mu_type1"]->at(vtx);
	// type=0=MUON
	// else, type=TPa
	
	float m3body = vf["vtx_mass"]->at(vtx);;
	int   code   = vi["vtx_code"]->at(vtx);;
	
	float muonMassMeV = 105.658367; // MeV
	TLorentzVector p1, p2, p3, p, p12, p13, p23;
	p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
	p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
	p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
	p12 = p1+p2;
	p13 = p1+p3;
	p23 = p2+p3;
	float q1 = (vf["mu_trkqoverp1"]->at(vtx)<0) ? -1. : +1;
	float q2 = (vf["mu_trkqoverp2"]->at(vtx)<0) ? -1. : +1;
	float q3 = (vf["mu_trkqoverp3"]->at(vtx)<0) ? -1. : +1;
	float q12 = fabs(q1+q2);
	float q13 = fabs(q1+q3);
	float q23 = fabs(q2+q3);
	
	float mRho   = 770;
	float mOmega = 782;
	float mPhi   = 1020;
	float mDs    = 1968;
	float range  = 30;
	
	if(q12==0.)
	{
		float mOS = p12.M();
		if(fabs(mOS-mRho)<range   && type3!=0) return true;
		if(fabs(mOS-mOmega)<range && type3!=0) return true;
		if(fabs(mOS-mPhi)<range   && type3!=0) return true;
	}
	if(q13==0.)
	{
		float mOS = p13.M();
		if(fabs(mOS-mRho)<range   && type2!=0) return true;
		if(fabs(mOS-mOmega)<range && type2!=0) return true;
		if(fabs(mOS-mPhi)<range   && type2!=0) return true;
	}
	if(q23==0.)
	{
		float mOS = p23.M();
		if(fabs(mOS-mRho)<range   && type1!=0) return true;
		if(fabs(mOS-mOmega)<range && type1!=0) return true;
		if(fabs(mOS-mPhi)<range   && type1!=0) return true;
	}
	
	// if(fabs(m3body-mDs)<range && code!=0) return true;
	
	
	return false;
}

bool passBjetVeto(unsigned int vtx, TMapTSP2vf& vf, bool isAntiKt4LCJet, TString tightness)
{
	float GeV2MeV = 1000.;
	float ptj1 = -1; float mj1  = -1; float fwj1  = -99;
	float ptj2 = -1; float mj2  = -1; float fwj2  = -99;
	float ptj3 = -1; float mj3  = -1; float fwj3  = -99;
	float ptj4 = -1; float mj4  = -1; float fwj4  = -99;

	float btaggingCut = (isAntiKt4LCJet) ? 0.3511 : 0.39;
	ptj1 = vf["jet_pt1"]->at(vtx); mj1 = vf["jet_m1"]->at(vtx); fwj1 = vf["jet_MV1w1"]->at(vtx);
	ptj2 = vf["jet_pt2"]->at(vtx); mj2 = vf["jet_m2"]->at(vtx); fwj2 = vf["jet_MV1w2"]->at(vtx);
	ptj3 = vf["jet_pt3"]->at(vtx); mj3 = vf["jet_m3"]->at(vtx); fwj3 = vf["jet_MV1w3"]->at(vtx);
	ptj4 = vf["jet_pt4"]->at(vtx); mj4 = vf["jet_m4"]->at(vtx); fwj4 = vf["jet_MV1w4"]->at(vtx);	

	bool bjetveto1 = (tightness=="loose") ? (ptj1>(minJetPtGeV+10.)*GeV2MeV && fwj1>=btaggingCut)  :  (ptj1>minJetPtGeV*GeV2MeV && fwj1>=btaggingCut);
	bool bjetveto2 = (tightness=="loose") ? (ptj2>(minJetPtGeV+10.)*GeV2MeV && fwj2>=btaggingCut)  :  (ptj2>minJetPtGeV*GeV2MeV && fwj2>=btaggingCut);
	bool bjetveto3 = (tightness=="loose") ? (ptj3>(minJetPtGeV+10.)*GeV2MeV && fwj3>=btaggingCut)  :  (ptj3>minJetPtGeV*GeV2MeV && fwj3>=btaggingCut);
	bool bjetveto4 = (tightness=="loose") ? (ptj4>(minJetPtGeV+10.)*GeV2MeV && fwj4>=btaggingCut)  :  (ptj4>minJetPtGeV*GeV2MeV && fwj4>=btaggingCut);

	if(bjetveto1+bjetveto2+bjetveto3+bjetveto4>0) return false;
	return true;
}

bool passCollinearJetVeto(unsigned int vtx, TMapTSP2vf& vf, TString tightness)
{
	float GeV2MeV = 1000.;
	float ptj1 = -1;
	float dphi = -1;
	ptj1 = vf["jet_pt1"]->at(vtx);
	dphi = vf["jet_dphi3muJ1"]->at(vtx);
	// float dphilow  = (tightness=="loose") ? 0.15        : 0.2;
	float dphihigh = (tightness=="loose") ? 2.95        : 2.7;
	float ptcut    = (tightness=="loose") ? 50.*GeV2MeV : 45.*GeV2MeV;
	// if(ptj1>ptcut && (dphi<dphilow || dphi>dphihigh)) return false;
	if(ptj1>ptcut && dphi>dphihigh) return false;
	return true;
}

bool passDiJetVeto(unsigned int vtx, TMapTSP2vf& vf, TString tightness)
{
	float GeV2MeV = 1000.;
	float ptj1 = -1;
	float ptj2 = -1;
	float dphi = -1;
	float sumpt = -1;
	ptj1 = vf["jet_pt1"]->at(vtx);
	ptj2 = vf["jet_pt2"]->at(vtx);
	sumpt = vf["jet_sumpt12"]->at(vtx);
	dphi  = vf["jet_dphiJ1J2"]->at(vtx);
	float sumptcut = (tightness=="loose") ? 80.*GeV2MeV : 70.*GeV2MeV;// 80  : 60;
	float dphicut  = (tightness=="loose") ? 2.7         : 2.5;
	if(ptj1>minJetPtGeV*GeV2MeV && ptj2>minJetPtGeV*GeV2MeV && sumpt>sumptcut && dphi>dphicut) return false;
	return true;
}

bool passMET(TMapTSf& f, TString metType, TString tightness)
{
	float GeV2MeV = 1000.;
	float met_calo  = f["met_"+metType+"_et"];
	float met_track = f["met_track_et"];
	float metcut_calo  = (tightness=="loose") ? 10.*GeV2MeV : 20.*GeV2MeV;
	float metcut_track = (tightness=="loose") ? 10.*GeV2MeV : 20.*GeV2MeV;
	if(met_calo<metcut_calo || met_track<metcut_track || met_calo>250.*GeV2MeV || met_track>250.*GeV2MeV) return false;
	return true;
}
bool passMT(unsigned int vtx, TMapTSP2vf& vf, TString metType, TString tightness)
{
	float GeV2MeV = 1000.;
	float mt_calo  = vf["met_"+metType+"_mT"]->at(vtx);
	float mt_track = vf["met_track_mT"]->at(vtx);
	float mtcut_calo  = (tightness=="loose") ? 20.*GeV2MeV : 60.*GeV2MeV;
	float mtcut_track = (tightness=="loose") ? 20.*GeV2MeV : 60.*GeV2MeV;
	if(mt_calo<mtcut_calo ||  mt_track<mtcut_track) return false;
	return true;
}
bool passDphi3bodyMET(unsigned int vtx, TMapTSP2vf& vf, TString metType, TString tightness)
{
	float dphi_calo = vf["met_"+metType+"_dPhi3mu"]->at(vtx);
	float dphi_track = vf["met_track_dPhi3mu"]->at(vtx);
	float dphicut_calo  = (tightness=="loose") ? 1.5 : 2.0;
	float dphicut_track = (tightness=="loose") ? 1.5 : 2.0;
	if(dphi_calo<dphicut_calo || dphi_track<dphicut_track) return false;
	// if(dphi_calo<dphicut_calo) return false;
	return true;
}
bool passDphiMETs(unsigned int vtx, TMapTSP2vf& vfdec, TString tightness)
{
	float dphi  = vfdec["mets_dphi"]->at(vtx);
	float dphicut = (tightness=="loose") ? 2.0 : 1.5;
	if(dphi>dphicut) return false;
	return true;
}

bool passHT(unsigned int vtx, TMapTSP2vf& vfdec, TString tightness)
{
	float GeV2MeV = 1000.;
	float ht = vfdec["ht_pt"]->at(vtx);
	float htcut = (tightness=="loose") ? 20.*GeV2MeV : 25.*GeV2MeV;
	if(ht<htcut) return false;
	return true;
}
bool passMHT(unsigned int vtx, TMapTSP2vf& vfdec, TString tightness)
{
	float GeV2MeV = 1000.;
	float mht_cal = vfdec["ht_mT"]->at(vtx);
	float mht_trk = vfdec["ht_mT_mettrk"]->at(vtx);
	float mhtcut_cal = (tightness=="loose") ? 20.*GeV2MeV : 60.*GeV2MeV;
	float mhtcut_trk = (tightness=="loose") ? 30.*GeV2MeV : 60.*GeV2MeV;
	if(mht_cal<mhtcut_cal || mht_trk<mhtcut_trk) return false;
	return true;
}
bool passDphihtMET(unsigned int vtx, TMapTSP2vf& vfdec, TString metType, TString tightness)
{
	float dphi    = vfdec["ht_dphimet_"+metType]->at(vtx);
	float dphitrk = vfdec["ht_dphimet_track"]->at(vtx);
	float dphicut    = (tightness=="loose") ? 1.5 : 2.0;
	float dphitrkcut = (tightness=="loose") ? 1.5 : 2.0;
	if(dphi<dphicut || dphitrk<dphitrkcut) return false;
	// if(dphi<dphicut) return false;
	return true;
}
bool passDrht3body(unsigned int vtx, TMapTSP2vf& vfdec, TString tightness)
{
	float dr = vfdec["ht_dr3body"]->at(vtx);
	float drcut = (tightness=="loose") ? 3.0 : 2.5;
	if(dr>drcut) return false;
	return true;
}



#endif
