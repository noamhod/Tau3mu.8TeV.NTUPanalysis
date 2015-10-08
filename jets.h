#ifndef JETS_H
#define JETS_H

#include "std.h"

static float minJetPtGeV = 30.;

enum jetmetMode
{
	JETS_UNCALIB,
	JETS_JESUP, JETS_JESDWN,
	JETS_JERUP, JETS_JERDWN,
	JETS_CALIB
};

enum trkmetMode
{
	TRK_UNCALIB,
	TRK_SOFTUP, TRK_SOFTDWN,
	TRK_PARA,   TRK_PERP,
	TRK_JETUP,  TRK_JETDWN,
	TRK_CALIB
};

enum trkjetmetMode
{
	ALL_UNCALIB,
	ALL_JESUP,  ALL_JESDWN,
	ALL_JERUP,  ALL_JERDWN,
	ALL_SOFTUP, ALL_SOFTDWN,
	ALL_PARA,   ALL_PERP,
	ALL_JETUP,  ALL_JETDWN,
	ALL_CALIB
};

TString getJetMETword(unsigned int mode)
{
	TString jetmetword = "";
	switch(mode)
	{
		case JETS_CALIB:   jetmetword="";          break;
		case JETS_UNCALIB: jetmetword="_uncalib";  break;
		case JETS_JESUP:   jetmetword="_jes_up";   break;
		case JETS_JESDWN:  jetmetword="_jes_dwn";  break;
		case JETS_JERUP:   jetmetword="_jer_up";   break;
		case JETS_JERDWN:  jetmetword="_jer_dwn";  break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return jetmetword;
}
TString getJetMETname(unsigned int mode)
{
	TString jetmetname = "";
	switch(mode)
	{
		case JETS_CALIB:   jetmetname="Calibrated";    break;
		case JETS_UNCALIB: jetmetname="Uncalibrated";  break;
		case JETS_JESUP:   jetmetname="JESUP       ";  break;
		case JETS_JESDWN:  jetmetname="JESDWN      ";  break;
		case JETS_JERUP:   jetmetname="JERUP       ";  break;
		case JETS_JERDWN:  jetmetname="JERDWN      ";  break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return jetmetname;
}

TString getTrkMETword(unsigned int mode)
{
	TString jetmetword = "";
	switch(mode)
	{
		case TRK_CALIB:   jetmetword="";                 break;
		case TRK_UNCALIB: jetmetword="_uncalib";         break;
		case TRK_SOFTUP:  jetmetword="_softtrk_up";      break;
		case TRK_SOFTDWN: jetmetword="_softtrk_dwn";     break;
		case TRK_PARA:    jetmetword="_softtrkres_para"; break;
		case TRK_PERP:    jetmetword="_softtrkres_perp"; break;
		case TRK_JETUP:   jetmetword="_jettrk_up";       break;
		case TRK_JETDWN:  jetmetword="_jettrk_dwn";      break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return jetmetword;
}
TString getTrkMETname(unsigned int mode)
{
	TString jetmetname = "";
	switch(mode)
	{
		case TRK_CALIB:   jetmetname="Calibrated";    break;
		case TRK_UNCALIB: jetmetname="Uncalibrated";  break;
		case TRK_SOFTUP:  jetmetname="SOFTUP";        break;
		case TRK_SOFTDWN: jetmetname="SOFTDWN";       break;
		case TRK_PARA:    jetmetname="PARA";          break;
		case TRK_PERP:    jetmetname="PERP";          break;
		case TRK_JETUP:   jetmetname="JETUP";         break;
		case TRK_JETDWN:  jetmetname="JETDWN";        break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return jetmetname;
}

TString getTrkJetMETword(unsigned int mode)
{
	TString word = "";
	switch(mode)
	{
		case ALL_CALIB:   word="";                 break;
		case ALL_UNCALIB: word="_uncalib";         break;
		case ALL_JESUP:   word="_jes_up";          break;
		case ALL_JESDWN:  word="_jes_dwn";         break;
		case ALL_JERUP:   word="_jer_up";          break;
		case ALL_JERDWN:  word="_jer_dwn";         break;
		case ALL_SOFTUP:  word="_softtrk_up";      break;
		case ALL_SOFTDWN: word="_softtrk_dwn";     break;
		case ALL_PARA:    word="_softtrkres_para"; break;
		case ALL_PERP:    word="_softtrkres_perp"; break;
		case ALL_JETUP:   word="_jettrk_up";       break;
		case ALL_JETDWN:  word="_jettrk_dwn";      break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return word;
}
TString getTrkJetMETname(unsigned int mode)
{
	TString name = "";
	switch(mode)
	{
		case ALL_CALIB:   name="Calibrated";    break;
		case ALL_UNCALIB: name="Uncalibrated";  break;
		case ALL_JESUP:   name="JESUP";         break;
		case ALL_JESDWN:  name="JESDWN";        break;
		case ALL_JERUP:   name="JERUP";         break;
		case ALL_JERDWN:  name="JERDWN";        break;
		case ALL_SOFTUP:  name="SOFTUP";        break;
		case ALL_SOFTDWN: name="SOFTDWN";       break;
		case ALL_PARA:    name="PARA";          break;
		case ALL_PERP:    name="PERP";          break;
		case ALL_JETUP:   name="JETUP";         break;
		case ALL_JETDWN:  name="JETDWN";        break;
		default: _FAT("mode:"<<mode<<" is unknown");
	}
	return name;
}


unsigned int getNjets0to2(unsigned int vtx, TMapTSP2vf& vf)
{
	float GeV2MeV = 1000.;
	unsigned int njets = 0;
	if     (vf["jet_pt1"]->at(vtx)<=minJetPtGeV*GeV2MeV && vf["jet_pt2"]->at(vtx)<=minJetPtGeV*GeV2MeV) njets = 0;
	else if(vf["jet_pt1"]->at(vtx)>minJetPtGeV*GeV2MeV  && vf["jet_pt2"]->at(vtx)<=minJetPtGeV*GeV2MeV) njets = 1;
	else if(vf["jet_pt1"]->at(vtx)>minJetPtGeV*GeV2MeV  && vf["jet_pt2"]->at(vtx)>minJetPtGeV*GeV2MeV)  njets = 2;
	else _FAT("njets="<<njets<<" is not supported");
	return njets;
}
unsigned int getNjetsAll(unsigned int vtx, TMapTSP2vf& vf)
{
	float GeV2MeV = 1000.;
	int njets = (vf["jet_pt1"]->at(vtx)>minJetPtGeV*GeV2MeV) +
				(vf["jet_pt2"]->at(vtx)>minJetPtGeV*GeV2MeV) +
				(vf["jet_pt3"]->at(vtx)>minJetPtGeV*GeV2MeV) +
				(vf["jet_pt4"]->at(vtx)>minJetPtGeV*GeV2MeV);
	return njets;
}
unsigned int getNbjetsAll(unsigned int vtx, TMapTSP2vf& vf, float MV1min)
{
	float GeV2MeV = 1000.;
	int nbjets = (vf["jet_pt1"]->at(vtx)>minJetPtGeV*GeV2MeV && vf["jet_MV1w1"]->at(vtx)>MV1min) +
				 (vf["jet_pt2"]->at(vtx)>minJetPtGeV*GeV2MeV && vf["jet_MV1w2"]->at(vtx)>MV1min) +
				 (vf["jet_pt3"]->at(vtx)>minJetPtGeV*GeV2MeV && vf["jet_MV1w3"]->at(vtx)>MV1min) +
				 (vf["jet_pt4"]->at(vtx)>minJetPtGeV*GeV2MeV && vf["jet_MV1w4"]->at(vtx)>MV1min);
	return nbjets;
}

TLorentzVector getP3body(unsigned int vtx, TMapTSP2vf& vf)
{
	float muonMassMeV = 105.658367;
	TLorentzVector p1, p2, p3;
	p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
	p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
	p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
	TLorentzVector psum = p1+p2+p3;
	return psum;
}
TLorentzVector getPsum(unsigned int vtx, TMapTSP2vf& vf, int njets=4, TString jetMode="")
{
	float GeV2MeV = 1000.;
	TLorentzVector p3body = getP3body(vtx,vf);
	TLorentzVector pSum = p3body;
	njets = (njets>4) ? 4 : njets;
	for(int i=1 ; i<=njets ; ++i)
	{
		stringstream strm;                      string str;
		strm << setprecision(0) << fixed << i; strm >> str;
		TString si = str;
		if(vf["jet_pt"+si+jetMode]->at(vtx)<=minJetPtGeV*GeV2MeV) continue;
		
		TLorentzVector Jet;
		Jet.SetPtEtaPhiE(vf["jet_pt"+si+jetMode]->at(vtx),vf["jet_eta"+si+jetMode]->at(vtx),vf["jet_phi"+si+jetMode]->at(vtx),vf["jet_E"+si+jetMode]->at(vtx));
		pSum += Jet;
	}
	return pSum;
}

void setMETmode(unsigned int vtx, TMapTSf& f, TMapTSP2vf& vf, TString metMode)
{
	TString metmodeword;
	if     (metMode=="calib")    return; // this is the default content anyway...
	else if(metMode=="uncalib")  metmodeword="_uncalib";
	else if(metMode=="JESUP")    metmodeword="_jes_up";
	else if(metMode=="JESDWN")   metmodeword="_jes_dwn";
	else if(metMode=="JERUP")    metmodeword="_jer_up";
	else if(metMode=="JERDWN")   metmodeword="_jer_dwn";
	else _FAT("metMode:"<<metMode<<" is unknown");

	//// set the met variables
	f["met_reffinal_et"]                = f["met_reffinal_et"+metmodeword];
	f["met_reffinal_phi"]               = f["met_reffinal_phi"+metmodeword];
	vf["met_reffinal_mT"]->at(vtx)      = vf["met_reffinal_mT"+metmodeword]->at(vtx);
	vf["met_reffinal_dPhi3mu"]->at(vtx) = vf["met_reffinal_dPhi3mu"+metmodeword]->at(vtx);

	f["met_muons_et"]                = f["met_muons_et"+metmodeword];
	f["met_muons_phi"]               = f["met_muons_phi"+metmodeword];
	vf["met_muons_mT"]->at(vtx)      = vf["met_muons_mT"+metmodeword]->at(vtx);
	vf["met_muons_dPhi3mu"]->at(vtx) = vf["met_muons_dPhi3mu"+metmodeword]->at(vtx);
	
	// f["met_track_et"]                = f["met_track_et"+metmodeword];
	// f["met_track_phi"]               = f["met_track_phi"+metmodeword];
	// vf["met_track_mT"]->at(vtx)      = vf["met_track_mT"+metmodeword]->at(vtx);
	// vf["met_track_dPhi3mu"]->at(vtx) = vf["met_track_dPhi3mu"+metmodeword]->at(vtx);
}

void setTrkMETmode(unsigned int vtx, TMapTSf& f, TMapTSP2vf& vf, TString trkmetMode)
{
	TString metmodeword;
	if     (trkmetMode=="calib")    return; // this is the default content anyway...
	else if(trkmetMode=="uncalib")  metmodeword="_uncalib";
	else if(trkmetMode=="SOFTUP")  metmodeword="_softtrk_up";     
	else if(trkmetMode=="SOFTDWN") metmodeword="_softtrk_dwn";    
	else if(trkmetMode=="PARA")    metmodeword="_softtrkres_para";
	else if(trkmetMode=="PERP")    metmodeword="_softtrkres_perp";
	else if(trkmetMode=="JETUP")   metmodeword="_jettrk_up";      
	else if(trkmetMode=="JETDWN")  metmodeword="_jettrk_dwn";     
	else _FAT("metMode:"<<trkmetMode<<" is unknown");

	//// set the met variables
	f["met_track_et"]                = f["met_track_et"+metmodeword];
	f["met_track_phi"]               = f["met_track_phi"+metmodeword];
	vf["met_track_mT"]->at(vtx)      = vf["met_track_mT"+metmodeword]->at(vtx);
	vf["met_track_dPhi3mu"]->at(vtx) = vf["met_track_dPhi3mu"+metmodeword]->at(vtx);
}

void setJETmode(unsigned int vtx, TMapTSP2vi& vi, TMapTSP2vf& vf, TString jetMode)
{
	TString jetmodeword;
	if     (jetMode=="calib")    return; // this is the default content anyway...
	else if(jetMode=="uncalib")  jetmodeword="_uncalib";
	else if(jetMode=="JESUP")    jetmodeword="_jes_up";
	else if(jetMode=="JESDWN")   jetmodeword="_jes_dwn";
	else if(jetMode=="JERUP")    jetmodeword="_jer_up";
	else if(jetMode=="JERDWN")   jetmodeword="_jer_dwn";
	else _FAT("jetMode:"<<jetMode<<" is unknown");

	//// set the jet variables
	for(int i=1 ; i<=4 ; i++)
	{
		TString si;
		stringstream strm;
		strm << setprecision(0) << fixed << i;
		strm >> si;
		
		vf["jet_pt"+si]->at(vtx)   = vf["jet_pt"+si+jetmodeword]->at(vtx);
		vf["jet_eta"+si]->at(vtx)  = vf["jet_eta"+si+jetmodeword]->at(vtx);
		vf["jet_phi"+si]->at(vtx)  = vf["jet_phi"+si+jetmodeword]->at(vtx);
		vf["jet_E"+si]->at(vtx)    = vf["jet_E"+si+jetmodeword]->at(vtx);
		vf["jet_m"+si]->at(vtx)    = vf["jet_m"+si+jetmodeword]->at(vtx);
		vf["jet_MV1w"+si]->at(vtx) = vf["jet_MV1w"+si+jetmodeword]->at(vtx);
		vf["jet_vtxf"+si]->at(vtx) = vf["jet_vtxf"+si+jetmodeword]->at(vtx);
		vi["jet_ntrk"+si]->at(vtx) = vi["jet_ntrk"+si+jetmodeword]->at(vtx);
	}

	//// set the jet composite variables
	vf["jet_sumpt12"]->at(vtx)   = vf["jet_sumpt12"+jetmodeword]->at(vtx);
	vf["jet_dphi3muJ1"]->at(vtx) = vf["jet_dphi3muJ1"+jetmodeword]->at(vtx);
	vf["jet_dR3muJ1"]->at(vtx)   = vf["jet_dR3muJ1"+jetmodeword]->at(vtx);
	vf["jet_dphiJ1J2"]->at(vtx)  = vf["jet_dphiJ1J2"+jetmodeword]->at(vtx);
	vf["jet_dRJ1J2"]->at(vtx)    = vf["jet_dRJ1J2"+jetmodeword]->at(vtx);
}

void fixJets(unsigned int vtx, TMapTSP2vi& vi, TMapTSP2vf& vf)
{	
	float GeV2MeV = 1000.;
	float defaultvalue  = -1.;
	float defaultpt     = 0.;
	
	///////////////////////////////////////////////////////////////////////
	//// Fix the calibrated jet info - need to fill jet histos before this
	if(vf["jet_pt1"]->at(vtx)<minJetPtGeV*GeV2MeV)
	{
		vf["jet_pt1"]->at(vtx)       = defaultpt;
		vf["jet_dphi3muJ1"]->at(vtx) = defaultvalue;
		vf["jet_dR3muJ1"]->at(vtx)   = defaultvalue;
		vf["jet_dphiJ1J2"]->at(vtx)  = defaultvalue;
		vf["jet_dRJ1J2"]->at(vtx)    = defaultvalue;
		vf["jet_vtxf1"]->at(vtx)     = defaultvalue;
		vi["jet_ntrk1"]->at(vtx)     = defaultvalue;
	}
	if(vf["jet_pt2"]->at(vtx)<minJetPtGeV*GeV2MeV)
	{
		vf["jet_pt2"]->at(vtx)      = defaultpt;
		vf["jet_dphiJ1J2"]->at(vtx) = defaultvalue;
		vf["jet_dRJ1J2"]->at(vtx)   = defaultvalue;
		vf["jet_vtxf2"]->at(vtx)    = defaultvalue;
		vi["jet_ntrk2"]->at(vtx)    = defaultvalue;
	}
	if(vf["jet_pt3"]->at(vtx)<minJetPtGeV*GeV2MeV) vf["jet_pt3"]->at(vtx) = defaultpt;
	if(vf["jet_pt4"]->at(vtx)<minJetPtGeV*GeV2MeV) vf["jet_pt4"]->at(vtx) = defaultpt;
	vf["jet_sumpt12"]->at(vtx) = vf["jet_pt1"]->at(vtx)+vf["jet_pt2"]->at(vtx); // that will automatically set it to the right value
	
	/////////////////////////////////////////////////////////////////////////
	//// Fix the uncalibrated jet info - need to fill jet histos before this
	if(vf["jet_pt1_uncalib"]->at(vtx)<minJetPtGeV*GeV2MeV)
	{
		vf["jet_pt1_uncalib"]->at(vtx)       = defaultpt;
		vf["jet_dphi3muJ1_uncalib"]->at(vtx) = defaultvalue;
		vf["jet_dR3muJ1_uncalib"]->at(vtx)   = defaultvalue;
		vf["jet_dphiJ1J2_uncalib"]->at(vtx)  = defaultvalue;
		vf["jet_dRJ1J2_uncalib"]->at(vtx)    = defaultvalue;
		vf["jet_vtxf1_uncalib"]->at(vtx)     = defaultvalue;
		vi["jet_ntrk1_uncalib"]->at(vtx)     = defaultvalue;
	}
	if(vf["jet_pt2_uncalib"]->at(vtx)<minJetPtGeV*GeV2MeV)
	{
		vf["jet_pt2_uncalib"]->at(vtx)      = defaultpt;
		vf["jet_dphiJ1J2_uncalib"]->at(vtx) = defaultvalue;
		vf["jet_dRJ1J2_uncalib"]->at(vtx)   = defaultvalue;
		vf["jet_vtxf2_uncalib"]->at(vtx)    = defaultvalue;
		vi["jet_ntrk2_uncalib"]->at(vtx)    = defaultvalue;
	}
	if(vf["jet_pt3_uncalib"]->at(vtx)<minJetPtGeV*GeV2MeV) vf["jet_pt3_uncalib"]->at(vtx) = defaultpt;
	if(vf["jet_pt4_uncalib"]->at(vtx)<minJetPtGeV*GeV2MeV) vf["jet_pt4_uncalib"]->at(vtx) = defaultpt;
	vf["jet_sumpt12_uncalib"]->at(vtx) = vf["jet_pt1_uncalib"]->at(vtx)+vf["jet_pt2_uncalib"]->at(vtx); // that will automatically set it to the right value
}


#endif
