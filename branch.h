#ifndef BRANCH_H
#define BRANCH_H


#include "std.h"
#include "type.h"

enum objecttype
{
	INT,FLT,DBL,STR,
	VINT,VFLT,VDBL,VSTR,
	VVINT,VVFLT,VVDBL,VVSTR
};

TMapTSi branchtypes;
TMapTSi ints;
TMapTSf floats;
TMapTSs strings;
TMapTSP2vi vints;
TMapTSP2vf vfloats;
TMapTSP2vs vstrings;
TMapTSP2vi vints_decorations;
TMapTSP2vf vfloats_decorations;

void setBranch(TString name, int type)
{
	if     (type==INT)  ints.insert(make_pair(name, -9999));
	else if(type==FLT)  floats.insert(make_pair(name, -9999.));
	else if(type==STR)  strings.insert(make_pair(name, ""));
	else if(type==VINT) vints.insert(make_pair(name, (vector<int>*)NULL));
	else if(type==VFLT) vfloats.insert(make_pair(name, (vector<float>*)NULL));
	else if(type==VSTR) vstrings.insert(make_pair(name, (vector<string>*)NULL));
	else _FAT("unsupported type of branch");
	branchtypes.insert(make_pair(name,type));
}
void setBranches(TTree* itree)
{
	// event info
	setBranch("evt_RunNumber",         INT);
	setBranch("evt_lbn",               INT);
	setBranch("evt_EventNumber",       INT);
	setBranch("evt_actualIntPerXing",  INT);
	setBranch("evt_averageIntPerXing", INT);
	
	// weights info
	setBranch("wgt_shapeFONLL", FLT);
	setBranch("wgt_normFONLL",  FLT);
	setBranch("wgt_luminosity", FLT);
	setBranch("wgt_kfactor",    FLT);
	setBranch("wgt_dijets",     FLT);
	setBranch("wgt_pileup",     FLT);
	setBranch("wgt_total",      FLT);
	
	// trigger info
	setBranch("EF_3mu4T",                INT);
	setBranch("EF_3mu6",                 INT);
	setBranch("EF_3mu6_MSonly",          INT);
	setBranch("EF_2mu13",                INT);
	setBranch("EF_mu18_tight_mu8_EFFS",  INT);
	setBranch("EF_mu18_tight_2mu4_EFFS", INT);
	setBranch("EF_2mu8_EFxe30_tclcw",    INT);
	setBranch("EF_mu24_tight_EFxe40",    INT);
	setBranch("EF_mu24i_tight",          INT);
	setBranch("EF_mu36_tight",           INT);
	
	// vertex info
	setBranch("vtx_n",         INT);
	setBranch("vtx_type",      VSTR);
	setBranch("vtx_code",      VINT);
	setBranch("vtx_charge",    VFLT);
	setBranch("vtx_mass",      VFLT);
	setBranch("vtx_mOS1",      VFLT);
	setBranch("vtx_mOS2",      VFLT);
	setBranch("vtx_mSS",       VFLT);
	setBranch("vtx_pt",        VFLT);
	setBranch("vtx_rapidity",  VFLT);
	setBranch("vtx_chi2",      VFLT);
	setBranch("vtx_ndf",       VFLT);
	setBranch("vtx_chi2ndf",   VFLT);
	setBranch("vtx_pval",      VFLT);
	setBranch("vtx_lxy",       VFLT);
	setBranch("vtx_lxyErr",    VFLT);
	setBranch("vtx_a0",        VFLT);
	setBranch("vtx_a0Err",     VFLT);
	setBranch("vtx_a0xy",      VFLT);
	setBranch("vtx_a0xyErr",   VFLT);
	setBranch("vtx_cosT",      VFLT);
	setBranch("vtx_cosTxy",    VFLT);
	setBranch("vtx_tau",       VFLT);
	setBranch("vtx_tauErr",    VFLT);
	setBranch("vtx_ptfrac12",  VFLT);
	setBranch("vtx_ptfrac23",  VFLT);
	setBranch("vtx_ptfrac13",  VFLT);
	setBranch("vtx_dpt12",     VFLT);
	setBranch("vtx_dpt23",     VFLT);
	setBranch("vtx_dpt13",     VFLT);
	setBranch("vtx_dRmax",     VFLT);
	setBranch("vtx_dRmin",     VFLT);
	setBranch("vtx_pvNtrk",    VINT);
	setBranch("vtx_nElectrons",VINT);
	setBranch("vtx_isolation000", VFLT);
	setBranch("vtx_isolation001", VFLT);
	setBranch("vtx_isolation002", VFLT);
	setBranch("vtx_isolation003", VFLT);
	setBranch("vtx_isolation004", VFLT);
	setBranch("vtx_isolation005", VFLT);
	setBranch("vtx_isolation006", VFLT);
	setBranch("vtx_isolation007", VFLT);
	setBranch("vtx_isolation008", VFLT);
	setBranch("vtx_isolation009", VFLT);
	setBranch("vtx_isolation010", VFLT);
	setBranch("vtx_isolation012", VFLT);
	setBranch("vtx_isolation014", VFLT);
	setBranch("vtx_isolation016", VFLT);
	setBranch("vtx_isolation018", VFLT);
	setBranch("vtx_isolation020", VFLT);
	setBranch("vtx_isolation022", VFLT);
	setBranch("vtx_isolation024", VFLT);
	setBranch("vtx_isolation026", VFLT);
	setBranch("vtx_isolation028", VFLT);
	setBranch("vtx_isolation030", VFLT);
	
	// MET_RefFinal-based kinematics
	setBranch("met_reffinal_et_uncalib",      FLT);
	setBranch("met_reffinal_phi_uncalib",     FLT);
	setBranch("met_reffinal_mT_uncalib",      VFLT);
	setBranch("met_reffinal_dPhi3mu_uncalib", VFLT);
	
	setBranch("met_reffinal_et",      FLT);
	setBranch("met_reffinal_phi",     FLT);
	setBranch("met_reffinal_mT",      VFLT);
	setBranch("met_reffinal_dPhi3mu", VFLT);
	
	setBranch("met_reffinal_et_jes_up",       FLT);
	setBranch("met_reffinal_phi_jes_up",      FLT);
	setBranch("met_reffinal_mT_jes_up",       VFLT);
	setBranch("met_reffinal_dPhi3mu_jes_up",  VFLT);
	
	setBranch("met_reffinal_et_jes_dwn",      FLT);
	setBranch("met_reffinal_phi_jes_dwn",     FLT);
	setBranch("met_reffinal_mT_jes_dwn",      VFLT);
	setBranch("met_reffinal_dPhi3mu_jes_dwn", VFLT);
	
	setBranch("met_reffinal_et_jer_up",       FLT);
	setBranch("met_reffinal_phi_jer_up",      FLT);
	setBranch("met_reffinal_mT_jer_up",       VFLT);
	setBranch("met_reffinal_dPhi3mu_jer_up",  VFLT);
	
	setBranch("met_reffinal_et_jer_dwn",      FLT);
	setBranch("met_reffinal_phi_jer_dwn",     FLT);
	setBranch("met_reffinal_mT_jer_dwn",      VFLT);
	setBranch("met_reffinal_dPhi3mu_jer_dwn", VFLT);
	
	// MET_Muons-based kinematics
	setBranch("met_muons_et_uncalib",      FLT);
	setBranch("met_muons_phi_uncalib",     FLT);
	setBranch("met_muons_mT_uncalib",      VFLT);
	setBranch("met_muons_dPhi3mu_uncalib", VFLT);
	
	setBranch("met_muons_et",      FLT);
	setBranch("met_muons_phi",     FLT);
	setBranch("met_muons_mT",      VFLT);
	setBranch("met_muons_dPhi3mu", VFLT);
	
	setBranch("met_muons_et_jes_up",       FLT);
	setBranch("met_muons_phi_jes_up",      FLT);
	setBranch("met_muons_mT_jes_up",       VFLT);
	setBranch("met_muons_dPhi3mu_jes_up",  VFLT);
	
	setBranch("met_muons_et_jes_dwn",      FLT);
	setBranch("met_muons_phi_jes_dwn",     FLT);
	setBranch("met_muons_mT_jes_dwn",      VFLT);
	setBranch("met_muons_dPhi3mu_jes_dwn", VFLT);
	
	setBranch("met_muons_et_jer_up",       FLT);
	setBranch("met_muons_phi_jer_up",      FLT);
	setBranch("met_muons_mT_jer_up",       VFLT);
	setBranch("met_muons_dPhi3mu_jer_up",  VFLT);
	
	setBranch("met_muons_et_jer_dwn",      FLT);
	setBranch("met_muons_phi_jer_dwn",     FLT);
	setBranch("met_muons_mT_jer_dwn",      VFLT);
	setBranch("met_muons_dPhi3mu_jer_dwn", VFLT);
	
	// MET_Track-based kinematics
	setBranch("met_track_et_uncalib",      FLT);
	setBranch("met_track_phi_uncalib",     FLT);
	setBranch("met_track_mT_uncalib",      VFLT);
	setBranch("met_track_dPhi3mu_uncalib", VFLT);
	
	setBranch("met_track_et",      FLT);
	setBranch("met_track_phi",     FLT);
	setBranch("met_track_mT",      VFLT);
	setBranch("met_track_dPhi3mu", VFLT);
	
	setBranch("met_track_et_jes_up",        FLT);
	setBranch("met_track_phi_jes_up",       FLT);
	setBranch("met_track_mT_jes_up",        VFLT);
	setBranch("met_track_dPhi3mu_jes_up",   VFLT);
	setBranch("met_track_et_jes_dwn",      FLT);
	setBranch("met_track_phi_jes_dwn",     FLT);
	setBranch("met_track_mT_jes_dwn",      VFLT);
	setBranch("met_track_dPhi3mu_jes_dwn", VFLT);
	
	setBranch("met_track_et_jer_up",        FLT);
	setBranch("met_track_phi_jer_up",       FLT);
	setBranch("met_track_mT_jer_up",        VFLT);
	setBranch("met_track_dPhi3mu_jer_up",   VFLT);
	setBranch("met_track_et_jer_dwn",      FLT);
	setBranch("met_track_phi_jer_dwn",     FLT);
	setBranch("met_track_mT_jer_dwn",      VFLT);
	setBranch("met_track_dPhi3mu_jer_dwn", VFLT);
	
	setBranch("met_track_et_jettrk_nom",        FLT);
	setBranch("met_track_phi_jettrk_nom",       FLT);
	setBranch("met_track_mT_jettrk_nom",        VFLT);
	setBranch("met_track_dPhi3mu_jettrk_nom",   VFLT);
	setBranch("met_track_et_jettrk_up",        FLT);
	setBranch("met_track_phi_jettrk_up",       FLT);
	setBranch("met_track_mT_jettrk_up",        VFLT);
	setBranch("met_track_dPhi3mu_jettrk_up",   VFLT);
	setBranch("met_track_et_jettrk_dwn",        FLT);
	setBranch("met_track_phi_jettrk_dwn",       FLT);
	setBranch("met_track_mT_jettrk_dwn",        VFLT);
	setBranch("met_track_dPhi3mu_jettrk_dwn",   VFLT);
	
	setBranch("met_track_et_softtrk_nom",        FLT);
	setBranch("met_track_phi_softtrk_nom",       FLT);
	setBranch("met_track_mT_softtrk_nom",        VFLT);
	setBranch("met_track_dPhi3mu_softtrk_nom",   VFLT);
	setBranch("met_track_et_softtrk_up",        FLT);
	setBranch("met_track_phi_softtrk_up",       FLT);
	setBranch("met_track_mT_softtrk_up",        VFLT);
	setBranch("met_track_dPhi3mu_softtrk_up",   VFLT);
	setBranch("met_track_et_softtrk_dwn",        FLT);
	setBranch("met_track_phi_softtrk_dwn",       FLT);
	setBranch("met_track_mT_softtrk_dwn",        VFLT);
	setBranch("met_track_dPhi3mu_softtrk_dwn",   VFLT);
	
	setBranch("met_track_et_softtrkres_para",        FLT);
	setBranch("met_track_phi_softtrkres_para",       FLT);
	setBranch("met_track_mT_softtrkres_para",        VFLT);
	setBranch("met_track_dPhi3mu_softtrkres_para",   VFLT);
	setBranch("met_track_et_softtrkres_perp",        FLT);
	setBranch("met_track_phi_softtrkres_perp",       FLT);
	setBranch("met_track_mT_softtrkres_perp",        VFLT);
	setBranch("met_track_dPhi3mu_softtrkres_perp",   VFLT);
	setBranch("met_track_et_softtrkres_corr",        FLT);
	setBranch("met_track_phi_softtrkres_corr",       FLT);
	setBranch("met_track_mT_softtrkres_corr",        VFLT);
	setBranch("met_track_dPhi3mu_softtrkres_corr",   VFLT);
	
	// Muons
	setBranch("mu_order1",          VINT);
	setBranch("mu_order2",          VINT);
	setBranch("mu_order3",          VINT);
	setBranch("mu_type1",           VINT);
	setBranch("mu_type2",           VINT);
	setBranch("mu_type3",           VINT);
	setBranch("mu_pt1",             VFLT);
	setBranch("mu_pt2",             VFLT);
	setBranch("mu_pt3",             VFLT);
	setBranch("mu_eta1",            VFLT);
	setBranch("mu_eta2",            VFLT);
	setBranch("mu_eta3",            VFLT);
	setBranch("mu_phi1",            VFLT);
	setBranch("mu_phi2",            VFLT);
	setBranch("mu_phi3",            VFLT);
	setBranch("mu_sctangsig1",      VFLT);
	setBranch("mu_sctangsig2",      VFLT);
	setBranch("mu_sctangsig3",      VFLT);
	setBranch("mu_sctngbsig1",      VFLT);
	setBranch("mu_sctngbsig2",      VFLT);
	setBranch("mu_sctngbsig3",      VFLT);
	setBranch("mu_pbalsig1",        VFLT);
	setBranch("mu_pbalsig2",        VFLT);
	setBranch("mu_pbalsig3",        VFLT);
	setBranch("mu_chi2trkfit1",     VFLT);
	setBranch("mu_chi2trkfit2",     VFLT);
	setBranch("mu_chi2trkfit3",     VFLT);
	setBranch("mu_ndftrkfit1",      VFLT);
	setBranch("mu_ndftrkfit2",      VFLT);
	setBranch("mu_ndftrkfit3",      VFLT);
	setBranch("mu_chi2ndftrkfit1",  VFLT);
	setBranch("mu_chi2ndftrkfit2",  VFLT);
	setBranch("mu_chi2ndftrkfit3",  VFLT);
	setBranch("mu_pvaltrkfit1",     VFLT);
	setBranch("mu_pvaltrkfit2",     VFLT);
	setBranch("mu_pvaltrkfit3",     VFLT);
	setBranch("mu_srcqoverp1",      VFLT);
	setBranch("mu_srcqoverp2",      VFLT);
	setBranch("mu_srcqoverp3",      VFLT);
	setBranch("mu_trkqoverp1",      VFLT);
	setBranch("mu_trkqoverp2",      VFLT);
	setBranch("mu_trkqoverp3",      VFLT);
	setBranch("mu_ptfrac1",         VFLT);
	setBranch("mu_ptfrac2",         VFLT);
	setBranch("mu_ptfrac3",         VFLT);
	setBranch("mu_pixeldEdx1",      VFLT);
	setBranch("mu_pixeldEdx2",      VFLT);
	setBranch("mu_pixeldEdx3",      VFLT);
	setBranch("mu_isMedium1",       VINT);
	setBranch("mu_isMedium2",       VINT);
	setBranch("mu_isMedium3",       VINT);
	setBranch("mu_nPIXhits1",       VINT);
	setBranch("mu_nPIXhits2",       VINT);
	setBranch("mu_nPIXhits3",       VINT);
	setBranch("mu_nDeadPIX1",       VINT);
	setBranch("mu_nDeadPIX2",       VINT);
	setBranch("mu_nDeadPIX3",       VINT);
	setBranch("mu_nPIXholes1",      VINT);
	setBranch("mu_nPIXholes2",      VINT);
	setBranch("mu_nPIXholes3",      VINT);
	setBranch("mu_nSCThits1",       VINT);
	setBranch("mu_nSCThits2",       VINT);
	setBranch("mu_nSCThits3",       VINT);
	setBranch("mu_nDeadSCT1",       VINT);
	setBranch("mu_nDeadSCT2",       VINT);
	setBranch("mu_nDeadSCT3",       VINT);
	setBranch("mu_nSCTholes1",      VINT);
	setBranch("mu_nSCTholes2",      VINT);
	setBranch("mu_nSCTholes3",      VINT);
	setBranch("mu_nTRThits1",       VINT);
	setBranch("mu_nTRThits2",       VINT);
	setBranch("mu_nTRThits3",       VINT);
	setBranch("mu_nTRToutliers1",   VINT);
	setBranch("mu_nTRToutliers2",   VINT);
	setBranch("mu_nTRToutliers3",   VINT);
	setBranch("mu_htTRThits1",      VINT);
	setBranch("mu_htTRThits2",      VINT);
	setBranch("mu_htTRThits3",      VINT);
	setBranch("mu_nUsedHitsdEdx1",  VINT);	
	setBranch("mu_nUsedHitsdEdx2",  VINT);	
	setBranch("mu_nUsedHitsdEdx3",  VINT);               
	setBranch("mu_nMDThits1",                  VINT);
	setBranch("mu_nMDThits2",                  VINT);
	setBranch("mu_nMDThits3",                  VINT);
	setBranch("mu_nTGCPhiHits1",               VINT);
	setBranch("mu_nTGCPhiHits2",               VINT);
	setBranch("mu_nTGCPhiHits3",               VINT);
	setBranch("mu_nTGCEtaHits1",               VINT);
	setBranch("mu_nTGCEtaHits2",               VINT);
	setBranch("mu_nTGCEtaHits3",               VINT);
	setBranch("mu_nCSCPhiHits1",               VINT);
	setBranch("mu_nCSCPhiHits2",               VINT);
	setBranch("mu_nCSCPhiHits3",               VINT);
	setBranch("mu_nCSCEtaHits1",               VINT);
	setBranch("mu_nCSCEtaHits2",               VINT);
	setBranch("mu_nCSCEtaHits3",               VINT);
	setBranch("mu_nRPCPhiHits1",               VINT);
	setBranch("mu_nRPCPhiHits2",               VINT);
	setBranch("mu_nRPCPhiHits3",               VINT);
	setBranch("mu_nRPCEtaHits1",               VINT);
	setBranch("mu_nRPCEtaHits2",               VINT);
	setBranch("mu_nRPCEtaHits3",               VINT);
	setBranch("mu_nCSCEtaHoles1",              VINT);
	setBranch("mu_nCSCEtaHoles2",              VINT);
	setBranch("mu_nCSCEtaHoles3",              VINT);
	setBranch("mu_nCSCPhiHoles1",              VINT);
	setBranch("mu_nCSCPhiHoles2",              VINT);
	setBranch("mu_nCSCPhiHoles3",              VINT);
	setBranch("mu_nRPCEtaHoles1",              VINT);
	setBranch("mu_nRPCEtaHoles2",              VINT);
	setBranch("mu_nRPCEtaHoles3",              VINT);
	setBranch("mu_nRPCPhiHoles1",              VINT);
	setBranch("mu_nRPCPhiHoles2",              VINT);
	setBranch("mu_nRPCPhiHoles3",              VINT);
	setBranch("mu_nMDTHoles1",                 VINT);
	setBranch("mu_nMDTHoles2",                 VINT);
	setBranch("mu_nMDTHoles3",                 VINT);
	setBranch("mu_nTGCEtaHoles1",              VINT);
	setBranch("mu_nTGCEtaHoles2",              VINT);
	setBranch("mu_nTGCEtaHoles3",              VINT);
	setBranch("mu_nTGCPhiHoles1",              VINT);
	setBranch("mu_nTGCPhiHoles2",              VINT);
	setBranch("mu_nTGCPhiHoles3",              VINT);
	setBranch("mu_nOutliersOnTrack1",          VINT);
	setBranch("mu_nOutliersOnTrack2",          VINT);
	setBranch("mu_nOutliersOnTrack3",          VINT);
	setBranch("mu_standardDeviationOfChi2OS1", VINT);
	setBranch("mu_standardDeviationOfChi2OS2", VINT);
	setBranch("mu_standardDeviationOfChi2OS3", VINT);

	setBranch("mu_nPrecisionHits1",        VINT);
	setBranch("mu_nPrecisionHits2",        VINT);
	setBranch("mu_nPrecisionHits3",        VINT);
	setBranch("mu_nPhiLayers1",            VINT);
	setBranch("mu_nPhiLayers2",            VINT);
	setBranch("mu_nPhiLayers3",            VINT);
	setBranch("mu_nEtaPhiLayers1",         VINT);
	setBranch("mu_nEtaPhiLayers2",         VINT);
	setBranch("mu_nEtaPhiLayers3",         VINT);
	setBranch("mu_nPrecisionHoles1",       VINT);
	setBranch("mu_nPrecisionHoles2",       VINT);
	setBranch("mu_nPrecisionHoles3",       VINT);
	setBranch("mu_nEtaTriggerHoleLayers1", VINT);
	setBranch("mu_nEtaTriggerHoleLayers2", VINT);
	setBranch("mu_nEtaTriggerHoleLayers3", VINT);
	setBranch("mu_nPhiHoleLayers1",        VINT);
	setBranch("mu_nPhiHoleLayers2",        VINT);
	setBranch("mu_nPhiHoleLayers3",        VINT);
	setBranch("mu_nPrecisionOutliers1",    VINT);
	setBranch("mu_nPrecisionOutliers2",    VINT);
	setBranch("mu_nPrecisionOutliers3",    VINT);
	
	// Jets
	setBranch("jet_pt1_uncalib", VFLT);
	setBranch("jet_pt2_uncalib", VFLT);
	setBranch("jet_pt3_uncalib", VFLT);
	setBranch("jet_pt4_uncalib", VFLT);
	setBranch("jet_pt1",         VFLT);
	setBranch("jet_pt2",         VFLT);
	setBranch("jet_pt3",         VFLT);
	setBranch("jet_pt4",         VFLT);
	setBranch("jet_pt1_jes_up",  VFLT);
	setBranch("jet_pt2_jes_up",  VFLT);
	setBranch("jet_pt3_jes_up",  VFLT);
	setBranch("jet_pt4_jes_up",  VFLT);
	setBranch("jet_pt1_jes_dwn", VFLT);
	setBranch("jet_pt2_jes_dwn", VFLT);
	setBranch("jet_pt3_jes_dwn", VFLT);
	setBranch("jet_pt4_jes_dwn", VFLT);
	setBranch("jet_pt1_jer_up",  VFLT);
	setBranch("jet_pt2_jer_up",  VFLT);
	setBranch("jet_pt3_jer_up",  VFLT);
	setBranch("jet_pt4_jer_up",  VFLT);
	setBranch("jet_pt1_jer_dwn", VFLT);
	setBranch("jet_pt2_jer_dwn", VFLT);
	setBranch("jet_pt3_jer_dwn", VFLT);
	setBranch("jet_pt4_jer_dwn", VFLT);
	
	setBranch("jet_E1_uncalib", VFLT);
	setBranch("jet_E2_uncalib", VFLT);
	setBranch("jet_E3_uncalib", VFLT);
	setBranch("jet_E4_uncalib", VFLT);
	setBranch("jet_E1",         VFLT);
	setBranch("jet_E2",         VFLT);
	setBranch("jet_E3",         VFLT);
	setBranch("jet_E4",         VFLT);
	setBranch("jet_E1_jes_up",  VFLT);
	setBranch("jet_E2_jes_up",  VFLT);
	setBranch("jet_E3_jes_up",  VFLT);
	setBranch("jet_E4_jes_up",  VFLT);
	setBranch("jet_E1_jes_dwn", VFLT);
	setBranch("jet_E2_jes_dwn", VFLT);
	setBranch("jet_E3_jes_dwn", VFLT);
	setBranch("jet_E4_jes_dwn", VFLT);
	setBranch("jet_E1_jer_up",  VFLT);
	setBranch("jet_E2_jer_up",  VFLT);
	setBranch("jet_E3_jer_up",  VFLT);
	setBranch("jet_E4_jer_up",  VFLT);
	setBranch("jet_E1_jer_dwn", VFLT);
	setBranch("jet_E2_jer_dwn", VFLT);
	setBranch("jet_E3_jer_dwn", VFLT);
	setBranch("jet_E4_jer_dwn", VFLT);
	
	setBranch("jet_m1_uncalib", VFLT);
	setBranch("jet_m2_uncalib", VFLT);
	setBranch("jet_m3_uncalib", VFLT);
	setBranch("jet_m4_uncalib", VFLT);
	setBranch("jet_m1",         VFLT);
	setBranch("jet_m2",         VFLT);
	setBranch("jet_m3",         VFLT);
	setBranch("jet_m4",         VFLT);
	setBranch("jet_m1_jes_up",  VFLT);
	setBranch("jet_m2_jes_up",  VFLT);
	setBranch("jet_m3_jes_up",  VFLT);
	setBranch("jet_m4_jes_up",  VFLT);
	setBranch("jet_m1_jes_dwn", VFLT);
	setBranch("jet_m2_jes_dwn", VFLT);
	setBranch("jet_m3_jes_dwn", VFLT);
	setBranch("jet_m4_jes_dwn", VFLT);
	setBranch("jet_m1_jer_up",  VFLT);
	setBranch("jet_m2_jer_up",  VFLT);
	setBranch("jet_m3_jer_up",  VFLT);
	setBranch("jet_m4_jer_up",  VFLT);
	setBranch("jet_m1_jer_dwn", VFLT);
	setBranch("jet_m2_jer_dwn", VFLT);
	setBranch("jet_m3_jer_dwn", VFLT);
	setBranch("jet_m4_jer_dwn", VFLT);
	
	setBranch("jet_eta1_uncalib", VFLT);
	setBranch("jet_eta2_uncalib", VFLT);
	setBranch("jet_eta3_uncalib", VFLT);
	setBranch("jet_eta4_uncalib", VFLT);
	setBranch("jet_eta1",         VFLT);
	setBranch("jet_eta2",         VFLT);
	setBranch("jet_eta3",         VFLT);
	setBranch("jet_eta4",         VFLT);
	setBranch("jet_eta1_jes_up",  VFLT);
	setBranch("jet_eta2_jes_up",  VFLT);
	setBranch("jet_eta3_jes_up",  VFLT);
	setBranch("jet_eta4_jes_up",  VFLT);
	setBranch("jet_eta1_jes_dwn", VFLT);
	setBranch("jet_eta2_jes_dwn", VFLT);
	setBranch("jet_eta3_jes_dwn", VFLT);
	setBranch("jet_eta4_jes_dwn", VFLT);
	setBranch("jet_eta1_jer_up",  VFLT);
	setBranch("jet_eta2_jer_up",  VFLT);
	setBranch("jet_eta3_jer_up",  VFLT);
	setBranch("jet_eta4_jer_up",  VFLT);
	setBranch("jet_eta1_jer_dwn", VFLT);
	setBranch("jet_eta2_jer_dwn", VFLT);
	setBranch("jet_eta3_jer_dwn", VFLT);
	setBranch("jet_eta4_jer_dwn", VFLT);
	
	setBranch("jet_phi1_uncalib", VFLT);
	setBranch("jet_phi2_uncalib", VFLT);
	setBranch("jet_phi3_uncalib", VFLT);
	setBranch("jet_phi4_uncalib", VFLT);
	setBranch("jet_phi1",         VFLT);
	setBranch("jet_phi2",         VFLT);
	setBranch("jet_phi3",         VFLT);
	setBranch("jet_phi4",         VFLT);
	setBranch("jet_phi1_jes_up",  VFLT);
	setBranch("jet_phi2_jes_up",  VFLT);
	setBranch("jet_phi3_jes_up",  VFLT);
	setBranch("jet_phi4_jes_up",  VFLT);
	setBranch("jet_phi1_jes_dwn", VFLT);
	setBranch("jet_phi2_jes_dwn", VFLT);
	setBranch("jet_phi3_jes_dwn", VFLT);
	setBranch("jet_phi4_jes_dwn", VFLT);
	setBranch("jet_phi1_jer_up",  VFLT);
	setBranch("jet_phi2_jer_up",  VFLT);
	setBranch("jet_phi3_jer_up",  VFLT);
	setBranch("jet_phi4_jer_up",  VFLT);
	setBranch("jet_phi1_jer_dwn", VFLT);
	setBranch("jet_phi2_jer_dwn", VFLT);
	setBranch("jet_phi3_jer_dwn", VFLT);
	setBranch("jet_phi4_jer_dwn", VFLT);
	
	setBranch("jet_MV1w1_uncalib", VFLT);
	setBranch("jet_MV1w2_uncalib", VFLT);
	setBranch("jet_MV1w3_uncalib", VFLT);
	setBranch("jet_MV1w4_uncalib", VFLT);
	setBranch("jet_MV1w1",         VFLT);
	setBranch("jet_MV1w2",         VFLT);
	setBranch("jet_MV1w3",         VFLT);
	setBranch("jet_MV1w4",         VFLT);
	setBranch("jet_MV1w1_jes_up",  VFLT);
	setBranch("jet_MV1w2_jes_up",  VFLT);
	setBranch("jet_MV1w3_jes_up",  VFLT);
	setBranch("jet_MV1w4_jes_up",  VFLT);
	setBranch("jet_MV1w1_jes_dwn", VFLT);
	setBranch("jet_MV1w2_jes_dwn", VFLT);
	setBranch("jet_MV1w3_jes_dwn", VFLT);
	setBranch("jet_MV1w4_jes_dwn", VFLT);
	setBranch("jet_MV1w1_jer_up",  VFLT);
	setBranch("jet_MV1w2_jer_up",  VFLT);
	setBranch("jet_MV1w3_jer_up",  VFLT);
	setBranch("jet_MV1w4_jer_up",  VFLT);
	setBranch("jet_MV1w1_jer_dwn", VFLT);
	setBranch("jet_MV1w2_jer_dwn", VFLT);
	setBranch("jet_MV1w3_jer_dwn", VFLT);
	setBranch("jet_MV1w4_jer_dwn", VFLT);
	
	setBranch("jet_vtxf1_uncalib", VFLT);
	setBranch("jet_vtxf2_uncalib", VFLT);
	setBranch("jet_vtxf3_uncalib", VFLT);
	setBranch("jet_vtxf4_uncalib", VFLT);
	setBranch("jet_vtxf1",         VFLT);
	setBranch("jet_vtxf2",         VFLT);
	setBranch("jet_vtxf3",         VFLT);
	setBranch("jet_vtxf4",         VFLT);
	setBranch("jet_vtxf1_jes_up",  VFLT);
	setBranch("jet_vtxf2_jes_up",  VFLT);
	setBranch("jet_vtxf3_jes_up",  VFLT);
	setBranch("jet_vtxf4_jes_up",  VFLT);
	setBranch("jet_vtxf1_jes_dwn", VFLT);
	setBranch("jet_vtxf2_jes_dwn", VFLT);
	setBranch("jet_vtxf3_jes_dwn", VFLT);
	setBranch("jet_vtxf4_jes_dwn", VFLT);
	setBranch("jet_vtxf1_jer_up",  VFLT);
	setBranch("jet_vtxf2_jer_up",  VFLT);
	setBranch("jet_vtxf3_jer_up",  VFLT);
	setBranch("jet_vtxf4_jer_up",  VFLT);
	setBranch("jet_vtxf1_jer_dwn", VFLT);
	setBranch("jet_vtxf2_jer_dwn", VFLT);
	setBranch("jet_vtxf3_jer_dwn", VFLT);
	setBranch("jet_vtxf4_jer_dwn", VFLT);
	
	setBranch("jet_ntrk1_uncalib", VINT);
	setBranch("jet_ntrk2_uncalib", VINT);
	setBranch("jet_ntrk3_uncalib", VINT);
	setBranch("jet_ntrk4_uncalib", VINT);
	setBranch("jet_ntrk1",         VINT);
	setBranch("jet_ntrk2",         VINT);
	setBranch("jet_ntrk3",         VINT);
	setBranch("jet_ntrk4",         VINT);
	setBranch("jet_ntrk1_jes_up",  VINT);
	setBranch("jet_ntrk2_jes_up",  VINT);
	setBranch("jet_ntrk3_jes_up",  VINT);
	setBranch("jet_ntrk4_jes_up",  VINT);
	setBranch("jet_ntrk1_jes_dwn", VINT);
	setBranch("jet_ntrk2_jes_dwn", VINT);
	setBranch("jet_ntrk3_jes_dwn", VINT);
	setBranch("jet_ntrk4_jes_dwn", VINT);
	setBranch("jet_ntrk1_jer_up",  VINT);
	setBranch("jet_ntrk2_jer_up",  VINT);
	setBranch("jet_ntrk3_jer_up",  VINT);
	setBranch("jet_ntrk4_jer_up",  VINT);
	setBranch("jet_ntrk1_jer_dwn", VINT);
	setBranch("jet_ntrk2_jer_dwn", VINT);
	setBranch("jet_ntrk3_jer_dwn", VINT);
	setBranch("jet_ntrk4_jer_dwn", VINT);
	
	setBranch("jet_JES_shift1", VFLT);
	setBranch("jet_JES_shift2", VFLT);
	setBranch("jet_JES_shift3", VFLT);
	setBranch("jet_JES_shift4", VFLT);
	setBranch("jet_JER_shift1", VFLT);
	setBranch("jet_JER_shift2", VFLT);
	setBranch("jet_JER_shift3", VFLT);
	setBranch("jet_JER_shift4", VFLT);
	
	setBranch("jet_sumpt12_uncalib",   VFLT);
	setBranch("jet_dphi3muJ1_uncalib", VFLT);
	setBranch("jet_dR3muJ1_uncalib",   VFLT);
	setBranch("jet_dphiJ1J2_uncalib",  VFLT);
	setBranch("jet_dRJ1J2_uncalib",    VFLT);
	setBranch("jet_sumpt12",           VFLT);
	setBranch("jet_dphi3muJ1",         VFLT);
	setBranch("jet_dR3muJ1",           VFLT);
	setBranch("jet_dphiJ1J2",          VFLT);
	setBranch("jet_dRJ1J2",            VFLT);
	setBranch("jet_sumpt12_jes_up",    VFLT);
	setBranch("jet_dphi3muJ1_jes_up",  VFLT);
	setBranch("jet_dR3muJ1_jes_up",    VFLT);
	setBranch("jet_dphiJ1J2_jes_up",   VFLT);
	setBranch("jet_dRJ1J2_jes_up",     VFLT);
	setBranch("jet_sumpt12_jes_dwn",   VFLT);
	setBranch("jet_dphi3muJ1_jes_dwn", VFLT);
	setBranch("jet_dR3muJ1_jes_dwn",   VFLT);
	setBranch("jet_dphiJ1J2_jes_dwn",  VFLT);
	setBranch("jet_dRJ1J2_jes_dwn",    VFLT);
	setBranch("jet_sumpt12_jer_up",    VFLT);
	setBranch("jet_dphi3muJ1_jer_up",  VFLT);
	setBranch("jet_dR3muJ1_jer_up",    VFLT);
	setBranch("jet_dphiJ1J2_jer_up",   VFLT);
	setBranch("jet_dRJ1J2_jer_up",     VFLT);
	setBranch("jet_sumpt12_jer_dwn",   VFLT);
	setBranch("jet_dphi3muJ1_jer_dwn", VFLT);
	setBranch("jet_dR3muJ1_jer_dwn",   VFLT);
	setBranch("jet_dphiJ1J2_jer_dwn",  VFLT);
	setBranch("jet_dRJ1J2_jer_dwn",    VFLT);
	
	// Attach the branches
	for(TMapTSi::iterator it=branchtypes.begin() ; it!=branchtypes.end() ; ++it)
	{
		TString name = it->first;
		int     type = it->second;
		if     (type==INT)  itree->SetBranchAddress(name, &ints[name]);
		else if(type==FLT)  itree->SetBranchAddress(name, &floats[name]);
		else if(type==STR)  itree->SetBranchAddress(name, &strings[name]);
		else if(type==VINT) itree->SetBranchAddress(name, &vints[name]);
		else if(type==VFLT) itree->SetBranchAddress(name, &vfloats[name]);
		else if(type==VSTR) itree->SetBranchAddress(name, &vstrings[name]);
		else _FAT("unsupported type of branch");
	}
}

#endif