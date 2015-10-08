//////////////////////////////////////////////
//// run ./execute here and follow ///////////
//// the instructions ////////////////////////
//////////////////////////////////////////////

#include "../PileupReweighting/PileupReweighting/TPileupReweighting.h"
#include "../ApplyJetCalibration/ApplyJetCalibration/ApplyJetCalibration.h"
#include "../JetUncertainties/JetUncertainties/MultijetJESUncertaintyProvider.h"
#include "../ApplyJetResolutionSmearing/ApplyJetResolutionSmearing/ApplyJetSmearing.h"
#include "../JetUncertainties/JetUncertainties/JESUncertaintyProvider.h"
#include "../MissingETUtility/MissingETUtility/METUtility.h"
#include "../METSystematics/METSystematics/METSystTool.h"
#include "../METTrackUtil/METTrackUtil/TrackMETMaker.h"
#include "../TileTripReader/TileTripReader/TTileTripReader.h"
#include "../BCHCleaningTool/BCHCleaningTool/BCHCleaningToolRoot.h"


enum systshift
{
	SHIFTDWN=-1,
	NOSHIFT,
	SHIFTUP
};

enum uncertainties
{
	NOJES,
	NOMINAL,
	JESUP,JESDWN,
	JERUP,JERDWN,
	JETTRKNOM, JETTRKUP, JETTRKDWN,
	SOFTTRKNOM, SOFTTRKUP, SOFTTRKDWN,
	SOFTTRKRESPARA, SOFTTRKRESPERP, SOFTTRKRESCORR
};

enum met_types
{
	METSTACO,METMUONS,METTRACK
};


TString JESconfig  = "Full2012"; // "Moriond" or "Full2012"
TString JetQuality = "VeryLooseBad";

vector<TLorentzVector> calibJets;
vector<bool>           badBCHJets;
vector<TLorentzVector> smearedJets;
vector<double>         calibJetsUnc;
vector<double>         smearedJetsUnc;
vector<unsigned int>   calibJetsIndex;

vector<float>* AntiKt4LCTopoJets_calibrated_pt;
vector<float>* AntiKt4LCTopoJets_calibrated_E;
vector<float>* AntiKt4LCTopoJets_calibrated_m;
vector<float>* AntiKt4LCTopoJets_calibrated_eta;
vector<float>* AntiKt4LCTopoJets_calibrated_phi;

vector<float>* akt4lc_jet_pt;
vector<float>* akt4lc_jet_E;
vector<float>* akt4lc_jet_eta;
vector<float>* akt4lc_jet_phi;

METUtil::METObject uncalibMET;
METUtil::METObject calibMET_nominal;
METUtil::METObject calibMET_jes_up;
METUtil::METObject calibMET_jes_dwn;
METUtil::METObject calibMET_jer_up;
METUtil::METObject calibMET_jer_dwn;
METUtil::METObject uncalibMUMET;
METUtil::METObject calibMUMET_nominal;
METUtil::METObject calibMUMET_jes_up;
METUtil::METObject calibMUMET_jes_dwn;
METUtil::METObject calibMUMET_jer_up;
METUtil::METObject calibMUMET_jer_dwn;
METUtil::METObject uncalibMETTRK;
METUtil::METObject calibMETTRK_nominal;

METUtil::METObject calibMETTRK_jettrk_up;
METUtil::METObject calibMETTRK_jettrk_dwn;
METUtil::METObject calibMETTRK_softtrk_up;
METUtil::METObject calibMETTRK_softtrk_dwn;
METUtil::METObject calibMETTRK_softtrkres_para;
METUtil::METObject calibMETTRK_softtrkres_perp;



int met_RefFinal_source;
double met_RefFinal_etX;
double met_RefFinal_etY;
double met_RefFinal_sumet;
double met_RefFinal_et;
double met_RefFinal_phi;

int met_Track_source;
double met_Track_etX;
double met_Track_etY;
double met_Track_sumet;
double met_Track_et;
double met_Track_phi;

vector<vector<int> >*    muons_cal_subCalo_subCaloId;
vector<vector<double> >* muons_cal_subCalo_energyDeposited;
vector<vector<double> >* muons_cal_subCalo_muonEnergyLoss;
vector<int>*             muons_cal_energyLossType;
vector<int>*             muons_cal_caloMuonIdTag;
vector<double>*          muons_cal_caloLRLikelihood;
vector<double>*          muons_cal_fsrCandidateEnergy;
vector<double>*          muons_cal_etCore;


int            AntiKt4LCTopoJets_n;
vector<float>* AntiKt4LCTopoJets_E;
vector<float>* AntiKt4LCTopoJets_pt;
vector<float>* AntiKt4LCTopoJets_m;
vector<float>* AntiKt4LCTopoJets_eta;
vector<float>* AntiKt4LCTopoJets_phi;
vector<float>* AntiKt4LCTopoJets_EtaOrigin;
vector<float>* AntiKt4LCTopoJets_PhiOrigin;
vector<float>* AntiKt4LCTopoJets_MOrigin;
vector<float>* AntiKt4LCTopoJets_WIDTH;
vector<float>* AntiKt4LCTopoJets_n90;
vector<float>* AntiKt4LCTopoJets_Timing;
vector<float>* AntiKt4LCTopoJets_LArQuality;
vector<float>* AntiKt4LCTopoJets_nTrk;
vector<float>* AntiKt4LCTopoJets_sumPtTrk;
vector<float>* AntiKt4LCTopoJets_OriginIndex;
vector<float>* AntiKt4LCTopoJets_HECQuality;
vector<float>* AntiKt4LCTopoJets_NegativeE;
vector<float>* AntiKt4LCTopoJets_AverageLArQF;
vector<float>* AntiKt4LCTopoJets_BCH_CORR_CELL;
vector<float>* AntiKt4LCTopoJets_BCH_CORR_DOTX;
vector<float>* AntiKt4LCTopoJets_BCH_CORR_JET;
vector<float>* AntiKt4LCTopoJets_BCH_CORR_JET_FORCELL;
vector<float>* AntiKt4LCTopoJets_ENG_BAD_CELLS;
vector<float>* AntiKt4LCTopoJets_N_BAD_CELLS;
vector<float>* AntiKt4LCTopoJets_N_BAD_CELLS_CORR;
vector<float>* AntiKt4LCTopoJets_BAD_CELLS_CORR_E;
vector<float>* AntiKt4LCTopoJets_NumTowers;
vector<float>* AntiKt4LCTopoJets_ootFracCells5;
vector<float>* AntiKt4LCTopoJets_ootFracCells10;
vector<float>* AntiKt4LCTopoJets_ootFracClusters5;
vector<float>* AntiKt4LCTopoJets_ootFracClusters10;
vector<int>*   AntiKt4LCTopoJets_SamplingMax;
vector<float>* AntiKt4LCTopoJets_fracSamplingMax;
vector<float>* AntiKt4LCTopoJets_hecf;
vector<float>* AntiKt4LCTopoJets_tgap3f;
vector<int>*   AntiKt4LCTopoJets_isUgly;
vector<int>*   AntiKt4LCTopoJets_isBadLooseMinus;
vector<int>*   AntiKt4LCTopoJets_isBadLoose;
vector<int>*   AntiKt4LCTopoJets_isBadMedium;
vector<int>*   AntiKt4LCTopoJets_isBadTight;
vector<float>* AntiKt4LCTopoJets_emfrac;
vector<float>* AntiKt4LCTopoJets_Offset;
vector<float>* AntiKt4LCTopoJets_EMJES;
vector<float>* AntiKt4LCTopoJets_EMJES_EtaCorr;
vector<float>* AntiKt4LCTopoJets_EMJESnooffset;
vector<float>* AntiKt4LCTopoJets_LCJES;
vector<float>* AntiKt4LCTopoJets_LCJES_EtaCorr;
vector<float>* AntiKt4LCTopoJets_constscale_E;
vector<float>* AntiKt4LCTopoJets_constscale_pt;
vector<float>* AntiKt4LCTopoJets_constscale_m;
vector<float>* AntiKt4LCTopoJets_constscale_eta;
vector<float>* AntiKt4LCTopoJets_constscale_phi;
vector<float>* AntiKt4LCTopoJets_jvtx_x;
vector<float>* AntiKt4LCTopoJets_jvtx_y;
vector<float>* AntiKt4LCTopoJets_jvtx_z;
vector<float>* AntiKt4LCTopoJets_jvtxf;
vector<float>* AntiKt4LCTopoJets_GSCFactorF;
vector<float>* AntiKt4LCTopoJets_WidthFraction;
vector<float>* AntiKt4LCTopoJets_e_PreSamplerB;
vector<float>* AntiKt4LCTopoJets_e_EMB1;
vector<float>* AntiKt4LCTopoJets_e_EMB2;
vector<float>* AntiKt4LCTopoJets_e_EMB3;
vector<float>* AntiKt4LCTopoJets_e_PreSamplerE;
vector<float>* AntiKt4LCTopoJets_e_EME1;
vector<float>* AntiKt4LCTopoJets_e_EME2;
vector<float>* AntiKt4LCTopoJets_e_EME3;
vector<float>* AntiKt4LCTopoJets_e_HEC0;
vector<float>* AntiKt4LCTopoJets_e_HEC1;
vector<float>* AntiKt4LCTopoJets_e_HEC2;
vector<float>* AntiKt4LCTopoJets_e_HEC3;
vector<float>* AntiKt4LCTopoJets_e_TileBar0;
vector<float>* AntiKt4LCTopoJets_e_TileBar1;
vector<float>* AntiKt4LCTopoJets_e_TileBar2;
vector<float>* AntiKt4LCTopoJets_e_TileGap1;
vector<float>* AntiKt4LCTopoJets_e_TileGap2;
vector<float>* AntiKt4LCTopoJets_e_TileGap3;
vector<float>* AntiKt4LCTopoJets_e_TileExt0;
vector<float>* AntiKt4LCTopoJets_e_TileExt1;
vector<float>* AntiKt4LCTopoJets_e_TileExt2;
vector<float>* AntiKt4LCTopoJets_e_FCAL0;
vector<float>* AntiKt4LCTopoJets_e_FCAL1;
vector<float>* AntiKt4LCTopoJets_e_FCAL2;
vector<vector<float> >* AntiKt4LCTopoJets_shapeBins;
vector<int>*  AntiKt4LCTopoJets_Nconst;
vector<vector<float> >* AntiKt4LCTopoJets_ptconst_default;
vector<vector<float> >* AntiKt4LCTopoJets_econst_default;
vector<vector<float> >* AntiKt4LCTopoJets_etaconst_default;
vector<vector<float> >* AntiKt4LCTopoJets_phiconst_default;
vector<vector<float> >* AntiKt4LCTopoJets_weightconst_default;
vector<float>* AntiKt4LCTopoJets_emscale_E;
vector<float>* AntiKt4LCTopoJets_emscale_pt;
vector<float>* AntiKt4LCTopoJets_emscale_m;
vector<float>* AntiKt4LCTopoJets_emscale_eta;
vector<float>* AntiKt4LCTopoJets_emscale_phi;
vector<float>* AntiKt4LCTopoJets_LArBadHVEnergy;
vector<float>* AntiKt4LCTopoJets_LArBadHVRatio;
vector<float>* AntiKt4LCTopoJets_flavor_weight_Comb;
vector<float>* AntiKt4LCTopoJets_flavor_weight_IP2D;
vector<float>* AntiKt4LCTopoJets_flavor_weight_IP3D;
vector<float>* AntiKt4LCTopoJets_flavor_weight_SV0;
vector<float>* AntiKt4LCTopoJets_flavor_weight_SV1;
vector<float>* AntiKt4LCTopoJets_flavor_weight_SV2;
vector<float>* AntiKt4LCTopoJets_flavor_weight_SoftMuonTagChi2;
vector<float>* AntiKt4LCTopoJets_flavor_weight_SecondSoftMuonTagChi2;
vector<float>* AntiKt4LCTopoJets_flavor_weight_JetFitterTagNN;
vector<float>* AntiKt4LCTopoJets_flavor_weight_JetFitterCOMBNN;
vector<float>* AntiKt4LCTopoJets_flavor_weight_MV1;
vector<float>* AntiKt4LCTopoJets_flavor_weight_MV2;
vector<float>* AntiKt4LCTopoJets_flavor_weight_GbbNN;
vector<int>*   AntiKt4LCTopoJets_flavor_truth_label;
vector<float>* AntiKt4LCTopoJets_flavor_truth_dRminToB;
vector<float>* AntiKt4LCTopoJets_flavor_truth_dRminToC;
vector<float>* AntiKt4LCTopoJets_flavor_truth_dRminToT;
vector<int>*   AntiKt4LCTopoJets_flavor_truth_BHadronpdg;
vector<float>* AntiKt4LCTopoJets_flavor_truth_vx_x;
vector<float>* AntiKt4LCTopoJets_flavor_truth_vx_y;
vector<float>* AntiKt4LCTopoJets_flavor_truth_vx_z;
vector<int>*   AntiKt4LCTopoJets_flavor_putruth_label;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_dRminToB;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_dRminToC;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_dRminToT;
vector<int>*   AntiKt4LCTopoJets_flavor_putruth_BHadronpdg;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_vx_x;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_vx_y;
vector<float>* AntiKt4LCTopoJets_flavor_putruth_vx_z;
vector<float>* AntiKt4LCTopoJets_flavor_component_ip2d_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_ip2d_pb;
vector<int>*   AntiKt4LCTopoJets_flavor_component_ip2d_isValid;
vector<int>*   AntiKt4LCTopoJets_flavor_component_ip2d_ntrk;
vector<float>* AntiKt4LCTopoJets_flavor_component_ip3d_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_ip3d_pb;
vector<int>  * AntiKt4LCTopoJets_flavor_component_ip3d_isValid;
vector<int>  * AntiKt4LCTopoJets_flavor_component_ip3d_ntrk;
vector<int>  * AntiKt4LCTopoJets_flavor_component_jetprob_ntrk;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv1_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv1_pb;
vector<int>*   AntiKt4LCTopoJets_flavor_component_sv1_isValid;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv2_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv2_pb;
vector<int>*   AntiKt4LCTopoJets_flavor_component_sv2_isValid;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_pb;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_pc;
vector<int>*   AntiKt4LCTopoJets_flavor_component_jfit_isValid;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfitcomb_pu;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfitcomb_pb;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfitcomb_pc;
vector<int>*   AntiKt4LCTopoJets_flavor_component_jfitcomb_isValid;
vector<int>*   AntiKt4LCTopoJets_flavor_component_gbbnn_nMatchingTracks;
vector<double>* AntiKt4LCTopoJets_flavor_component_gbbnn_trkJetWidth;
vector<double>* AntiKt4LCTopoJets_flavor_component_gbbnn_trkJetMaxDeltaR;
vector<int>*   AntiKt4LCTopoJets_flavor_component_jfit_nvtx;
vector<int>*   AntiKt4LCTopoJets_flavor_component_jfit_nvtx1t;
vector<int>*   AntiKt4LCTopoJets_flavor_component_jfit_ntrkAtVx;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_efrc;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_mass;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_sig3d;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_deltaPhi;
vector<float>* AntiKt4LCTopoJets_flavor_component_jfit_deltaEta;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_d0val;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_d0sig;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_z0val;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_z0sig;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_w2D;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_w3D;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_pJP;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_pJPneg;
vector<vector<int> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_grade;
vector<vector<int> >* AntiKt4LCTopoJets_flavor_component_ipplus_trk_isFromV0;
vector<int>* AntiKt4LCTopoJets_flavor_component_svp_isValid;
vector<int>* AntiKt4LCTopoJets_flavor_component_svp_ntrkv;
vector<int>* AntiKt4LCTopoJets_flavor_component_svp_ntrkj;
vector<int>* AntiKt4LCTopoJets_flavor_component_svp_n2t;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_mass;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_efrc;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_x;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_y;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_z;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_err_x;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_err_y;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_err_z;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_cov_xy;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_cov_xz;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_cov_yz;
vector<float>* AntiKt4LCTopoJets_flavor_component_svp_chi2;
vector<int>*  AntiKt4LCTopoJets_flavor_component_svp_ndof;
vector<int>*  AntiKt4LCTopoJets_flavor_component_svp_ntrk;
vector<int>*  AntiKt4LCTopoJets_flavor_component_sv0p_isValid;
vector<int>*  AntiKt4LCTopoJets_flavor_component_sv0p_ntrkv;
vector<int>*  AntiKt4LCTopoJets_flavor_component_sv0p_ntrkj;
vector<int>*  AntiKt4LCTopoJets_flavor_component_sv0p_n2t;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_mass;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_efrc;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_x;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_y;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_z;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_err_x;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_err_y;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_err_z;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_cov_xy;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_cov_xz;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_cov_yz;
vector<float>* AntiKt4LCTopoJets_flavor_component_sv0p_chi2;
vector<int>*   AntiKt4LCTopoJets_flavor_component_sv0p_ndof;
vector<int>*   AntiKt4LCTopoJets_flavor_component_sv0p_ntrk;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuoninfo_muon_w;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuoninfo_muon_pTRel;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuoninfo_muon_dRJet;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuonchi2info_muon_w;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuonchi2info_muon_pTRel;
vector<vector<float> >* AntiKt4LCTopoJets_flavor_component_softmuonchi2info_muon_dRJet;
vector<float>* AntiKt4LCTopoJets_el_dr;
vector<int>*   AntiKt4LCTopoJets_el_matched;
vector<float>* AntiKt4LCTopoJets_mu_dr;
vector<int>*   AntiKt4LCTopoJets_mu_matched;
vector<float>* AntiKt4LCTopoJets_L1_dr;
vector<int>*   AntiKt4LCTopoJets_L1_matched;
vector<float>* AntiKt4LCTopoJets_L2_dr;
vector<int>*   AntiKt4LCTopoJets_L2_matched;
vector<float>* AntiKt4LCTopoJets_EF_dr;
vector<int>*   AntiKt4LCTopoJets_EF_matched;
vector<float>* AntiKt4LCTopoJets_LikeLihood_0;
vector<float>* AntiKt4LCTopoJets_ActiveArea;
vector<float>* AntiKt4LCTopoJets_ActiveAreaPx;
vector<float>* AntiKt4LCTopoJets_ActiveAreaPy;
vector<float>* AntiKt4LCTopoJets_ActiveAreaPz;
vector<float>* AntiKt4LCTopoJets_ActiveAreaE;
vector<float>* AntiKt4LCTopoJets_VoronoiArea;
vector<float>* AntiKt4LCTopoJets_VoronoiAreaPx;
vector<float>* AntiKt4LCTopoJets_VoronoiAreaPy;
vector<float>* AntiKt4LCTopoJets_VoronoiAreaPz;
vector<float>* AntiKt4LCTopoJets_VoronoiAreaE;
vector<float>* AntiKt4LCTopoJets_LowEtConstituentsFrac;
vector<float>* AntiKt4LCTopoJets_pt_truth;
vector<float>* AntiKt4LCTopoJets_IsoKR20Perp;
vector<float>* AntiKt4LCTopoJets_IsoKR20Par;
vector<float>* AntiKt4LCTopoJets_IsoKR20SumPt;
vector<float>* AntiKt4LCTopoJets_IsoDelta2Perp;
vector<float>* AntiKt4LCTopoJets_IsoDelta2Par;
vector<float>* AntiKt4LCTopoJets_IsoDelta2SumPt;
vector<float>* AntiKt4LCTopoJets_IsoFixedCone8Perp;
vector<float>* AntiKt4LCTopoJets_IsoFixedCone8Par;
vector<float>* AntiKt4LCTopoJets_IsoFixedCone8SumPt;
vector<float>* AntiKt4LCTopoJets_IsoFixedArea13Perp;
vector<float>* AntiKt4LCTopoJets_IsoFixedArea13Par;
vector<float>* AntiKt4LCTopoJets_IsoFixedArea13SumPt;
vector<float>* AntiKt4LCTopoJets_Iso6To88Perp;
vector<float>* AntiKt4LCTopoJets_Iso6To88Par;
vector<float>* AntiKt4LCTopoJets_Iso6To88SumPt;
vector<float>* AntiKt4LCTopoJets_KtDr;
vector<float>* AntiKt4LCTopoJets_Centroid_r;
vector<float>* AntiKt4LCTopoJets_nTrk_pv0_1GeV;
vector<float>* AntiKt4LCTopoJets_sumPtTrk_pv0_1GeV;
vector<float>* AntiKt4LCTopoJets_nTrk_allpv_1GeV;
vector<float>* AntiKt4LCTopoJets_sumPtTrk_allpv_1GeV;
vector<float>* AntiKt4LCTopoJets_nTrk_pv0_500MeV;
vector<float>* AntiKt4LCTopoJets_sumPtTrk_pv0_500MeV;
vector<float>* AntiKt4LCTopoJets_trackWIDTH_pv0_1GeV;
vector<float>* AntiKt4LCTopoJets_trackWIDTH_allpv_1GeV;
vector<float>* AntiKt4LCTopoJets_TrackMFindex;
vector<float>* AntiKt4LCTopoJets_TrackMFPt;
vector<float>* AntiKt4LCTopoJets_TruthMFindex;
vector<float>* AntiKt4LCTopoJets_TruthMFPt;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float MET_RefFinal_etx;
float MET_RefFinal_ety;
float MET_RefFinal_phi;
float MET_RefFinal_et;
float MET_RefFinal_sumet;

float MET_RefMuon_etx;
float MET_RefMuon_ety;
float MET_RefMuon_phi;
float MET_RefMuon_et;
float MET_RefMuon_sumet;

float MET_RefTau_etx;
float MET_RefTau_ety;
float MET_RefTau_phi;
float MET_RefTau_et;
float MET_RefTau_sumet;

float MET_CellOut_Eflow_etx;
float MET_CellOut_Eflow_ety;
float MET_CellOut_Eflow_phi;
float MET_CellOut_Eflow_et;
float MET_CellOut_Eflow_sumet;

float MET_SoftJets_etx;
float MET_SoftJets_ety;
float MET_SoftJets_phi;
float MET_SoftJets_et;
float MET_SoftJets_sumet;

float MET_Track_etx;
float MET_Track_ety;
float MET_Track_phi;
float MET_Track_et;
float MET_Track_sumet;

float MET_MuonBoy_etx;
float MET_MuonBoy_ety;
float MET_MuonBoy_phi;
float MET_MuonBoy_et;
float MET_MuonBoy_sumet;

float MET_Muon_etx;
float MET_Muon_ety;
float MET_Muon_phi;
float MET_Muon_et;
float MET_Muon_sumet;

float MET_Muons_etx;
float MET_Muons_ety;
float MET_Muons_phi;
float MET_Muons_et;
float MET_Muons_sumet;

float MET_Muid_etx;
float MET_Muid_ety;
float MET_Muid_phi;
float MET_Muid_et;
float MET_Muid_sumet;

float MET_RefGamma_etx;
float MET_RefGamma_ety;
float MET_RefGamma_phi;
float MET_RefGamma_et;
float MET_RefGamma_sumet;

float MET_RefEle_etx;
float MET_RefEle_ety;
float MET_RefEle_phi;
float MET_RefEle_et;
float MET_RefEle_sumet;

float MET_RefJet_etx;
float MET_RefJet_ety;
float MET_RefJet_phi;
float MET_RefJet_et;
float MET_RefJet_sumet;

float MET_Truth_NonInt_etx;
float MET_Truth_NonInt_ety;
float MET_Truth_NonInt_phi;
float MET_Truth_NonInt_et;
float MET_Truth_NonInt_sumet;

Int_t          trk_n;
vector<float>* trk_pt;
vector<float>* trk_eta;
vector<float>* trk_phi_wrtPV;


int  el_MET_RefFinal_comp_n;
vector<vector<float> > *el_MET_RefFinal_comp_wpx;
vector<vector<float> > *el_MET_RefFinal_comp_wpy;
vector<vector<float> > *el_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *el_MET_RefFinal_comp_statusWord;
int           ph_MET_RefFinal_comp_n;
vector<vector<float> > *ph_MET_RefFinal_comp_wpx;
vector<vector<float> > *ph_MET_RefFinal_comp_wpy;
vector<vector<float> > *ph_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *ph_MET_RefFinal_comp_statusWord;
int           mu_staco_MET_RefFinal_comp_n;
vector<vector<float> > *mu_staco_MET_RefFinal_comp_wpx;
vector<vector<float> > *mu_staco_MET_RefFinal_comp_wpy;
vector<vector<float> > *mu_staco_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *mu_staco_MET_RefFinal_comp_statusWord;
int           mu_muid_MET_RefFinal_comp_n;
vector<vector<float> > *mu_muid_MET_RefFinal_comp_wpx;
vector<vector<float> > *mu_muid_MET_RefFinal_comp_wpy;
vector<vector<float> > *mu_muid_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *mu_muid_MET_RefFinal_comp_statusWord;
int           mu_MET_RefFinal_comp_n;
vector<vector<float> > *mu_MET_RefFinal_comp_wpx;
vector<vector<float> > *mu_MET_RefFinal_comp_wpy;
vector<vector<float> > *mu_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *mu_MET_RefFinal_comp_statusWord;
int           tau_MET_RefFinal_comp_n;
vector<vector<float> > *tau_MET_RefFinal_comp_wpx;
vector<vector<float> > *tau_MET_RefFinal_comp_wpy;
vector<vector<float> > *tau_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *tau_MET_RefFinal_comp_statusWord;
int           jet_antikt4LCtopo_MET_RefFinal_comp_n;
vector<vector<float> > *jet_antikt4LCtopo_MET_RefFinal_comp_wpx;
vector<vector<float> > *jet_antikt4LCtopo_MET_RefFinal_comp_wpy;
vector<vector<float> > *jet_antikt4LCtopo_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *jet_antikt4LCtopo_MET_RefFinal_comp_statusWord;
int           cl_MET_RefFinal_comp_n;
vector<vector<float> > *cl_MET_RefFinal_comp_wpx;
vector<vector<float> > *cl_MET_RefFinal_comp_wpy;
vector<vector<float> > *cl_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *cl_MET_RefFinal_comp_statusWord;
int           trk_MET_RefFinal_comp_n;
vector<vector<float> > *trk_MET_RefFinal_comp_wpx;
vector<vector<float> > *trk_MET_RefFinal_comp_wpy;
vector<vector<float> > *trk_MET_RefFinal_comp_wet;
vector<vector<unsigned int> > *trk_MET_RefFinal_comp_statusWord;


int  el_MET_Muons_comp_n;
vector<vector<float> > *el_MET_Muons_comp_wpx;
vector<vector<float> > *el_MET_Muons_comp_wpy;
vector<vector<float> > *el_MET_Muons_comp_wet;
vector<vector<unsigned int> > *el_MET_Muons_comp_statusWord;
int           ph_MET_Muons_comp_n;
vector<vector<float> > *ph_MET_Muons_comp_wpx;
vector<vector<float> > *ph_MET_Muons_comp_wpy;
vector<vector<float> > *ph_MET_Muons_comp_wet;
vector<vector<unsigned int> > *ph_MET_Muons_comp_statusWord;
int           mu_staco_MET_Muons_comp_n;
vector<vector<float> > *mu_staco_MET_Muons_comp_wpx;
vector<vector<float> > *mu_staco_MET_Muons_comp_wpy;
vector<vector<float> > *mu_staco_MET_Muons_comp_wet;
vector<vector<unsigned int> > *mu_staco_MET_Muons_comp_statusWord;
int           mu_muid_MET_Muons_comp_n;
vector<vector<float> > *mu_muid_MET_Muons_comp_wpx;
vector<vector<float> > *mu_muid_MET_Muons_comp_wpy;
vector<vector<float> > *mu_muid_MET_Muons_comp_wet;
vector<vector<unsigned int> > *mu_muid_MET_Muons_comp_statusWord;
int           mu_MET_Muons_comp_n;
vector<vector<float> > *mu_MET_Muons_comp_wpx;
vector<vector<float> > *mu_MET_Muons_comp_wpy;
vector<vector<float> > *mu_MET_Muons_comp_wet;
vector<vector<unsigned int> > *mu_MET_Muons_comp_statusWord;
int           tau_MET_Muons_comp_n;
vector<vector<float> > *tau_MET_Muons_comp_wpx;
vector<vector<float> > *tau_MET_Muons_comp_wpy;
vector<vector<float> > *tau_MET_Muons_comp_wet;
vector<vector<unsigned int> > *tau_MET_Muons_comp_statusWord;
int           jet_antikt4LCtopo_MET_Muons_comp_n;
vector<vector<float> > *jet_antikt4LCtopo_MET_Muons_comp_wpx;
vector<vector<float> > *jet_antikt4LCtopo_MET_Muons_comp_wpy;
vector<vector<float> > *jet_antikt4LCtopo_MET_Muons_comp_wet;
vector<vector<unsigned int> > *jet_antikt4LCtopo_MET_Muons_comp_statusWord;
int           cl_MET_Muons_comp_n;
vector<vector<float> > *cl_MET_Muons_comp_wpx;
vector<vector<float> > *cl_MET_Muons_comp_wpy;
vector<vector<float> > *cl_MET_Muons_comp_wet;
vector<vector<unsigned int> > *cl_MET_Muons_comp_statusWord;
int           trk_MET_Muons_comp_n;
vector<vector<float> > *trk_MET_Muons_comp_wpx;
vector<vector<float> > *trk_MET_Muons_comp_wpy;
vector<vector<float> > *trk_MET_Muons_comp_wet;
vector<vector<unsigned int> > *trk_MET_Muons_comp_statusWord;


float rhorhoKt3EM;
float rhorhoKt4EM;
float rhorhoKt3LC;
float rhorhoKt4LC;

Int_t                   musp_n;
vector<float>*          musp_eta;
vector<float>*          musp_phi;
vector<unsigned short>* musp_trigHits;
vector<unsigned short>* musp_innerHits;
vector<unsigned short>* musp_middleHits;
vector<unsigned short>* musp_outerHits;
vector<unsigned short>* musp_innerSegments;
vector<unsigned short>* musp_middleSegments;
vector<unsigned short>* musp_outerSegments;

Int_t          vxp_n;
vector<float>* vxp_x;
vector<float>* vxp_y;
vector<float>* vxp_z;
vector<int>*   vxp_type;
vector<float>* vxp_chi2;
vector<int>*   vxp_ndof;
vector<float>* vxp_px;
vector<float>* vxp_py;
vector<float>* vxp_pz;
vector<float>* vxp_E;
vector<float>* vxp_m;
vector<int>*   vxp_nTracks;
vector<float>* vxp_sumPt;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

Root::TPileupReweighting* pileupTool;
void initializePileup()
{
	pileupTool = new Root::TPileupReweighting("pileuptool");
	if(makepufile)
	{
		_INFO("making pileup file");
		pileupTool->UsePeriodConfig("MC12b");
	}
	else
	{
		_INFO("reading pileup file");
		pileupTool->AddConfigFile("Wtaunu_3mu.prw.root");
		pileupTool->AddLumiCalcFile("ilumicalc_histograms_None_200842-215643.root");
	}
	pileupTool->Initialize();
}
float getPileupWeight()
{
	/////////////////////////////
	float pileup_weight = 1.; ///
	/////////////////////////////

	// NOTE (23/01/2013): A bug has been found in the d3pd making code,
	// causing all MC12 samples to have a few of the averageIntPerXing
	// values incorrectly set (some should be 0 but are set to 1).
	// The bug does not affect data. To resolve this, when reading this branch,
	// for both prw file generating and for when retrieving pileup weights,
	// you should amend the value with the following line of code:

	float averageIntPerXing_fixed = (lbn==1 && int(phys_averageIntPerXing+0.5)==1) ? 0. : phys_averageIntPerXing;
	if(makepufile) pileupTool->Fill(phys_RunNumber,phys_mc_channel_number,phys_mc_event_weight,averageIntPerXing_fixed);
	else           pileup_weight = pileupTool->GetCombinedWeight(phys_RunNumber,phys_mc_channel_number,averageIntPerXing_fixed);
	return pileup_weight;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

BCHTool::BCHCleaningToolRoot* BCH;
Root::TTileTripReader* TTR;
void initBCH(bool isData)
{
	if(TTR) delete TTR;
	if(BCH) delete BCH;

	TTR = new Root::TTileTripReader("TileTripReader");	
	TTR->setTripFile("/afs/cern.ch/user/h/hod/TileTripReader/data/CompleteTripList_2011-2012.root" );
	
	// initialize the tool, giving the path to the file to emulate masked modules (used in MC only)
	BCH = new BCHTool::BCHCleaningToolRoot("BCH");
	BCH->InitializeTool(isData, TTR, "/afs/cern.ch/user/h/hod/BCHCleaningTool/share/FractionsRejectedJetsMC.root");
}

bool isBadCalibJetBCH(unsigned int jet)
{
	int randomRunNumber = pileupTool->GetRandomRunNumber(phys_RunNumber,phys_averageIntPerXing);
	int randomLumiBlockNumber = 0;
	if(randomRunNumber>0)
	{
		randomLumiBlockNumber = pileupTool->GetRandomLumiBlockNumber(randomRunNumber);
		double calibrated_calo_jet_pt  = AntiKt4LCTopoJets_calibrated_pt->at(jet);
		double calibrated_calo_jet_eta = AntiKt4LCTopoJets_calibrated_eta->at(jet);
		double calibrated_calo_jet_phi = AntiKt4LCTopoJets_calibrated_phi->at(jet);
		double BCH_CORR_CELL           = AntiKt4LCTopoJets_BCH_CORR_CELL->at(jet);
		double emfrac                  = AntiKt4LCTopoJets_emfrac->at(jet);
		int isUp                       = 0; // isUp (MC only, default=0): variation due to flavor dependence, can be either -1 (sample = 100% light quark jets), 0 (default QCD composition), 1 (sample = 100% gluon jets)
		
		 bool isBad = false;
		
		// bool isInMaskedRegion = BCH->IsInMaskedRegion(randomRunNumber, randomLumiBlockNumber, calibrated_calo_jet_eta, calibrated_calo_jet_phi);
		
		// OR: The medium cut, which comes with associated systematics in MC
		BCH->SetSeed(phys_RunNumber); // For setting the random seed used within IsBadMediumBCH(). Make sure to use this when comparing to others to make sure the same events are selected
		isBad = BCH->IsBadMediumBCH(randomRunNumber, randomLumiBlockNumber,
									calibrated_calo_jet_eta, calibrated_calo_jet_phi,
									BCH_CORR_CELL, emfrac, calibrated_calo_jet_pt, isUp);
		// OR: The tight cut
		// BCH->IsBadTightBCH(const int run, const int lbn, const double eta, const double phi, 
		//                    const double BCH_CORR_CELL, const double emfrac, const double pt);
		if(isBad) return true;
	}
	return false;
}

bool isGoodJet(unsigned int jet, TString configuration = "VeryLooseBad")
{
	bool isGood = true;
	
	// first do BCH cleaning
	if(badBCHJets.size()>0 && badBCHJets[jet]) return false;
	
	// then the rest of usual jet cleaning which can be taken from the tool if you are *not* using RootCore.
	double emf          = AntiKt4LCTopoJets_emfrac->at(jet);
	double hecf         = AntiKt4LCTopoJets_hecf->at(jet);
	double larq         = AntiKt4LCTopoJets_LArQuality->at(jet);
	double hecq         = AntiKt4LCTopoJets_HECQuality->at(jet);
	double time         = AntiKt4LCTopoJets_Timing->at(jet); //in ns
	double sumpttrk     = AntiKt4LCTopoJets_sumPtTrk_pv0_500MeV->at(jet); //in MeV, same as sumpttrk
	double eta          = AntiKt4LCTopoJets_constscale_eta->at(jet); // constscale/emscale Eta
	double pt           = (!skim && glob_doJetCalib) ? AntiKt4LCTopoJets_calibrated_pt->at(jet) : AntiKt4LCTopoJets_pt->at(jet); // should be calibrated jet pT in the analysis mode (!skim)
	double fmax         = AntiKt4LCTopoJets_fracSamplingMax->at(jet);
	double negE         = AntiKt4LCTopoJets_NegativeE->at(jet); //in MeV
	double AverageLArQF = AntiKt4LCTopoJets_AverageLArQF->at(jet);
	
	// -----------------------------------------------------------
	// Do the actual selection
	double chf=sumpttrk/pt;
	
	//=============================================================
	//VeryLoose cuts
	//=============================================================
	//Non-collision background & cosmics
	if(emf<0.05 && chf<0.05 && fabs(eta)<2)                                   isGood = false;
	if(emf<0.05 && fabs(eta)>=2)                                              isGood = false;
	if(fmax>0.99 && fabs(eta)<2)                                              isGood = false;
	//HEC spike
	if(fabs(negE/1000.)>60)                                                   isGood = false;
	if(hecf>0.5 && fabs(hecq)>0.5 && AverageLArQF/65535>0.8)                  isGood = false;
	//EM calo noise
	if(emf>0.95 && fabs(larq)>0.8 && fabs(eta)<2.8 && AverageLArQF/65535>0.8) isGood = false;
	if ("VeryLooseBad"==configuration) return isGood;
	
	//=============================================================
	//Loose cuts
	//=============================================================
	//Non-collision background & cosmics
	if(fabs(time)>25)                                  isGood = false;
	//HEC spike
	if(hecf>0.5 && fabs(hecq)>0.5)                     isGood = false;
	//EM calo noise
	if(emf>0.95 && fabs(larq)>0.8 && fabs(eta)<2.8)    isGood = false;
	if ("LooseBad"==configuration) return isGood;

	//=============================================================
	//Additionnal medium cuts
	//=============================================================
	//Non-collision background & cosmics
	if(fabs(time)>10)                                isGood = false;
	if(emf<0.05 && chf<0.1  && fabs(eta)<2)          isGood = false;
	if(emf>0.95 && chf<0.05 && fabs(eta)<2)          isGood = false;
	//HEC spike
	if(hecf>1-fabs(hecq))                            isGood = false;
	//EM calo noise
	if(emf>0.9 && fabs(larq)>0.8 && fabs(eta)<2.8)   isGood = false;
	if ("MediumBad"==configuration) return isGood;
	
	//=============================================================
	//Additionnal tight cuts
	//=============================================================
	//Non-collision background & cosmics
	if(emf<0.1 && chf<0.2 && fabs(eta)<2.5)          isGood = false;
	if(emf<0.1 && fabs(eta)>2.5 )                    isGood = false;
	if(emf>0.9 && chf<0.1 && fabs(eta)<2.5)          isGood = false;
	//EM calo noise
	if(fabs(larq)>0.95)                              isGood = false;
	if(emf>0.98 && fabs(larq)>0.05)                  isGood = false;
	if(chf<0.01 && fabs(eta)<2.5 )                   isGood = false;
	if ("TightBad"==configuration) return isGood;
	
	// We should never arrive here!
	_FATAL("Unknown configuration: "+(string)configuration+" in jet cleaning");

	return isGood;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

JetAnalysisCalib::JetCalibrationTool* JES;
void initJES(bool isdata, TString tag="Moriond") 
{
	if(JES) delete JES;
	
	if(AntiKt4LCTopoJets_calibrated_pt)  delete AntiKt4LCTopoJets_calibrated_pt;
	if(AntiKt4LCTopoJets_calibrated_E)   delete AntiKt4LCTopoJets_calibrated_E;
	if(AntiKt4LCTopoJets_calibrated_m)   delete AntiKt4LCTopoJets_calibrated_m;
	if(AntiKt4LCTopoJets_calibrated_eta) delete AntiKt4LCTopoJets_calibrated_eta;
	if(AntiKt4LCTopoJets_calibrated_phi) delete AntiKt4LCTopoJets_calibrated_phi;
	AntiKt4LCTopoJets_calibrated_pt  = new vector<float>;
	AntiKt4LCTopoJets_calibrated_E   = new vector<float>;
	AntiKt4LCTopoJets_calibrated_m   = new vector<float>;
	AntiKt4LCTopoJets_calibrated_eta = new vector<float>;
	AntiKt4LCTopoJets_calibrated_phi = new vector<float>;

	TString JES_config_file = "../ApplyJetCalibration/data/CalibrationConfigs/";
	TString jetAlgo="AntiKt4LCTopo";
	
	if(tag=="Moriond") JES_config_file += "JES_Full2012dataset_Preliminary_Jan13.config"; // Data/MC12 Configuration file
	else
	{
		if(isdata) JES_config_file += "JES_Full2012dataset_May2014.config";       // Data Configuration file
		else       JES_config_file += "JES_Full2012dataset_MC12b_May2014.config"; // MC12b Configuration file (Fullsim)
	}
	JES = new JetAnalysisCalib::JetCalibrationTool(jetAlgo,JES_config_file,isdata);
}

int findMuspContainer(float jet_eta, float jet_phi)
{
	double delR = 0;
	double delta_phi;
	double delta_eta;
	// Match mspn container to jet delR<0.4
	delR = 100;
	int index_musp=-1;
	for(unsigned int i=0 ; i<musp_phi->size(); i++)
	{
		delta_phi = fabs(jet_phi-musp_phi->at(i)); //calculate the distance in phi
		if(delta_phi>TMath::Pi()) delta_phi = (2*TMath::Pi()) - delta_phi; // always take the smaller angle (below 180°)
		delta_eta = jet_eta - musp_eta->at(i); // distance in eta
		if(sqrt(delta_phi*delta_phi + delta_eta*delta_eta) < delR)
		{
			delR = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
			index_musp = i;
		}
	}
	if(delR<0.4) return index_musp;
	return -1;
}
double getNsegments(unsigned int jet)
{
	double eta_det   = AntiKt4LCTopoJets_constscale_eta->at(jet);
	double phi       = AntiKt4LCTopoJets_constscale_phi->at(jet);
	double Nsegments = 0;  //see below
	int musp_index = findMuspContainer(eta_det, phi); //Detector jet eta and phi, obtained as above.
	if(musp_index>=0) Nsegments = musp_innerSegments->at(musp_index) + musp_outerSegments->at(musp_index) + musp_middleSegments->at(musp_index);
	return Nsegments;
}
int getNPV()
{
	int NPV=0; // count the number of vertices with 2 or more tracks
	// for(unsigned tracki=0 ; tracki<vxp_nTracks->size(); tracki++)
	for(unsigned tracki=0 ; tracki<pv_ntrk->size() ; tracki++)
	{
		// if(vxp_nTracks->at(tracki)>=2) NPV++;
		if(pv_ntrk->at(tracki)>=2) NPV++;
	}
	return NPV;
}
void clearCalibratedJets()
{
	if(!AntiKt4LCTopoJets_calibrated_pt) return;
	
	AntiKt4LCTopoJets_calibrated_pt->clear();
	AntiKt4LCTopoJets_calibrated_E->clear();
	AntiKt4LCTopoJets_calibrated_m->clear();
	AntiKt4LCTopoJets_calibrated_eta->clear();
	AntiKt4LCTopoJets_calibrated_phi->clear();
}
void addCalibratedJet(TLorentzVector Jet)
{
	if(!AntiKt4LCTopoJets_calibrated_pt) return;
	
	AntiKt4LCTopoJets_calibrated_pt->push_back(Jet.Pt());
	AntiKt4LCTopoJets_calibrated_E->push_back(Jet.E());
	AntiKt4LCTopoJets_calibrated_m->push_back(Jet.M());
	AntiKt4LCTopoJets_calibrated_eta->push_back(Jet.Eta());
	AntiKt4LCTopoJets_calibrated_phi->push_back(Jet.Phi());
}

TLorentzVector getJES(unsigned int jet, TString tag="Moriond")
{
	TLorentzVector Jet;
	
	if(tag=="Moriond")
	{
		_DEBUG("");
		
		double Eraw = AntiKt4LCTopoJets_constscale_E->at(jet);
		double eta  = AntiKt4LCTopoJets_constscale_eta->at(jet);
		double phi  = AntiKt4LCTopoJets_constscale_phi->at(jet);
		double m    = AntiKt4LCTopoJets_constscale_m->at(jet);
		double Ax   = AntiKt4LCTopoJets_ActiveAreaPx->at(jet);
		double Ay   = AntiKt4LCTopoJets_ActiveAreaPy->at(jet);
		double Az   = AntiKt4LCTopoJets_ActiveAreaPz->at(jet);
		double Ae   = AntiKt4LCTopoJets_ActiveAreaE->at(jet);
		double rho  = rhorhoKt4LC; // or rhorhoKt4EM;
		
		// For the pile-up correction, we need mu and NPV(2+ tracks)
		double mu = phys_averageIntPerXing;
		int NPV = getNPV(); // count the number of vertices with 2 or more tracks
		
		Jet = JES->ApplyJetAreaOffsetEtaJES(Eraw,eta,phi,m,Ax,Ay,Az,Ae,rho,mu,NPV);
		
		_DEBUG("");
	}
	else
	{
		_DEBUG("");
		
		double Eraw       = AntiKt4LCTopoJets_constscale_E->at(jet);
		double eta_det    = AntiKt4LCTopoJets_constscale_eta->at(jet);
		double phi        = AntiKt4LCTopoJets_constscale_phi->at(jet);
		double m          = AntiKt4LCTopoJets_constscale_m->at(jet);
		double eta_origin = AntiKt4LCTopoJets_EtaOrigin->at(jet);
		double phi_origin = AntiKt4LCTopoJets_PhiOrigin->at(jet);
		double m_origin   = AntiKt4LCTopoJets_MOrigin->at(jet);
		double Ax         = AntiKt4LCTopoJets_ActiveAreaPx->at(jet);
		double Ay         = AntiKt4LCTopoJets_ActiveAreaPy->at(jet);
		double Az         = AntiKt4LCTopoJets_ActiveAreaPz->at(jet);
		double Ae         = AntiKt4LCTopoJets_ActiveAreaE->at(jet);
		double rho        = rhorhoKt4LC;       
		double fEM3       = (AntiKt4LCTopoJets_e_EMB3->at(jet)+AntiKt4LCTopoJets_e_EME3->at(jet))/AntiKt4LCTopoJets_constscale_E->at(jet);
		double fTile0     = (AntiKt4LCTopoJets_e_TileBar0->at(jet)+AntiKt4LCTopoJets_e_TileExt0->at(jet))/AntiKt4LCTopoJets_constscale_E->at(jet);
		double nTrk       = AntiKt4LCTopoJets_nTrk_pv0_1GeV->at(jet); 
		double trackWIDTH = AntiKt4LCTopoJets_trackWIDTH_pv0_1GeV->at(jet); //This variable may also be called AntiKt4LCTopoJets_trackWIDTH_pv0_1GeV
		
		_DEBUG("");
		
		double Nsegments = getNsegments(jet);
		
		_DEBUG("");		
		// For the pile-up correction, we need mu and NPV(2+ tracks)
		double mu = phys_averageIntPerXing;
		int NPV = getNPV(); // count the number of vertices with 2 or more tracks
		
		_DEBUG("");
		
		Jet = JES->ApplyJetAreaOffsetOriginEtaJESGSC(Eraw,eta_det,phi,m,eta_origin,phi_origin,m_origin,Ax,Ay,Az,Ae,rho,trackWIDTH,nTrk,fTile0,fEM3,Nsegments,mu,NPV);
		
		_DEBUG("");
	}
	
	_DEBUG("");
	
	return Jet;
}

MultijetJESUncertaintyProvider* JUNCP;
void initJUN(TString tag="Moriond")
{
	if(JUNCP) delete JUNCP;
	
	if(tag=="Moriond")
	{
		TString multijetConfig = "JES_2012/Moriond2013/MultijetJES_2012.config";
		TString jesprovConfig  = "JES_2012/Moriond2013/InsituJES2012_14NP.config"; // Two other nominal options
		TString jetAlgorithm   = "AntiKt4LCTopo"; // One of {AntiKt4TopoEM,AntiKt4LCTopo,AntiKt6TopoEM,AntiKt6LCTopo}
		TString mcType         = "MC12a"; // Also accepts {Pythia8,MC12c,AFII} (note: MC12a==Pythia8)
		TString path           = "../JetUncertainties/share/"; // Path to the share directory
		JUNCP = new MultijetJESUncertaintyProvider(multijetConfig,jesprovConfig,jetAlgorithm,mcType,path);
	}
	else
	{
		TString base_path_for_config_files = "../JetUncertainties/share/";
		TString flavour_and_topology_uncertainties_config_file = "JES_2012/Final/MultijetJES_2012.config";
		TString baseline_uncertainties_uncertainties_config_file = "JES_2012/Final/InsituJES2012_AllNuisanceParameters.config";
		TString collection_name = "AntiKt4LCTopo"; //  AntiKt4TopoEM, AntiKt6TopoEM, AntiKt4TopoLC and AntiKt6TopoLC.
		TString MC_youre_running_on = "MC12a"; // "Pythia8";
		JUNCP = new MultijetJESUncertaintyProvider(flavour_and_topology_uncertainties_config_file,
											  	   baseline_uncertainties_uncertainties_config_file,
											  	   collection_name,
											  	   MC_youre_running_on,
											  	   base_path_for_config_files);
	}
}


double getJUNC(unsigned int jet, TLorentzVector Jet, TMapVL& JetShiftsUp, TMapVL& JetShiftsDwn, TString tag="Moriond")
{
	double total_shift = 0;
	JetShiftsUp.clear();
	JetShiftsDwn.clear();
	TLorentzVector JetShiftUp, JetShiftDwn;
	
	/////////////////////////
	// Baseline uncertainties
	int Ncomp = JUNCP->getNUncertaintyComponents();
	for(int icomp=0 ; icomp<Ncomp ; ++icomp)
	{
		TString compName = JUNCP->getComponentNames().at(icomp); // returns the name of the nuisance parameter
		TString compDesc = JUNCP->getComponentDescriptions().at(icomp); // returns a one-sentence description of the nuisance parameter
		compCategory categoryEnum = JUNCP->getComponentCategories().at(icomp); // returns a enumeration of the category for the nuisance parameter
		TString categoryName = JUNCP->getCategoryStringFromEnum(categoryEnum); // returns a string for the category name from the enumeration
		double shift = JUNCP->getRelUncertComponent(icomp, Jet.Pt(), Jet.Eta()); // nuisance paramter amplitude (with sign) ("relative uncertainty" of component) 
		// can also use  jUNCP->getRelUncertComponent(compName, Jet.Pt(), Jet.Eta());
		
		JetShiftUp  = Jet;
		JetShiftDwn = Jet;
		
		// now scale your jet (the full 4-vector)
		// Upward variation: Jet *= 1+shift
		// Downward variation: Jet *= 1-shift
		JetShiftUp *= (1.0+shift);
		JetShiftDwn *= (1.0-shift);
		
		total_shift += shift*shift;
		
		JetShiftsUp.insert(make_pair(compName,JetShiftUp));
		JetShiftsDwn.insert(make_pair(compName,JetShiftDwn));
		
		// print the uncertainty and information
		// printf("JES uncertainty source %2d: %s\n",icomp,compName.Data());  
		// printf("Description : %s\n",compDesc.Data());
		// printf("Category: %s\n",categoryName.Data());
		// printf("The uncertainty for (pT,eta) = (%.1f,%5.1f) is %.1f%%\n",Jet.Pt(),Jet.Eta(),unc*100);
	}
	
	///////////////////////
	// Pileup uncertainties
	double mu = phys_averageIntPerXing;
	int NPV = getNPV(); // count the number of vertices with 2 or more tracks

	double shift = 0;
	
	shift = JUNCP->getRelNPVOffsetTerm(Jet.Pt(),Jet.Eta(),NPV); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("NPVOffsetTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("NPVOffsetTerm",JetShiftDwn));
	
	shift = JUNCP->getRelMuOffsetTerm(Jet.Pt(),Jet.Eta(),mu); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("MuOffsetTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("MuOffsetTerm",JetShiftDwn));
	
	shift = JUNCP->getRelPileupPtTerm(Jet.Pt(),Jet.Eta(),NPV,mu); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("PileupPtTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("PileupPtTerm",JetShiftDwn));
	
	shift = JUNCP->getRelPileupRhoTopology(Jet.Pt(),Jet.Eta()); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("PileupRhoTopology", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("PileupRhoTopology",JetShiftDwn));
	
	
	///////////////////////////////////////
	//// b-tagging / flavor uncertainty
	shift = JUNCP->getRelFlavorCompUncert(Jet.Pt(),Jet.Eta(),true);  total_shift += shift*shift; // ONLY FOR LIGHT JETS
	shift = JUNCP->getRelFlavorResponseUncert(Jet.Pt(),Jet.Eta());   total_shift += shift*shift; // ONLY FOR LIGHT JETS
	shift = JUNCP->getRelBJESUncert(Jet.Pt(),Jet.Eta());             total_shift += shift*shift; // ONLY FOR b-JETS
	
	
	if(tag=="Moriond")
	{
		// ????
	}
	else
	{
		double Nsegments = getNsegments(jet);
		shift = JUNCP->getRelPunchThroughUncert(Jet.Pt(),Jet.Eta(),Nsegments);  total_shift += shift*shift; // Punch throgh
	}
	
	return sqrt(total_shift);
}


JetSmearingTool* JER;
void initJER()
{
	if(JER) delete JER;

	TString JERinputFile = "/afs/cern.ch/user/h/hod/JetResolution/share/JERProviderPlots_2012.root";
	TString jetAlgo = "AntiKt4LCTopo";
	JER = new JetSmearingTool(jetAlgo,JERinputFile);
	JER->init();
}
void seedJER()
{
	JER->SetSeed(phys_RunNumber);
}
TLorentzVector getJER(TLorentzVector Jet, double& shift)
{
	TLorentzVector smearedJet = Jet;
	// The below is systematic evaluation, and ONLY for MC
	// Smear the jet to match the MC resolution+1 sigma!
	// JER->SmearJet_Syst(smearedJet); // FOR 7 TeV 2011 and 8 TeV 2012 MC12 Fullsim
	double factor = JER->GetRandomSmearingFactorSyst(Jet);
	smearedJet *= factor;
	shift = fabs(factor-1.);
	return smearedJet;
}


METUtility *METU;
void initMET()
{
	if(METU) delete METU;
	if(METUfortrk) delete METUfortrk;

	METU = new METUtility;
	METUfortrk = new METUtility;
	
	METU->configMissingET(true,false);
	METUfortrk->configMissingET(true,false);
	
	METUfortrk->configMETSyst("METTrack_2012.config");
	
	// METU->setVerbosity(isVerbose);
	// METU->setSoftJetCut(20); // soft jet cut
	
	// // set pre-defined regions as above
	// m_util->setCaloRegion(region);
	// // region = METUtil::Central, METUtil::EndCap or METUtil::Forward
	// // or set eta cuts explicitly
	// m_util->setObjectEtaCut(etaLow, etaHigh);

	if(akt4lc_jet_pt)  akt4lc_jet_pt->clear();  else akt4lc_jet_pt  = new vector<float>;
	if(akt4lc_jet_E)   akt4lc_jet_E->clear();   else akt4lc_jet_E   = new vector<float>;
	if(akt4lc_jet_eta) akt4lc_jet_eta->clear(); else akt4lc_jet_eta = new vector<float>;
	if(akt4lc_jet_phi) akt4lc_jet_phi->clear(); else akt4lc_jet_phi = new vector<float>;
}
void setJetVectorPointers(vector<TLorentzVector>& nominalJets, vector<double>& shifts, int shift_code)
{
	akt4lc_jet_pt->clear();
	akt4lc_jet_E->clear();
	akt4lc_jet_eta->clear();
	akt4lc_jet_phi->clear();
	
	unsigned int nJets = nominalJets.size();
	for(unsigned int jet=0 ; jet<nJets ; ++jet)
	{
		TLorentzVector Jet = nominalJets[jet];
		double shift = shifts[jet];
		if     (shift_code==SHIFTUP)  Jet *= (1.+shift);
		else if(shift_code==SHIFTDWN) Jet *= (1.-shift);
		else if(shift_code==NOSHIFT)  Jet *= 1.;
		else _FATAL("unsupported shift code = "+_s(shift_code));
	
		// cout << "jet=" << jet << ", shift_code=" << shift_code << ", shift=" << shift << ", pT:" << nominalJets[jet].Pt() << " --> " << Jet.Pt() << endl;
		
		akt4lc_jet_pt->push_back(Jet.Pt());
		akt4lc_jet_E->push_back(Jet.E());
		akt4lc_jet_eta->push_back(Jet.Eta());
		akt4lc_jet_phi->push_back(Jet.Phi());
	}
}
void setJetVectorPointers(vector<TLorentzVector>& Jets)
{
	akt4lc_jet_pt->clear();
	akt4lc_jet_E->clear();
	akt4lc_jet_eta->clear();
	akt4lc_jet_phi->clear();
	
	unsigned int nJets = Jets.size();
	for(unsigned int jet=0 ; jet<nJets ; ++jet)
	{
		TLorentzVector Jet = Jets[jet];
		
		akt4lc_jet_pt->push_back(Jet.Pt());
		akt4lc_jet_E->push_back(Jet.E());
		akt4lc_jet_eta->push_back(Jet.Eta());
		akt4lc_jet_phi->push_back(Jet.Phi());
	}
}
METUtil::METObject getMETU(int mettype)
{
	METU->reset();
	
	METU->defineMissingET(true,true,true,true,true,true,true);
	// METU->setIsMuid(true);
	if(mettype==METMUONS) METU->setIsMuons(true);
	
	//// soft muons and other leptons
	METU->setMETTerm(METUtil::RefMuon,   MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet); // must be done in addition to the staco/Muons code below
	METU->setMETTerm(METUtil::RefTau,    MET_RefTau_etx,  MET_RefTau_ety,  MET_RefTau_sumet);
	// METU->setElectronParameters(el_pt, el_eta, el_phi, el_MET_RefFinal_comp_wet, el_MET_RefFinal_comp_wpx, el_MET_RefFinal_comp_wpy, el_MET_RefFinal_comp_statusWord);
	METU->setMETTerm(METUtil::RefEle, MET_RefEle_etx, MET_RefEle_ety, MET_RefEle_sumet);
	// METU->setPhotonParameters(ph_pt, ph_eta, ph_phi, ph_MET_RefFinal_comp_wet, ph_MET_RefFinal_comp_wpx, ph_MET_RefFinal_comp_wpy, ph_MET_RefFinal_comp_statusWord);
	METU->setMETTerm(METUtil::RefGamma, MET_RefGamma_etx, MET_RefGamma_ety, MET_RefGamma_sumet);
	
	//// soft terms
	METU->setMETTerm(METUtil::SoftTerms, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet);
	// METU->setMETTerm(METUtil::SoftTerms, MET_CellOut_Eflow_etx+MET_SoftJets_etx, MET_CellOut_Eflow_ety+MET_SoftJets_ety, MET_CellOut_Eflow_sumet+MET_SoftJets_sumet);
	
	//// Jets
	METU->setJetParameters(akt4lc_jet_pt, akt4lc_jet_eta, akt4lc_jet_phi,akt4lc_jet_E, jet_antikt4LCtopo_MET_RefFinal_comp_wet, jet_antikt4LCtopo_MET_RefFinal_comp_wpx, jet_antikt4LCtopo_MET_RefFinal_comp_wpy, jet_antikt4LCtopo_MET_RefFinal_comp_statusWord);
	// 	METU->setMETTerm(METUtil::RefJet, MET_RefJet_etx, MET_RefJet_ety, MET_RefJet_sumet);
	

	//// "Hard" muons
	if(mettype==METSTACO)
	{
		// METU->setMuonParameters(mu_staco_pt, mu_staco_eta, mu_staco_phi, mu_staco_MET_RefFinal_comp_wet, mu_staco_MET_RefFinal_comp_wpx, mu_staco_MET_RefFinal_comp_wpy, mu_staco_MET_RefFinal_comp_statusWord);
		// The default pt for muons is from the combined ID/MS track.
		// // The default pt for muons is from the combined ID/MS track.
		// // Spectro-only muons need to have the MS momentum set separately.
		// // This version of the method uses commonly available D3PD branches.
		// METU->setExtraMuonParameters(mu_staco_ms_qoverp, mu_staco_ms_theta, mu_staco_ms_phi, mu_staco_charge);
		METU->setMETTerm(METUtil::MuonTotal, MET_MuonBoy_etx, MET_MuonBoy_ety, MET_MuonBoy_sumet);
	}
	else if(mettype==METMUONS)
	{
		// METU->setMuonParameters(mu_muons_pt, mu_muons_eta, mu_muons_phi, mu_muons_MET_RefFinal_comp_wet, mu_MET_RefFinal_comp_wpx, mu_MET_RefFinal_comp_wpy, mu_MET_RefFinal_comp_statusWord);
		// The default pt for muons is from the combined ID/MS track.
		// // The default pt for muons is from the combined ID/MS track.
		// // Spectro-only muons need to have the MS momentum set separately.
		// // This version of the method uses commonly available D3PD branches.
		// METU->setExtraMuonParameters(mu_muons_ms_qoverp, mu_muons_ms_theta, mu_muons_ms_phi, mu_muons_charge);
		METU->setMETTerm(METUtil::MuonTotal, MET_Muons_etx, MET_Muons_ety, MET_Muons_sumet);
	}
	else _FATAL("Unsupported MET type: "+_s(mettype));


	METU->setAverageIntPerXing(phys_averageIntPerXing);
	if(glob_isMC) METU->setMETTerm(METUtil::Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet);

	
	
	//// Now get the MET (RefFinal)
	METUtil::METObject refFinal = METU->getMissingET(METUtil::RefFinal);
	// METUtil::METObject refFinal_test(MET_RefFinal_etx, MET_RefFinal_ety, MET_RefFinal_sumet);
	// bool check_refFinal = METU->checkConsistency(refFinal_test,METUtil::RefFinal);
	// if(check_refFinal) cout << "RefFinal checks out!" << endl;
	// else               cout << "RefFinal doesn't check out!" << endl;
	
	// if(MET_RefEle_et>0 || MET_RefGamma_et>0)
	// {
	// 	cout << "MET_RefJet_et="   << MET_RefJet_et   << ", refJet_et="   << METU->getMissingET(METUtil::RefJet).et() << endl;
	// 	cout << "MET_RefEle_et="   << MET_RefEle_et   << ", refEle_et="   << METU->getMissingET(METUtil::RefEle).et() << endl;
	// 	cout << "MET_RefGamma_et=" << MET_RefGamma_et << ", refGamma_et=" << METU->getMissingET(METUtil::RefGamma).et() << endl;
	// }
	
 	return refFinal;
}


METTrack::TrackMETMaker *METTRK;
void initMETTRK()
{
	if(METTRK) delete METTRK;
	
	METTRK = new METTrack::TrackMETMaker;
	
	METTRK->init("TrkMetNom.config", "/afs/cern.ch/user/h/hod/METTrackUtil/share/");
}

void setTracksFloat()
{
	
	/// THE TRACK BRANCHES SHOULD HAVE DIFFERENT NEAMES IN THE NTUP_COMMON THAN WHAT YOU SEE HERE !!!
	ftrk_pt->clear();
	ftrk_eta->clear();
	ftrk_phi0->clear();
	ftrk_d0->clear();
	ftrk_z0->clear();
	ftrk_qoverp->clear();
	ftrk_qoverpCov->clear();
	
	for(unsigned int i=0 ; i<trks_pt->size() ; ++i)
	{
		ftrk_pt->push_back(trks_pt->at(i));
		ftrk_eta->push_back(trks_eta->at(i));
		ftrk_phi0->push_back(trks_phi0->at(i));
		ftrk_d0->push_back(trks_extrapD0->at(i) /*trks_d0->at(i)*/);
		ftrk_z0->push_back(trks_extrapZ0->at(i) /*trks_z0->at(i)*/);
		ftrk_qoverp->push_back(trks_qoverp->at(i));
		ftrk_qoverpCov->push_back(trks_qoverpErr->at(i)*trks_qoverpErr->at(i));
	}
}

void setMETTRK()
{
	METTRK->reset();
	METTRK->fillTracks(ftrk_pt,ftrk_eta,ftrk_phi0,ftrk_d0,ftrk_z0,trks_nPix,trks_nSCT/*trkpt4_pt,trkpt4_phi*/); // default to 0 pointer
	METTRK->filterTracks(ftrk_qoverpCov,ftrk_qoverp,cl_lc_pt,cl_lc_eta,cl_lc_phi);
	METTRK->setMuons(mu_muons_pt,mu_muons_eta,mu_muons_phi,mu_muons_isCombinedMuon,mu_muons_id_qoverp_exPV,mu_muons_id_theta_exPV,mu_muons_id_phi_exPV);
	METTRK->setElectrons(el_cl_E,el_eta,el_mediumPP,el_author,el_tracketa,el_trackphi,el_Unrefittedtrack_pt,el_Unrefittedtrack_eta,el_Unrefittedtrack_phi/*el_trk_index*/); // optional -- if left 0 will match via dR and relative pt
	METTRK->setJets(akt4lc_jet_pt,akt4lc_jet_eta,akt4lc_jet_phi,akt4lc_jet_E,AntiKt4LCTopoJets_jvtxf);
					// 0, /*jet_AntiKt4LCTopo_trackAssoc_index,*/ // optional -- if left 0, will select using dR
					// AntiKt4LCTopoJets_isBadLooseMinus, // optional -- if left 0, will not select on isBad 
					// AntiKt4LCTopoJets_isUgly // optional -- if left 0, will not select on isUgly 
					// );
	METTRK->assignTracks();
}
METUtil::METObject getMETTRK()
{
	return METTRK->getTrackMET();
}

// METSyst::METSystTool* METTRKU;
// void initMETTRKU()
// {
// 	METTRKU = new METSyst::METSystTool();
// 	// METTRKU->initialise("/afs/cern.ch/user/h/hod/METSystematics/data/METTrack_2012.config");	
// 	METTRKU->initialise("METTrack_2012.config");	
// }
METUtil::METObject METTRKU;
void setMETTRKU()
{
	METUfortrk->reset();

	METUfortrk->setMETTerm(METUtil::RefJet,METTRK->getTrackMET(METTrack::TrkJet).etx(),METTRK->getTrackMET(METTrack::TrkJet).ety(),METTRK->getTrackMET(METTrack::TrkJet).sumet());
	METUfortrk->setMETTerm(METUtil::RefEle,METTRK->getTrackMET(METTrack::TrkEle).etx(),METTRK->getTrackMET(METTrack::TrkEle).ety(),METTRK->getTrackMET(METTrack::TrkEle).sumet());
	METUfortrk->setMETTerm(METUtil::MuonTotal,METTRK->getTrackMET(METTrack::TrkMu).etx(),METTRK->getTrackMET(METTrack::TrkMu).ety(),METTRK->getTrackMET(METTrack::TrkMu).sumet());
	METUfortrk->setMETTerm(METUtil::SoftTerms,METTRK->getTrackMET(METTrack::TrkSoft).etx(),METTRK->getTrackMET(METTrack::TrkSoft).ety(),METTRK->getTrackMET(METTrack::TrkSoft).sumet());
	
	vector<const METTrack::Track*> jettracks = METTRK->getTermTracks(METTrack::TrkJet);
	METUfortrk->setJetTracks(akt4lc_jet_pt,akt4lc_jet_eta,jettracks);
	
	METUfortrk->setAverageIntPerXing(phys_averageIntPerXing);
	if(glob_isMC) METUfortrk->setMETTerm(METUtil::Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet);
	
	METTRKU = METUfortrk->getMissingET(METUtil::RefFinal);
	
}
METUtil::METObject getMETTRKU(int syst)
{
	_DEBUG("");

	METUtil::METObject mettrk_syst;
	METUtil::METObject jettrkUP = METUfortrk->getJetTrackMET(METSyst::JetTrk_Up);
	METUtil::METObject jettrkDN = METUfortrk->getJetTrackMET(METSyst::JetTrk_Down);
	
	if     (syst==-1)             mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal);
	else if(syst==SOFTTRKUP)      mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal, METUtil::SoftTrackScaleUp);
	else if(syst==SOFTTRKDWN)     mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal, METUtil::SoftTrackScaleDown);
	else if(syst==SOFTTRKRESPARA) mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal, METUtil::SoftTrackResoPara);
	else if(syst==SOFTTRKRESPERP) mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal, METUtil::SoftTrackResoPerp);
	else if(syst==JETTRKUP)
	{
	        METUfortrk->setMETTerm(METUtil::RefJet,jettrkUP.etx(),jettrkUP.ety(),jettrkUP.sumet());
	        mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal);
	}
	else if(syst==JETTRKDWN)
	{
	        METUfortrk->setMETTerm(METUtil::RefJet,jettrkDN.etx(),jettrkDN.ety(),jettrkDN.sumet());
	        mettrk_syst = METUfortrk->getMissingET(METUtil::RefFinal);
	}
	else _FATAL("Unknown enum for systematic in TrackMET");
	return mettrk_syst;

}

void analysis()
{
	initializePileup();        // PU tool is needed for BCH cleaning !
	initBCH(isdata);           // BCH cleaning should be used either for data or signal
	initJES(isdata,JESconfig); // Jet Energy Scale
	initJER();                 // Jet Energy Resolution
	initJUN(JESconfig);        // Jet Energy Scale uncertainty
	initMET();                 // MET Utility
	initMETTRK();              // MET Track tool
	// initMETTRKU();             // MET Systematics (for Track MET only)
	
	for(/* loop on the tree entires */)
	{
		seedJER();
		calibJets.clear();
		smearedJets.clear();
		badBCHJets.clear();
		calibJetsUnc.clear();
		smearedJetsUnc.clear();
		calibJetsIndex.clear();
		clearCalibratedJets();
		
		vector<TLorentzVector> jets_uncalibrated;
		vector<TLorentzVector> jets_jes_nominal;
		vector<TLorentzVector> jets_jes_up;
		vector<TLorentzVector> jets_jes_dwn;
		vector<TLorentzVector> jets_jer_up;
		vector<TLorentzVector> jets_jer_dwn;
		
		unsigned int nJets = AntiKt4LCTopoJets_pt->size();
		TMapVL JetShiftsUp, JetShiftsDwn;
		
		// cout << "\n\n--------- New Event, nJets=" << nJets << " ---------" << endl;
		for(unsigned int jet=0 ; jet<nJets ; jet++)
		{
			TLorentzVector ucJet;
			ucJet.SetPtEtaPhiE(AntiKt4LCTopoJets_pt->at(jet),AntiKt4LCTopoJets_eta->at(jet),AntiKt4LCTopoJets_phi->at(jet),AntiKt4LCTopoJets_E->at(jet));
			jets_uncalibrated.push_back(ucJet);
		
			TLorentzVector Jet = getJES(jet,JESconfig);
			calibJets.push_back(Jet);
			addCalibratedJet(Jet);
			jets_jes_nominal.push_back(Jet);
			
			bool badBCHjet = isBadCalibJetBCH(jet);
			badBCHJets.push_back(badBCHjet);
			
			double quadJES  = getJUNC(jet,Jet,JetShiftsUp,JetShiftsDwn,JESconfig);
			calibJetsUnc.push_back(quadJES);
			calibJetsIndex.push_back(jet);
			jets_jes_up.push_back(Jet*(1+quadJES));
			jets_jes_dwn.push_back(Jet*(1-quadJES));
			
			double quadJER = 0;
			TLorentzVector smearedJet = getJER(Jet,quadJER);
			smearedJets.push_back(smearedJet);
			smearedJetsUnc.push_back(quadJER);
			jets_jer_up.push_back(Jet*(1+quadJER));
			jets_jer_dwn.push_back(Jet*(1-quadJER));
		}
		
		
		// double METx = MET_RefMuon_etx+MET_RefTau_etx+MET_CellOut_Eflow_etx+MET_RefEle_etx+MET_RefGamma_etx+MET_RefJet_etx+MET_MuonBoy_etx;
		// double METy = MET_RefMuon_ety+MET_RefTau_ety+MET_CellOut_Eflow_ety+MET_RefEle_ety+MET_RefGamma_ety+MET_RefJet_ety+MET_MuonBoy_ety;
		// cout << "X: RefMuon=" << MET_RefMuon_etx << ", RefTau=" << MET_RefTau_etx << ", CellOut_Eflow=" << MET_CellOut_Eflow_etx << ", RefEle=" << MET_RefEle_etx << ", RefGamma=" << MET_RefGamma_etx << ", RefJet=" << MET_RefJet_etx << ", MuonBoy=" << MET_MuonBoy_etx << " --> Sum=" << METx << "  :  RefFinal=" << MET_RefFinal_etx << endl;
		// cout << "Y: RefMuon=" << MET_RefMuon_ety << ", RefTau=" << MET_RefTau_ety << ", CellOut_Eflow=" << MET_CellOut_Eflow_ety << ", RefEle=" << MET_RefEle_ety << ", RefGamma=" << MET_RefGamma_ety << ", RefJet=" << MET_RefJet_ety << ", MuonBoy=" << MET_MuonBoy_ety << " --> Sum=" << METy << "  :  RefFinal=" << MET_RefFinal_ety << endl;
		
		// cout << "MET_RefFinal_et=" << MET_RefFinal_et << ", met_reffinal_et=" << met_RefFinal_et << endl;
		setJetVectorPointers(jets_uncalibrated); uncalibMET       = getMETU(METSTACO); uncalibMUMET       = getMETU(METMUONS); setMETTRK(); uncalibMETTRK       = getMETTRK();
		setJetVectorPointers(jets_jes_nominal);  calibMET_nominal = getMETU(METSTACO); calibMUMET_nominal = getMETU(METMUONS); setMETTRK(); calibMETTRK_nominal = getMETTRK();
		setJetVectorPointers(jets_jes_up);       calibMET_jes_up  = getMETU(METSTACO); calibMUMET_jes_up  = getMETU(METMUONS); setMETTRK();
		setJetVectorPointers(jets_jes_dwn);      calibMET_jes_dwn = getMETU(METSTACO); calibMUMET_jes_dwn = getMETU(METMUONS); setMETTRK();
		setJetVectorPointers(jets_jer_up);       calibMET_jer_up  = getMETU(METSTACO); calibMUMET_jer_up  = getMETU(METMUONS); setMETTRK();
		setJetVectorPointers(jets_jer_dwn);      calibMET_jer_dwn = getMETU(METSTACO); calibMUMET_jer_dwn = getMETU(METMUONS); setMETTRK();
		
		setJetVectorPointers(jets_jes_nominal); setMETTRK(); setMETTRKU(); // setup Track MET systematics with the nominally calibrated jets !!!
		
		calibMETTRK_softtrk_up      = getMETTRKU(SOFTTRKUP);
		calibMETTRK_softtrk_dwn     = getMETTRKU(SOFTTRKDWN);
		calibMETTRK_softtrkres_para = getMETTRKU(SOFTTRKRESPARA);
		calibMETTRK_softtrkres_perp = getMETTRKU(SOFTTRKRESPERP);
		calibMETTRK_jettrk_up       = getMETTRKU(JETTRKUP);
		calibMETTRK_jettrk_dwn      = getMETTRKU(JETTRKDWN);
	}
}