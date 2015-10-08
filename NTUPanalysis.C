//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// root -b -l -q  compile.C  --chnl=all  --master=muons  --method=cuts  --jetmode=calib  --trkmode=calib --metmode=calib  --mettype=muons --blinded=yes  [--ndata=1000] ////
//// see also compile.C //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "std.h"
#include "type.h"
// #include "enums.h"
#include "branch.h"
#include "select.h"
#include "tmvaRead.h"
// #include "count.h"
#include "draw.h"
#include "mvaout.h"
#include "cutsout.h"
// #include "clonetree.h"
#include "runinfo.h"
#include "postBDTcuts.h"

static const int vis = 0;
ofstream* ofstr  = new ofstream("events.txt");
ofstream* ofstr1 = new ofstream("eventrecord.txt");
bool isAntiKt4LCJet = true;
int nData = -1;

void correctLumiWeight(TString name)
{
	if(name.Contains("Data")) return; // don't touch the data
	if(name.Contains("NpX") || name.Contains("JZxW")) // ignore the binned dijet / EW samples
	{
		// for the dijet it will not be correct
		// because it was not assumed to be normalized to 1fb-1
		// while making the bulk of the MCs with the data in lxplus
		floats["wgt_luminosity"] *= luminosity;
		return;
	}
	
	float sigma,Nevents,genEff;
	float sigmaSM=-1.;
	if     (name.Contains("Wtaunu_3mu"))      { sigma=9.5753E+00*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=9.1081E+00*nb2fb; }
	else if(name.Contains("Wtaunu_200k_3mu")) { sigma=9.5753E+00*nb2fb; Nevents=200000.;   genEff=1.0000E+00; sigmaSM=9.1081E+00*nb2fb; }
	else if(name.Contains("bbTotau10_3mu"))   { sigma=3.5388E+02*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=3.5388E+02*nb2fb; }
	else if(name.Contains("ccTotau10_3mu"))   { sigma=2.9674E+02*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=2.9674E+02*nb2fb; }
	else if(name.Contains("bb_mu4mu4"))       { sigma=1.1464E+02*nb2fb; Nevents=19980978.; genEff=1.0000E+00; }
	else if(name.Contains("bbTomu15"))        { sigma=1.9898E+02*nb2fb; Nevents=4998088.;  genEff=1.0000E+00; }
	else if(name.Contains("ccTomu15"))        { sigma=8.0088E+01*nb2fb; Nevents=4998690.;  genEff=1.0000E+00; }
	else if(name.Contains("bb_Jpsimu4mu4"))   { sigma=2.0874E+02*nb2fb; Nevents=9969994.;  genEff=1.0000E+00; } // cross section is wrong !
                                              
	else if(name.Contains("JZ0W"))            { sigma=7.2850E+07*nb2fb; Nevents=3998693.;  genEff=4.2413E-04;}
	else if(name.Contains("JZ1W"))            { sigma=4.1440E+06*nb2fb; Nevents=1999694.;  genEff=3.4217E-05;}
	else if(name.Contains("JZ2W"))            { sigma=5.0147E+03*nb2fb; Nevents=2499692.;  genEff=7.0769E-04; }
	else if(name.Contains("JZ3W"))            { sigma=5.4418E+02*nb2fb; Nevents=998688.;   genEff=6.7835E-05; }
	                                          
	else if(name.Contains("WmunuNp0"))        { sigma=8.03135E+00*nb2fb; Nevents=3469591.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp1"))        { sigma=1.58358E+00*nb2fb; Nevents=2499893.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp2"))        { sigma=4.80668E-01*nb2fb; Nevents=3769890.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp3"))        { sigma=1.33657E-01*nb2fb; Nevents=1009896.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp4"))        { sigma=3.57460E-02*nb2fb; Nevents=255000.;   genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp5"))        { sigma=1.05518E-02*nb2fb; Nevents=20000.;    genEff=1.0000E+00;}
	                                          
	else if(name.Contains("ZmumuNp0"))        { sigma=7.09421E-01*nb2fb; Nevents=6619489.;  genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp1"))        { sigma=1.53982E-01*nb2fb; Nevents=1334706.;  genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp2"))        { sigma=4.89215E-02*nb2fb; Nevents=404997.;   genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp3"))        { sigma=1.41532E-02*nb2fb; Nevents=110000.;   genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp4"))        { sigma=3.80883E-03*nb2fb; Nevents=30000.;    genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp5"))        { sigma=1.10938E-03*nb2fb; Nevents=10000.;    genEff=1.0000E+00;}
	
	else if(name.Contains("WtaunuNp0"))       { sigma=8.0362E+00*nb2fb; Nevents=3364789.;  genEff=1.0000E+00;}
	else if(name.Contains("WtaunuNp1"))       { sigma=1.5787E+00*nb2fb; Nevents=2449991.;  genEff=1.0000E+00;}
	else if(name.Contains("WtaunuNp2"))       { sigma=4.7780E-01*nb2fb; Nevents=3719888.;  genEff=1.0000E+00;}
	else if(name.Contains("WtaunuNp3"))       { sigma=1.3401E-01*nb2fb; Nevents=1009993.;  genEff=1.0000E+00;}
	else if(name.Contains("WtaunuNp4"))       { sigma=3.5258E-02*nb2fb; Nevents=249898.;   genEff=1.0000E+00;}
	else if(name.Contains("WtaunuNp5"))       { sigma=1.0646E-02*nb2fb; Nevents=65000.;    genEff=1.0000E+00;}
	                                          
	else if(name.Contains("ZtautauNp0"))      { sigma=7.102515E-01*nb2fb; Nevents=6619683.;  genEff=1.0000E+00;}
	else if(name.Contains("ZtautauNp1"))      { sigma=1.558205E-01*nb2fb; Nevents=1334996.;  genEff=1.0000E+00;}
	else if(name.Contains("ZtautauNp2"))      { sigma=4.890875E-02*nb2fb; Nevents=404997.;   genEff=1.0000E+00;}
	else if(name.Contains("ZtautauNp3"))      { sigma=1.409320E-02*nb2fb; Nevents=110000.;   genEff=1.0000E+00;}
	else if(name.Contains("ZtautauNp4"))      { sigma=3.769315E-03*nb2fb; Nevents=30000.;    genEff=1.0000E+00;}
	else if(name.Contains("ZtautauNp5"))      { sigma=1.112190E-03*nb2fb; Nevents=10000.;    genEff=1.0000E+00;}
	
	else _FAT("sample name: "<<name<<" is not supported");
	float Lmc = (sigmaSM>0.) ? (Nevents/genEff)/(sigmaSM*(BRbelle*BRfactor)) : (Nevents/genEff)/sigma;
	floats["wgt_luminosity"] = luminosity/Lmc; // correct the value that is there by default
}

void loop(TFile* ifile, TFile* ofile_loose, TFile* ofile_tight, TString name, TString jetMode, TString trkMode, TString metMode, TString metType, bool doMVA=false, bool doMVAout=false, bool isBlinded=true)
{	
	///////////////
	//// input tree
	TTree* itree = (TTree*)ifile->Get("flatout_"+name);
	setBranches(itree);
	_INF(vis,"");
	
	
	////////////////
	//// output tree
	// TTree *otree_loose = NULL;
	// if(!doMVA)
	// {
	// 	ofile_loose->cd();
	// 	otree_loose = itree->CloneTree(0);
	// 	
	// 	////////////////
	// 	//// decorations
	// 	for(TMapTSP2vi::iterator it=vints_decorations.begin()   ; it!=vints_decorations.end()   ; ++it) { it->second->clear(); otree_loose->Branch(it->first,&it->second); }
	// 	for(TMapTSP2vf::iterator it=vfloats_decorations.begin() ; it!=vfloats_decorations.end() ; ++it) { it->second->clear(); otree_loose->Branch(it->first,&it->second); }
	// }
	ofile_loose->cd();
	TTree *otree_loose = itree->CloneTree(0);
	for(TMapTSP2vi::iterator it=vints_decorations.begin()   ; it!=vints_decorations.end()   ; ++it) { it->second->clear(); otree_loose->Branch(it->first,&it->second); }
	for(TMapTSP2vf::iterator it=vfloats_decorations.begin() ; it!=vfloats_decorations.end() ; ++it) { it->second->clear(); otree_loose->Branch(it->first,&it->second); }
	_INF(vis,"");
	ofile_tight->cd();
	TTree *otree_tight = itree->CloneTree(0);
	for(TMapTSP2vi::iterator it=vints_decorations.begin()   ; it!=vints_decorations.end()   ; ++it) { it->second->clear(); otree_tight->Branch(it->first,&it->second); }
	for(TMapTSP2vf::iterator it=vfloats_decorations.begin() ; it!=vfloats_decorations.end() ; ++it) { it->second->clear(); otree_tight->Branch(it->first,&it->second); }
	_INF(vis,"");
	
	
	/////////////////////////////////
	//// simple output with MVA score
	if(doMVAout) setMVAout(name);
	if(!doMVA)   setCUTSout(name);
	_INF(vis,"");
	
	
	//// TMVA initialization
	if(doMVA)
	{
		clearMVAvars();
		bookMVAvars(metType);
		initTMVA("weights/","BDTG");
		_INF(vis,"");
	}
	_INF(vis,"");
	
	
	/////////////
	//// counters
	clearCounters();
	initCounters(doMVA,doMVAout,isBlinded);
	_INF(vis,"");
	
	
	/////////////////////////////////////////////////
	//// clear RunInfo vector (for duplicated events)
	clearRunInfo();
	_INF(vis,"");
	
	
	/////////////////////////
	//// Start the event loop
	Long64_t nentries = itree->GetEntries();
	if(name=="Data" && nData>0) nentries=nData;
	for(Long64_t i=0;i<nentries; i++)
	{
		itree->GetEntry(i);
		// int allPassing = getCounter("nPassing_evt_all");
		_INF(vis,"");
		
		
		///////////////////////////
		//// clear MVA out branches
		if(doMVAout) clearMVAout();
		if(!doMVA)   clearCUTSout();
		_INF(vis,"");
		
		
		////////////
		//// weights
		// float wgt = floats["wgt_total"];
		correctLumiWeight(name);
		float wgt = floats["wgt_luminosity"]*
					floats["wgt_kfactor"]*
					floats["wgt_dijets"]*
					floats["wgt_pileup"];
		_INF(vis,"");
		
		
		/////////////
		//// counters
		resetCounterFlags();
		_INF(vis,"");

		//////////////////
		//// clear decorations
		for(TMapTSP2vi::iterator it=vints_decorations.begin()   ; it!=vints_decorations.end()   ; ++it) it->second->clear();
		for(TMapTSP2vf::iterator it=vfloats_decorations.begin() ; it!=vfloats_decorations.end() ; ++it) it->second->clear();
		_INF(vis,"");
		
		
		// //////////////////////////////////
		// //// find if RunInfo is duplicated
		// bool isDup = isDuplicatedEvent(ints["evt_EventNumber"],ints["evt_lbn"],ints["evt_RunNumber"]);
		// if(isDup) continue;
		// _INF(vis,"");
		
		
		int nstoredvertices = ints["vtx_n"];
		if(isCounter("nPassing_evt_all") && nstoredvertices<1) continue;
		incrementCounter("nPassing_evt_all",wgt);
		_INF(vis,"");
		
		
		//////////////////////
		//// event level cuts
		bool trigger = false;
		unsigned int ntriggers = triggers.size();
		for(unsigned int t=0 ; t<ntriggers ; ++t) { if(ints[triggers[t]]==1) { trigger = true ; break; } }
		if(isCounter("nPassing_evt_trigger") && !trigger) continue;
		incrementCounter("nPassing_evt_trigger",wgt);
		_INF(vis,"");
		

		//////////////////////
		//// vertex level cuts
		vector<unsigned int> vertices;
		unsigned int nvtx = vfloats["vtx_mass"]->size();
		_INF(vis,"");
		
		bool oneCand = (nvtx==1);
		if(isCounter("nPassing_evt_oneCandidate") && !oneCand) continue;
		incrementCounter("nPassing_evt_oneCandidate",wgt);
		_INF(vis,"");
		
		for(unsigned int vtx=0 ; vtx<nvtx ; ++vtx)
		{
			///////////////////////////////////////////////
			//// set the default jets and fix the jet info
			setJETmode(vtx,vints,vfloats,jetMode);
			setMETmode(vtx,floats,vfloats,metMode);
			setTrkMETmode(vtx,floats,vfloats,trkMode);
			fillJetCalibrationHistos(vtx,name,"_all",floats,vfloats,metType,wgt);
			fixJets(vtx,vints,vfloats);
			_INF(vis,"");
			
			////////////////////////////////////////
			//// fill the decorations before any cut
			TLorentzVector p3body = getP3body(vtx,vfloats);
			float trksprob    = vfloats["mu_pvaltrkfit1"]->at(vtx)*vfloats["mu_pvaltrkfit2"]->at(vtx)*vfloats["mu_pvaltrkfit3"]->at(vtx);
			float trksdEdx    = vfloats["mu_pixeldEdx1"]->at(vtx)+vfloats["mu_pixeldEdx2"]->at(vtx)+vfloats["mu_pixeldEdx3"]->at(vtx);
			float lxySig      = vfloats["vtx_lxy"]->at(vtx)/vfloats["vtx_lxyErr"]->at(vtx);
			float a0xySig     = vfloats["vtx_a0xy"]->at(vtx)/vfloats["vtx_a0xyErr"]->at(vtx);	
			float minPbalSig  = +1e20; float minSignedPbalSig  = +1e20;
			float maxPbalSig  = -1e20; float maxSignedPbalSig  = -1e20;
			if(fabs(vfloats["mu_pbalsig1"]->at(vtx))<minPbalSig && vints["mu_type1"]->at(vtx)==MUON) { minPbalSig = fabs(vfloats["mu_pbalsig1"]->at(vtx)); minSignedPbalSig = vfloats["mu_pbalsig1"]->at(vtx); }
			if(fabs(vfloats["mu_pbalsig2"]->at(vtx))<minPbalSig && vints["mu_type2"]->at(vtx)==MUON) { minPbalSig = fabs(vfloats["mu_pbalsig2"]->at(vtx)); minSignedPbalSig = vfloats["mu_pbalsig2"]->at(vtx); }
			if(fabs(vfloats["mu_pbalsig3"]->at(vtx))<minPbalSig && vints["mu_type3"]->at(vtx)==MUON) { minPbalSig = fabs(vfloats["mu_pbalsig3"]->at(vtx)); minSignedPbalSig = vfloats["mu_pbalsig3"]->at(vtx); }
			if(fabs(vfloats["mu_pbalsig1"]->at(vtx))>maxPbalSig && vints["mu_type1"]->at(vtx)==MUON) { maxPbalSig = fabs(vfloats["mu_pbalsig1"]->at(vtx)); maxSignedPbalSig = vfloats["mu_pbalsig1"]->at(vtx); }
			if(fabs(vfloats["mu_pbalsig2"]->at(vtx))>maxPbalSig && vints["mu_type2"]->at(vtx)==MUON) { maxPbalSig = fabs(vfloats["mu_pbalsig2"]->at(vtx)); maxSignedPbalSig = vfloats["mu_pbalsig2"]->at(vtx); }
			if(fabs(vfloats["mu_pbalsig3"]->at(vtx))>maxPbalSig && vints["mu_type3"]->at(vtx)==MUON) { maxPbalSig = fabs(vfloats["mu_pbalsig3"]->at(vtx)); maxSignedPbalSig = vfloats["mu_pbalsig3"]->at(vtx); }
			float minScatSig  = +1e20; float minSignedScatSig  = +1e20;
			float maxScatSig  = -1e20; float maxSignedScatSig  = -1e20;
			if(fabs(vfloats["mu_sctangsig1"]->at(vtx))<minScatSig && vints["mu_type1"]->at(vtx)==MUON) { minScatSig = fabs(vfloats["mu_sctangsig1"]->at(vtx)); minSignedScatSig = vfloats["mu_sctangsig1"]->at(vtx); }
			if(fabs(vfloats["mu_sctangsig2"]->at(vtx))<minScatSig && vints["mu_type2"]->at(vtx)==MUON) { minScatSig = fabs(vfloats["mu_sctangsig2"]->at(vtx)); minSignedScatSig = vfloats["mu_sctangsig2"]->at(vtx); }
			if(fabs(vfloats["mu_sctangsig3"]->at(vtx))<minScatSig && vints["mu_type3"]->at(vtx)==MUON) { minScatSig = fabs(vfloats["mu_sctangsig3"]->at(vtx)); minSignedScatSig = vfloats["mu_sctangsig3"]->at(vtx); }
			if(fabs(vfloats["mu_sctangsig1"]->at(vtx))>maxScatSig && vints["mu_type1"]->at(vtx)==MUON) { maxScatSig = fabs(vfloats["mu_sctangsig1"]->at(vtx)); maxSignedScatSig = vfloats["mu_sctangsig1"]->at(vtx); }
			if(fabs(vfloats["mu_sctangsig2"]->at(vtx))>maxScatSig && vints["mu_type2"]->at(vtx)==MUON) { maxScatSig = fabs(vfloats["mu_sctangsig2"]->at(vtx)); maxSignedScatSig = vfloats["mu_sctangsig2"]->at(vtx); }
			if(fabs(vfloats["mu_sctangsig3"]->at(vtx))>maxScatSig && vints["mu_type3"]->at(vtx)==MUON) { maxScatSig = fabs(vfloats["mu_sctangsig3"]->at(vtx)); maxSignedScatSig = vfloats["mu_sctangsig3"]->at(vtx); }			
			float minNegbSig  = +1e20; float minSignedNegbSig  = +1e20;
			float maxNegbSig  = -1e20; float maxSignedNegbSig  = -1e20;
			if(fabs(vfloats["mu_sctngbsig1"]->at(vtx))<minNegbSig && vints["mu_type1"]->at(vtx)==MUON) { minNegbSig = fabs(vfloats["mu_sctngbsig1"]->at(vtx)); minSignedNegbSig = vfloats["mu_sctngbsig1"]->at(vtx); }
			if(fabs(vfloats["mu_sctngbsig2"]->at(vtx))<minNegbSig && vints["mu_type2"]->at(vtx)==MUON) { minNegbSig = fabs(vfloats["mu_sctngbsig2"]->at(vtx)); minSignedNegbSig = vfloats["mu_sctngbsig2"]->at(vtx); }
			if(fabs(vfloats["mu_sctngbsig3"]->at(vtx))<minNegbSig && vints["mu_type3"]->at(vtx)==MUON) { minNegbSig = fabs(vfloats["mu_sctngbsig3"]->at(vtx)); minSignedNegbSig = vfloats["mu_sctngbsig3"]->at(vtx); }
			if(fabs(vfloats["mu_sctngbsig1"]->at(vtx))>maxNegbSig && vints["mu_type1"]->at(vtx)==MUON) { maxNegbSig = fabs(vfloats["mu_sctngbsig1"]->at(vtx)); maxSignedNegbSig = vfloats["mu_sctngbsig1"]->at(vtx); }
			if(fabs(vfloats["mu_sctngbsig2"]->at(vtx))>maxNegbSig && vints["mu_type2"]->at(vtx)==MUON) { maxNegbSig = fabs(vfloats["mu_sctngbsig2"]->at(vtx)); maxSignedNegbSig = vfloats["mu_sctngbsig2"]->at(vtx); }
			if(fabs(vfloats["mu_sctngbsig3"]->at(vtx))>maxNegbSig && vints["mu_type3"]->at(vtx)==MUON) { maxNegbSig = fabs(vfloats["mu_sctngbsig3"]->at(vtx)); maxSignedNegbSig = vfloats["mu_sctngbsig3"]->at(vtx); }

			_INF(vis,"");

			vfloats_decorations["mva_score"]->push_back(-2);
			vfloats_decorations["trks_fitprob"]->push_back(trksprob);
			vfloats_decorations["trks_pixdEdx"]->push_back(trksdEdx);
			vfloats_decorations["geo_lxySig"]->push_back(lxySig);
			vfloats_decorations["geo_a0xySig"]->push_back(a0xySig);
			vfloats_decorations["muons_minpbalsig"]->push_back(minSignedPbalSig);
			vfloats_decorations["muons_maxpbalsig"]->push_back(maxSignedPbalSig);
			vfloats_decorations["muons_minscatsig"]->push_back(minSignedScatSig);
			vfloats_decorations["muons_maxscatsig"]->push_back(maxSignedScatSig);
			vfloats_decorations["muons_minnegbsig"]->push_back(minSignedNegbSig);
			vfloats_decorations["muons_maxnegbsig"]->push_back(maxSignedNegbSig);

			_INF(vis,"");

			//////////////////////////
			//// Jet/MET variables 
			//// start from the end and finish with the calibrated one
			float bTagCut       = (isAntiKt4LCJet) ? 0.3511 : 0.39;
			TLorentzVector pSum         ; 
			int   ngoodjets     = -9999 ; 
			int   ngoodbjets    = -9999 ;
			float dphiHTMETref  = -9999.; 
			float dphiHTMETmu   = -9999.; 
			float dphiHTMETtrk  = -9999.; 
			float dphiHT3body   = -9999.; 
			float dRHT3body     = -9999.; 
			float mtHT          = -9999.; 
			float mtHTtrk       = -9999.; 
			float METsig        = -9999.; 
			float metsDphi      = -9999.; 
			float metsdpTrel    = -9999.; 
			float metsdpTrelTrk = -9999.; 
			float metsdpTrelCal = -9999.; 
			float metsdHTrelTrk = -9999.; 
			float metsdHTrelCal = -9999.;
			for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
			{
				TString word = getTrkJetMETword(mode);
				_INF(vis,"mode="<<mode<<", word="<<word);
				
				TString trkmetword = word;
				TString jetmetword = word;
				if(mode>=ALL_JESUP  && mode<=ALL_JERDWN) trkmetword="";
				if(mode>=ALL_SOFTUP && mode<=ALL_JETDWN) jetmetword="";
				TString jointword = (mode==ALL_UNCALIB) ? "_uncalib" : jetmetword+trkmetword;
				_INF(vis,"jetmetword="<<jetmetword<<", trkmetword="<<trkmetword);
			
			    pSum          = getPsum(vtx,vfloats,1,jetmetword);
				ngoodjets     = (vfloats["jet_pt1"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV)+(vfloats["jet_pt2"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV)+(vfloats["jet_pt3"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV)+(vfloats["jet_pt4"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV);
				ngoodbjets    = (vfloats["jet_pt1"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV && vfloats["jet_MV1w1"+jetmetword]->at(vtx)>bTagCut)+(vfloats["jet_pt2"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV && vfloats["jet_MV1w2"+jetmetword]->at(vtx)>bTagCut)+(vfloats["jet_pt3"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV && vfloats["jet_MV1w3"+jetmetword]->at(vtx)>bTagCut)+(vfloats["jet_pt4"+jetmetword]->at(vtx)>minJetPtGeV*GeV2MeV && vfloats["jet_MV1w4"+jetmetword]->at(vtx)>bTagCut);
				dphiHTMETref  = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-floats["met_reffinal_phi"+jetmetword]));
				dphiHTMETmu   = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-floats["met_muons_phi"+jetmetword]));
				dphiHTMETtrk  = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-floats["met_track_phi"+trkmetword]));
				dphiHT3body   = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-p3body.Phi()));
				dRHT3body     = fabs(pSum.DeltaR(p3body));
				mtHT          = sqrt(2*pSum.Pt()*floats["met_"+metType+"_et"+jetmetword]*(1-TMath::Cos(pSum.Phi()-floats["met_"+metType+"_phi"+jetmetword])));
				mtHTtrk       = sqrt(2*pSum.Pt()*floats["met_track_et"+trkmetword]*(1-TMath::Cos(pSum.Phi()-floats["met_track_phi"+trkmetword])));
				METsig        = floats["met_"+metType+"_et"+jetmetword]/sqrt(pSum.Pt());
				metsDphi      = fabs(TVector2::Phi_mpi_pi(floats["met_"+metType+"_phi"+jetmetword]-floats["met_track_phi"+trkmetword]));
				metsdpTrel    = (vfloats["vtx_pt"]->at(vtx)-(floats["met_"+metType+"_et"+jetmetword]+floats["met_track_et"+trkmetword])/2.)/vfloats["vtx_pt"]->at(vtx);
				metsdpTrelTrk = vfloats["vtx_pt"]->at(vtx)/floats["met_track_et"+trkmetword]-1.;
				metsdpTrelCal = vfloats["vtx_pt"]->at(vtx)/floats["met_"+metType+"_et"+jetmetword]-1.;
				metsdHTrelTrk = pSum.Pt()/floats["met_track_et"+trkmetword]-1.;
				metsdHTrelCal = pSum.Pt()/floats["met_"+metType+"_et"+jetmetword]-1.;
			
				vints_decorations["jets_n"+jetmetword]->push_back(ngoodjets);
				vints_decorations["jets_b"+jetmetword]->push_back(ngoodbjets);
				vfloats_decorations["ht_pt"+jetmetword]->push_back(pSum.Pt());
				vfloats_decorations["ht_mT"+jetmetword]->push_back(mtHT);
				vfloats_decorations["ht_mT_mettrk"+trkmetword]->push_back(mtHTtrk);
				vfloats_decorations["ht_metsig"+jetmetword]->push_back(METsig);
				vfloats_decorations["ht_dphimet_reffinal"+jetmetword]->push_back(dphiHTMETref);
				vfloats_decorations["ht_dphimet_muons"+jetmetword]->push_back(dphiHTMETmu);
				vfloats_decorations["ht_dphimet_track"+jointword]->push_back(dphiHTMETtrk);
				vfloats_decorations["ht_dr3body"+jetmetword]->push_back(dRHT3body);
				vfloats_decorations["ht_dphi3body"+jetmetword]->push_back(dphiHT3body);
				vfloats_decorations["mets_dphi"+jointword]->push_back(metsDphi);
				vfloats_decorations["mets_dptrelavg"+jointword]->push_back(metsdpTrel);
				vfloats_decorations["mets_dptreltrk"+trkmetword]->push_back(metsdpTrelTrk);
				vfloats_decorations["mets_dptrelcal"+jetmetword]->push_back(metsdpTrelCal);
				vfloats_decorations["mets_dhtreltrk"+trkmetword]->push_back(metsdHTrelTrk);
				vfloats_decorations["mets_dhtrelcal"+jetmetword]->push_back(metsdHTrelCal);
			}
			pSum = getPsum(vtx,vfloats,1,""); // just to make sure the "calibrated" is taken elsewhere (blanc suffix)
			_INF(vis,"");
			
			//////////////////////////////////////////////////////////////////
			//// initialize the flag with 0 and fill it later with 1 if passes
			vints_decorations["pass_loose"]->push_back(0);      //// this is very important to do before any vertex cut!!!
			vints_decorations["pass_tight"]->push_back(0);      //// this is very important to do before any vertex cut!!!
			vints_decorations["pass_tight_scat"]->push_back(0); //// this is very important to do before any vertex cut!!!
			vints_decorations["pass_tight_pbal"]->push_back(0); //// this is very important to do before any vertex cut!!!
			vints_decorations["pass_tight_fitq"]->push_back(0); //// this is very important to do before any vertex cut!!!
			_INF(vis,"");
			
			
			
			/////////////////////////////////////////////////////////
			//// the category, hit and trigger histos before any cut
			int category = vints["vtx_code"]->at(vtx);
			fillPVHistos(vtx,name,"_after_vertexing",vints,wgt);
			fillCategories(category,name,"tripletCategories_afterVertexing");
			fillCategories(category,name,"tripletCategories_norm_afterVertexing");
			fillHitHistos(vtx,name,vints,wgt);
			fillTriggerHistos(name,"after_vertexing",vints,wgt);
			fillRhoOmegaPhiHistos(vtx,name,"_beforeObjectcuts",floats,vfloats,metType,wgt);
			_INF(vis,"");
			
			
			
			////////////////////
			//// begin cuts ////
			////////////////////
			
			
	
			
			bool vtxcat = ((vints["mu_type1"]->at(vtx)==MUON || vints["mu_type1"]->at(vtx)==TPA) &&
						   (vints["mu_type2"]->at(vtx)==MUON || vints["mu_type2"]->at(vtx)==TPA) &&
						   (vints["mu_type3"]->at(vtx)==MUON || vints["mu_type3"]->at(vtx)==TPA) &&
						  category!=MUONS1TPA2);
			if(isCounter("nPassing_evt_goodtriplets") && !vtxcat) continue;
			incrementCounter("nPassing_evt_goodtriplets",wgt);
			_INF(vis,"");
	
			
			fillAllObjectHistos(vtx,name,vints,vfloats,wgt);
			_INF(vis,"");


			bool medium = (isMedium(vtx,1,vints)              && isMedium(vtx,2,vints)              && isMedium(vtx,3,vints));
			bool mcpid  = (MCP(vtx,1,vints,vfloats)           && MCP(vtx,2,vints,vfloats)           && MCP(vtx,3,vints,vfloats));
			bool trt    = (highThresholdTRTRhits(vtx,1,vints) && highThresholdTRTRhits(vtx,2,vints) && highThresholdTRTRhits(vtx,3,vints));
			bool layers = (TPaEtaPhiMSlayers(vtx,1,vints)     && TPaEtaPhiMSlayers(vtx,2,vints)     && TPaEtaPhiMSlayers(vtx,3,vints));
			bool tgcrpc = (TGCRPChits(vtx,1,vints)            && TGCRPChits(vtx,2,vints)            && TGCRPChits(vtx,3,vints));
			bool hits   = (trt && layers && tgcrpc);
			if(isCounter("nPassing_mu_mcp") && !(mcpid && medium)) continue;
			incrementCounter("nPassing_mu_mcp",wgt);
			_INF(vis,"");
			
			
			fillObjectValidationHistos(vtx,name,"_onJpsi",vints,vfloats,wgt);
			fillObjectValidationHistos(vtx,name,"_offJpsi",vints,vfloats,wgt);
			_INF(vis,"");
			
			
			if(isCounter("nPassing_mu_hits") && !hits) continue;
			incrementCounter("nPassing_mu_hits",wgt);
			_INF(vis,"");
			
			
			bool mupt1 = (vfloats["mu_pt1"]->at(vtx)>5.5*GeV2MeV);
			bool mupt2 = (vfloats["mu_pt2"]->at(vtx)>3.5*GeV2MeV);
			bool mupt3 = (vfloats["mu_pt3"]->at(vtx)>2.5*GeV2MeV);
			mupt1 = (vints["mu_type1"]->at(vtx)==TPB) ? (mupt1 && vfloats["mu_pt1"]->at(vtx)>3.0*GeV2MeV) : mupt1;
			mupt2 = (vints["mu_type2"]->at(vtx)==TPB) ? (mupt2 && vfloats["mu_pt2"]->at(vtx)>3.0*GeV2MeV) : mupt2;
			mupt3 = (vints["mu_type3"]->at(vtx)==TPB) ? (mupt3 && vfloats["mu_pt3"]->at(vtx)>3.0*GeV2MeV) : mupt3;
			if(isCounter("nPassing_mu_pt") && mupt1+mupt2+mupt3!=3) continue;
			incrementCounter("nPassing_mu_pt",wgt);
			_INF(vis,"");
			
			
			bool mueta1 = (fabs(vfloats["mu_eta1"]->at(vtx))<2.5);
			bool mueta2 = (fabs(vfloats["mu_eta2"]->at(vtx))<2.5);
			bool mueta3 = (fabs(vfloats["mu_eta3"]->at(vtx))<2.5);
			if(isCounter("nPassing_mu_eta") && mueta1+mueta2+mueta3!=3) continue;
			incrementCounter("nPassing_mu_eta",wgt);
			_INF(vis,"");
			
			
			bool loosesctangsig = ((vints["mu_type1"]->at(vtx)!=MUON || (vints["mu_type1"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig1"]->at(vtx))<4.0 && fabs(vfloats["mu_sctngbsig1"]->at(vtx))<4.0)) &&
								   (vints["mu_type2"]->at(vtx)!=MUON || (vints["mu_type2"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig2"]->at(vtx))<4.0 && fabs(vfloats["mu_sctngbsig2"]->at(vtx))<4.0)) &&
								   (vints["mu_type3"]->at(vtx)!=MUON || (vints["mu_type3"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig3"]->at(vtx))<4.0 && fabs(vfloats["mu_sctngbsig3"]->at(vtx))<4.0)));
			bool tightsctangsig = ((vints["mu_type1"]->at(vtx)!=MUON || (vints["mu_type1"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig1"]->at(vtx))<2.5 && fabs(vfloats["mu_sctngbsig1"]->at(vtx))<2.5)) &&
			/*for phi rejection*/  (vints["mu_type2"]->at(vtx)!=MUON || (vints["mu_type2"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig2"]->at(vtx))<2.5 && fabs(vfloats["mu_sctngbsig2"]->at(vtx))<2.5)) &&
								   (vints["mu_type3"]->at(vtx)!=MUON || (vints["mu_type3"]->at(vtx)==MUON && fabs(vfloats["mu_sctangsig3"]->at(vtx))<2.5 && fabs(vfloats["mu_sctngbsig3"]->at(vtx))<2.5)));
			vints_decorations["pass_tight_scat"]->at(vtx) = tightsctangsig;
			if(isCounter("nPassing_mu_sctangsig") && !loosesctangsig) continue;                                                                       
			incrementCounter("nPassing_mu_sctangsig",wgt);
			_INF(vis,"");
			
																																		/*was: -3 */						/*was: +3*/	
			bool loosepbalsig = ((vints["mu_type1"]->at(vtx)!=MUON || (vints["mu_type1"]->at(vtx)==MUON && (vfloats["mu_pbalsig1"]->at(vtx)>-3. && vfloats["mu_pbalsig1"]->at(vtx)<4.))) &&
							     (vints["mu_type2"]->at(vtx)!=MUON || (vints["mu_type2"]->at(vtx)==MUON && (vfloats["mu_pbalsig2"]->at(vtx)>-3. && vfloats["mu_pbalsig2"]->at(vtx)<4.))) &&
							     (vints["mu_type3"]->at(vtx)!=MUON || (vints["mu_type3"]->at(vtx)==MUON && (vfloats["mu_pbalsig3"]->at(vtx)>-3. && vfloats["mu_pbalsig3"]->at(vtx)<4.))));
			bool tightpbalsig = ((vints["mu_type1"]->at(vtx)!=MUON || (vints["mu_type1"]->at(vtx)==MUON && (vfloats["mu_pbalsig1"]->at(vtx)>-2. && vfloats["mu_pbalsig1"]->at(vtx)<3.))) &&
			/*for phi rejection*/(vints["mu_type2"]->at(vtx)!=MUON || (vints["mu_type2"]->at(vtx)==MUON && (vfloats["mu_pbalsig2"]->at(vtx)>-2. && vfloats["mu_pbalsig2"]->at(vtx)<3.))) &&
							     (vints["mu_type3"]->at(vtx)!=MUON || (vints["mu_type3"]->at(vtx)==MUON && (vfloats["mu_pbalsig3"]->at(vtx)>-2. && vfloats["mu_pbalsig3"]->at(vtx)<3.))));
			vints_decorations["pass_tight_pbal"]->at(vtx) = tightpbalsig;
			if(isCounter("nPassing_mu_pbalsig") && !loosepbalsig) continue;
			incrementCounter("nPassing_mu_pbalsig",wgt);
			_INF(vis,"");
							
																																			/* was: 0.01*/
			bool loosetrkquality = (((vints["mu_type1"]->at(vtx)==MUON) || (vints["mu_type1"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit1"]->at(vtx)>0.01) || (vints["mu_type1"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit1"]->at(vtx)>0.1)) &&
			                        ((vints["mu_type2"]->at(vtx)==MUON) || (vints["mu_type2"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit2"]->at(vtx)>0.01) || (vints["mu_type2"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit2"]->at(vtx)>0.1)) &&
							        ((vints["mu_type3"]->at(vtx)==MUON) || (vints["mu_type3"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit3"]->at(vtx)>0.01) || (vints["mu_type3"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit3"]->at(vtx)>0.1)));
			bool tighttrkquality = (((vints["mu_type1"]->at(vtx)==MUON) || (vints["mu_type1"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit1"]->at(vtx)>0.03) || (vints["mu_type1"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit1"]->at(vtx)>0.2)) &&
			/*for phi rejection*/   ((vints["mu_type2"]->at(vtx)==MUON) || (vints["mu_type2"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit2"]->at(vtx)>0.03) || (vints["mu_type2"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit2"]->at(vtx)>0.2)) &&
							        ((vints["mu_type3"]->at(vtx)==MUON) || (vints["mu_type3"]->at(vtx)==TPA && vfloats["mu_pvaltrkfit3"]->at(vtx)>0.03) || (vints["mu_type3"]->at(vtx)==TPB && vfloats["mu_pvaltrkfit3"]->at(vtx)>0.2)));
			//// loosetrkquality = (loosetrkquality && vfloats_decorations["trks_fitprob"]->at(vtx)>1.e-9); // (1e-9)^(1/3) = 0.001
			vints_decorations["pass_tight_fitq"]->at(vtx) = tighttrkquality;
			if(isCounter("nPassing_mu_trkquality") && !loosetrkquality) continue;
			incrementCounter("nPassing_mu_trkquality",wgt);
			_INF(vis,"");
			
			
			fillObjectHistos(vtx,name,"_after_muons",vints,vfloats,wgt); 
			fillCategories(category,name,"tripletCategories_after_muons");
			fillCategories(category,name,"tripletCategories_norm_after_muons");
			fillMETHistos(vtx,name,"_after_muons",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_muons",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_muons",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_muons",floats,vfloats,metType,wgt);
			fillVertexHistos(vtx,name,vfloats,wgt);
			_INF(vis,"");
			
			
			bool trksProbLoose = (vfloats_decorations["trks_fitprob"]->at(vtx)>1.e-9); // (1e-9)^(1/3) = 0.001
			if(isCounter("nPassing_3body_trksprobLoose") && !trksProbLoose) continue;
			incrementCounter("nPassing_3body_trksprobLoose",wgt);
			_INF(vis,"");
			
			
			bool pvalueLoose = (vfloats["vtx_pval"]->at(vtx)>1.e-7);
			if(isCounter("nPassing_3body_pvalueLoose") && !pvalueLoose) continue;
			incrementCounter("nPassing_3body_pvalueLoose",wgt);
			_INF(vis,"");
			
			
			//bool vtxclean = (vfloats["vtx_lxy"]->at(vtx)>-1. &&  vfloats["vtx_lxy"]->at(vtx)<12.  &&  fabs(vfloats["vtx_a0xy"]->at(vtx))<0.1);
			// bool vtxclean = (vfloats["vtx_lxy"]->at(vtx)>-5. &&  vfloats["vtx_lxy"]->at(vtx)<15.  &&  fabs(vfloats["vtx_a0xy"]->at(vtx))<0.5);
			bool vtxclean  = (vfloats_decorations["geo_lxySig"]->at(vtx)>-10 && vfloats_decorations["geo_lxySig"]->at(vtx)<50 && vfloats_decorations["geo_a0xySig"]->at(vtx)<25);
			if(isCounter("nPassing_3body_vtxcleanLoose") && !vtxclean) continue;
			incrementCounter("nPassing_3body_vtxcleanLoose",wgt);
			_INF(vis,"");
			
			
			bool massLoose = (vfloats["vtx_mass"]->at(vtx)<4.*GeV2MeV);
			if(isCounter("nPassing_3body_m") && !massLoose) continue;
			incrementCounter("nPassing_3body_m",wgt);
			_INF(vis,"");
			
			
			bool ptLoose = (vfloats["vtx_pt"]->at(vtx)>10.*GeV2MeV);
			if(isCounter("nPassing_3body_ptLoose") && !ptLoose) continue;
			incrementCounter("nPassing_3body_ptLoose",wgt);
			_INF(vis,"");
			
			
			bool charge = (fabs(vfloats["vtx_charge"]->at(vtx))==1);
			if(isCounter("nPassing_3body_charge") && !charge) continue;
			incrementCounter("nPassing_3body_charge",wgt);
			_INF(vis,"");
			
			
			bool passm2body = (vfloats["vtx_mOS1"]->at(vtx)>220 && vfloats["vtx_mOS2"]->at(vtx)>220 && vfloats["vtx_mSS"]->at(vtx)>220);
			if(isCounter("nPassing_3body_m2body") && !passm2body) continue;
			incrementCounter("nPassing_3body_m2body",wgt);
			_INF(vis,"");
			

			float m2SS = vfloats["vtx_mSS"]->at(vtx)*vfloats["vtx_mSS"]->at(vtx);
			float m2OS1 = vfloats["vtx_mOS1"]->at(vtx)*vfloats["vtx_mOS1"]->at(vtx);
			float m2OS2 = vfloats["vtx_mOS2"]->at(vtx)*vfloats["vtx_mOS2"]->at(vtx);
			bool allFound = (m2SS>0. && m2OS1>0. && m2OS2>0.);
			bool mSSOSDalitzLow1  = ((m2SS+m2OS1)<500.e3);
			bool mSSOSDalitzLow2  = ((m2SS+m2OS2)<500.e3);
			bool mSSOSDalitzHigh1 = ((m2SS+m2OS1)>7300.e3);
			bool mSSOSDalitzHigh2 = ((m2SS+m2OS2)>7300.e3);
			bool mOSOSDalitzLow   = ((m2OS2+m2OS1)<500.e3);
			bool mOSOSDalitzHigh  = ((m2OS2+m2OS1)>7300.e3);
			bool failDalitz  = (allFound && (mSSOSDalitzHigh1 || mSSOSDalitzHigh2 || mSSOSDalitzLow1 || mSSOSDalitzLow2 || mOSOSDalitzLow || mOSOSDalitzHigh));
			if(isCounter("nPassing_3body_dalitz") && failDalitz) continue;
			incrementCounter("nPassing_3body_dalitz",wgt);
			_INF(vis,"");
			

			// // float range = SigmaRhoOmegaPhiMeV;
			// // bool onPhi   = (fabs(vfloats["vtx_mOS1"]->at(vtx)-1020.)<range || fabs(vfloats["vtx_mOS2"]->at(vtx)-1020.)<range); // || fabs(vfloats["vtx_mSS"]->at(vtx)-1020.)<range);
			// // bool onOmega = (fabs(vfloats["vtx_mOS1"]->at(vtx)-782.)<range  || fabs(vfloats["vtx_mOS2"]->at(vtx)-782.)<range ); // || fabs(vfloats["vtx_mSS"]->at(vtx)-782.)<range);
			// // bool onRho   = (fabs(vfloats["vtx_mOS1"]->at(vtx)-770.)<range  || fabs(vfloats["vtx_mOS2"]->at(vtx)-770.)<range  || fabs(vfloats["vtx_mSS"]->at(vtx)-770.)<range);
			// // bool failRhoOmegaPhi = (onPhi||onOmega||onRho);
			// // bool failRhoOmegaPhi = (onPhi||onOmega);
			// // if(isCounter("nPassing_3body_rhoomegaphi") && failRhoOmegaPhi) continue;
			// // bool failmasses     = (vfloats["vtx_mOS1"]->at(vtx)<=250 || vfloats["vtx_mOS2"]->at(vtx)<=250 || vfloats["vtx_mSS"]->at(vtx)<=250  ||  vfloats["vtx_mOS1"]->at(vtx)>=1600 || vfloats["vtx_mOS2"]->at(vtx)>=1600 || vfloats["vtx_mSS"]->at(vtx)>=1600);
			// bool failmasses     = (vfloats["vtx_mOS1"]->at(vtx)<=220 || vfloats["vtx_mOS2"]->at(vtx)<=220 || vfloats["vtx_mSS"]->at(vtx)<=220);
			// bool faildistances  = (vfloats["geo_lxySig"]->at(vtx)<=-10 || vfloats["geo_lxySig"]->at(vtx)>=50 || vfloats["geo_a0xySig"]->at(vtx)>=25);
			// bool failpvalues    = (vfloats_decorations["trks_fitprob"]->at(vtx)<=1e-9);
			// bool failkinematics = (vfloats["vtx_pt"]->at(vtx)<=10000 || (floats["met_"+metType+"_et"]<=10000 || floats["met_"+metType+"_et"]>=200000) || (floats["met_track_et"]<=10000 || floats["met_track_et"]>=200000) || vfloats["met_"+metType+"_mT"]->at(vtx)<=20000 || vfloats["met_track_mT"]->at(vtx)<=20000);
			// bool failisol       = (vfloats["vtx_isolation020"]->at(vtx)>=0.3 || vfloats["vtx_isolation003"]->at(vtx)>=0.1);
			// bool failHT         = (vfloats_decorations["ht_pt"]->at(vtx)<=20000 || vfloats_decorations["ht_dphimet_track"]->at(vtx)<=2);
			// bool failonecand    = (vfloats["vtx_mass"]->size()!=1);
			// // if(isCounter("nPassing_3body_rhoomegaphi") && (failRhoOmegaPhi || failmasses || faildistances || failpvalues || failkinematics || failisol || failonecand)) continue;
			// if(isCounter("nPassing_3body_rhoomegaphi") && (failpvalues || failonecand || failisol || failmasses || failkinematics || failHT)) continue;
			// incrementCounter("nPassing_3body_rhoomegaphi",wgt);
			// _INF(vis,"");
			

			float maxJpT = 30.*GeV2MeV;
			int nj = (vfloats["jet_pt1"]->at(vtx)>maxJpT)+(vfloats["jet_pt2"]->at(vtx)>maxJpT)+(vfloats["jet_pt3"]->at(vtx)>maxJpT)+(vfloats["jet_pt4"]->at(vtx)>maxJpT);
			bool failJets = (nj>0);
			if(isCounter("nPassing_3body_njets") && failJets) continue;
			incrementCounter("nPassing_3body_njets",wgt);
			_INF(vis,"");
			
			
			fillCategories(category,name,"tripletCategories_after_triplet");
			fillCategories(category,name,"tripletCategories_norm_after_triplet");
			fillMETHistos(vtx,name,"_after_triplet",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_triplet",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_triplet",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_triplet",floats,vfloats,metType,wgt);
			fillObjectHistos(vtx,name,"_after_triplet",vints,vfloats,wgt);
			_INF(vis,"");
			
			

			
			
			
			
			fillIsolationHistos(vtx,name,vfloats,wgt); // only before applying the loose isolation cut
			fillPVHistos(vtx,name,"_before_isolation",vints,wgt);
			_INF(vis,"");
		
			
			// bool isolationLoose = (vfloats["vtx_isolation003"]->at(vtx)<0.1 && vfloats["vtx_isolation020"]->at(vtx)<0.3);
			bool isolationLoose = (vfloats["vtx_isolation020"]->at(vtx)<0.3 && vfloats["vtx_isolation030"]->at(vtx)<1);
			if(isCounter("nPassing_3body_isolationLoose") && !isolationLoose) continue;
			incrementCounter("nPassing_3body_isolationLoose",wgt);
			_INF(vis,"");
			
			
			fillPVHistos(vtx,name,"_after_isolation",vints,wgt);
			fillCategories(category,name,"tripletCategories_after_hadclean");
			fillCategories(category,name,"tripletCategories_norm_after_hadclean");
			fillMETHistos(vtx,name,"_after_hadclean",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_hadclean",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_hadclean",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_hadclean",floats,vfloats,metType,wgt);
			fillObjectHistos(vtx,name,"_after_hadclean",vints,vfloats,wgt);
			fillRhoOmegaPhiHistos(vtx,name,"_beforeWcuts",floats,vfloats,metType,wgt);
			_INF(vis,"");

			
			
			
			if(isCounter("nPassing_met_metLoose") && !passMET(floats,metType,"loose")) continue;
			incrementCounter("nPassing_met_metLoose",wgt);
			_INF(vis,"");

			if(isCounter("nPassing_met_mTLoose") && !passMT(vtx,vfloats,metType,"loose")) continue;
			incrementCounter("nPassing_met_mTLoose",wgt);
			_INF(vis,"");
			
	
			fillCategories(category,name,"tripletCategories_after_met");
			fillCategories(category,name,"tripletCategories_norm_after_met");
			fillMETHistos(vtx,name,"_after_met",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_met",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_met",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_met",floats,vfloats,metType,wgt);
			fillObjectHistos(vtx,name,"_after_met",vints,vfloats,wgt);
			_INF(vis,"");
			
			

			if(isCounter("nPassing_ht_mhTLoose") && !passMHT(vtx,vfloats_decorations,"loose")) continue;
			incrementCounter("nPassing_ht_mhTLoose",wgt);
			_INF(vis,"");
			

			fillPVHistos(vtx,name,"_after_ht",vints,wgt);
			fillCategories(category,name,"tripletCategories_after_ht");
			fillCategories(category,name,"tripletCategories_norm_after_ht");
			fillMETHistos(vtx,name,"_after_ht",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_ht",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_ht",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_ht",floats,vfloats,metType,wgt);
			fillObjectHistos(vtx,name,"_after_ht",vints,vfloats,wgt);
			fillRhoOmegaPhiHistos(vtx,name,"_afterLoosecuts",floats,vfloats,metType,wgt);
			_INF(vis,"");
			
			
			
			/////////////////////
			//// fill decorations
			vints_decorations["pass_loose"]->at(vtx) = 1; // passed loose cuts ==> change flag from 0 to 1
			_INF(vis,"");
			
			
			//////////////////////////////////////////
			//// get MVA score and fill the decoration
			if(doMVA || doMVAout)
			{
				setMVAvars(vtx,ints,floats,vints,vfloats,vints_decorations,vfloats_decorations,minJetPtGeV*GeV2MeV);
				setMVAspect(vtx,ints,floats,vints,vfloats);
				float score = getMVAscore();
				vfloats_decorations["mva_score"]->at(vtx) = score;
				_INF(vis,"");
			}
			_INF(vis,"");
			
			////////////////////////////
			//// OK vertex for writeout
			vertices.push_back(vtx);
			_INF(vis,"");
		}
		
		
		_INF(vis,"");
		
		
		//////////////////////////////////////////
		//// fill the tree for TMVA classification
		ofile_loose->cd();
		if(vertices.size()>0) otree_loose->Fill();
		if(i%10000==0 && i!=0)
		{
			otree_loose->FlushBaskets();
			// otree_loose->Write("", TObject::kOverwrite);
		}
		olddir->cd();
		_INF(vis,"");
		
		
		//////////////////////////////////
		//// no need to continue this even
		//// if there are no loose vertices
		if(vertices.size()<1) continue;
		_INF(vis,"");
		/////////////////////////////////
		//// Only one candidate per event
		// unsigned int nloosevertices = tightvertices.size();
		// if(isCounter("nPassing_evt_onecandidate") && nloosevertices!=1) continue;
		// incrementCounter("nPassing_evt_onecandidate",wgt);
		// _INF(vis,"");
	
	
		
		vector<unsigned int> tightvertices;
		for(unsigned int v=0 ; v<vertices.size() ; ++v)
		{
			unsigned int vtx = vertices[v];
			int category = vints["vtx_code"]->at(vtx);
			_INF(vis,"");
			
			
			if(doMVA)
			{
				// setMVAvars(vtx,ints,floats,vints,vfloats,vints_decorations,vfloats_decorations,minJetPtGeV*GeV2MeV);
				// setMVAspect(vtx,ints,floats,vints,vfloats);
				// float score = getMVAscore();
				// fillMVAHistos(vtx,name,vfloats,score,wgt);
				// _INF(vis,"");
				// 
				// 
				// /////////////////////
				// //// fill decorations
				// vfloats_decorations["mva_score"]->at(vtx) = score;
				// _INF(vis,"");
				float score = vfloats_decorations["mva_score"]->at(vtx);
				_INF(vis,"");
				
				
				if(doMVAout)
				{
					//// fill the simple MVA out branches
					fillMVAoutVecVars(vtx,score,floats,vints,vfloats,vints_decorations,vfloats_decorations,metType);
					_INF(vis,"");
				}
				else
				{
					if(score>-0.95) fillMVAevoHistos(vtx,name,"_neg_veryloose",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>-0.80) fillMVAevoHistos(vtx,name,"_neg_loose",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>-0.50) fillMVAevoHistos(vtx,name,"_neg_medium",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>-0.30) fillMVAevoHistos(vtx,name,"_neg_tight",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>0.000) fillMVAevoHistos(vtx,name,"_zero",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>+0.30) fillMVAevoHistos(vtx,name,"_pos_loose",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>+0.50) fillMVAevoHistos(vtx,name,"_pos_medium",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>+0.80) fillMVAevoHistos(vtx,name,"_pos_tight",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					if(score>+0.95) fillMVAevoHistos(vtx,name,"_pos_verytight",floats,vfloats,vints_decorations,vfloats_decorations,metType,score,wgt);
					
					
					//// do not cut on the MVA score since
					//// we want to write all candidates
					bool failMVA = (score<optBDTcut);
					if(isCounter("nPassing_optMVA") && failMVA) continue;
					incrementCounter("nPassing_optMVA",wgt);
					_INF(vis,"");
				}
				
				_INF(vis,"");
			}
			else
			{
				// if(isWsignal(name))
				// {
				// 	(*ofstr1) << ints["evt_EventNumber"] << " " << ints["evt_lbn"] << endl;
				// 	if(vfloats["jet_pt1"]->at(vtx)>0) (*ofstr1) << "  pt1= " << vfloats["jet_pt1"]->at(vtx) << ", dphi3muj1=" << vfloats["jet_dphi3muJ1"]->at(vtx) << endl;
				// 	if(vfloats["jet_pt2"]->at(vtx)>0) (*ofstr1) << "  pt1= " << vfloats["jet_pt1"]->at(vtx) << ", pt2= "<< vfloats["jet_pt2"]->at(vtx) << ", dphij1j2=" << vfloats["jet_dphiJ1J2"]->at(vtx) << endl;
				// 	(*ofstr1) << "-------------------------------------------" << endl;
				// }
				
				// if(isCounter("nPassing_jets_bjetveto_std") && !passBjetVeto(vtx,vfloats,isAntiKt4LCJet,"tight")) continue; // loose is the same as tight
				// incrementCounter("nPassing_jets_bjetveto_std",wgt);
				// _INF(vis,"");
				// 
				// 
				// if(isCounter("nPassing_jets_coljetveto") && !passCollinearJetVeto(vtx,vfloats,"tight")) continue;
				// incrementCounter("nPassing_jets_coljetveto",wgt);
				// _INF(vis,"");
				// 
				// if(isWsignal(name)) (*ofstr) << ints["evt_EventNumber"] << " " << ints["evt_lbn"] << endl;
				// 
				// if(isCounter("nPassing_jets_dijetveto") && !passDiJetVeto(vtx,vfloats,"tight")) continue;
				// incrementCounter("nPassing_jets_dijetveto",wgt);
				// _INF(vis,"");
				
				
				fillCategories(category,name,"tripletCategories_after_hadcleanTight");
				fillCategories(category,name,"tripletCategories_norm_after_hadcleanTight");
				_INF(vis,"");
				
				
				bool ptTight = (vfloats["vtx_pt"]->at(vtx)>25.*GeV2MeV);
				if(isCounter("nPassing_3body_pt") && !ptTight) continue;
				incrementCounter("nPassing_3body_pt",wgt);
				_INF(vis,"");
				
				
				bool pvalueTight = (vfloats["vtx_pval"]->at(vtx)>0.1);
				if(isCounter("nPassing_3body_pvalue") && !pvalueTight) continue;
				incrementCounter("nPassing_3body_pvalue",wgt);
				_INF(vis,"");
				
				
				bool vtxclean = (vfloats["vtx_lxy"]->at(vtx)>0.  &&  fabs(vfloats["vtx_a0xy"]->at(vtx))<0.05);
				if(isCounter("nPassing_3body_vtxclean") && !vtxclean) continue;
				incrementCounter("nPassing_3body_vtxclean",wgt);
				_INF(vis,"");
				
				
				fillCategories(category,name,"tripletCategories_after_tripletTight");
				fillCategories(category,name,"tripletCategories_norm_after_tripletTight");
				_INF(vis,"");
				
				
				bool isolation = (vfloats["vtx_isolation003"]->at(vtx)<=0.08);
				if(isCounter("nPassing_3body_isolation") && !isolation) continue;
				incrementCounter("nPassing_3body_isolation",wgt);
				_INF(vis,"");
				
				
				fillCategories(category,name,"tripletCategories_after_isolation");
				fillCategories(category,name,"tripletCategories_norm_after_isolation");
				_INF(vis,"");
				
				

				if(isCounter("nPassing_met_MET") && !passMET(floats,metType,"tight")) continue;
				incrementCounter("nPassing_met_MET",wgt);
				_INF(vis,"");
				
				if(isCounter("nPassing_met_dphi3bodyMET") && !passDphi3bodyMET(vtx,vfloats,metType,"tight")) continue;
				incrementCounter("nPassing_met_dphi3bodyMET",wgt);
				_INF(vis,"");
				
				if(isCounter("nPassing_met_mT") && !passMT(vtx,vfloats,metType,"tight")) continue;
				incrementCounter("nPassing_met_mT",wgt);
				_INF(vis,"");
				
				if(isCounter("nPassing_met_dphiCalTrk") && !passDphiMETs(vtx,vfloats_decorations,"tight")) continue;
				incrementCounter("nPassing_met_dphiCalTrk",wgt);
				_INF(vis,"");
				
				
				
				if(isCounter("nPassing_ht_ht") && !passHT(vtx,vfloats_decorations,"tight")) continue;
				incrementCounter("nPassing_ht_ht",wgt);
				_INF(vis,"");
				
				if(isCounter("nPassing_ht_dphihtMET") && !passDphihtMET(vtx,vfloats_decorations,metType,"tight")) continue;
				incrementCounter("nPassing_ht_dphihtMET",wgt);
				_INF(vis,"");

				if(isCounter("nPassing_ht_mhT") && !passMHT(vtx,vfloats_decorations,"tight")) continue;
				incrementCounter("nPassing_ht_mhT",wgt);
				_INF(vis,"");
				
				if(isCounter("nPassing_ht_dr3body") && !passDrht3body(vtx,vfloats_decorations,"tight")) continue;
				incrementCounter("nPassing_ht_dr3body",wgt);
				_INF(vis,"");
				
				
				
				//// fill the simple CUTS out branches
				fillCUTSoutVecVars(vtx,vints,vfloats,vints_decorations,vfloats_decorations);
				_INF(vis,"");
			}
			
			
			
			
			//// everything should be in the sidebands
			bool inmasssidebands = (vfloats["vtx_mass"]->at(vtx)>mSideBandLeftLowerMeVGlob && vfloats["vtx_mass"]->at(vtx)<mSideBandRightUpperMeVGlob);
			if(isCounter("nPassing_3body_sidebands") && !inmasssidebands) continue;
			incrementCounter("nPassing_3body_sidebands",wgt);
			_INF(vis,"");
			
			
			//// everything should be blinded except for the MC'S
			// bool inblindregion = (!name.Contains("_3mu") && (vfloats["vtx_mass"]->at(vtx)>mBlindMinGlob && vfloats["vtx_mass"]->at(vtx)<mBlindMaxGlob));
			bool inblindregion = (name.Contains("Data") && (vfloats["vtx_mass"]->at(vtx)>mBlindMinGlob && vfloats["vtx_mass"]->at(vtx)<mBlindMaxGlob));
			if(isCounter("nPassing_3body_blinded") && inblindregion) continue;
			incrementCounter("nPassing_3body_blinded",wgt);
			_INF(vis,"");
			
			vector<float>* metcal = new vector<float>; metcal->push_back(floats["met_"+metType+"_et"]);
			vector<float>* mettrk = new vector<float>; mettrk->push_back(floats["met_track_et"]);
			TMapTSP2vi vi;
			TMapTSP2vf vf;
			vi.insert(make_pair("pass_loose", vints_decorations["pass_loose"]));
			vi.insert(make_pair("njets"     , vints_decorations["jets_n"]));
			vf.insert(make_pair("pt3body"   , vfloats["vtx_pt"]));
			vf.insert(make_pair("m3body"    , vfloats["vtx_mass"]));
			vf.insert(make_pair("mSS"       , vfloats["vtx_mSS"]));
			vf.insert(make_pair("mOS1"      , vfloats["vtx_mOS1"]));
			vf.insert(make_pair("mOS2"      , vfloats["vtx_mOS2"]));
			vf.insert(make_pair("score"     , vfloats_decorations["mva_score"]));
			vf.insert(make_pair("pval"      , vfloats["vtx_pval"]));
			vf.insert(make_pair("iso003"    , vfloats["vtx_isolation003"]));
			vf.insert(make_pair("iso020"    , vfloats["vtx_isolation020"]));
			vf.insert(make_pair("iso030"    , vfloats["vtx_isolation030"]));
			vf.insert(make_pair("trkspval"  , vfloats_decorations["trks_fitprob"]));
			vf.insert(make_pair("lxysig"    , vfloats_decorations["geo_lxySig"]));
			vf.insert(make_pair("a0xysig"   , vfloats_decorations["geo_a0xySig"]));
			vf.insert(make_pair("mettrk"    , mettrk));
			vf.insert(make_pair("metcal"    , metcal));
			vf.insert(make_pair("mttrk"     , vfloats["met_track_mT"]));
			vf.insert(make_pair("mtcal"     , vfloats["met_"+metType+"_mT"]));
			vf.insert(make_pair("dphihttrk" , vfloats_decorations["ht_dphimet_track"]));
			vf.insert(make_pair("dphihtcal" , vfloats_decorations["ht_dphimet_muons"]));
			TString chan = (name.EndsWith("_3mu")) ? "sig" : "bkg";
			bool passPostTraining = passPostBDTcut(vtx,vf,vi,mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob, minBDTcut,optBDTcut, chan, mSRminMeVGlob,mSRmaxMeVGlob,isBlinded);			
			delete metcal;
			delete mettrk;
			// bool passRes = passPostBDTcut(vfloats["vtx_mass"]->at(vtx),-1,vfloats["vtx_mOS1"]->at(vtx),vfloats["vtx_mOS2"]->at(vtx),vfloats["vtx_mSS"]->at(vtx), vints_decorations["jets_n"]->at(vtx), vfloats_decorations["ht_pt"]->at(vtx), vfloats_decorations["ht_dr3body"]->at(vtx),
			// 							  mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob,
			// 							  -1.,-1.,SigmaRhoOmegaPhiMeV,chan,"CR1",isBlinded);
			if(isCounter("nPassing_evt_rhoomegaphi") && !passPostTraining) continue;
			incrementCounter("nPassing_evt_rhoomegaphi",wgt);
			_INF(vis,"");
			
			
			fillCategories(category,name,"tripletCategories_after_lowmassres");
			fillCategories(category,name,"tripletCategories_norm_after_lowmassres");
			fillMETHistos(vtx,name,"_after_lowmassres",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"_after_lowmassres",vfloats,wgt);
			fillJetHistos(vtx,name,"_after_lowmassres",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillJetCalibrationHistos(vtx,name,"_after_lowmassres",floats,vfloats,metType,wgt);
			fillObjectHistos(vtx,name,"_after_lowmassres",vints,vfloats,wgt);
			_INF(vis,"");
			
			
			////////////////////////////////
			//// vertices survivnig all cuts
			tightvertices.push_back(vtx);
			_INF(vis,"");
		}

		
		//////////////////////////
		//// Fill the MVA out tree
		if(doMVAout) fillMVAoutTree(ints,floats);
		if(!doMVA)   fillCUTSoutTree(ints,floats);
		_INF(vis,"");
		
		
		
		//////////////////////////
		//// End of full selection
		unsigned int ntightvertices = tightvertices.size();
		if(isCounter("nPassing_evt_onecandidate") && ntightvertices!=1) continue;
		incrementCounter("nPassing_evt_onecandidate",wgt);
		_INF(vis,"");
		
		
		
		
		for(unsigned int v=0 ; v<ntightvertices ; ++v)
		{
			unsigned int vtx = tightvertices[v];
			int category = vints["vtx_code"]->at(vtx);
		
			/////////////////////
			//// fill decorations
			vints_decorations["pass_tight"]->at(vtx) = 1; // passed tight cuts ==> change flag from 0 to 1
			_INF(vis,"");
		
			/////////////////////
			//// fill some histos
			fillCategories(category,name,"tripletCategories_endOfSelection");
			fillCategories(category,name,"tripletCategories_norm_endOfSelection");
			fillMETHistos(vtx,name,"",floats,vfloats,metType,wgt);
			fill3BodyHistos(vtx,name,"",vfloats,wgt);
			fillJetHistos(vtx,name,"",vfloats,floats,minJetPtGeV*GeV2MeV,metType,wgt);
			fillObjectHistos(vtx,name,"",vints,vfloats,wgt);
			_INF(vis,"");
		}
		if(ntightvertices>0) fillTriggerHistos(name,"endOfSelection",vints,wgt);
		_INF(vis,"");
	
		////////////////////////
		//// fill the tight file
		ofile_tight->cd();
		otree_tight->Fill();
		if(i%10000==0 && i!=0)
		{
			otree_tight->FlushBaskets();
			// otree_tight->Write("", TObject::kOverwrite);
		}
		olddir->cd();
		
	}
	
	
	_INF(vis,"");
	
	_INF(1,"Processed "<<getCounter("nPassing_evt_all")<<" entries for tree: "<<name);
	// if(!doMVA) otree_loose->AutoSave();
	ofile_loose->cd();
	otree_loose->AutoSave();
	ofile_tight->cd();
	otree_tight->AutoSave();
	
	string scounters = printCounters("cutflow",name);
	fillCutFlowHisto(name,histos1);
	_INF(vis,"");
	
	if(isWsignal(name)) ofstr->close();
	if(isWsignal(name)) ofstr1->close();
	_INF(vis,"");
}

void NTUPanalysis(TString channel, TString master, TString method, TString jetMode, TString trkMode, TString metMode, TString metType, TString blinded, TString nDataTS="-1")
{	
	////////////////////
	//// input variables
	_INF(1,"channel = "<<channel);
	_INF(1,"master  = "<<master);
	_INF(1,"method  = "<<method);
	_INF(1,"jetMode = "<<jetMode);
	_INF(1,"trkMode = "<<trkMode);
	_INF(1,"metMode = "<<metMode);
	_INF(1,"metType = "<<metType);
	_INF(1,"blinded = "<<blinded);
	_INF(1,"nData   = "<<nDataTS);
	bool doMVA     = (method.Contains("mva"));
	bool doMVAout  = (method.Contains("mvaout"));
	bool isBlinded = (method.Contains("yes"));
	nData = stoi((string)nDataTS);
	
	////////////
	//// styles
	setStyle();
	_INF(vis,"");
	
	
	////////////////
	//// decorations
	vints_decorations.insert(make_pair("pass_loose",            new vector<int>));
	vints_decorations.insert(make_pair("pass_tight",            new vector<int>));
	vints_decorations.insert(make_pair("pass_tight_scat",       new vector<int>));
	vints_decorations.insert(make_pair("pass_tight_pbal",       new vector<int>));
	vints_decorations.insert(make_pair("pass_tight_fitq",       new vector<int>));
	
	vfloats_decorations.insert(make_pair("mva_score",           new vector<float>));
	vfloats_decorations.insert(make_pair("trks_fitprob",        new vector<float>));
	vfloats_decorations.insert(make_pair("trks_pixdEdx",        new vector<float>));
	vfloats_decorations.insert(make_pair("geo_lxySig",          new vector<float>));
	vfloats_decorations.insert(make_pair("geo_a0xySig",         new vector<float>));
	vfloats_decorations.insert(make_pair("muons_minpbalsig",    new vector<float>));
	vfloats_decorations.insert(make_pair("muons_maxpbalsig",    new vector<float>));
	vfloats_decorations.insert(make_pair("muons_minscatsig",    new vector<float>));
	vfloats_decorations.insert(make_pair("muons_maxscatsig",    new vector<float>));
	vfloats_decorations.insert(make_pair("muons_minnegbsig",    new vector<float>));
	vfloats_decorations.insert(make_pair("muons_maxnegbsig",    new vector<float>));
	
	for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
	{
		TString word = getTrkJetMETword(mode);
		_INF(vis,"mode="<<mode<<", word="<<word);
		
		TString trkmetword = word;
		TString jetmetword = word;
		if(mode>=ALL_JESUP  && mode<=ALL_JERDWN) trkmetword="";
		if(mode>=ALL_SOFTUP && mode<=ALL_JETDWN) jetmetword="";		
		TString jointword = (mode==ALL_UNCALIB) ? "_uncalib" : jetmetword+trkmetword;
		_INF(vis,"jetmetword="<<jetmetword<<", trkmetword="<<trkmetword);
		
		vints_decorations.insert(make_pair("jets_n"+jetmetword,                new vector<int>));
		vints_decorations.insert(make_pair("jets_b"+jetmetword,                new vector<int>));
		vfloats_decorations.insert(make_pair("ht_pt"+jetmetword,               new vector<float>));
		vfloats_decorations.insert(make_pair("ht_mT"+jetmetword,               new vector<float>));
		vfloats_decorations.insert(make_pair("ht_mT_mettrk"+trkmetword,        new vector<float>));
		vfloats_decorations.insert(make_pair("ht_metsig"+jetmetword,           new vector<float>));
		vfloats_decorations.insert(make_pair("ht_dphimet_reffinal"+jetmetword, new vector<float>));
		vfloats_decorations.insert(make_pair("ht_dphimet_muons"+jetmetword,    new vector<float>));
		vfloats_decorations.insert(make_pair("ht_dphimet_track"+jointword,     new vector<float>));
		vfloats_decorations.insert(make_pair("ht_dr3body"+jetmetword,          new vector<float>));
		vfloats_decorations.insert(make_pair("ht_dphi3body"+jetmetword,        new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dphi"+jointword,            new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dptrelavg"+jointword,       new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dptreltrk"+trkmetword,      new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dptrelcal"+jetmetword,      new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dhtreltrk"+trkmetword,      new vector<float>));
		vfloats_decorations.insert(make_pair("mets_dhtrelcal"+jetmetword,      new vector<float>));
	}
	_INF(vis,"");
	
	
	////////////////
	//// Input file
	TString ifbasenae = "flatout.periodall+MC.muons.cuts.analysis.n0.j0";
	TString ifname = (doMVAout) ? ifbasenae+".loose.root" : ifbasenae+".root";
	TFile *ifile = new TFile(ifname);
	_INF(vis,"");


	////////////////
	//// Output files 
	TString ofname_loose = ifname;
	if(!doMVAout) ofname_loose.ReplaceAll(".root",".loose.root");
	if(method=="mva")    ofname_loose.ReplaceAll("cuts","mva");
	if(method=="mvaout") ofname_loose.ReplaceAll("cuts","mvaout");
	if(channel!="all")   ofname_loose.ReplaceAll("periodall+MC",channel);
	// TFile *ofile_loose = NULL;
	// if(!doMVA) ofile_loose = new TFile(ofname,"recreate");
	TFile *ofile_loose = new TFile(ofname_loose,"recreate");
	TString ofname_tight = ofname_loose;
	ofname_tight.ReplaceAll("loose","tight");
	TFile *ofile_tight = new TFile(ofname_tight,"recreate");
	_INF(1,"ifname       = "<<ifname);
	_INF(1,"ofname_loose = "<<ofname_loose);
	_INF(1,"ofname_tight = "<<ofname_tight);
	_INF(vis,"");
	
	
	//////////////////////
	//// simple MVA output
	TString mvaoutName = "mvaout."+master+".root";
	if(channel!="all") mvaoutName = channel+"."+mvaoutName;
	if(doMVAout) initMVAout(mvaoutName);
	_INF(vis,"");
	
	//////////////////////
	//// simple CUTS output
	TString cutsoutName = "cutsout."+master+".root";
	if(channel!="all") cutsoutName = channel+"."+cutsoutName;
	if(!doMVA) initCUTSout(cutsoutName);
	_INF(vis,"");
	
	
	//////////////////////
	//// add the channels
	TString chnl = "";
	int counter = 0;
	chnl = "Wtaunu_3mu";      if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "Wtaunu_200k_3mu"; if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "bbTotau10_3mu";   if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "ccTotau10_3mu";   if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "bb_mu4mu4";       if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	// chnl = "bbTomu15";     if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "ccTomu15";        if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "JZxW";            if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "WmunuNpX";        if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "ZmumuNpX";        if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "WtaunuNpX";       if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "ZtautauNpX";      if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "bb_Jpsimu4mu4";   if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	chnl = "Data";            if(channel=="all" || channel==chnl) channels.insert(make_pair(increment(counter),chnl));
	_INF(vis,"");

	
	//////////////////////////
	//// decorate the channels
	// properties("Wtaunu_3mu",     "#it{W#rightarrow#nu#tau#rightarrow3#mu}#times"+sBRf,  kAzure-3,  1001, "f",   1); // 1001
	// properties("bbTotau10_3mu",  "#it{bb#rightarrow#tau10#rightarrow3#mu}#times"+sBRf,  kGreen+2,  1001, "f",   0);
	// properties("ccTotau10_3mu",  "#it{cc#rightarrow#tau10#rightarrow3#mu}#times"+sBRf,  kGreen-6,  1001, "f",   0);
	properties("Wtaunu_3mu",     "#it{W#rightarrow#nu#tau#rightarrow3#mu}",  kAzure-3,  1001, "f",   1); // 1001
	properties("Wtaunu_200k_3mu","#it{W#rightarrow#nu#tau#rightarrow3#mu}",  kAzure-3,  1001, "f",   1); // 1001
	properties("bbTotau10_3mu",  "#it{bb#rightarrow#tau10#rightarrow3#mu}",  kGreen+2,  1001, "f",   0);
	properties("ccTotau10_3mu",  "#it{cc#rightarrow#tau10#rightarrow3#mu}",  kGreen-6,  1001, "f",   0);
	properties("bb_Jpsimu4mu4",  "#it{J/#psi#rightarrow#mu4#mu4}",           kRed-9,    3015, "f",   0);
	properties("bbTomu15",       "#it{bb#rightarrow#mu15}",                  kRed-3,    3005, "f",   0);
	properties("bb_mu4mu4",      "#it{bb#rightarrow#mu4#mu4}",               kRed+3,    3005, "f",   1);
	properties("ccTomu15",       "#it{cc#rightarrow#mu15}",                  kGreen+3,  3004, "f",   1);
	properties("JZxW",           "JZxW (dijet)",                             kGray+2,   3022, "f",   0);
	properties("WmunuNpX",       "#it{W#rightarrow#mu#nu}",                  kViolet+4, 3144, "f",   0);
	properties("ZmumuNpX",       "#it{Z#rightarrow#mu#mu}",                  kViolet+3, 3144, "f",   0);
	properties("WtaunuNpX",      "#it{W#rightarrow#tau#nu}",                 kViolet+2, 3144, "f",   0);
	properties("ZtautauNpX",     "#it{Z#rightarrow#tau#tau}",                kViolet+1, 3144, "f",   0);
	properties("Data",           "Data",                                     kBlack,    3004, "ple", 1);
	_INF(vis,"");


	///////////////
	//// categories
	makeCategories();
	_INF(vis,"");
	
	
	/////////////
	//// triggers
	makeTriggers();
	_INF(vis,"");
	
	

	/////////////
	//// analysis
	for(TMapuiTS::iterator it=channels.begin() ; it!=channels.end() ; it++)
	{
		TString name = it->second;
		
		////////////////////
		//// book all histos
		clearCounters();
		initCounters(doMVA,doMVAout,isBlinded); // have to build the counters
		bookHistos(name);
		clearCounters(); // have to clear the counters
		_INF(vis,"");
		
		////////////////////////////////
		//// Loop and copy the data tree
		loop(ifile,ofile_loose,ofile_tight,name,jetMode,trkMode,metMode,metType,doMVA,doMVAout,isBlinded);
		_INF(vis,"");
	}
	_INF(vis,"");
	
	
	///////////////////
	//// Finalize trees
	delete ifile;
	// if(!doMVA) delete ofile_loose;
	olddir->cd();
	ofile_loose->cd();
	ofile_loose->Write();
	ofile_loose->Close();
	olddir->cd();
	ofile_tight->cd();
	ofile_tight->Write();
	ofile_tight->Close();
	// delete ofile_loose;
	// delete ofile_tight;
	_INF(1,"Written: "<<ofname_loose);
	_INF(1,"Written: "<<ofname_tight);
	_INF(vis,"");
	
	////////////////////
	//// Finalize histos
	olddir->cd();
	finalizeHistos(master,isBlinded,channel,doMVA);
	_INF(vis,"");

	//////////////////
	//// Write MVA out
	olddir->cd();
	if(doMVAout) finalizeMVAout();
	_INF(vis,"");
	
	//////////////////
	//// Write CUTS out
	olddir->cd();
	if(!doMVA) finalizeCUTSout();
	_INF(vis,"");

	/////////////////////////////
	//// exit if only one channel
	if(channel!="all") return;
	_INF(vis,"");
	
	/////////////////////
	//// exit if doMVAout
	if(doMVAout) return;
	_INF(vis,"");
	
	/////////////////
	//// Draw pdf/eps
	olddir->cd();
	finalizeFigures(master,doMVA);
	_INF(vis,"");
}
