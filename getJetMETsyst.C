/////////////////////////////////////////
//// root -b -l -q getJetMETsyst.C++ ////
/////////////////////////////////////////

#include "postBDTcuts.h"
#include "type.h"
#include "const.h"
#include "jets.h"

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}

TPaveText* ptxt;
void makeAtlasLabel()
{
	// ptxt = new TPaveText(0.62,0.70,0.87,0.87,"NDC");
	ptxt = new TPaveText(0.15,0.70,0.40,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

TLegend* leg;
TLegend* legR;
void makeLegend()
{
	leg = new TLegend(0.2,0.3,0.5,0.45,NULL,"brNDC");
	// leg = new TLegend(0.6,0.75,0.89,0.88,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	
	legR = new TLegend(0.6,0.74,0.89,0.87,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
}

float maxdistance(vector<float>& v, float reference, bool print=false)
{
	float maxd = -1.e20;
	for(unsigned int i=0 ; i<v.size() ; ++i)
	{
		float d = fabs(v[i]-reference);
		if(print) cout << "["<<i<<"] " << d << endl;
		maxd = (d>maxd) ? d : maxd;
	}
	if(print) cout << "max distance=" << maxd <<endl;
	return maxd;
}

float quadrature(vector<float>& v, float reference, bool print=false)
{
	float quad = 0;
	for(unsigned int i=0 ; i<v.size() ; ++i) quad += (v[i]-reference)*(v[i]-reference);
	quad = sqrt(quad);
	if(print) cout << "quadrature=" << quad <<endl;
	return quad;
}


void setLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
}


int readData(TTree* t, TString channel, float xFullMin,float xBlindMin,float xBlindMax,float xFullMax, TString mode)
{
	vector<int>*   loose           = 0;
	vector<int>*   tight           = 0;
	vector<int>*   pvntrk          = 0;
	vector<int>*   njets           = 0;
	vector<int>*   muonauthor      = 0;
	vector<float>* drmax           = 0;
	vector<float>* mass            = 0;
	vector<float>* pt              = 0;
	vector<float>* mSS             = 0;
	vector<float>* mOS1            = 0;
	vector<float>* mOS2            = 0;
	vector<float>* trkprob         = 0;
	vector<float>* iso003          = 0;
	vector<float>* iso010          = 0;
	vector<float>* iso014          = 0;
	vector<float>* iso020          = 0;
	vector<float>* iso030          = 0;
	vector<float>* maxpbalsig      = 0;
	vector<float>* pvalue          = 0;
	vector<float>* lxy             = 0;
	vector<float>* dlxy            = 0;
	vector<float>* a0xy            = 0;
	vector<float>* da0xy           = 0;
	float          calo_met        = 0;
	float          trk_met         = 0;
	vector<float>* calo_mt         = 0;
	vector<float>* trk_mt          = 0;
	vector<float>* calo_dphimet    = 0;
	vector<float>* trk_dphimet     = 0;
	vector<float>* dphimets        = 0;
	vector<float>* ht              = 0;
	vector<float>* calo_mht        = 0;
	vector<float>* trk_mht         = 0;
	vector<float>* calo_dphimet_ht = 0;
	vector<float>* trk_dphimet_ht  = 0;
	vector<float>* dr_ht           = 0;
	vector<float>* score           = 0;
	vector<float>* lxySig          = 0;
	vector<float>* a0xySig         = 0;
	vector<float>* dptreltrk       = 0;
	vector<float>* dptrelcal       = 0;
	
	t->SetBranchAddress("pass_loose",       &loose);
	t->SetBranchAddress("pass_tight",       &tight);
	t->SetBranchAddress("vtx_pvNtrk",       &pvntrk);
	t->SetBranchAddress("vtx_dRmax",        &drmax);
	t->SetBranchAddress("vtx_pt",           &pt);
	t->SetBranchAddress("vtx_mass",         &mass);
	t->SetBranchAddress("vtx_mSS",          &mSS);
	t->SetBranchAddress("vtx_mOS1",         &mOS1);
	t->SetBranchAddress("vtx_mOS2",         &mOS2);
	t->SetBranchAddress("trks_fitprob",     &trkprob);
	t->SetBranchAddress("vtx_isolation003", &iso003);
	t->SetBranchAddress("vtx_isolation010", &iso010);
	t->SetBranchAddress("vtx_isolation014", &iso014);
	t->SetBranchAddress("vtx_isolation020", &iso020);
	t->SetBranchAddress("vtx_isolation030", &iso030);
	t->SetBranchAddress("muons_maxpbalsig", &maxpbalsig);
	t->SetBranchAddress("vtx_pval",         &pvalue);
	t->SetBranchAddress("vtx_lxy",          &lxy);
	t->SetBranchAddress("vtx_lxyErr",       &dlxy);
	t->SetBranchAddress("geo_lxySig",       &lxySig);
	t->SetBranchAddress("vtx_a0xy",         &a0xy);
	t->SetBranchAddress("vtx_a0xyErr",      &da0xy);
	t->SetBranchAddress("geo_a0xySig",      &a0xySig);
	t->SetBranchAddress("met_muons_et",     &calo_met);
	t->SetBranchAddress("met_track_et",     &trk_met);
	t->SetBranchAddress("met_muons_mT",     &calo_mt);
	t->SetBranchAddress("met_track_mT",     &trk_mt);
	t->SetBranchAddress("met_muons_dPhi3mu",&calo_dphimet);
	t->SetBranchAddress("met_track_dPhi3mu",&trk_dphimet);
	t->SetBranchAddress("mets_dphi",        &dphimets);
	t->SetBranchAddress("mets_dptreltrk",   &dptreltrk);
	t->SetBranchAddress("mets_dptrelcal",   &dptrelcal);
	t->SetBranchAddress("ht_pt",            &ht);
	t->SetBranchAddress("ht_mT",            &calo_mht);
	t->SetBranchAddress("ht_mT_mettrk",     &trk_mht);
	t->SetBranchAddress("ht_dphimet_muons", &calo_dphimet_ht);
	t->SetBranchAddress("ht_dphimet_track", &trk_dphimet_ht);
	t->SetBranchAddress("ht_dr3body",       &dr_ht);
	t->SetBranchAddress("jets_n",           &njets);
	t->SetBranchAddress("vtx_code",         &muonauthor);
	t->SetBranchAddress("mva_score",        &score);

	vector<bool> rangeOSflags;
	int ncandidates = 0;
	for(Int_t entry=1 ; entry<=t->GetEntries() ; entry++)
	{
		t->GetEntry(entry);
		
		if(mass->size()!=1) continue;
		
		bool isFilled = false;
		
		TString type   = (channel.Contains("_3mu")) ? "sig" : "bkg";
		bool doBDT     = true;
		bool isBlinded = true;
		vector<float>* metcal = new vector<float>; metcal->push_back(calo_met);
		vector<float>* mettrk = new vector<float>; mettrk->push_back(trk_met);
		TMapTSP2vi vi;
		TMapTSP2vf vf;
		vi.insert(make_pair("njets"     , njets));
		vi.insert(make_pair("pass_loose", loose));
		vf.insert(make_pair("pt3body"   , pt));
		vf.insert(make_pair("m3body"    , mass));
		vf.insert(make_pair("mSS"       , mSS));
		vf.insert(make_pair("mOS1"      , mOS1));
		vf.insert(make_pair("mOS2"      , mOS2));
		vf.insert(make_pair("score"     , score));
		vf.insert(make_pair("pval"      , pvalue));
		vf.insert(make_pair("iso003"    , iso003));
		vf.insert(make_pair("iso020"    , iso020));
		vf.insert(make_pair("iso030"    , iso030));
		vf.insert(make_pair("trkspval"  , trkprob));
		vf.insert(make_pair("lxysig"    , lxySig));
		vf.insert(make_pair("a0xysig"   , a0xySig));
		vf.insert(make_pair("mettrk"    , mettrk));
		vf.insert(make_pair("metcal"    , metcal));
		vf.insert(make_pair("mttrk"     , trk_mt));
		vf.insert(make_pair("mtcal"     , calo_mt));
		vf.insert(make_pair("dphihttrk" , trk_dphimet_ht));
		vf.insert(make_pair("dphihtcal" , calo_dphimet_ht));
		bool pass = passPostBDTcut(0,vf,vi,xFullMin,xBlindMin,xBlindMax,xFullMax, minBDTcut,minBDTcut, type,mSRminMeVGlob,mSRmaxMeVGlob,isBlinded,doBDT,mode);
		delete metcal;
		delete mettrk;
		if(!pass) continue;
		
		ncandidates++;
	}
	return ncandidates;
}

Double_t getMaxErr(TH1* h)
{
	Double_t max = -1e20;
	for(Int_t b=1 ; b<=h->GetNbinsX() ; ++b) 
	{
		Double_t height = h->GetBinContent(b)+h->GetBinError(b);
		max = (height>max) ? height : max;
	}
	return max;
}

void plot(TMapTSP2TH1& histos1, TString name, TString title, unsigned int rangestart, unsigned int rangeend, bool doQuad=true)
{
	for(Int_t b=1 ; b<histos1[name]->GetNbinsX() ; ++b) 
	{	
		vector<float> values;
		for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
		{
			if(mode==ALL_UNCALIB) continue;
			if(mode==ALL_CALIB)   continue;
			TString word = getTrkJetMETword(mode);
			if(mode>=rangestart && mode<=rangeend) values.push_back(histos1[name+word]->GetBinContent(b));
		}
		float maxdistbin = maxdistance(values,histos1[name]->GetBinContent(b));
		float quadbin    = quadrature(values,histos1[name]->GetBinContent(b));
		float binerror = (doQuad) ? quadbin : maxdistbin;
		histos1[name]->SetBinError(b,binerror);
	}
	
	TH1* hLine;
	
	TCanvas* cnv = new TCanvas("cnv","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	histos1[name]->SetMinimum(0); histos1[name]->SetMaximum(getMaxErr(histos1[name])*1.35);
	histos1[name]->Draw("e2"); hLine=(TH1*)histos1[name]->Clone("hHT_line"); hLine->SetFillStyle(0); hLine->Draw("hist same");
	histos1[name+"_uncalib"]->Draw("same");
	ptxt->Draw("same");
	legR->Clear();
	legR->AddEntry(histos1[name+"_uncalib"],"Uncalibrated","F");
	legR->AddEntry(histos1[name],"Calibrated","fl");
	legR->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/JetMetTrk.syst.pdf");
	cnv->SaveAs("figures/JetMetTrk.syst."+title+".pdf");
	cnv->SaveAs("figures/JetMetTrk.syst."+title+".png");
	cnv->SaveAs("figures/JetMetTrk.syst."+title+".eps");
}


void getJetMETsyst()
{
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);
	// gStyle->SetPadTopMargin(0.13);
	// gStyle->SetPadRightMargin(0.15);
	// gStyle->SetPadBottomMargin(0.14);
	// gStyle->SetPadLeftMargin(0.12);
	Int_t font=42;
	Double_t tsize=0.04;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	gStyle->SetLabelSize(tsize,"x");
	gStyle->SetTitleSize(tsize,"x");
	gStyle->SetLabelSize(tsize,"y");
	gStyle->SetTitleSize(tsize,"y");
	gStyle->SetLabelSize(tsize,"z");
	gStyle->SetTitleSize(tsize,"z");
	gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatX(0);
	gStyle->SetStatY(0);
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0);
	gStyle->SetOptStat(0);
	gStyle->SetStatW(0);
	gStyle->SetStatH(0);
	gStyle->SetTitleX(0.55); //title X location 
	gStyle->SetTitleY(0.96); //title Y location 
	gStyle->SetTitleW(0.5); //title width 
	gStyle->SetTitleH(0.05); //title height
	gStyle->SetTitleBorderSize(0);
	
	makeAtlasLabel();
	makeLegend();
	
	//// calculate the signal efficiency
	// TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	// TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	
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
	
	TMapTSP2TH1 histos1;
	for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
	{
		TString word = getTrkJetMETword(mode);
		TString name = "";
		name = "hJet1"+word      ; histos1.insert(make_pair(name, new TH1F(name, ";Leading-jet #it{p}_{T} [MeV];Events",35,30000,100000)));                                        
		name = "hHT"+word        ; histos1.insert(make_pair(name, new TH1F(name, ";#it{H}_{T} [MeV];Events",50,0,100000)));                                        
		name = "hMETcal"+word    ; histos1.insert(make_pair(name, new TH1F(name, ";#it{E}_{T,cal}^{miss} [MeV];Events",45,10000,100000)));                     
		name = "hMETtrk"+word    ; histos1.insert(make_pair(name, new TH1F(name, ";#it{E}_{T,trk}^{miss} [MeV];Events",45,10000,100000)));                     
		name = "hMTcal"+word     ; histos1.insert(make_pair(name, new TH1F(name, ";#it{m}_{T}^{cal} [MeV];Events",65,20000,150000)));                           
		name = "hMTtrk"+word     ; histos1.insert(make_pair(name, new TH1F(name, ";#it{m}_{T}^{trk} [MeV];Events",60,30000,150000)));                           
		name = "hDphical"+word   ; histos1.insert(make_pair(name, new TH1F(name, ";#Delta#phi(3body,#it{E}_{T,cal}^{miss});Events",32,0,TMath::Pi())));       
		name = "hDphitrk"+word   ; histos1.insert(make_pair(name, new TH1F(name, ";#Delta#phi(3body,#it{E}_{T,trk}^{miss});Events",32,0,TMath::Pi())));       
		name = "hDphiHTcal"+word ; histos1.insert(make_pair(name, new TH1F(name, ";#Delta#phi(#it{H}_{T},#it{E}_{T,cal}^{miss});Events",32,0,TMath::Pi())));
		name = "hDphiHTtrk"+word ; histos1.insert(make_pair(name, new TH1F(name, ";#Delta#phi(#it{H}_{T},#it{E}_{T,trk}^{miss});Events",32,0,TMath::Pi())));
		name = "hDphicaltrk"+word; histos1.insert(make_pair(name, new TH1F(name, ";#Delta#it{#phi}_{trk}^{cal};Events",32,0,TMath::Pi())));                
		name = "hDpTcal"+word    ; histos1.insert(make_pair(name, new TH1F(name, ";p_{T}^{3body}/#it{E}_{T,cal}^{miss}-1;Events",60,-1,+5)));                  
		name = "hDpTtrk"+word    ; histos1.insert(make_pair(name, new TH1F(name, ";p_{T}^{3body}/#it{E}_{T,trk}^{miss}-1;Events",60,-1,+5)));                  
	}
	for(TMapTSP2TH1::iterator it=histos1.begin() ; it!=histos1.end() ; ++it)
	{
		if(it->first.Contains("_uncalib"))
		{
			it->second->SetLineColor(kBlack);
			it->second->SetLineWidth(2);
			it->second->SetLineStyle(2);
			continue;
		}
		it->second->SetLineColor(kBlack);
		it->second->SetLineWidth(2);
		it->second->SetFillColor(kGray+1);
	}
	
	cout << endl;
	
	float acceffNom = -1;
	vector<float> vacceff;
	for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
	{
		TString word = getTrkJetMETword(mode);
		TString name = getTrkJetMETname(mode);
		TString trkmetword = word;
		TString jetmetword = word;
		if(mode>=ALL_JESUP  && mode<=ALL_JERDWN) trkmetword="";
		if(mode>=ALL_SOFTUP && mode<=ALL_JETDWN) jetmetword="";
		
		TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"sig", sxSRmin,sxSRmax, false,true,"postTraining",ntuple,jetmetword,trkmetword);
		Float_t npassedS = tS->GetEntries(cuts_sig);
		Float_t ninitS   = (((TString)tS->GetName()).Contains("Wtaunu_200k_3mu")) ? 200000 : 99900;
		float AccEffSig  = npassedS/ninitS*100;
		cout << name << " signal Acc*Eff=" << AccEffSig << "\%" << " (npassedS=" << npassedS << ")" << endl;
		
		TString cuts_sig_noBDT = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", "","", false,false,"preTraining",ntuple,jetmetword,trkmetword);


		if(mode==ALL_CALIB)
		{
			tS->Draw("jet_pt1>>hJet1", cuts_sig_noBDT);
			tS->Draw("ht_pt>>hHT", cuts_sig_noBDT);
			tS->Draw("met_muons_et>>hMETcal", cuts_sig_noBDT);
			tS->Draw("met_muons_mT>>hMTcal", cuts_sig_noBDT);
			tS->Draw("met_muons_dPhi3mu>>hDphical", cuts_sig_noBDT);
			tS->Draw("ht_dphimet_muons>>hDphiHTcal", cuts_sig_noBDT);
			tS->Draw("mets_dptrelcal>>hDpTcal", cuts_sig_noBDT);
			
			tS->Draw("met_track_et>>hMETtrk", cuts_sig_noBDT);
			tS->Draw("met_track_mT>>hMTtrk", cuts_sig_noBDT);
			tS->Draw("met_track_dPhi3mu>>hDphitrk", cuts_sig_noBDT);
			tS->Draw("mets_dptreltrk>>hDpTtrk", cuts_sig_noBDT);
			
			tS->Draw("ht_dphimet_track>>hDphiHTtrk", cuts_sig_noBDT);
			tS->Draw("mets_dphi>>hDphicaltrk", cuts_sig_noBDT);
		}
		else if(mode==ALL_UNCALIB)
		{
			tS->Draw("jet_pt1_uncalib>>hJet1_uncalib", cuts_sig_noBDT);
			tS->Draw("ht_pt_uncalib>>hHT_uncalib", cuts_sig_noBDT);
			tS->Draw("met_muons_et_uncalib>>hMETcal_uncalib", cuts_sig_noBDT);
			tS->Draw("met_muons_mT_uncalib>>hMTcal_uncalib", cuts_sig_noBDT);
			tS->Draw("met_muons_dPhi3mu_uncalib>>hDphical_uncalib", cuts_sig_noBDT);
			tS->Draw("ht_dphimet_muons_uncalib>>hDphiHTcal_uncalib", cuts_sig_noBDT);
			tS->Draw("mets_dptrelcal_uncalib>>hDpTcal_uncalib", cuts_sig_noBDT);
			
			tS->Draw("met_track_et_uncalib>>hMETtrk_uncalib", cuts_sig_noBDT);
			tS->Draw("met_track_mT_uncalib>>hMTtrk_uncalib", cuts_sig_noBDT);
			tS->Draw("met_track_dPhi3mu_uncalib>>hDphitrk_uncalib", cuts_sig_noBDT);
			tS->Draw("mets_dptreltrk_uncalib>>hDpTtrk_uncalib", cuts_sig_noBDT);
			
			tS->Draw("ht_dphimet_track_uncalib>>hDphiHTtrk_uncalib", cuts_sig_noBDT);
			tS->Draw("mets_dphi_uncalib>>hDphicaltrk_uncalib", cuts_sig_noBDT);
		}
		else if(mode<=ALL_JERDWN)
		{
			tS->Draw("jet_pt1"+jetmetword+">>hJet1"+jetmetword, cuts_sig_noBDT);
			tS->Draw("ht_pt"+jetmetword+">>hHT"+jetmetword, cuts_sig_noBDT);
			tS->Draw("met_muons_et"+jetmetword+">>hMETcal"+jetmetword, cuts_sig_noBDT);
			tS->Draw("met_muons_mT"+jetmetword+">>hMTcal"+jetmetword, cuts_sig_noBDT);
			tS->Draw("met_muons_dPhi3mu"+jetmetword+">>hDphical"+jetmetword, cuts_sig_noBDT);
			tS->Draw("ht_dphimet_muons"+jetmetword+">>hDphiHTcal"+jetmetword, cuts_sig_noBDT);
			tS->Draw("mets_dptrelcal"+jetmetword+">>hDpTcal"+jetmetword, cuts_sig_noBDT);
			
			tS->Draw("ht_dphimet_track"+jetmetword+">>hDphiHTtrk"+jetmetword, cuts_sig_noBDT);
			tS->Draw("mets_dphi"+jetmetword+">>hDphicaltrk"+jetmetword, cuts_sig_noBDT);
		}
		else
		{
			tS->Draw("met_track_et"+trkmetword+">>hMETtrk"+trkmetword, cuts_sig_noBDT);
			tS->Draw("met_track_mT"+trkmetword+">>hMTtrk"+trkmetword, cuts_sig_noBDT);
			tS->Draw("met_track_dPhi3mu"+trkmetword+">>hDphitrk"+trkmetword, cuts_sig_noBDT);
			tS->Draw("mets_dptreltrk"+trkmetword+">>hDpTtrk"+trkmetword, cuts_sig_noBDT);
			
			tS->Draw("ht_dphimet_track"+trkmetword+">>hDphiHTtrk"+trkmetword, cuts_sig_noBDT);
			tS->Draw("mets_dphi"+trkmetword+">>hDphicaltrk"+trkmetword, cuts_sig_noBDT);
		}

		if     (mode==ALL_UNCALIB) continue;
		else if(mode==ALL_CALIB)   acceffNom = AccEffSig;
		else                       vacceff.push_back(AccEffSig);
	}
	float maxdist = maxdistance(vacceff,acceffNom);
	float quad = quadrature(vacceff,acceffNom);
	cout << "-------------------------------------------------" << endl;
	cout << "Uncertainty on Acc*Eff at x=x1 is quad:" << quad/acceffNom*100 << "\%, and maxdist: " << maxdist/acceffNom*100 << "\%"  << endl;
	cout << endl;
	
	TCanvas* cnv = new TCanvas("cnv","",600,400); cnv->Draw(); cnv->SaveAs("figures/JetMetTrk.syst.pdf(");
	
	
	plot(histos1,"hJet1","Jet1",           ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hHT","HT",               ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hMETcal","METcal",       ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hMTcal","MTcal",         ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hDphical","Dphical",     ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hDphiHTcal","DphiHTcal", ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	plot(histos1,"hDpTcal","DpTcal",       ALL_JESUP,ALL_JERDWN); // JES/JER-only syst range
	
	plot(histos1,"hMETtrk","METtrk",       ALL_SOFTUP,ALL_JETDWN); // Trk-only syst range
	plot(histos1,"hMTtrk","MTtrk",         ALL_SOFTUP,ALL_JETDWN); // Trk-only syst range
	plot(histos1,"hDphitrk","Dphitrk",     ALL_SOFTUP,ALL_JETDWN); // Trk-only syst range
	plot(histos1,"hDpTtrk","DpTtrk",       ALL_SOFTUP,ALL_JETDWN); // Trk-only syst range
	
	plot(histos1,"hDphiHTtrk","DphiHTtrk",    ALL_JESUP,ALL_JETDWN); // JES/JER+Trk syst range
	plot(histos1,"hDphicaltrk","hDphicaltrk", ALL_JESUP,ALL_JETDWN); // JES/JER+Trk syst range


	// cnv = new TCanvas("cnv","",600,400); cnv->Draw(); cnv->SaveAs("figures/JetMetTrk.syst.pdf)");
	// return;
	_INF(1,"Done kinematics plots");
	
	

	TString cuts_sig = "";
	TString cuts_bkg = "";
	Float_t ninitS   = -1;
	Float_t npassedS = -1;
	Float_t npassedD = -1;
	bool blinded = true; bool unblinded = false;
	bool doBDT   = true; bool noBDT     = false;
		
	// const Int_t    nlogxbins = 100;
	// const Double_t logxmin   = 0.2;
	// const Double_t logxmax   = 1.0;
	// Double_t logxbins[nlogxbins+1];
	// setLogBins(nlogxbins,logxmin,logxmax,logxbins);
	TH1F* hEff0 = new TH1F("hEff0",";BDT score cut (x_{1});Signal Acceptance#timesEfficiency [\%]",160,0.2,1.0); hEff0->SetLineColor(kBlack); hEff0->SetLineWidth(2); hEff0->SetLineStyle(2);
	TH1F* hEff1 = new TH1F("hEff1",";BDT score cut (x_{1});Signal Acceptance#timesEfficiency [\%]",160,0.2,1.0); hEff1->SetLineColor(kBlack); hEff1->SetLineWidth(2); hEff1->SetLineStyle(1); hEff1->SetFillColor(kGray+1);
	for(int b=1 ; b<=hEff0->GetNbinsX() ; ++b)
	{
		float x = hEff0->GetBinCenter(b);
		TString scurrentBDTcut = tstr(x,5);
		
		TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,scurrentBDTcut,"sig", sxSRmin,sxSRmax, false,true,"postTraining",ntuple);

		Float_t npassedS0 = tS0->GetEntries(cuts_sig);
		Float_t ninitS0   = 99900;
		float AccEffSig0  = npassedS0/ninitS0*100;
		
		Float_t npassedS1 = tS->GetEntries(cuts_sig);
		Float_t ninitS1   = 200000;
		float AccEffSig1  = npassedS1/ninitS1*100;
		
		vector<float> vacceff;
		for(unsigned int mode=ALL_UNCALIB ; mode<=ALL_CALIB ; ++mode)
		{
			if(mode==ALL_UNCALIB) continue;
			if(mode==ALL_CALIB)   continue;
			
			TString word = getTrkJetMETword(mode);
			TString name = getTrkJetMETname(mode);
			TString trkmetword = word;
			TString jetmetword = word;
			if(mode>=ALL_JESUP  && mode<=ALL_JERDWN) trkmetword="";
			if(mode>=ALL_SOFTUP && mode<=ALL_JETDWN) jetmetword="";

			TString cuts_sig_x = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,scurrentBDTcut,"sig", sxSRmin,sxSRmax, false,true,"postTraining",ntuple,jetmetword,trkmetword);
			npassedS = tS->GetEntries(cuts_sig_x);
			float AccEffSigx  = npassedS/200000*100;
			vacceff.push_back(AccEffSigx);
		}
		// float maxdist = maxdistance(vacceff,AccEffSig1);
		float quad = quadrature(vacceff,AccEffSig1);
		hEff0->SetBinContent(b,AccEffSig0);
		hEff1->SetBinContent(b,AccEffSig1);
		hEff1->SetBinError(b,quad);
	}
	cnv = new TCanvas("cnv","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	hEff1->SetMinimum(1);   hEff0->SetMinimum(1);
	hEff1->SetMaximum(3.5); hEff0->SetMaximum(3.5);
	hEff1->Draw("e2");  TH1* hLine=(TH1*)hEff1->Clone("hLine"); hLine->SetFillStyle(0); hLine->Draw("hist same");
	hEff0->Draw("hist same");
	ptxt->Draw("same");
	leg->AddEntry(hEff0,"Training sample","F");
	leg->AddEntry(hEff1,"Independent sample","fl");
	leg->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/JetMetTrk.syst.pdf");
	cnv->SaveAs("figures/Efficiencies.pdf");
	cnv->SaveAs("figures/Efficiencies.png");
	cnv->SaveAs("figures/Efficiencies.eps");
	_INF(1,"Done efficiency plots");
	
	
	cnv = new TCanvas("cnv","",600,400); cnv->Draw(); cnv->SaveAs("figures/JetMetTrk.syst.pdf)");	
	
	
	
	
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
	
	//// loose selection unblinded
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", sxSRmin,sxSRmax, unblinded,noBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkgsr", sxSRmin,sxSRmax, unblinded,noBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose only (SR data): " << npassedD << endl;
	
	//// loose+x>x0 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose+x>x0 only (SB data): " << npassedD << endl;
	
	//// loose+x>x0 selection unblinded
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"preTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkgsr", sxSRmin,sxSRmax,unblinded,doBDT,"preTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Loose+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Loose+x>x0 only (SR data): " << npassedD << endl;
	
	//// tight selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", sxSRmin,sxSRmax, unblinded,noBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkg", sxSRmin,sxSRmax, blinded,noBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight only (SB data): " << npassedD << endl;
	
	//// tight selection unblinded
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig", sxSRmin,sxSRmax, unblinded,noBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkgsr", sxSRmin,sxSRmax, unblinded,noBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight only (SR data): " << npassedD << endl;
	
	//// tight+x>x0 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x0 only (SB data): " << npassedD << endl;
	
	//// tight+x>x0 selection unblinded
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkgsr", sxSRmin,sxSRmax,unblinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x0 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x0 only (SR data): " << npassedD << endl;
	
	//// tight+x>x1 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"bkg", sxSRmin,sxSRmax, blinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x1 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x1 only (SB data): " << npassedD << endl;
	
	//// tight+x>x1 selection in the SR
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,soptBDTcut,"bkgsr", sxSRmin,sxSRmax, unblinded,doBDT,"postTraining",ntuple);
	npassedS = tS->GetEntries(cuts_sig); cout << "Tight+x>x1 only (signal) : " << npassedS << endl;
	npassedD = tD->GetEntries(cuts_bkg); cout << "Tight+x>x1 only (SR): " << npassedD << endl;
	
	cout << endl;

	npassedS = readData(tS,"Wtaunu_3mu", mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob,"preTraining");
	npassedD = readData(tD,"Data",       mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob,"preTraining");
	cout << "Loose+x>x0 only (signal) : " << npassedS << endl;
	cout << "Loose+x>x0 only (SB data): " << npassedD << endl;
	
	npassedS = readData(tS,"Wtaunu_3mu", mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob,"postTraining");
	npassedD = readData(tD,"Data",       mSideBandLeftLowerMeVGlob,mBlindMinGlob,mBlindMaxGlob,mSideBandRightUpperMeVGlob,"postTraining");
	cout << "Tight+x>x0 only (signal) : " << npassedS << endl;
	cout << "Tight+x>x0 only (SB data): " << npassedD << endl;
}
