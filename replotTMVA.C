////////////////////////////////////////
//// root -b -l -q replotTMVA.C++ //////
////////////////////////////////////////
#include "std.h"
#include "const.h"
#include "type.h"

int vis = 1;

enum trees
{
	TRAIN,TEST
};
enum channels
{
	SIG,BKG
};

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
	ptxt = new TPaveText(0.56,0.77,0.90,0.93,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} Internal");
	ptxt->AddText("#sqrt{s}=8 TeV, "+slumi);
}

void plotAtlasLabel()
{
	Double_t x = 0.60;
	Double_t y = 0.88;
	
	TLatex t1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
	t1.SetNDC();
	t1.SetTextFont(72);
	t1.SetTextColor(kBlack);
	double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
	// double dely = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
	t1.DrawLatex(x,y,"ATLAS");
	
	TLatex t2; 
	t2.SetNDC();
	t2.SetTextFont(42);
	t2.SetTextColor(kBlack);
	t2.DrawLatex(x+delx,y,"Internal");
	
	TLatex t3; 
	t3.SetNDC();
	t3.SetTextFont(42);
	t3.SetTextColor(kBlack);
	t3.DrawLatex(x,y-0.7*delx,"#sqrt{s}=8 TeV, "+slumi);
}

TLegend* legR;
TLegend* legL;
TLegend* legTL;
TLegend* legTM;
void makeLegend()
{	
	legR = new TLegend(0.56,0.57,0.91,0.775,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
	
	legTM = new TLegend(0.2,0.72,0.55,0.925,NULL,"brNDC");
	legTM->SetFillStyle(4000); //will be transparent
	legTM->SetFillColor(0);
	legTM->SetTextFont(42);
	legTM->SetBorderSize(0);
	
	legTL = new TLegend(0.15,0.78,0.45,0.93,NULL,"brNDC");
	legTL->SetFillStyle(4000); //will be transparent
	legTL->SetFillColor(0);
	legTL->SetTextFont(42);
	legTL->SetBorderSize(0);
	
	legL  = new TLegend(0.15,0.78,0.45,0.93,NULL,"brNDC");
	legL->SetFillStyle(4000); //will be transparent
	legL->SetFillColor(0);
	legL->SetTextFont(42);
	legL->SetBorderSize(0);
}

void setLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
}

void addHist(TMapTSP2TH1& histos, int id, int channel, TString name, TString titles, Int_t nbins, Double_t* xpoints)
{
	if(nbins<=0) _FAT("nbins="<<nbins<<" for histogram: "<<name);
	
	TString channelname = (channel==SIG) ? "signal" : "background";
	TString idname      = (id==TRAIN)    ? "train"  : "test";
	TString hname       = name+"_"+channelname+"_"+idname;
	
	histos.insert(make_pair(hname, new TH1F(hname,titles,nbins,xpoints)));
	histos.insert(make_pair(hname+"_frame", new TH1F(hname+"_frame",titles,nbins,xpoints)));
	histos[hname]->Sumw2();
	
	if(channel==BKG)
	{
		if(id==TRAIN)
		{
			histos[hname]->SetMarkerStyle(20);
			histos[hname]->SetMarkerSize(1.2);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetMarkerColor(kBlack);
			histos[hname]->SetLineColor(kBlack);
		}
		else
		{
			histos[hname]->SetLineColor(kBlack);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetFillColor(kBlack);
			histos[hname]->SetFillStyle(3354);
		}
	}
	else
	{
		if(id==TRAIN)
		{
			histos[hname]->SetMarkerStyle(24);
			histos[hname]->SetMarkerSize(1.2);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetMarkerColor(kBlack);
			histos[hname]->SetLineColor(kBlack);
		}
		else
		{
			histos[hname]->SetFillColor(kGray);
			histos[hname]->SetLineColor(kGray);
		}
	}
	
	_DBG(vis,"Added histo: "<<hname);
}


void addHist(TMapTSP2TH1& histos, int id, int channel, TString name, TString titles, int nbins, float xmin, float xmax)
{
	if(nbins<=0) _FAT("nbins="<<nbins<<" for histogram: "<<name);
	
	TString channelname = (channel==SIG) ? "signal" : "background";
	TString idname      = (id==TRAIN)    ? "train"  : "test";
	TString hname       = name+"_"+channelname+"_"+idname;
	
	histos.insert(make_pair(hname, new TH1F(hname,titles,nbins,xmin,xmax)));
	histos.insert(make_pair(hname+"_frame", new TH1F(hname+"_frame",titles,nbins,xmin,xmax)));
	histos[hname]->Sumw2();
	
	if(channel==BKG)
	{
		if(id==TRAIN)
		{
			histos[hname]->SetMarkerStyle(20);
			histos[hname]->SetMarkerSize(1.2);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetMarkerColor(kBlack);
			histos[hname]->SetLineColor(kBlack);
		}
		else
		{
			histos[hname]->SetLineColor(kBlack);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetFillColor(kBlack);
			histos[hname]->SetFillStyle(3354);
		}
	}
	else
	{
		if(id==TRAIN)
		{
			histos[hname]->SetMarkerStyle(24);
			histos[hname]->SetMarkerSize(1.2);
			histos[hname]->SetLineWidth(1);
			histos[hname]->SetMarkerColor(kBlack);
			histos[hname]->SetLineColor(kBlack);
		}
		else
		{
			histos[hname]->SetFillColor(kGray);
			histos[hname]->SetLineColor(kGray);
		}
	}
	
	_DBG(vis,"Added histo: "<<hname);
}

float Sum(TH1* h, bool addUunderFlow=false, bool addOverFlow=false)
{
	float I=0.;
	for(int i=1 ; i<=h->GetNbinsX() ; i++)// without under- and over-flow
	{
		I += h->GetBinContent(i);
	}
	if(addUunderFlow) I+=h->GetBinContent(0);
	if(addOverFlow)   I+=h->GetBinContent(h->GetNbinsX()+1);
	return I;
}

void plot(TString prefix, TString var, TMapTSP2TH1& histos1, TLegend* leg, bool doLogy=false, bool doLogx=false)
{
	TCanvas* cnv = new TCanvas("cnv","",800,600);
	gPad->SetTicks(1,1);
	if(doLogy) gPad->SetLogy();
	if(doLogx) gPad->SetLogx();
	
	float dx = histos1[var+"_background_test"]->GetBinWidth(1);
	
	float scaleBkgTest  = 1./dx/Sum(histos1[var+"_background_test"]);
	float scaleBkgTrain = 1./dx/Sum(histos1[var+"_background_train"]);
	float scaleSigTest  = 1./dx/Sum(histos1[var+"_signal_test"]);
	float scaleSigTrain = 1./dx/Sum(histos1[var+"_signal_train"]);
	
	histos1[var+"_background_test"]->Scale(scaleBkgTest);
	histos1[var+"_background_train"]->Scale(scaleBkgTrain);
	histos1[var+"_signal_test"]->Scale(scaleSigTest);
	histos1[var+"_signal_train"]->Scale(scaleSigTrain);
	
	histos1[var+"_background_test"]->GetYaxis()->SetTitle("(1/N)dN/dx");
	histos1[var+"_background_train"]->GetYaxis()->SetTitle("(1/N)dN/dx");
	histos1[var+"_signal_test"]->GetYaxis()->SetTitle("(1/N)dN/dx");
	histos1[var+"_signal_train"]->GetYaxis()->SetTitle("(1/N)dN/dx");
	
	histos1[var+"_signal_train_frame"]->GetYaxis()->SetTitle("(1/N)dN/dx");
	
	Int_t maxbinSig = -1;
	Int_t maxbinDat = -1;
	float maxSig = -1;
	float maxDat = -1;
	
	maxbinSig = histos1[var+"_signal_train"]->GetMaximumBin();
	maxbinDat = histos1[var+"_background_train"]->GetMaximumBin();
	maxSig = histos1[var+"_signal_train"]->GetBinContent(maxbinSig)+histos1[var+"_signal_train"]->GetBinError(maxbinSig);
	maxDat = histos1[var+"_background_train"]->GetBinContent(maxbinDat)+histos1[var+"_background_train"]->GetBinError(maxbinDat);
	float maxTrain = (maxSig>maxDat) ? maxSig : maxDat;
	
	maxbinSig = histos1[var+"_signal_test"]->GetMaximumBin();
	maxbinDat = histos1[var+"_background_test"]->GetMaximumBin();
	maxSig = histos1[var+"_signal_test"]->GetBinContent(maxbinSig)+histos1[var+"_signal_test"]->GetBinError(maxbinSig);
	maxDat = histos1[var+"_background_test"]->GetBinContent(maxbinDat)+histos1[var+"_background_test"]->GetBinError(maxbinDat);
	float maxTest = (maxSig>maxDat) ? maxSig : maxDat;
	
	float max = (maxTest>maxTrain) ? maxTest : maxTrain;
	histos1[var+"_signal_train_frame"]->SetMaximum((doLogy) ? max*20 : max*1.1);
	histos1[var+"_signal_train_frame"]->SetMinimum((doLogy) ? 5.e-3  : 0);
	
	histos1[var+"_signal_train_frame"]->Draw();
	histos1[var+"_signal_test"]->Draw("hist same");
	histos1[var+"_background_test"]->Draw("hist same");
	histos1[var+"_signal_train"]->Draw("p1x1 same");
	histos1[var+"_background_train"]->Draw("p1x1 same");	
	
	plotAtlasLabel();
	
	TH1F* hTmp = new TH1F("hTmp","",1,0,1);
	hTmp->SetMarkerStyle(4);
	hTmp->SetMarkerSize(1.2);
	hTmp->SetLineWidth(1);
	hTmp->SetMarkerColor(kBlack);
	hTmp->SetLineColor(kBlack);
	
	leg->Clear();
	leg->AddEntry(histos1[var+"_background_train"],"Background (training sample)","elp");
	// leg->AddEntry(histos1[var+"_signal_train"],"Signal (training sample)","elp");
	leg->AddEntry(hTmp,"Signal (training sample)","elp");
	leg->AddEntry(histos1[var+"_background_test"],"Background (test sample)","f");
	leg->AddEntry(histos1[var+"_signal_test"],"Signal (test sample)","f");
	leg->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	
	TString suffix = (doLogy) ? ".logy" : "";
	
	cnv->SaveAs("figures/"+prefix+".pdf");
	cnv->SaveAs("figures/"+prefix+"."+var+suffix+".png");
	cnv->SaveAs("figures/"+prefix+"."+var+suffix+".eps");
	cnv->SaveAs("figures/"+prefix+"."+var+suffix+".pdf");
	
	delete hTmp;
}

void plot8(TString prefix, TString name, vector<TString>& vars, vector<TLegend*>& legs, TMapTSP2TH1& histos1, bool doLogy=false)
{
	TCanvas* cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,4);
	
	for(unsigned int pad=0 ; pad<vars.size() ; ++pad)
	{
		TString var = vars[pad];
		cnv->cd(pad+1);
		gPad->SetTicks(1,1);
		if(doLogy) gPad->SetLogy();
		
		Int_t maxbinSig = -1;
		Int_t maxbinDat = -1;
		float maxSig = -1;
		float maxDat = -1;
	
		maxbinSig = histos1[var+"_signal_train"]->GetMaximumBin();
		maxbinDat = histos1[var+"_background_train"]->GetMaximumBin();
		maxSig = histos1[var+"_signal_train"]->GetBinContent(maxbinSig)+histos1[var+"_signal_train"]->GetBinError(maxbinSig);
		maxDat = histos1[var+"_background_train"]->GetBinContent(maxbinDat)+histos1[var+"_background_train"]->GetBinError(maxbinDat);
		float maxTrain = (maxSig>maxDat) ? maxSig : maxDat;

		maxbinSig = histos1[var+"_signal_test"]->GetMaximumBin();
		maxbinDat = histos1[var+"_background_test"]->GetMaximumBin();
		maxSig = histos1[var+"_signal_test"]->GetBinContent(maxbinSig)+histos1[var+"_signal_test"]->GetBinError(maxbinSig);
		maxDat = histos1[var+"_background_test"]->GetBinContent(maxbinDat)+histos1[var+"_background_test"]->GetBinError(maxbinDat);
		float maxTest = (maxSig>maxDat) ? maxSig : maxDat;

		float max = (maxTest>maxTrain) ? maxTest : maxTrain;
		histos1[var+"_signal_train_frame"]->SetMaximum(max*1.1);

		histos1[var+"_signal_train_frame"]->Draw();
		histos1[var+"_signal_test"]->Draw("hist same");
		histos1[var+"_background_test"]->Draw("hist same");
		histos1[var+"_signal_train"]->Draw("p1x1 same");
		histos1[var+"_background_train"]->Draw("p1x1 same");	

		// plotAtlasLabel();
		ptxt->Draw("same");
		
		legs[pad]->Clear();
		legs[pad]->AddEntry(histos1[var+"_background_train"],"Background (training sample)","ple");
		legs[pad]->AddEntry(histos1[var+"_signal_train"],"Signal (training sample)","ple");
		legs[pad]->AddEntry(histos1[var+"_background_test"],"Background (test sample)","f");
		legs[pad]->AddEntry(histos1[var+"_signal_test"],"Signal (test sample)","f");
		legs[pad]->Draw("same");
		
		gPad->Update();
		gPad->RedrawAxis();
	}
	
	cnv->Update();
	cnv->SaveAs("figures/"+prefix+"."+name+".png");
	cnv->SaveAs("figures/"+prefix+"."+name+".eps");
	cnv->SaveAs("figures/"+prefix+"."+name+".pdf");
	cnv->SaveAs("figures/"+prefix+".pdf");
	delete cnv;
}


void replotTMVA()
{
	// use plain black on white colors
	Int_t icol=0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects
	
	// set the paper & margin sizes
	gStyle->SetPaperSize(20,26);
	
	// set margin sizes
	// gStyle->SetPadTopMargin(0.03);
	// gStyle->SetPadRightMargin(0.07);
	// gStyle->SetPadBottomMargin(0.125);
	// gStyle->SetPadLeftMargin(0.115);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.12);
	
	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.05);
	gStyle->SetTitleYOffset(0.95);
	
	// use large fonts
	//Int_t font=72; // Helvetica italics
	Int_t font=42; // Helvetica
	Double_t tsize=0.05;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	
	gStyle->SetLabelSize(tsize*0.85,"x");
	gStyle->SetTitleSize(tsize*1.10,"x");
	gStyle->SetLabelSize(tsize*0.85,"y");
	gStyle->SetTitleSize(tsize*1.10,"y");
	gStyle->SetLabelSize(tsize*0.85,"z");
	gStyle->SetTitleSize(tsize*1.10,"z");
	
	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2.);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
	// get rid of X error bars 
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);
	
	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
	
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	
	
	
	makeAtlasLabel();
	makeLegend();
	
	_DBG(vis,"");
	
	
	int      nPbins = 25;
	Double_t Pbins[nPbins+1];
	setLogBins(nPbins,1e-9,1.,Pbins);
	
	TMapTSP2TH1 histos1;
	for(unsigned int j=SIG ; j<=BKG ; ++j)
	{
		for(unsigned int i=TRAIN ; i<=TEST ; ++i)
		{
			addHist(histos1,i,j,"score",         ";BDT score;Events",40,-1,+1);
			addHist(histos1,i,j,"PVNtrk",        ";#it{N}_{trk}^{PV};Events",50,0,200);
			addHist(histos1,i,j,"pT3body",       ";#it{p}_{T}^{3body} [GeV];Events", 40,0.,100.);
			addHist(histos1,i,j,"isolation020",  ";#Sigma#it{p}_{T}^{trk}(cone #Delta#it{R}_{max}+0.20)/#it{p}_{T}^{3body};Events", 60,0.,0.3);
			addHist(histos1,i,j,"trksfitprob",   ";#it{P}_{trks};Events",50,0.,1.);
			// addHist(histos1,i,j,"trksfitprob",   ";#it{P}_{trks};Events",nPbins,Pbins);
			addHist(histos1,i,j,"pvalue",        ";#it{p}-value (three-body vertex);Events", 50,0.,1.);				
			addHist(histos1,i,j,"SLxy",          ";#it{S}(#it{L}_{xy});Events", 60,-10.,+50.);
			addHist(histos1,i,j,"Sa0xy",         ";#it{S}(#it{a}_{0}^{xy});Events", 50,0.,25.);
			addHist(histos1,i,j,"calo_met",      ";#it{E}_{T,cal}^{miss} [GeV];Events",25,10.,110.);
			addHist(histos1,i,j,"calo_mt",       ";#it{m}_{T}^{cal} [GeV];Events",40,20.,140.);
			addHist(histos1,i,j,"calo_dphi3mu",  ";#Delta#it{#phi}_{3body}^{cal};Events",32,0.,TMath::Pi());
			addHist(histos1,i,j,"trk_met",       ";#it{E}_{T,trk}^{miss} [GeV];Events",25,10.,110.);
			addHist(histos1,i,j,"trk_mt",        ";#it{m}_{T}^{trk} [GeV];Events",40,20.,140.);
			addHist(histos1,i,j,"calo_trk_dphi", ";#Delta#it{#phi}_{trk}^{cal};Events",32,0.,TMath::Pi());
			addHist(histos1,i,j,"dptreltrk",     ";p_{T}^{3body}/#it{E}_{T,trk}^{miss}-1;Events",60,-1.,+5.);
			addHist(histos1,i,j,"dptrelcal",     ";p_{T}^{3body}/#it{E}_{T,cal}^{miss}-1;Events",60,-1.,+5.);
			addHist(histos1,i,j,"ht",            ";#it{#Sigma}_{T} [GeV];Events",40,0.,100.);
		}
	}
	
	_DBG(vis,"");
	
	
	
	TFile* f = new TFile("TMVA.root","READ");
	for(unsigned int tree=TRAIN ; tree<=TEST ; ++tree)
	{
		TTree* t = (tree==TRAIN) ? (TTree*)f->Get("TrainTree") : (TTree*)f->Get("TestTree");
		_DBG(vis,"Tree is: "<<tree);
		
		Int_t           classID;
		Float_t         vtx_isolation020;
		Float_t         vtx_pval;
		Float_t         trks_fitprob;
		Float_t         vtx_pt;
		Float_t         ht_pt;
		Float_t         mets_dphi;
		Float_t         met_muons_et;
		Float_t         met_muons_mT;
		Float_t         mets_dptrelcal;
		Float_t         met_track_et;
		Float_t         met_track_mT;
		Float_t         mets_dptreltrk;
		Float_t         met_muons_dPhi3mu;
		Float_t         geo_lxySig;
		Float_t         geo_a0xySig;
		Float_t         vtx_pvNtrk;
		Float_t         weight;
		Float_t         BDTG;
		
		t->SetBranchAddress("classID", &classID);
		t->SetBranchAddress("vtx_isolation020", &vtx_isolation020);
		t->SetBranchAddress("vtx_pval", &vtx_pval);
		t->SetBranchAddress("trks_fitprob", &trks_fitprob);
		t->SetBranchAddress("vtx_pt", &vtx_pt);
		t->SetBranchAddress("ht_pt", &ht_pt);
		t->SetBranchAddress("mets_dphi", &mets_dphi);
		t->SetBranchAddress("met_muons_et", &met_muons_et);
		t->SetBranchAddress("met_muons_mT", &met_muons_mT);
		t->SetBranchAddress("mets_dptrelcal", &mets_dptrelcal);
		t->SetBranchAddress("met_track_et", &met_track_et);
		t->SetBranchAddress("met_track_mT", &met_track_mT);
		t->SetBranchAddress("mets_dptreltrk", &mets_dptreltrk);
		t->SetBranchAddress("met_muons_dPhi3mu", &met_muons_dPhi3mu);
		t->SetBranchAddress("geo_lxySig", &geo_lxySig);
		t->SetBranchAddress("geo_a0xySig", &geo_a0xySig);
		t->SetBranchAddress("vtx_pvNtrk", &vtx_pvNtrk);
		t->SetBranchAddress("weight", &weight);
		t->SetBranchAddress("BDTG", &BDTG);
		
		Long64_t nentries = t->GetEntriesFast();
		for (Long64_t jentry=0 ; jentry<nentries ; jentry++)
		{
			t->GetEntry(jentry);
			TString channel = (classID==SIG) ? "_signal" : "_background";
			TString id      = (tree==TRAIN)  ? "_train"  : "_test";
			
			//// fill histos
			histos1["score"+channel+id]->Fill(BDTG,weight);
			histos1["PVNtrk"+channel+id]->Fill(vtx_pvNtrk,weight);
			histos1["pT3body"+channel+id]->Fill(vtx_pt*MeV2GeV,weight);
			histos1["isolation020"+channel+id]->Fill(vtx_isolation020,weight);
			histos1["trksfitprob"+channel+id]->Fill(trks_fitprob,weight);
			histos1["pvalue"+channel+id]->Fill(vtx_pval,weight);		
			histos1["SLxy"+channel+id]->Fill(geo_lxySig,weight);
			histos1["Sa0xy"+channel+id]->Fill(geo_a0xySig,weight);
			histos1["calo_met"+channel+id]->Fill(met_muons_et*MeV2GeV,weight);
			histos1["calo_mt"+channel+id]->Fill(met_muons_mT*MeV2GeV,weight);
			histos1["calo_dphi3mu"+channel+id]->Fill(met_muons_dPhi3mu,weight);
			histos1["trk_met"+channel+id]->Fill(met_track_et*MeV2GeV,weight);
			histos1["trk_mt"+channel+id]->Fill(met_track_mT*MeV2GeV,weight);
			histos1["calo_trk_dphi"+channel+id]->Fill(mets_dphi,weight);
			histos1["dptreltrk"+channel+id]->Fill(mets_dptreltrk,weight);
			histos1["dptrelcal"+channel+id]->Fill(mets_dptrelcal,weight);
			histos1["ht"+channel+id]->Fill(ht_pt*MeV2GeV,weight);
		}
   }

	_DBG(vis,"");

	TCanvas* cnv = new TCanvas("cnv","",800,600);
	TString pdffilename = "tmvaplots";
	cnv->SaveAs("figures/"+pdffilename+".pdf(");

	plot(pdffilename,"score",         histos1,legTM);
	plot(pdffilename,"calo_mt",       histos1,legR);
	plot(pdffilename,"trk_met",       histos1,legR);
	plot(pdffilename,"isolation020",  histos1,legR);
	plot(pdffilename,"ht",            histos1,legR);
	plot(pdffilename,"trk_mt",        histos1,legR);
	plot(pdffilename,"calo_trk_dphi", histos1,legR);
	plot(pdffilename,"calo_met",      histos1,legR);
	plot(pdffilename,"dptreltrk",     histos1,legR);
	plot(pdffilename,"calo_dphi3mu",  histos1,legTM);
	plot(pdffilename,"pvalue",        histos1,legR);
	plot(pdffilename,"Sa0xy",         histos1,legR);
	// plot(pdffilename,"trksfitprob",   histos1,legR,false,true);
	plot(pdffilename,"trksfitprob",   histos1,legR);
	plot(pdffilename,"pT3body",       histos1,legR);
	plot(pdffilename,"PVNtrk",        histos1,legR);
	plot(pdffilename,"SLxy",          histos1,legR);
	plot(pdffilename,"dptrelcal",     histos1,legR);

	plot(pdffilename,"score",         histos1,legTM,true);
	plot(pdffilename,"calo_mt",       histos1,legR,true);
	plot(pdffilename,"trk_met",       histos1,legR,true);
	plot(pdffilename,"isolation020",  histos1,legR,true);
	plot(pdffilename,"ht",            histos1,legR,true);
	plot(pdffilename,"trk_mt",        histos1,legR,true);
	plot(pdffilename,"calo_trk_dphi", histos1,legR,true);
	plot(pdffilename,"calo_met",      histos1,legR,true);
	plot(pdffilename,"dptreltrk",     histos1,legR,true);
	plot(pdffilename,"calo_dphi3mu",  histos1,legTM,true);
	plot(pdffilename,"pvalue",        histos1,legR,true);
	plot(pdffilename,"Sa0xy",         histos1,legR,true);
	// plot(pdffilename,"trksfitprob",   histos1,legR,true,true);
	plot(pdffilename,"trksfitprob",   histos1,legR,true);
	plot(pdffilename,"pT3body",       histos1,legR,true);
	plot(pdffilename,"PVNtrk",        histos1,legR,true);
	plot(pdffilename,"SLxy",          histos1,legR,true);
	plot(pdffilename,"dptrelcal",     histos1,legR,true);
	
	_DBG(vis,"");
	
	
	vector<TString>  vars;
	vector<TLegend*> legs;
	
	vars.clear(); legs.clear();
	vars.push_back("calo_mt");      vars.push_back("trk_met");
	vars.push_back("isolation020"); vars.push_back("ht");
	vars.push_back("trk_mt");       vars.push_back("calo_trk_dphi");
	vars.push_back("calo_met");     vars.push_back("dptreltrk");
	legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR);
	plot8(pdffilename,"TrainingBDTinputs01to08",vars,legs,histos1);
	
	vars.clear(); legs.clear();
	vars.push_back("calo_dphi3mu"); vars.push_back("pvalue");
	vars.push_back("Sa0xy");        vars.push_back("trksfitprob");
	vars.push_back("pT3body");      vars.push_back("PVNtrk");
	vars.push_back("SLxy");         vars.push_back("dptrelcal");
	legs.push_back(legTM); legs.push_back(legR);
	legs.push_back(legR);  legs.push_back(legR);
	legs.push_back(legR);  legs.push_back(legR);
	legs.push_back(legR);  legs.push_back(legR);
	plot8(pdffilename,"TrainingBDTinputs09to16",vars,legs,histos1);
	
	
	cnv = new TCanvas("cnv","",800,600);
	cnv->SaveAs("figures/"+pdffilename+".pdf)");
	
	_DBG(vis,"Done");
}
