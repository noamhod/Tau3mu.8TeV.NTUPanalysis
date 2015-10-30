#include "std.h"
#include "type.h"
#include "postBDTcuts.h"
#include "const.h"

void rgbPalette(Double_t r, Double_t g, Double_t b, Int_t nb=50)
{
	const UInt_t Number = 3;
	Double_t Red[Number]    = { r,   0.0, 0.0 };
	Double_t Green[Number]  = { g,   0.0, 0.0 };
	Double_t Blue[Number]   = { b,   0.0, 0.0 };
	Double_t Length[Number] = { 0.1, 0.5, 1.0 };
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}
void gryPalette(Int_t nb=50)
{
	const UInt_t Number = 2;
	Double_t Red[Number]    = { 1.0, 0.0 };
	Double_t Green[Number]  = { 1.0, 0.0 };
	Double_t Blue[Number]   = { 1.0, 0.0 };
	Double_t Length[Number] = { 0.0, 1.0,};
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}
TExec* exeRed   = new TExec("exeRed",   "rgbPalette(0.7,0,0,50);");
TExec* exeGreen = new TExec("exeGreen", "rgbPalette(0,0.5,0,50);");
TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0,0.8,10);");
TExec* exeGray  = new TExec("exeGray",  "rgbPalette(1,1,1,50);");

void Dalitz()
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
	gStyle->SetPadRightMargin(0.16);
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
	// gStyle->SetMarkerStyle(20);
	// gStyle->SetMarkerSize(1.2);
	// gStyle->SetHistLineWidth(2.);
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
	
	TFile* fPass = new TFile("mvaout.muons.root","READ");
	TTree* tPass = (TTree*)fPass->Get("fltmva_Wtaunu_200k_3mu");
	
	Int_t nbins = 350000+1;
	Double_t xmin = 0-0.5;
	Double_t xmax = 350000+0.5;
	TH1D* hEvents = new TH1D("hEvents","",nbins,xmin,xmax);
	TString cuts_sig = postBDTcut("1450","1690","1870","2110","-0.9","0.933","sig","1713","1841",false,true,"postTraining","mvaout");
	Double_t npassedS = tPass->Draw("EventNumber>>hEvents",cuts_sig);
	cout << "Tight+x>x1: " << npassedS << endl;
	
	vector<float>* mOS1 = 0;
	vector<float>* mOS2 = 0;
	vector<float>* mSS  = 0;
	vector<float>* mass = 0;
	int            evt;
	TFile* fTruth = new TFile("fM2body.Wtaunu_200k_3mu.root","READ");
	TTree* tTruth = (TTree*)fTruth->Get("tM2body");
	tTruth->SetBranchAddress("mOS1", &mOS1);
	tTruth->SetBranchAddress("mOS2", &mOS2);
	tTruth->SetBranchAddress("mSS",  &mSS);
	tTruth->SetBranchAddress("mass", &mass);
	tTruth->SetBranchAddress("evt",  &evt);
	
	TH2D* hAllDalitz0 = new TH2D("Dalitz0.all", ";#it{m}_{OS1}^{2} [GeV^{2}];#it{m}_{OS2}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]",100,0.,1.75*1.75, 100,0.,1.75*1.75); 
	TH2D* hAllDalitz1 = new TH2D("Dalitz1.all", ";#it{m}_{OS1}^{2} [GeV^{2}];#it{m}_{SS}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]", 100,0.,1.75*1.75, 100,0.,1.75*1.75); 
	TH2D* hAllDalitz2 = new TH2D("Dalitz2.all", ";#it{m}_{OS2}^{2} [GeV^{2}];#it{m}_{SS}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]", 100,0.,1.75*1.75, 100,0.,1.75*1.75);
	
	TH2D* hPasDalitz0 = new TH2D("Dalitz0.pas", ";#it{m}_{OS1}^{2} [GeV^{2}];#it{m}_{OS2}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]",100,0.,1.75*1.75, 100,0.,1.75*1.75); 
	TH2D* hPasDalitz1 = new TH2D("Dalitz1.pas", ";#it{m}_{OS1}^{2} [GeV^{2}];#it{m}_{SS}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]", 100,0.,1.75*1.75, 100,0.,1.75*1.75); 
	TH2D* hPasDalitz2 = new TH2D("Dalitz2.pas", ";#it{m}_{OS2}^{2} [GeV^{2}];#it{m}_{SS}^{2} [GeV^{2}];Acceptance#timesEfficiency [%]", 100,0.,1.75*1.75, 100,0.,1.75*1.75);
	
	
	hAllDalitz0->Sumw2();
	hAllDalitz1->Sumw2();
	hAllDalitz2->Sumw2();
	
	hPasDalitz0->Sumw2();
	hPasDalitz1->Sumw2();
	hPasDalitz2->Sumw2();
	
	Int_t nall = 0;
	Int_t npas = 0;
	for(Int_t entry=1 ; entry<tTruth->GetEntries() ; ++entry)
	{
		tTruth->GetEntry(entry);
		if(mOS1->size()!=1) continue;
		nall++;
		
		float mOS1sq = mOS1->at(0)*MeV2GeV*mOS1->at(0)*MeV2GeV;
		float mOS2sq = mOS2->at(0)*MeV2GeV*mOS2->at(0)*MeV2GeV;
		float mSSsq  = mSS->at(0)*MeV2GeV*mSS->at(0)*MeV2GeV;
		
		hAllDalitz0->Fill(mOS1sq,mOS2sq);
		hAllDalitz1->Fill(mOS1sq,mSSsq);
		hAllDalitz2->Fill(mOS2sq,mSSsq);
		
		Int_t bin = hEvents->FindBin(evt);
		if(hEvents->GetBinContent(bin)<1) continue;
		npas++;
		
		hPasDalitz0->Fill(mOS1sq,mOS2sq);
		hPasDalitz1->Fill(mOS1sq,mSSsq);
		hPasDalitz2->Fill(mOS2sq,mSSsq);
	}
	
	hPasDalitz0->Divide(hAllDalitz0);
	hPasDalitz1->Divide(hAllDalitz1);
	hPasDalitz2->Divide(hAllDalitz2);
	
	hPasDalitz0->Scale(100);
	hPasDalitz1->Scale(100);
	hPasDalitz2->Scale(100);
	
	TCanvas* cnv = new TCanvas("cnv","",1200,1000);
	cnv->Draw();
	cnv->Divide(2,2);
	cnv->cd(1);
	hPasDalitz0->SetContour(20);
	hPasDalitz0->Draw();
	hPasDalitz0->Draw("colz");
	exeGray->Draw();
	hPasDalitz0->Draw("colz same");
	gPad->RedrawAxis();
	
	cnv->cd(2);
	hPasDalitz1->SetContour(20);
	hPasDalitz1->Draw();
	hPasDalitz1->Draw("colz");
	exeGray->Draw();
	hPasDalitz1->Draw("colz same");
	gPad->RedrawAxis();
	
	cnv->cd(3);
	hPasDalitz2->SetContour(20);
	hPasDalitz2->Draw();
	hPasDalitz2->Draw("colz");
	exeGray->Draw();
	hPasDalitz2->Draw("colz same");
	gPad->RedrawAxis();
	
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/Dalitz.png");
	cnv->SaveAs("figures/Dalitz.eps");

	cout << "nall=" << nall << ", npas=" << npas << endl;
}