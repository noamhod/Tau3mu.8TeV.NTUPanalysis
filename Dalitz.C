//////////////////////////////////
//// root -b -l -q Dalitz.C++ ////
//////////////////////////////////

#include "std.h"
#include "type.h"
#include "postBDTcuts.h"
#include "const.h"

TString drawopt = "colz text45"; //"colz_one_palette";
bool doLog      = false;

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}

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

void setLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
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

void plotAtlasLabel2D()
{
	Double_t x = 0.50;
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

void plotLabel(TString label, Double_t x, Double_t y)
{
	TLatex t; //l.SetTextAlign(12); l.SetTextSize(tsize); 
	t.SetNDC();
	t.SetTextFont(72);
	t.SetTextColor(kBlack);
	t.DrawLatex(x,y,label);
}


void printTable(TH2* h, TString xName, TString yName, int precision=2)
{
	cout << "%----------------- " << yName << " vs " << xName << " ------------------" << endl;
	for(Int_t x=0 ; x<=h->GetNbinsX() ; ++x)
	{
		if(x==0) { cout << "\\backslashbox{$"+yName+"$~[GeV$^2$]}{$"+xName+"$~[GeV$^2$]}"; continue; }
		TString xtitle = "["+tstr(h->GetXaxis()->GetBinLowEdge(x),2)+","+tstr(h->GetXaxis()->GetBinUpEdge(x),2)+"]";
		if(x>0) cout << " &" << xtitle << " ";
		if(x==h->GetNbinsX()) cout << "\\\\\\hline\\relax" << endl; 
	}
	for(Int_t y=1 ; y<=h->GetNbinsY() ; ++y)
	{
		for(Int_t x=1 ; x<=h->GetNbinsX() ; ++x)
		{
			TString ytitle = "["+tstr(h->GetYaxis()->GetBinLowEdge(y),2)+","+tstr(h->GetYaxis()->GetBinUpEdge(y),2)+"]";
			TString acceff = "";
			if(h->GetBinContent(x,y)>0)  acceff = tstr(h->GetBinContent(x,y),precision)+"\\pm "+tstr(h->GetBinError(x,y),precision);
			if(h->GetBinContent(x,y)==0) acceff = "0";
			if(h->GetBinContent(x,y)<0)  acceff = "--";
			
			if(h->GetBinContent(x,y)>=0) acceff = "$"+acceff+"$";
			
			if(x==1)                                   cout << ytitle << " &" << acceff << " "; 
			if(x>1 && x<h->GetNbinsX())                cout << " &" << acceff << " ";  
			if(x==h->GetNbinsX() && y<h->GetNbinsY())  cout << "\\\\\\relax" << endl;
			if(x==h->GetNbinsX() && y==h->GetNbinsY()) cout << "\\\\\\hline\\hline" << endl;
		}
	}
	
	cout << "%---------------------------------------------------" << endl;
}



void regularize(TH2* hEff, TH2* hAll, int nMin=100)
{
	for(Int_t y=1 ; y<=hEff->GetNbinsY() ; ++y)
	{
		for(Int_t x=1 ; x<=hEff->GetNbinsX() ; ++x)
		{
			if(hAll->GetBinContent(x,y)<nMin)
			{
				hEff->SetBinContent(x,y,-1);
				hEff->SetBinError(x,y,0);
			}
		}
	}
}


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
	if(drawopt.Contains("z")) gStyle->SetPadRightMargin(0.16);
	else                      gStyle->SetPadRightMargin(0.08);
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
	
	gStyle->SetPaintTextFormat("4.2f");
	
	
	
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
	
	
	int      nmbins1 = 56;
	int      nmbins  = 20;
	Double_t mbins[nmbins+1];
	Double_t xMmin1 = (doLog) ? 0.3 : 0.3;
	Double_t xMmin  = (doLog) ? 0.1*0.1 : 0.3*0.3;
	Double_t xMmax1 = 1.7;
	Double_t xMmax  = 1.8*1.8;
	setLogBins(nmbins,xMmin,xMmax,mbins);
	
	
	Double_t xmax0 = 3.0;
	Double_t ymax0 = 2.2;
	
	Double_t xmax1 = 2.8;
	Double_t ymax1 = 3.0;
	
	Double_t xmax2 = 2.2;
	Double_t ymax2 = 3.0;
	
	TH1D* h1AllDalitz0;
	TH1D* h1AllDalitz1;
	TH1D* h1AllDalitz2;
	
	TH1D* h1PasDalitz0;
	TH1D* h1PasDalitz1;
	TH1D* h1PasDalitz2;
	
	TH2D* h2AllDalitz0;
	TH2D* h2AllDalitz1;
	TH2D* h2AllDalitz2;
	
	TH2D* h2PasDalitz0;
	TH2D* h2PasDalitz1;
	TH2D* h2PasDalitz2;
		
	if(doLog)
	{
		h2AllDalitz0 = new TH2D("h2Dalitz0.all", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{OS2}^{2} [GeV^{2}];Events",nmbins,mbins, nmbins,mbins); 
		h2AllDalitz1 = new TH2D("h2Dalitz1.all", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,mbins, nmbins,mbins); 
		h2AllDalitz2 = new TH2D("h2Dalitz2.all", ";Truth #it{m}_{OS2}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,mbins, nmbins,mbins);

		h2PasDalitz0 = new TH2D("h2Dalitz0.pas", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{OS2}^{2} [GeV^{2}];Events",nmbins,mbins, nmbins,mbins); 
		h2PasDalitz1 = new TH2D("h2Dalitz1.pas", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,mbins, nmbins,mbins); 
		h2PasDalitz2 = new TH2D("h2Dalitz2.pas", ";Truth #it{m}_{OS2}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,mbins, nmbins,mbins);
	}
	else
	{
		h2AllDalitz0 = new TH2D("h2Dalitz0.all", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{OS2}^{2} [GeV^{2}];Events",nmbins,xMmin,xmax0, nmbins,xMmin,ymax0); 
		h2AllDalitz1 = new TH2D("h2Dalitz1.all", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,xMmin,xmax1, nmbins,xMmin,ymax1); 
		h2AllDalitz2 = new TH2D("h2Dalitz2.all", ";Truth #it{m}_{OS2}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,xMmin,xmax2, nmbins,xMmin,ymax2);

		h2PasDalitz0 = new TH2D("h2Dalitz0.pas", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{OS2}^{2} [GeV^{2}];Events",nmbins,xMmin,xmax0, nmbins,xMmin,ymax0); 
		h2PasDalitz1 = new TH2D("h2Dalitz1.pas", ";Truth #it{m}_{OS1}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,xMmin,xmax1, nmbins,xMmin,ymax1); 
		h2PasDalitz2 = new TH2D("h2Dalitz2.pas", ";Truth #it{m}_{OS2}^{2} [GeV^{2}];Truth #it{m}_{SS}^{2} [GeV^{2}];Events", nmbins,xMmin,xmax2, nmbins,xMmin,ymax2);
		
		h1AllDalitz0 = new TH1D("h1Dalitz0.all", ";Truth #it{m}_{SS} [GeV];Events",  nmbins1,xMmin1,xMmax1); 
		h1AllDalitz1 = new TH1D("h1Dalitz1.all", ";Truth #it{m}_{OS1} [GeV];Events", nmbins1,xMmin1,xMmax1); 
		h1AllDalitz2 = new TH1D("h1Dalitz2.all", ";Truth #it{m}_{OS2} [GeV];Events", nmbins1,xMmin1,xMmax1);

		h1PasDalitz0 = new TH1D("h1Dalitz0.pas", ";Truth #it{m}_{SS} [GeV];Events",  nmbins1,xMmin1,xMmax1); 
		h1PasDalitz1 = new TH1D("h1Dalitz1.pas", ";Truth #it{m}_{OS1} [GeV];Events", nmbins1,xMmin1,xMmax1); 
		h1PasDalitz2 = new TH1D("h1Dalitz2.pas", ";Truth #it{m}_{OS2} [GeV];Events", nmbins1,xMmin1,xMmax1);
	}
	
	h1AllDalitz0->SetMarkerStyle(20); h1AllDalitz0->SetMarkerSize(0.8); h1AllDalitz0->SetMarkerColor(kBlack); h1AllDalitz0->SetLineColor(kBlack); 
	h1AllDalitz1->SetMarkerStyle(20); h1AllDalitz1->SetMarkerSize(0.8); h1AllDalitz1->SetMarkerColor(kBlack); h1AllDalitz1->SetLineColor(kBlack); 
	h1AllDalitz2->SetMarkerStyle(20); h1AllDalitz2->SetMarkerSize(0.8); h1AllDalitz2->SetMarkerColor(kBlack); h1AllDalitz2->SetLineColor(kBlack); 

	h1PasDalitz0->SetMarkerStyle(20); h1PasDalitz0->SetMarkerSize(0.8); h1PasDalitz0->SetMarkerColor(kBlack); h1PasDalitz0->SetLineColor(kBlack); 
	h1PasDalitz1->SetMarkerStyle(20); h1PasDalitz1->SetMarkerSize(0.8); h1PasDalitz1->SetMarkerColor(kBlack); h1PasDalitz1->SetLineColor(kBlack); 
	h1PasDalitz2->SetMarkerStyle(20); h1PasDalitz2->SetMarkerSize(0.8); h1PasDalitz2->SetMarkerColor(kBlack); h1PasDalitz2->SetLineColor(kBlack); 
	
	
	
	
	

	h1AllDalitz0->Sumw2();
	h1AllDalitz1->Sumw2();
	h1AllDalitz2->Sumw2();

	h1PasDalitz0->Sumw2();
	h1PasDalitz1->Sumw2();
	h1PasDalitz2->Sumw2();
	
	h2AllDalitz0->Sumw2();
	h2AllDalitz1->Sumw2();
	h2AllDalitz2->Sumw2();
	
	h2PasDalitz0->Sumw2();
	h2PasDalitz1->Sumw2();
	h2PasDalitz2->Sumw2();
	
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
		
		h1AllDalitz0->Fill(mSS->at(0)*MeV2GeV);
		h1AllDalitz1->Fill(mOS1->at(0)*MeV2GeV);
		h1AllDalitz2->Fill(mOS2->at(0)*MeV2GeV);
		
		h2AllDalitz0->Fill(mOS1sq,mOS2sq);
		h2AllDalitz1->Fill(mOS1sq,mSSsq);
		h2AllDalitz2->Fill(mOS2sq,mSSsq);
		
		Int_t bin = hEvents->FindBin(evt);
		if(hEvents->GetBinContent(bin)<1) continue;
		npas++;
		
		h1PasDalitz0->Fill(mSS->at(0)*MeV2GeV);
		h1PasDalitz1->Fill(mOS1->at(0)*MeV2GeV);
		h1PasDalitz2->Fill(mOS2->at(0)*MeV2GeV);
		
		h2PasDalitz0->Fill(mOS1sq,mOS2sq);
		h2PasDalitz1->Fill(mOS1sq,mSSsq);
		h2PasDalitz2->Fill(mOS2sq,mSSsq);
	}
	
	TH1D* h1PasDalitz0orig = (TH1D*)h1PasDalitz0->Clone();
	TH1D* h1PasDalitz1orig = (TH1D*)h1PasDalitz1->Clone();
	TH1D* h1PasDalitz2orig = (TH1D*)h1PasDalitz2->Clone();
	
	TH2D* h2PasDalitz0orig = (TH2D*)h2PasDalitz0->Clone();
	TH2D* h2PasDalitz1orig = (TH2D*)h2PasDalitz1->Clone();
	TH2D* h2PasDalitz2orig = (TH2D*)h2PasDalitz2->Clone();
	
	
	h2PasDalitz0->Divide(h2PasDalitz0,h2AllDalitz0,1,1,"B");
	h2PasDalitz1->Divide(h2PasDalitz1,h2AllDalitz1,1,1,"B");
	h2PasDalitz2->Divide(h2PasDalitz2,h2AllDalitz2,1,1,"B");
	
	h2PasDalitz0->GetZaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	h2PasDalitz1->GetZaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	h2PasDalitz2->GetZaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	
	h2PasDalitz0->Scale(100);
	h2PasDalitz1->Scale(100);
	h2PasDalitz2->Scale(100);
	

	int nMin = 100;
	regularize(h2PasDalitz0,h2AllDalitz0,nMin);
	regularize(h2PasDalitz1,h2AllDalitz1,nMin);
	regularize(h2PasDalitz2,h2AllDalitz2,nMin);

	
	////////////////////////////////////////////////////
	
	h1PasDalitz0->Divide(h1PasDalitz0,h1AllDalitz0,1,1,"B");
	h1PasDalitz1->Divide(h1PasDalitz1,h1AllDalitz1,1,1,"B");
	h1PasDalitz2->Divide(h1PasDalitz2,h1AllDalitz2,1,1,"B");
	
	h1PasDalitz0->GetYaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	h1PasDalitz1->GetYaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	h1PasDalitz2->GetYaxis()->SetTitle("Acceptance#timesEfficinecy [%]");
	
	h1PasDalitz0->Scale(100);
	h1PasDalitz1->Scale(100);
	h1PasDalitz2->Scale(100);
	
	h1PasDalitz0->SetMaximum(10);
	h1PasDalitz1->SetMaximum(10);
	h1PasDalitz2->SetMaximum(10);
	
	h1PasDalitz0->SetMinimum(0);
	h1PasDalitz1->SetMinimum(0);
	h1PasDalitz2->SetMinimum(0);
	
	////////////////////////////////////////////////////
	
	
	Double_t ymax = 5*0.8;//cnv->GetUymax();
	Double_t y = ymax;
	Double_t xOmega1 = (782-50)*MeV2GeV;
	Double_t xOmega2 = (782+50)*MeV2GeV;
	Double_t xPhi1   = (1020-50)*MeV2GeV;
	Double_t xPhi2   = (1020+50)*MeV2GeV;
	
	TLine* lOmegaMin = new TLine(xOmega1,0,xOmega1,y);
	lOmegaMin->SetLineColor(kBlack);
	lOmegaMin->SetLineWidth(2);
	lOmegaMin->SetLineStyle(3);
	TLine* lOmegaMax = new TLine(xOmega2,0,xOmega2,y);
	lOmegaMax->SetLineColor(kBlack);
	lOmegaMax->SetLineWidth(2);
	lOmegaMax->SetLineStyle(3);
	TLine* lOmegaHor = new TLine(xOmega1,y,xOmega2,y);
	lOmegaHor->SetLineColor(kBlack);
	lOmegaHor->SetLineWidth(2);
	lOmegaHor->SetLineStyle(3);
	
	TLine* lPhiMin = new TLine(xPhi1,0,xPhi1,y);
	lPhiMin->SetLineColor(kBlack);
	lPhiMin->SetLineWidth(2);
	lPhiMin->SetLineStyle(3);
	TLine* lPhiMax = new TLine(xPhi2,0,xPhi2,y);
	lPhiMax->SetLineColor(kBlack);
	lPhiMax->SetLineWidth(2);
	lPhiMax->SetLineStyle(3);
	TLine* lPhiHor = new TLine(xPhi1,y,xPhi2,y);
	lPhiHor->SetLineColor(kBlack);
	lPhiHor->SetLineWidth(2);
	lPhiHor->SetLineStyle(3);
	
	/////////////////////////////////////////////////////////////////
	
	TCanvas* cnv;
	
	gStyle->SetPadRightMargin(0.08);
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h1PasDalitz0->Draw("p1x1");
	lOmegaMin->Draw("L same"); lOmegaMax->Draw("L same"); lOmegaHor->Draw("L same"); 
	lPhiMin->Draw("L same");   lPhiMax->Draw("L same");   lPhiHor->Draw("L same");   
	plotAtlasLabel(); plotLabel("#it{#rho/#omega}",0.37,0.5); plotLabel("#it{#phi}",0.525,0.5);
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf(");
	cnv->SaveAs("figures/DalitzEff.mSS.pdf");
	cnv->SaveAs("figures/DalitzEff.mSS.png");
	cnv->SaveAs("figures/DalitzEff.mSS.eps");
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h1PasDalitz1->Draw("p1x1");
	lOmegaMin->Draw("L same"); lOmegaMax->Draw("L same");  lOmegaHor->Draw("L same");
	lPhiMin->Draw("L same");   lPhiMax->Draw("L same");    lPhiHor->Draw("L same");  
	plotAtlasLabel(); plotLabel("#it{#rho/#omega}",0.37,0.5); plotLabel("#it{#phi}",0.525,0.5);
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS1.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS1.png");
	cnv->SaveAs("figures/DalitzEff.mOS1.eps");
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h1PasDalitz2->Draw("p1x1");
	lOmegaMin->Draw("L same"); lOmegaMax->Draw("L same");  lOmegaHor->Draw("L same");
	lPhiMin->Draw("L same");   lPhiMax->Draw("L same");    lPhiHor->Draw("L same");  
	plotAtlasLabel(); plotLabel("#it{#rho/#omega}",0.37,0.5); plotLabel("#it{#phi}",0.525,0.5);
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS2.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS2.png");
	cnv->SaveAs("figures/DalitzEff.mOS2.eps");
	
	
	
	
	
	
	////////////////////////////////////////////////////////////////////
	if(drawopt.Contains("z")) gStyle->SetPadRightMargin(0.16); /////////
	////////////////////////////////////////////////////////////////////
	
	
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h2PasDalitz0->Draw(drawopt);
	plotAtlasLabel2D();
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS2.vs.mOS1.pdf");
	cnv->SaveAs("figures/DalitzEff.mOS2.vs.mOS1.png");
	cnv->SaveAs("figures/DalitzEff.mOS2.vs.mOS1.eps");
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h2PasDalitz1->Draw(drawopt);
	plotAtlasLabel2D();
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS1.pdf");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS1.png");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS1.eps");
	
	cnv = new TCanvas("cnv","",800,600); cnv->cd(); cnv->Draw();
	h2PasDalitz2->Draw(drawopt);
	plotAtlasLabel2D();
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS2.pdf");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS2.png");
	cnv->SaveAs("figures/DalitzEff.mSS.vs.mOS2.eps");
	
	
	////////////////////////////////////////////////////////////////////
	
	
	
	cnv = new TCanvas("cnv","",1200,1000);
	cnv->Draw();
	cnv->Divide(2,2);
	cnv->cd(1);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz0->SetContour(20);
		h2PasDalitz0->Draw();
		h2PasDalitz0->Draw("colz");
		exeGray->Draw();
		h2PasDalitz0->Draw("colz same");
	}
	else h2PasDalitz0->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(2);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz1->SetContour(20);
		h2PasDalitz1->Draw();
		h2PasDalitz1->Draw("colz");
		exeGray->Draw();
		h2PasDalitz1->Draw("colz same");
	}
	else h2PasDalitz1->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(3);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz2->SetContour(20);
		h2PasDalitz2->Draw();
		h2PasDalitz2->Draw("colz");
		exeGray->Draw();
		h2PasDalitz2->Draw("colz same");
	}
	else h2PasDalitz2->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf(");
	cnv->SaveAs("figures/DalitzEff2.pdf");
	cnv->SaveAs("figures/DalitzEff2.png");
	cnv->SaveAs("figures/DalitzEff2.eps");

	
	
	
	////////////////////////////////////////
	gStyle->SetPaintTextFormat("4.0f"); ////
	////////////////////////////////////////
	
	
	
	cnv = new TCanvas("cnv","",1200,1000);
	cnv->Draw();
	cnv->Divide(2,2);
	cnv->cd(1);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz0orig->SetContour(20);
		h2PasDalitz0orig->Draw();
		h2PasDalitz0orig->Draw("colz");
		exeGray->Draw();
		h2PasDalitz0orig->Draw("colz same");
	}
	else h2PasDalitz0orig->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(2);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz1orig->SetContour(20);
		h2PasDalitz1orig->Draw();
		h2PasDalitz1orig->Draw("colz");
		exeGray->Draw();
		h2PasDalitz1orig->Draw("colz same");
	}
	else h2PasDalitz1orig->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(3);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2PasDalitz2orig->SetContour(20);
		h2PasDalitz2orig->Draw();
		h2PasDalitz2orig->Draw("colz");
		exeGray->Draw();
		h2PasDalitz2orig->Draw("colz same");
	}
	else h2PasDalitz2orig->Draw(drawopt);
	gPad->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf");
	cnv->SaveAs("figures/DalitzPass.pdf");
	cnv->SaveAs("figures/DalitzPass.png");
	cnv->SaveAs("figures/DalitzPass.eps");
	
	
	
	
	
	cnv = new TCanvas("cnv","",1200,1000);
	cnv->Draw();
	cnv->Divide(2,2);
	cnv->cd(1);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2AllDalitz0->SetContour(20);
		h2AllDalitz0->Draw();
		h2AllDalitz0->Draw("colz");
		exeGray->Draw();
		h2AllDalitz0->Draw("colz same");
	}
	else h2AllDalitz0->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(2);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2AllDalitz1->SetContour(20);
		h2AllDalitz1->Draw();
		h2AllDalitz1->Draw("colz");
		exeGray->Draw();
		h2AllDalitz1->Draw("colz same");
	}
	else h2AllDalitz1->Draw(drawopt);
	gPad->RedrawAxis();
	
	cnv->cd(3);
	// if(doLog) { gPad->SetLogx();  gPad->SetLogy(); }
	if(drawopt=="colz_one_palette")
	{
		h2AllDalitz2->SetContour(20);
		h2AllDalitz2->Draw();
		h2AllDalitz2->Draw("colz");
		exeGray->Draw();
		h2AllDalitz2->Draw("colz same");
	}
	else h2AllDalitz2->Draw(drawopt);
	gPad->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/Dalitz.pdf)");
	cnv->SaveAs("figures/DalitzAll.pdf");
	cnv->SaveAs("figures/DalitzAll.png");
	cnv->SaveAs("figures/DalitzAll.eps");

	cout << "nall=" << nall << ", npas=" << npas << endl;
	
	
	cout << endl; printTable(h2PasDalitz0, "\\mOSasq", "\\mOSbsq");
	cout << endl; printTable(h2PasDalitz1, "\\mOSasq", "\\mSSsq");
	cout << endl; printTable(h2PasDalitz2, "\\mOSbsq", "\\mSSsq");
}