///////////////////////////////////////////////////
//// root -l -b -q cutsroofit.C++\(1450,2290\) ////
///////////////////////////////////////////////////

#include "std.h"
#include "const.h"
#include "postBDTcuts.h"
#include "roofit.h"

Double_t m3bodyMin = mSideBandLeftLowerMeVGlob;
Double_t m3bodyMax = mSideBandRightUpperMeVGlob;
Double_t xbmin  = mBlindMinGlob;
Double_t xbmax  = mBlindMaxGlob;
Double_t m3bodyBinSize = mBinSize;

TRandom* randGen;

string str(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return str;
}

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
	ptxt = new TPaveText(0.62,0.70,0.87,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

TLegend* leg;
void makeLegend()
{
	leg = new TLegend(0.25,0.67,0.5,0.87,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
}

bool minuitStatus(TMinuit* m) 
{
	if(!m) return false;
	TString stat = gMinuit->fCstatu;
	_INF(1,"Minuit: "<<stat<<". ");
	if(stat.Contains("SUCCESSFUL") || stat.Contains("CONVERGED") || stat.Contains("OK") ) return true;
	return false;
}

double randomizeItialGuess(double min, double max)
{
	return min + (max-min)*randGen->Uniform(); // Uniform(x1=0) returns a uniform deviate on the interval [0,x1].
}

int readData(TTree* t, RooRealVar* m3body, RooAbsData* data_m3body,
			float xFullMin,float xBlindMin,float xBlindMax,float xFullMax, float sigma, TString type, TString CRi, float xSRmin=1715, float xSRmax=1842)
{
	data_m3body->Clear();
		
	vector<float>* dummy    = 0;
	vector<float>* mass     = 0;
	vector<float>* mOS1     = 0;
	vector<float>* mOS2     = 0;
	t->SetBranchAddress("m3body",&mass);
	t->SetBranchAddress("mOS1",&mOS1);
	t->SetBranchAddress("mOS2",&mOS2);
	int ncandidates = 0;
	for(Int_t entry=1 ; entry<=t->GetEntries() ; entry++)
	{
		t->GetEntry(entry);
		
		if(mass->size()!=1) continue;
		
		bool pass = passPostBDTcut(0,mass,dummy,mOS1,mOS2, xFullMin,xBlindMin,xBlindMax,xFullMax, -1,-1, sigma, type,CRi,true);
		if(type!="bkg") pass = (pass && (mass>xSRmin && mass<xSRmax));
	
		//// add this candidtate to the set
		if(pass)
		{
			ncandidates++;
			*m3body = mass->at(0);
			data_m3body->add(RooArgSet(*m3body));
		}
	}
	return ncandidates;
}

void cutsroofit(Double_t xSBmin=0, Double_t xSBmax=0)
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
	
	TFile* fmvaout = new TFile("cutsout.muons.root","READ");
	TTree* tD = (TTree*)fmvaout->Get("fltcuts_Data");
	TTree* tS = (TTree*)fmvaout->Get("fltcuts_Wtaunu_3mu");
	
	makeAtlasLabel();
	makeLegend();
	
	//// for systematics (SB range) - only with the nominal bdtcutoff
	if(xSBmin>0  && xSBmax>0)
	{
		_INF(1,"New sidebands: "<< xSBmin << "-->" << xSBmax);
		m3bodyMin = xSBmin;
		m3bodyMax = xSBmax;
	}
	TString sidebands = tstr(m3bodyMin,0)+"-"+tstr(xbmin,0)+"-"+tstr(xbmax,0)+"-"+tstr(m3bodyMax,0);
	
	
	float Sigma = SigmaRhoOmegaPhiMeV;
	Int_t nm3bodybins = (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize);
	
	randGen = new TRandom();
	randGen->SetSeed(0); // Note that the machine clock is returned with a precision of 1 second.
						 // If one calls SetSeed(0) within a loop and the loop time is less than 1s,
						 // all generated numbers will be identical!
	
	RooRealVar* m3body = new RooRealVar("m3body","m_{3body} [MeV]",m3bodyMin,m3bodyMax);
	RooRealVar* m3bodyS = new RooRealVar("m3bodyS","m_{3body} [MeV]",m3bodyMin,m3bodyMax);
	
	// Roo Data holder for the 3body mass
	RooAbsData*  UnbinnedDataSet_m3body = new RooDataSet("data_m3body","data_m3body",RooArgSet(*m3body));
	RooAbsData*  UnbinnedDataSet_m3bodyS = new RooDataSet("signal_m3body","signal_m3body",RooArgSet(*m3bodyS));
	RooDataHist* BinnedDataSet_m3body; // set later
	RooDataHist* BinnedDataSet_m3bodyS; // set later
	
	_INF(1,"");
	
	//// for the final estimation: xbdtmaxuser=optimal cut
	//// for the inputs for optimizations, xbdtmaxuser = +1
	//// This is relevant for the SB fit only and lot for the BDT fit.
	//// For the final estimation this has to be called agin later with xbdtmaxuser instead of xbdtmax !
	Int_t ncandidatesS = readData(tS,m3bodyS,UnbinnedDataSet_m3bodyS,m3bodyMin,xbmin,xbmax,m3bodyMax, Sigma, "full","");
	Int_t ncandidatesDfull = readData(tD,m3body,UnbinnedDataSet_m3body,m3bodyMin,xbmin,xbmax,m3bodyMax,Sigma, "bkg","CR0");
	
	
	_INF(1,"");
	
	
	TCanvas* cnv;
	
	
	
	
	
	
	
	//////////////////
	//// Fit the 3body mass for the signal
	RooDataSet* rds_m3bodyS = (RooDataSet*)UnbinnedDataSet_m3bodyS;
	BinnedDataSet_m3bodyS = rds_m3bodyS->binnedClone();
	
	m3bodyS->setRange("range_m3bodyS",m3bodyMin,m3bodyMax);
	m3bodyS->setBins(nm3bodybins);
	
	// Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
	// RooRealVar* mean   = new RooRealVar("mean","mean of gaussians",1777) ;
	// RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussians",30) ;
	// RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussians",30) ;
	// RooGaussian* sig1 = new RooGaussian("sig1","Signal component 1",*m3bodyS,*mean,*sigma1) ;  
	// RooGaussian* sig2 = new RooGaussian("sig2","Signal component 2",*m3bodyS,*mean,*sigma2) ;
	// // Sum the signal components into a composite signal p.d.f.
	// RooRealVar* sig1frac = new RooRealVar("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
	// RooAddPdf* sigm3bodypdf = new RooAddPdf("sig","Signal",RooArgList(*sig1,*sig2),*sig1frac) ;
	
	// RooRealVar* mean   = new RooRealVar("mean","mean of gaussians",1777);
	// RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussian1",30);
	// RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussian2",30);
	// RooGaussian* sig1 = new RooGaussian("sig1","Signal component 1",*m3bodyS,*mean,*sigma1);
	// RooGaussian* sig2 = new RooGaussian("sig2","Signal component 2",*m3bodyS,*mean,*sigma2);
	// RooAddPdf* sigm3bodypdf = new RooAddPdf("sigm3bodypdf","sigm3bodypdf",RooArgList(*sig1,*sig2)) ;	
	
	RooRealVar* mean   = new RooRealVar("Mean","Mean",1777,1750,1800);
	RooRealVar* sigma = new RooRealVar("#sigma","#sigma",30,0,100);
	RooGaussian* sigm3bodypdf = new RooGaussian("sigm3bodypdf","sigm3bodypdf",*m3bodyS,*mean,*sigma);
	
	// // Voigtian is convolution of Gaussian and Breit-Wigner (with numerical normalization integral)
	// RooRealVar* mass = new RooRealVar("M","M",1776.82);
	// RooRealVar* sigma = new RooRealVar("#sigma","#sigma",30,0.,100);
	// RooRealVar* width = new RooRealVar("#Gamma","#Gamma",0,0,100);
	// RooVoigtian* sigm3bodypdf = new RooVoigtian("sigm3bodypdf","sigm3bodypdf",*m3bodyS,*mass,*width,*sigma);
	
	RooFitResult* fitresult_m3bodyS;
	fitresult_m3bodyS = sigm3bodypdf->fitTo( *UnbinnedDataSet_m3bodyS,Minos(kTRUE),Range("range_m3bodyS"),NormRange("range_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	// fitresult_m3bodyS = sigm3bodypdf->fitTo( *UnbinnedDataSet_m3bodyS,Minos(kTRUE),Range("rrange_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	// fitresult_m3bodyS = sigm3bodypdf->fitTo( *BinnedDataSet_m3bodyS,Minos(kTRUE),Range("range_m3bodyS"),NormRange("range_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	TMinuit* gFit_m3bodyS = gMinuit;
	bool FitStatus_m3bodyS = minuitStatus(gFit_m3bodyS);
	fitresult_m3bodyS->Print("v");
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* m3bodySFrame = m3bodyS->frame(Name("m3bodySFrame"),Title("Signal 3body mass"));
	UnbinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("m3bodyS"),MarkerSize(1),Binning(nm3bodybins));
	// UnbinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("BDT m3bodyS"),XErrorSize(0),MarkerSize(1),Binning(nm3bodybins));
	// BinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("BDT m3bodyS"),XErrorSize(0),MarkerSize(0.3),Binning(nm3bodybins));
	
	// sigm3bodypdf->plotOn(m3bodySFrame,LineWidth(2),LineColor(kBlue),Range("range_SBleft,range_SBright"),NormRange("range_SBleft,range_SBright"));
	sigm3bodypdf->plotOn(m3bodySFrame,Name("Signal"),LineWidth(2),LineColor(kBlue),Range("range_m3bodyS"),NormRange("range_m3bodyS"));
	// sigm3bodypdf->plotOn(m3bodySFrame,LineWidth(2),LineColor(kBlue),Range("range_m3bodyS"));
	sigm3bodypdf->paramOn(m3bodySFrame,Layout(0.55,0.88,0.30), Format("NEU", AutoPrecision(2))); // X size of box is from 60% to 88% of Xaxis range, top of box is at 87% of Yaxis range)
	// m3bodySFrame->getAttText()->SetTextSize(0.03);
	
	Float_t xSRmin = mean->getVal()-2.*sigma->getVal(); // 1717.97;
	Float_t xSRmax = mean->getVal()+2.*sigma->getVal(); // 1838.06;
	m3body->setRange("range_m3body_SR",xSRmin,xSRmax);
	
	cnv->SetLeftMargin(0.2);
	m3bodySFrame->SetTitleOffset(2,"Y");
	m3bodySFrame->Draw();
	ptxt->Draw("same");
	cnv->SaveAs("figures/CutBasedFitResults.Gaus."+sidebands+".png");
	cnv->SaveAs("figures/CutBasedFitResults.Gaus."+sidebands+".eps");
	cnv->SaveAs("figures/CutBasedFitResults."+sidebands+".pdf(");

	
	Double_t chi2Signalpdf = m3bodySFrame->chiSquare("Signal","m3bodySFrame",2);
	cout << "chi2Signalpdf=" << chi2Signalpdf << endl;
	
	_INF(1,"------------------- Done fitting the signal shape -------------------");
	
	
	
	/////////////////////////////
	//// Fit the m3body sidebands
	//////////////////////////////
	
	
	//// for the final estimation: xbdtmaxuser=optimal cut
	//// for the inputs for optimizations, xbdtmaxuser = +1
	//// This is relevant for the SB fit only and lot for the BDT fit.
	//// For the final estimation this has to be called here with xbdtmaxuser (not xbdtmax !)
	Int_t ncandidatesD = readData(tD,m3body,UnbinnedDataSet_m3body, m3bodyMin,xbmin,xbmax,m3bodyMax, Sigma, "bkg","CR0");
	
	
	RooDataSet* rds_m3body = (RooDataSet*)UnbinnedDataSet_m3body;
	BinnedDataSet_m3body = rds_m3body->binnedClone();
	
	m3body->setRange("range_m3body",m3bodyMin,m3bodyMax);
	m3body->setBins(nm3bodybins);
	m3body->setRange("range_SBleft",  m3bodyMin,xbmin);
	m3body->setRange("range_SBright", xbmax,m3bodyMax);
	

	
	RooRealVar* a0 = new RooRealVar("a0","a0", randomizeItialGuess(-1.0,0.0)/*-1.42803e-01*/, -1.0, 1.0); a0->setError(0.00001);
	RooChebychev* bkgm3bodypdf0 = new RooChebychev("bkgm3bodypdf0","bkgm3bodypdf0",*m3body,RooArgSet(*a0));
	
	RooRealVar* b0 = new RooRealVar("b0","b0", randomizeItialGuess(-1.0,0.0)/*-1.42803e-01*/, -1.0, 1.0); b0->setError(0.00001);
	RooRealVar* b1 = new RooRealVar("b1","b1", randomizeItialGuess( 0.0,0.1)/* 2.87161e-02*/, -1.0, 1.0); b1->setError(0.00001);
	RooChebychev* bkgm3bodypdf1 = new RooChebychev("bkgm3bodypdf1","bkgm3bodypdf1",*m3body,RooArgSet(*b0,*b1));
	
	RooRealVar* c0 = new RooRealVar("c0","c0", randomizeItialGuess(-1.00,0.0)/*-1.42803e-01*/, -1.0, 1.0); c0->setError(0.00001);
	RooRealVar* c1 = new RooRealVar("c1","c1", randomizeItialGuess( 0.00,0.1)/* 2.87161e-02*/, -1.0, 1.0); c1->setError(0.00001);
	RooRealVar* c2 = new RooRealVar("c2","c2", randomizeItialGuess(-0.01,0.0)/*-1.45196e-03*/, -1.0, 1.0); c2->setError(0.00001);
	RooChebychev* bkgm3bodypdf2 = new RooChebychev("bkgm3bodypdf2","bkgm3bodypdf2",*m3body,RooArgSet(*c0,*c1,*c2));
	
	RooRealVar* d0 = new RooRealVar("d0","d0", randomizeItialGuess(-1.00,0.0)/*-1.42803e-01*/, -1.0, 1.0); d0->setError(0.00001);
	RooRealVar* d1 = new RooRealVar("d1","d1", randomizeItialGuess( 0.00,0.1)/* 2.87161e-02*/, -1.0, 1.0); d1->setError(0.00001);
	RooRealVar* d2 = new RooRealVar("d2","d2", randomizeItialGuess(-0.01,0.0)/*-1.45196e-03*/, -1.0, 1.0); d2->setError(0.00001);
	RooRealVar* d3 = new RooRealVar("d3","d3", randomizeItialGuess(-0.01,0.0)/*-6.00993e-04*/, -1.0, 1.0); d3->setError(0.00001);
	RooChebychev* bkgm3bodypdf3 = new RooChebychev("bkgm3bodypdf3","bkgm3bodypdf3",*m3body,RooArgSet(*d0,*d1,*d2,*d3));
	
	RooRealVar* e0 = new RooRealVar("e0","e0", randomizeItialGuess(-1.000,0.0)/*-1.42803e-01*/, -1.0, 1.0); e0->setError(0.00001);
	RooRealVar* e1 = new RooRealVar("e1","e1", randomizeItialGuess( 0.000,0.1)/* 2.87161e-02*/, -1.0, 1.0); e1->setError(0.00001);
	RooRealVar* e2 = new RooRealVar("e2","e2", randomizeItialGuess(-0.010,0.0)/*-1.45196e-03*/, -1.0, 1.0); e2->setError(0.00001);
	RooRealVar* e3 = new RooRealVar("e3","e3", randomizeItialGuess(-0.010,0.0)/*-6.00993e-04*/, -1.0, 1.0); e3->setError(0.00001);
	RooRealVar* e4 = new RooRealVar("e4","e4", randomizeItialGuess(-0.001,0.0)/*-6.00993e-04*/, -1.0, 1.0); e4->setError(0.00001);
	RooChebychev* bkgm3bodypdf4 = new RooChebychev("bkgm3bodypdf4","bkgm3bodypdf4",*m3body,RooArgSet(*e0,*e1,*e2,*e3,*e4));
	
	RooRealVar* f0 = new RooRealVar("f0","f0", randomizeItialGuess(-1.0000,0.0)/*-1.42803e-01*/, -1.0, 1.0); f0->setError(0.00001);
	RooRealVar* f1 = new RooRealVar("f1","f1", randomizeItialGuess( 0.0000,0.1)/* 2.87161e-02*/, -1.0, 1.0); f1->setError(0.00001);
	RooRealVar* f2 = new RooRealVar("f2","f2", randomizeItialGuess(-0.0100,0.0)/*-1.45196e-03*/, -1.0, 1.0); f2->setError(0.00001);
	RooRealVar* f3 = new RooRealVar("f3","f3", randomizeItialGuess(-0.0100,0.0)/*-6.00993e-04*/, -1.0, 1.0); f3->setError(0.00001);
	RooRealVar* f4 = new RooRealVar("f4","f4", randomizeItialGuess(-0.0010,0.0)/*-6.00993e-04*/, -1.0, 1.0); f4->setError(0.00001);
	RooRealVar* f5 = new RooRealVar("f5","f5", randomizeItialGuess(-0.0001,0.0)/*-6.00993e-04*/, -1.0, 1.0); f5->setError(0.00001);
	RooChebychev* bkgm3bodypdf5 = new RooChebychev("bkgm3bodypdf5","bkgm3bodypdf5",*m3body,RooArgSet(*f0,*f1,*f2,*f3,*f4,*f5));
	
	RooRealVar* g0 = new RooRealVar("g0","g0",1.,0.,100.); g0->setError(0.00001);
	RooRealVar* g1 = new RooRealVar("g1","g1",1.,0.,100.); g1->setError(0.00001);
	RooRealVar* g2 = new RooRealVar("g2","g2",1.,0.,100.); g2->setError(0.00001);
	RooRealVar* g3 = new RooRealVar("g3","g3",1.,0.,100.); g3->setError(0.00001);
	RooRealVar* g4 = new RooRealVar("g4","g4",1.,0.,100.); g4->setError(0.00001);
	RooRealVar* g5 = new RooRealVar("g5","g5",1.,0.,100.); g5->setError(0.00001);
	RooRealVar* g6 = new RooRealVar("g6","g6",1.,0.,100.); g6->setError(0.00001);
	RooBernstein* bkgm3bodypdf6 = new RooBernstein("bkgm3bodypdf6","bkgm3bodypdf6",*m3body,RooArgSet(*g0,*g1,*g2,*g3));
	
	
	// --- Build Argus+Bernstein background PDF ---
	Double_t argparmin = -50;
	if(m3bodyMin<1450 || m3bodyMax>2290) argparmin = -100;
	RooRealVar* argpar = new RooRealVar("argpar","argpar",-1.,argparmin,0.); argpar->setError(0.00001);
	RooArgusBG* argus = new RooArgusBG("argus","argus",*m3body,RooConst(m3bodyMax),*argpar);
	RooRealVar* p0 = new RooRealVar("p0","p0",1.,0.,100.); p0->setError(0.00001);
	RooRealVar* p1 = new RooRealVar("p1","p1",1.,0.,100.); p1->setError(0.00001);
	RooRealVar* p2 = new RooRealVar("p2","p2",1.,0.,100.); p2->setError(0.00001);
	RooBernstein* polyn = new RooBernstein("poly","poly",*m3body,RooArgSet(*p0,*p1,*p2));
	// RooChebychev* polyn = new RooChebychev("polyn","polyn",*m3body,RooArgSet(*p0,*p1));
	// RooPolynomial* polyn = new RooPolynomial("polyn","polyn",*m3body,RooArgSet(*p0,*p1));
	// --- Construct signal+background PDF ---
	RooRealVar* npolyn = new RooRealVar("npolyn","npolyn",randomizeItialGuess(1.,1000.),1.,1000.);
	RooRealVar* nargus = new RooRealVar("nargus","nargus",randomizeItialGuess(1.,1000.),1.,1000.);
	RooAddPdf* bkgm3bodypdf7 = new RooAddPdf("bkgm3bodypdf7","bkgm3bodypdf7",RooArgList(*polyn,*argus),RooArgList(*npolyn,*nargus));
	// RooAbsPdf* bkgm3bodypdf7 = argus;
	
	
	
		
	// // To construct a proper p.d.f, the formula expression is explicitly normalized internally by dividing 
	// // it by a numeric integral of the expresssion over x in the defined range
	// RooRealVar* a0 = new RooRealVar("a0","a0", randomizeItialGuess(-1.0,0.0)/*-1.42803e-01*/, -1.0, 1.0);
	// RooRealVar* a1 = new RooRealVar("a1","a1", randomizeItialGuess( 0.0,0.1)/* 2.87161e-02*/, -1.0, 1.0);
	// RooRealVar* a2 = new RooRealVar("a2","a2", randomizeItialGuess(-0.01,0.0)/*-1.45196e-03*/, -1.0, 1.0);
	// RooRealVar* a3 = new RooRealVar("a3","a3", randomizeItialGuess(-0.01,0.0)/*-6.00993e-04*/, -1.0, 1.0);
	// RooChebychev* cheb;
	// if(xbdtmin>-0.5) cheb = new RooChebychev("cheb","cheb",*m3body,RooArgSet(*a0));
	// else             cheb = new RooChebychev("cheb","cheb",*m3body,RooArgSet(*a0,*a1));
	// RooRealVar* meanPhi  = new RooRealVar("meanPhi","mean of phi"/*,randomizeItialGuess(1550.,1650.)*/,1550.,1650.);  meanPhi->setError(0.1);
	// RooRealVar* sigmaPhi = new RooRealVar("sigmaPhi","width of phi"/*,randomizeItialGuess(1.,50.)*/,1.,50.);  sigmaPhi->setError(0.1);
	// RooGaussian* phi = new RooGaussian("phi","phi",*m3body,*meanPhi,*sigmaPhi);
	// // Sum the background components into a composite signal p.d.f.
	// RooRealVar* Nphi = new RooRealVar("Nphi","yield of phi in background"/*,(int)randomizeItialGuess(0,50)*/,0,50); Nphi->setError(0.1);
	// RooRealVar* Ncheb = new RooRealVar("Ncheb","yield of Chebychev in background"/*,(int)randomizeItialGuess(0,2000)*/,0,2000); Ncheb->setError(1);
	// RooAddPdf* bkgm3bodypdf = new RooAddPdf("bkgm3bodypdf","bkgm3bodypdf",RooArgSet(*cheb,*phi),RooArgSet(*Ncheb,*Nphi));
		
		
	RooFitResult* fitresult_m3body0;
	fitresult_m3body0 = bkgm3bodypdf0->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body0 = gMinuit;
	bool FitStatus_m3body0 = minuitStatus(gFit_m3body0);
	fitresult_m3body0->Print("v");
	
	RooFitResult* fitresult_m3body1;
	fitresult_m3body1 = bkgm3bodypdf1->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body1 = gMinuit;
	bool FitStatus_m3body1 = minuitStatus(gFit_m3body1);
	fitresult_m3body1->Print("v");
	
	RooFitResult* fitresult_m3body2;
	fitresult_m3body2 = bkgm3bodypdf2->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body2 = gMinuit;
	bool FitStatus_m3body2 = minuitStatus(gFit_m3body2);
	fitresult_m3body2->Print("v");
	
	RooFitResult* fitresult_m3body3;
	fitresult_m3body3 = bkgm3bodypdf3->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body3 = gMinuit;
	bool FitStatus_m3body3 = minuitStatus(gFit_m3body3);
	fitresult_m3body3->Print("v");
	
	RooFitResult* fitresult_m3body4;
	fitresult_m3body4 = bkgm3bodypdf4->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body4 = gMinuit;
	bool FitStatus_m3body4 = minuitStatus(gFit_m3body4);
	fitresult_m3body4->Print("v");
	
	RooFitResult* fitresult_m3body5;
	fitresult_m3body5 = bkgm3bodypdf5->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body5 = gMinuit;
	bool FitStatus_m3body5 = minuitStatus(gFit_m3body5);
	fitresult_m3body5->Print("v");
	
	RooFitResult* fitresult_m3body6;
	fitresult_m3body6 = bkgm3bodypdf6->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body6 = gMinuit;
	bool FitStatus_m3body6 = minuitStatus(gFit_m3body6);
	fitresult_m3body6->Print("v");
	
	RooFitResult* fitresult_m3body7;
	fitresult_m3body7 = bkgm3bodypdf7->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body7 = gMinuit;
	bool FitStatus_m3body7 = minuitStatus(gFit_m3body7);
	fitresult_m3body7->Print("v");
	
	
	
	
	RooAbsPdf* bkgm3bodypdfNominal0 = bkgm3bodypdf7;
	RooAbsPdf* bkgm3bodypdfNominal1 = bkgm3bodypdf0;
	RooAbsPdf* bkgm3bodypdfNominal2 = bkgm3bodypdf1;
	RooAbsPdf* bkgm3bodypdfNominal3 = bkgm3bodypdf2; // bkgm3bodypdf3;
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* m3bodyFrame = m3body->frame(Name("m3bodyFrame"),Title("3body mass in the sidebands of CR_{0}"));
	UnbinnedDataSet_m3body->plotOn(m3bodyFrame,Name("m3body SB"),MarkerSize(1),Binning(nm3bodybins));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal1->plotOn(m3bodyFrame,Name("Alternative A"),LineWidth(2),LineColor(kRed),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal2->plotOn(m3bodyFrame,Name("Alternative B"),LineWidth(2),LineColor(kViolet),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal3->plotOn(m3bodyFrame,Name("Alternative C"),LineWidth(2),LineColor(kGreen+1),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	
	// bkgm3bodypdfNominal0->paramOn(m3bodyFrame,Layout(0.25,0.5,0.4), Format("NEU", AutoPrecision(3))); // X size of box is from 25% to 50% of Xaxis range, top of box is at 40% of Yaxis range)
	// m3bodyFrame->getAttText()->SetTextSize(0.018);
	
	leg->Clear();
	leg->AddEntry(m3bodyFrame->findObject("m3body SB"),"Data in SB_{1}","pl");
	leg->AddEntry(m3bodyFrame->findObject("Nominal"),"Nominal fit","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative A"),"Alternative A","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative B"),"Alternative B","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative C"),"Alternative C","l");
	
	cnv->SetLeftMargin(0.2);
	m3bodyFrame->SetMaximum(1.3*m3bodyFrame->GetMaximum());
	m3bodyFrame->SetTitleOffset(2,"Y");
	m3bodyFrame->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/CutBasedFitResults.SB."+sidebands+".png");
	cnv->SaveAs("figures/CutBasedFitResults.SB."+sidebands+".eps");
	cnv->SaveAs("figures/CutBasedFitResults."+sidebands+".pdf)");
	
	Double_t chi2SBpdf0 = m3bodyFrame->chiSquare("Nominal",      "m3body SB",6);
	Double_t chi2SBpdf1 = m3bodyFrame->chiSquare("Alternative A","m3body SB",2);
	Double_t chi2SBpdf2 = m3bodyFrame->chiSquare("Alternative B","m3body SB",3);
	Double_t chi2SBpdf3 = m3bodyFrame->chiSquare("Alternative C","m3body SB",4);
	cout << "chi2SBpdf0=" << chi2SBpdf0 << endl;
	cout << "chi2SBpdf1=" << chi2SBpdf1 << endl;
	cout << "chi2SBpdf2=" << chi2SBpdf2 << endl;
	cout << "chi2SBpdf3=" << chi2SBpdf3 << endl;
	
	
	
	
	
	

	
	
	RooAbsReal* integralSR0      = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB0left  = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB0right = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBright"));
	float nSB0 = integralSB0left->getVal()+integralSB0right->getVal();
	float nSR0 = integralSR0->getVal();
	float dnSR0 = (nSR0/nSB0)*sqrt(nSB0);
	float RooFitScale0 = ncandidatesD/nSB0;
	
	RooAbsReal* integralSR1      = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB1left  = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB1right = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_SBright"));
	float nSB1 = integralSB1left->getVal()+integralSB1right->getVal();
	float nSR1 = integralSR1->getVal();
	float dnSR1 = (nSR1/nSB1)*sqrt(nSB1);
	float RooFitScale1 = ncandidatesD/nSB1;
	
	RooAbsReal* integralSR2      = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB2left  = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB2right = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_SBright"));
	float nSB2 = integralSB2left->getVal()+integralSB2right->getVal();
	float nSR2 = integralSR2->getVal();
	float dnSR2 = (nSR2/nSB2)*sqrt(nSB2);
	float RooFitScale2 = ncandidatesD/nSB2;
	
	RooAbsReal* integralSR3      = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB3left  = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB3right = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_SBright"));
	float nSB3 = integralSB3left->getVal()+integralSB3right->getVal();
	float nSR3 = integralSR3->getVal();
	float dnSR3 = (nSR3/nSB3)*sqrt(nSB3);
	float RooFitScale3 = ncandidatesD/nSB3;
	
	
	int ncandidatesSsr1 = readData(tS,m3bodyS,UnbinnedDataSet_m3bodyS, m3bodyMin,xbmin,xbmax,m3bodyMax, Sigma, "sig","CR1",1742,1812);
	float efficiencyS = ((float)ncandidatesSsr1/99900.)*100.;
	
	int ncandidatesD_CR1 = readData(tD,m3body,UnbinnedDataSet_m3body, m3bodyMin,xbmin,xbmax,m3bodyMax, Sigma, "bkg","CR1");
	
	
	//// Write the TF1 object and more numbers to a root file
	TFile* fpdf = new TFile("CutBasedFitResults."+sidebands+".root","RECREATE");
	TVectorF NSB(1);
	NSB[0]=ncandidatesD;
	NSB.Write("nSB"); // have to see if this is "all" according to "absBDTmaxUser"
	
	TVectorF NSR(4);
	TVectorF chi2SB(4);
	NSR[0]=RooFitScale0*nSR0; chi2SB[0]=chi2SBpdf0;
	NSR[1]=RooFitScale1*nSR1; chi2SB[1]=chi2SBpdf1;
	NSR[2]=RooFitScale2*nSR2; chi2SB[2]=chi2SBpdf2;
	NSR[3]=RooFitScale3*nSR3; chi2SB[3]=chi2SBpdf3;
	NSR.Write("nSR"); // have to see if this is "all" according to "absBDTmaxUser"
	chi2SB.Write("chi2SB");
	// and to read the file do: TVectorF* NSB = (TVectorF*)fpdf->Get("nSB"); float nSB = ((*NSB))[0];
	fpdf->Write();
	fpdf->Close();
	delete fpdf;
	
	string filename = "CutBasedFitResults."+(string)sidebands+".txt";
	ofstream* ofstr = new ofstream(filename.c_str());
	stringstream strm; strm << "";
	string tmp = "";
	string summary = "";
	summary += "========================================================================================\n";
	summary += "m3body left sideband    ["+str(m3bodyMin,0)+","+str(xbmin,0)+"] MeV\n";
	summary += "m3body right sideband   ["+str(xbmax,0)+","+str(m3bodyMax,0)+"] MeV\n";
	summary += "m3body signal region    ["+str(xSRmin,0)+","+str(xSRmax,0)+"] MeV (for sigma="+str(sigma->getVal(),2)+" MeV)\n";
	summary += "========================================================================================\n";
	summary += "Ch2^2/DOF Signal        "+str(chi2Signalpdf,2)+"\n";
	summary += "Ch2^2/DOF SB            "+str(chi2SBpdf0,2)+" | "+str(chi2SBpdf1,2)+" | "+str(chi2SBpdf2,2)+" | "+str(chi2SBpdf3,2)+"\n";
	summary += "========================================================================================\n";
	summary += "Count:    nSB           "+str(ncandidatesD,0)+"\n";
	summary += "Mass fit: nSB           "+str(nSB0,3)+" | "+str(nSB1,3)+" | "+str(nSB2,3)+"\n";
	summary += "Mass fit: nSR           "+str(nSR0,3)+" | "+str(nSR2,3)+" | "+str(nSR2,3)+"\n";
	summary += "RooFit scale            "+str(RooFitScale0,5)+" | "+str(RooFitScale1,5)+" | "+str(RooFitScale2,5)+"\n";
	summary += "Scaled nSR              "+str(RooFitScale0*nSR0,3)+" +- "+str((nSR0/nSB0)*sqrt(RooFitScale0*nSB0),3)+" | "
										 +str(RooFitScale1*nSR1,3)+" +- "+str((nSR1/nSB1)*sqrt(RooFitScale1*nSB1),3)+" | "
										 +str(RooFitScale2*nSR2,3)+" +- "+str((nSR2/nSB2)*sqrt(RooFitScale2*nSB2),3)+" | "
										 +str(RooFitScale3*nSR3,3)+" +- "+str((nSR3/nSB3)*sqrt(RooFitScale3*nSB3),3)+"\n";
	float maxDiff = -1e20;
	maxDiff = (fabs(RooFitScale0*nSR0-RooFitScale1*nSR1)>maxDiff) ? fabs(RooFitScale0*nSR0-RooFitScale1*nSR1) : maxDiff;
	maxDiff = (fabs(RooFitScale0*nSR0-RooFitScale2*nSR2)>maxDiff) ? fabs(RooFitScale0*nSR0-RooFitScale2*nSR2) : maxDiff;
	maxDiff = (fabs(RooFitScale0*nSR0-RooFitScale3*nSR3)>maxDiff) ? fabs(RooFitScale0*nSR0-RooFitScale3*nSR3) : maxDiff;
	summary += "dSR syst.               "+str(maxDiff,3)+" ("+str(maxDiff/(RooFitScale0*nSR0)*100,1)+"\%)\n";
	summary += "Acc*Eff in SR           "+str(efficiencyS,3)+"\%\n";
	summary += "========================================================================================\n";
	cout << summary << endl;
	(*ofstr) << summary;
	ofstr->close();
}
