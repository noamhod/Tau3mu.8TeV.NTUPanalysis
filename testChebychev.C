/////////////////////////////////////////
//// root -l -b -q testChebychev.C++ ////
/////////////////////////////////////////

#include <iostream>
#include <stdlib.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMinuit.h"

#include <RooFit.h>
#include <RooRealVar.h>
#include <RooGenericPdf.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooVoigtian.h>
#include <RooExponential.h>
#include <RooBinning.h>
#include <RooExtendPdf.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooArgSet.h>
#include <RooMinuit.h>
#include <RooChebychev.h>
#include <RooBernstein.h>


using namespace RooFit;
using namespace std;

Double_t xFullMin   = 1450;
Double_t xFullMax   = 2290;
Double_t xBlindMin  = 1690;
Double_t xBlindMax  = 1870;
Double_t xSignalmin = 1720;
Double_t xSignalmax = 1840;
Double_t xBinSize   = 30;

void testChebychev()
{
	TCanvas* cnv = new TCanvas("cnv","",800,600);
	
	Int_t xNbins = (Int_t)((xFullMax-xFullMin)/xBinSize);
	RooRealVar* x = new RooRealVar("x","x [units]",xFullMin,xFullMax);
	x->setBins(xNbins);
	x->setRange("range_x",      xFullMin,xFullMax);
	x->setRange("range_left",   xFullMin,xBlindMin);
	x->setRange("range_right",  xBlindMax,xFullMax);
	x->setRange("range_blind",  xBlindMin,xBlindMax);
	x->setRange("range_signal", xSignalmin,xSignalmax);
	TString fitrange = "range_left,range_right"; // "range_x"	
	
	
	// RooRealVar* a0 = new RooRealVar("a0","a0", -1.42803e-01, -1.0, 1.0);
	// RooRealVar* a1 = new RooRealVar("a1","a1", +2.87161e-02, -1.0, 1.0);
	// RooRealVar* a2 = new RooRealVar("a2","a2", -1.45196e-03, -1.0, 1.0);
	// RooRealVar* a3 = new RooRealVar("a3","a3", -6.00993e-04, -1.0, 1.0);
	RooRealVar* a0 = new RooRealVar("a0","a0", -0.3497, -1.0, 1.0);
	RooRealVar* a1 = new RooRealVar("a1","a1", -0.1187, -1.0, 1.0);
	RooRealVar* a2 = new RooRealVar("a2","a2", -0.0885, -1.0, 1.0);
	RooRealVar* a3 = new RooRealVar("a3","a3", -0.0040, -1.0, 1.0);
	RooChebychev* pdfCheb = new RooChebychev("pdfCheb","pdfCheb",*x,RooArgSet(*a0,*a1,*a2,*a3));

	RooAbsData* xData = pdfCheb->generate(*x,1000,Range(fitrange));
	
	
	RooRealVar* b0 = new RooRealVar("b0","b0",1.,0.,100.);
	RooRealVar* b1 = new RooRealVar("b1","b1",1.,0.,100.);
	RooRealVar* b2 = new RooRealVar("b2","b2",1.,0.,100.);
	RooRealVar* b3 = new RooRealVar("b3","b3",1.,0.,100.);
	RooRealVar* b4 = new RooRealVar("b4","b4",1.,0.,100.);
	RooRealVar* b5 = new RooRealVar("b5","b5",1.,0.,100.);
	RooRealVar* b6 = new RooRealVar("b6","b6",1.,0.,100.);
	RooBernstein* pdfBern = new RooBernstein("pdfBern","pdfBern",*x,RooArgSet(*b0,*b1,*b2,*b3,*b4,*b5,*b6));
	
	
	
	// RooChebychev* pdf = pdfCheb;
	RooBernstein* pdf = pdfBern;
	
	
	
	
	RooFitResult* fitresult = pdf->fitTo(*xData,Minos(kTRUE),Range(fitrange),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TString stat = gMinuit->fCstatu;
	cout << "Minuit: " << stat << endl;
	fitresult->Print("v");
	
	RooPlot* xFrame = x->frame(Name("xFrame"),Title("Chebychev sideband fit test"));
	xData->plotOn(xFrame,Name("x"),MarkerSize(1),Binning(xNbins));
	pdf->plotOn(xFrame,LineWidth(2),LineColor(kBlue),Range(fitrange),NormRange(fitrange));
	pdf->plotOn(xFrame,LineWidth(2),LineColor(kRed),LineStyle(kDashed),Range("range_blind"),NormRange(fitrange));
	pdf->paramOn(xFrame,Layout(0.62,0.88,0.4), Format("NEU", AutoPrecision(3))); 
	xFrame->getAttText()->SetTextSize(0.03);
	cnv->SetLeftMargin(0.2);
	xFrame->SetTitleOffset(2,"Y");
	xFrame->Draw();
	cnv->SaveAs("testChebychev.pdf");
	
	RooAbsReal* integralFull    = pdf->createIntegral(*x,Range("range_x"));
	RooAbsReal* integralSR      = pdf->createIntegral(*x,Range("range_signal"));
	RooAbsReal* integralSBleft  = pdf->createIntegral(*x,Range("range_left"));
	RooAbsReal* integralSBright = pdf->createIntegral(*x,Range("range_right"));
	float nSB0 = integralSBleft->getVal()+integralSBright->getVal();
	float nSR0 = integralSR->getVal();
	float nAll = integralFull->getVal();
	cout << "nAll=" << nAll << endl;
	cout << "nSB0=" << nSB0 << endl;
	cout << "nSR0=" << nSR0 << endl;
}