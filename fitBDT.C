///////////////////////////////////
//// root -b -l -q fitBDT.C++ /////
///////////////////////////////////
#include "std.h"
#include "const.h"
#include "postBDTcuts.h"
#include <iterator>
#include <cstdlib>

Bool_t reject;
Double_t m3bodyBinSize = mBinSize;
Double_t m3bodyMin = mSideBandLeftLowerMeVGlob;
Double_t m3bodyMax = mSideBandRightUpperMeVGlob;
Double_t xbmin  = mBlindMinGlob;
Double_t xbmax  = mBlindMaxGlob;
Int_t    nparbdt1 = 5;
Int_t    nparbdt2 = 6;
Int_t    nparmass1 = 5;
Int_t    nparmass2 = 4;
Int_t    nparmass3 = 2;
Int_t    nparmass4 = 2;
Int_t    nparmass = 1;
Double_t xbdtmin  = minBDTcut;
Double_t xbdtmax  = +1.0;
Double_t xbdtopt  = optBDTcut;
Int_t    nbdtbins = nBDTbins;
Float_t  ninitS   = 99900; // generated events in the signal sample

// Int_t nparsig = 3;
// Double_t fM3bodySignal(Double_t* X, Double_t* par)
// {
// 	Double_t x = X[0];
// 	// return par[0]*exp(-pow((x-par[1])/par[2],2))+par[3]*exp(-pow((x-par[4])/par[5],2));
// 	return par[0]*exp(-pow((x-par[1])/par[2],2));
// }

Double_t fBkgBDT1(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	return par[0]+par[1]*exp(x*par[2])+par[3]*exp(x*par[4]);
}
Double_t fBkgBDT2(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	return par[0]+par[1]*exp(x*par[2])+par[3]*pow(fabs(x),x*par[4])*exp(x*par[5]);
}
Double_t fBkgMassCR0_1(Double_t* X, Double_t* par) // red
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0]+par[1]*x+par[2]*pow(x,2)+par[3]*pow(x,3);//+par[4]*pow(x,4)+par[5]*pow(x,5);
}
Double_t fBkgMassCR0_2(Double_t* X, Double_t* par) // blue
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0]*pow(x,par[1])/(par[2]-exp(x*par[3]));
}
Double_t fBkgMassCR0_3(Double_t* X, Double_t* par) // green
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	// return par[0]*exp(x*par[1])*pow(x,x*par[2]);
	return par[0]/pow(x,par[1]);
}
Double_t fBkgMassCR0_4(Double_t* X, Double_t* par) // orange
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0]*exp(x*par[1]);
}
Double_t fBkgMassCR1(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0];
}


Int_t npar0 = 3;//6;
Int_t npar1 = 3;
Double_t fBkg0(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	// return par[0] + par[1]*exp(x*par[2]) + par[3]*exp(-pow((x-par[4])/par[5],2));
	return par[0] + par[1]*exp(x*par[2]);
}
Double_t fBkg1(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*exp(x*par[2]);
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

void fitBDT()
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
	
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");

	vector<float>* bdtscore = 0;
	vector<float>* m3body   = 0;
	vector<float>* mOS1     = 0;
	vector<float>* mOS2     = 0;
	tD->SetBranchAddress("score",&bdtscore);
	tD->SetBranchAddress("m3body",&m3body);
	tD->SetBranchAddress("mOS1",&mOS1);
	tD->SetBranchAddress("mOS2",&mOS2);
	xbdtmax = -999.;
	Int_t ncandidates = 0;
	float resSigma = SigmaRhoOmegaPhiMeV;
	for(Int_t entry=1 ; entry<=tD->GetEntries() ; entry++)
	{
		tD->GetEntry(entry);
		if(bdtscore->size()>1) continue;
		bool pass = passPostBDTcut(0,m3body,bdtscore,mOS1,mOS2, m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, resSigma, "bkg","CR0");
		// xbdtmax = (pass && bdtscore->at(0)>xbdtmax) ? bdtscore->at(0) : xbdtmax;
		if(pass) ncandidates++;		
		// for(unsigned int i=0 ; i<bdtscore->size() ; ++i)
		// {
		// 	bool pass = (((m3body->at(i)>m3bodyMin && m3body->at(i)<xbmin) || (m3body->at(i)>xbmax && m3body->at(i)<m3bodyMax))  &&  bdtscore->at(i)>xbdtmin);
		// 	xbdtmax = (pass && bdtscore->at(i)>xbdtmax) ? bdtscore->at(i) : xbdtmax;
		// 	if(pass) ncandidates++;
		// }
	}
	// xbdtmax = getXBDTMax(tD,m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, resSigma, "bkg","CR0");
	// Int_t ncandidates = tD->GetEntries();
	cout << "xbdtmax = " << xbdtmax << endl;
	cout << "ncandidates = " << ncandidates << endl;

	TCanvas* cnv = new TCanvas("","",600,400);
	
	TString sm3bodyMin = tstr(m3bodyMin,0);
	TString sm3bodyMax = tstr(m3bodyMax,0);
	TString sm3bodyBlindMin = tstr(xbmin,0);
	TString sm3bodyBlindMax = tstr(xbmax,0);
	TString sxBDTMin   = tstr(xbdtmin,4);
	TString sxBDTOpt   = tstr(xbdtopt,4);
	TString Sigma = "30.";
	TString basecuts      = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "base","");
	TString sidebandscuts = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "sidebands","");
	TString cuts_bkg_SB0  = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "bkg","CR0");
	TString cuts_sig_SB0  = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "sig","CR0");
	TString cuts_bkg_SB1  = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "bkg","CR1");
	TString cuts_sig_SB1  = postBDTcut(sm3bodyMin,sm3bodyBlindMin,sm3bodyBlindMax,sm3bodyMax, sxBDTMin,sxBDTOpt, Sigma, "sig","CR1");
	cout << "sidebandscuts = " << sidebandscuts << endl;
	cout << "basecuts      = " << basecuts      << endl;
	cout << "cuts_bkg_SB0  = " << cuts_bkg_SB0  << endl;
	cout << "cuts_sig_SB0  = " << cuts_sig_SB0  << endl;
	cout << "cuts_bkg_SB1  = " << cuts_bkg_SB1  << endl;
	cout << "cuts_sig_SB1  = " << cuts_sig_SB1  << endl;
	
	Float_t ninitS = 99900; 
	Float_t nallS  = tS->GetEntries(cuts_sig_SB0);
	Float_t nallD  = tD->GetEntries(cuts_bkg_SB0);
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "Npassing rho/omega,phi: " << tD->GetEntries(basecuts) << " out of: " << tD->GetEntries(sidebandscuts) << endl;
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	
	
	Double_t chi2 = 0;
	Double_t pvalue = 0;
	Int_t ndf  = 0;
	
	_INF(1,"");
	
	TH1F* mS = new TH1F("mS", ";m_{3body} [MeV];Events", (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize),m3bodyMin,m3bodyMax);
	mS->Sumw2();
	tS->Draw("m3body>>mS",cuts_sig_SB0);
	// TF1* fSig = new TF1("Signal0",fM3bodySignal,m3bodyMin,m3bodyMax,nparsig);
	// fSig->SetParName(0,"Norm_Gaus1");                                    fSig->SetParLimits(0,0,1e4);
	// fSig->SetParName(1,"Mean_Gaus1");   fSig->SetParameter(1,1776.82);
	// fSig->SetParName(2,"Sigma_Gaus1");                                   fSig->SetParLimits(2,1,100);
	// // fSig->SetParName(3,"Norm_Gaus2");                                    fSig->SetParLimits(3,0,1e4);
	// // fSig->SetParName(4,"Mean_Gaus2");   fSig->SetParameter(4,1776.82); 
	// // fSig->SetParName(5,"Sigma_Gaus2");  fSig->SetParameter(5,10);        fSig->SetParLimits(5,0,100);
	// fSig->SetLineColor(kRed);
	// mS->Fit("Signal0","L0 EMR");
	// chi2   = fSig->GetChisquare();
	// ndf    = fSig->GetNDF();
	// pvalue = TMath::Prob(chi2,ndf);
	// cout << "Signal   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	// Double_t avg1  = fSig->GetParameter(1);
	// Double_t davg1 = fSig->GetParError(1);
	// Double_t sigma1  = fSig->GetParameter(2);
	// Double_t dsigma1 = fSig->GetParError(2);
	// // Double_t avg2  = fSig->GetParameter(4);
	// // Double_t davg2 = fSig->GetParError(4);
	// // Double_t sigma2  = fSig->GetParameter(5);
	// // Double_t dsigma2 = fSig->GetParError(5);
	// cout << "Gauss1 Average:  " << avg1   << " +- " << davg1 << endl;
	// // cout << "Gauss2 Average:  " << avg2   << " +- " << davg2 << endl;
	// cout << "Gauss1 Sigma:    " << sigma1 << " +- " << dsigma1 << endl;
	// // cout << "Gauss2 Sigma:    " << sigma2 << " +- " << dsigma2 << endl;
	// // Double_t ximin = (avg1+avg2)/2.-2.*(sigma1+sigma2)/2.;
	// // Double_t ximax = (avg1+avg2)/2.+2.*(sigma1+sigma2)/2.;
	// Double_t ximin = avg1-2.*sigma1;
	// Double_t ximax = avg1+2.*sigma1;
	
	mS->Fit("gaus","0");
	TF1* fGauss = mS->GetFunction("gaus");
	chi2 = fGauss->GetChisquare();
	ndf = fGauss->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	Double_t avg  = fGauss->GetParameter(1);
	Double_t davg = fGauss->GetParError(1);
	Double_t sigma  = fGauss->GetParameter(2);
	Double_t dsigma = fGauss->GetParError(2);
	cout << "Gauss chi2/ndf: " << chi2/ndf  << endl;
	cout << "Gauss pvalue:   " << pvalue  << endl;
	cout << "Gauss Average:  " << avg   << " +- " << davg << endl;
	cout << "Gauss Sigma:    " << sigma << " +- " << dsigma << endl;
	Double_t ximin = avg-2.*sigma;
	Double_t ximax = avg+2.*sigma;
	
	TPaveText *ptGaus = new TPaveText(0.6,0.5,0.85,0.75,"NDC");
	ptGaus->SetFillStyle(4000); //will be transparent
	ptGaus->SetFillColor(0);
	ptGaus->SetTextFont(42);
	ptGaus->SetBorderSize(0);
	ptGaus->AddText("#sigma: "+tstr(sigma,2)+"#pm"+tstr(dsigma,2)+" MeV");
	// ptGaus->AddText("#sigma_{1}: "+tstr(sigma1,2)+"#pm"+tstr(dsigma1,2)+" MeV");
	// ptGaus->AddText("#sigma_{2}: "+tstr(sigma1,2)+"#pm"+tstr(dsigma2,2)+" MeV");
	ptGaus->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptGaus->AddText("#it{p}-value: "+tstr(pvalue,5));
	ptGaus->AddText("Integral: m_{3body} #in["+tstr(ximin,1)+","+tstr(ximax,1)+"]");
	
	TLine* lmin = new TLine(ximin,mS->GetMaximum()/2.,ximin,0.);
	TLine* lmax = new TLine(ximax,mS->GetMaximum()/2.,ximax,0.);
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	mS->Draw("p");
	fGauss->Draw("same");
	// fSig->Draw("same");
	lmin->Draw("same");
	lmax->Draw("same");
	ptGaus->Draw("same");
	cnv->Update();
	cnv->SaveAs("figures/bdtfit.ROOT.Gaus.png");
	cnv->SaveAs("figures/bdtfit.ROOT.Gaus.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf(");
	
	_INF(1,"");
	
	
	TH1F* mD0 = new TH1F("mD0", ";m_{3body} [MeV];Events", (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize),m3bodyMin,m3bodyMax);
	mD0->SetLineColor(kBlack);
	mD0->SetLineWidth(1);
	mD0->SetLineStyle(1);
	mD0->SetMarkerStyle(20);
	mD0->SetMarkerSize(0.6);
	mD0->SetMarkerColor(kBlack);
	mD0->Sumw2();
	
	TH1F* mD1 = new TH1F("mD1", ";m_{3body} [MeV];Events", (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize),m3bodyMin,m3bodyMax);
	mD1->SetLineColor(kBlack);
	mD1->SetLineWidth(1);
	mD1->SetLineStyle(1);
	mD1->SetMarkerStyle(20);
	mD1->SetMarkerSize(0.6);
	mD1->SetMarkerColor(kBlack);
	mD1->Sumw2();
	
	TH1F* BDTb = new TH1F("BDTb",";BDT score;Events",nbdtbins,xbdtmin,xbdtmax);
	BDTb->SetLineColor(kBlack);
	BDTb->SetLineWidth(1);
	BDTb->SetLineStyle(1);
	BDTb->SetMarkerStyle(20);
	BDTb->SetMarkerSize(0.6);
	BDTb->SetMarkerColor(kBlack);
	BDTb->Sumw2();


	_INF(1,"");

	
	tS->SetEventList(0);
	tS->Draw(">>elistSopt",cuts_sig_SB1);
	TEventList* elistSopt = (TEventList*)gDirectory->Get("elistSopt");
	Float_t npassedSopt = elistSopt->GetN();

	_INF(1,"");

	tD->SetEventList(0);
	tD->Draw(">>elistDopt",cuts_bkg_SB1);
	TEventList* elistDopt = (TEventList*)gDirectory->Get("elistDopt");
	Float_t npassedDopt = elistDopt->GetN();
	
	_INF(1,"");
	
	tS->SetEventList(0);
	tS->Draw(">>elistS",cuts_sig_SB0);
	TEventList* elistS = (TEventList*)gDirectory->Get("elistS");
	Float_t npassedS = elistS->GetN();

	_INF(1,"");

	tD->SetEventList(0);
	tD->Draw(">>elistD",cuts_bkg_SB0);
	TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
	Float_t npassedD = elistD->GetN();
	
	_INF(1,"");
	
	
	tD->SetEventList(elistD);
	tD->Draw("score>>BDTb");
	BDTb->SetMinimum(0.);
	
	TF1* tf1BDTdummy = new TF1("fBDTdummy","([0]+[1]*TMath::Exp(x*[2])+[3]*TMath::Power(TMath::Abs(x),[4]*x)*TMath::Exp(x*[5]))",xbdtmin,xbdtmax);
	tf1BDTdummy->SetParameter(0,1.6001e-03);
	tf1BDTdummy->SetParameter(1,9.0413e-01);
	tf1BDTdummy->SetParameter(2,-7.0029e+00);
	tf1BDTdummy->SetParameter(3,4.7542e+01);
	tf1BDTdummy->SetParameter(4,5.5042e-01);
	tf1BDTdummy->SetParameter(5,-4.4837e-01);
	float dummyIntegral = tf1BDTdummy->Integral(xbdtmin,xbdtmax)/BDTb->GetBinWidth(1);
	TF1* tf1BDTub = new TF1("fBDTub","([0]+[1]*TMath::Exp(x*[2])+[3]*TMath::Power(TMath::Abs(x),[4]*x)*TMath::Exp(x*[5]))*[6]",xbdtmin,xbdtmax);
	tf1BDTub->SetParameter(0,1.6001e-03);
	tf1BDTub->SetParameter(1,9.0413e-01);
	tf1BDTub->SetParameter(2,-7.0029e+00);
	tf1BDTub->SetParameter(3,4.7542e+01);
	tf1BDTub->SetParameter(4,5.5042e-01);
	tf1BDTub->SetParameter(5,-4.4837e-01);
	tf1BDTub->SetParameter(6,269/dummyIntegral);
	tf1BDTub->SetLineColor(kGreen+1);
	
	TF1* tf1Bkg1 = new TF1("BackgroundBDT1",fBkgBDT1,xbdtmin,xbdtmax,nparbdt1);
	TF1* tf1Bkg2 = new TF1("BackgroundBDT2",fBkgBDT2,xbdtmin,xbdtmax,nparbdt2);
	BDTb->Fit("BackgroundBDT1","L0 EMR");
	// tD->SetEventList(0);
	// // tD->UnbinnedFit("BackgroundBDT1", "score", cuts_bkg_SB0, "EMV");
	// tD->UnbinnedFit("BackgroundBDT2", "score", cuts_bkg_SB0, "EMV");
	chi2   = tf1Bkg1->GetChisquare();
	ndf    = tf1Bkg1->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "BDT1:  fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	BDTb->Fit("BackgroundBDT2","L0 EMR");
	chi2   = tf1Bkg2->GetChisquare();
	ndf    = tf1Bkg2->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "BDT12  fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fBDTb1 = new TF1("fBDTb1",fBkgBDT1,xbdtmin,xbdtmax,nparbdt1); fBDTb1->SetLineColor(kRed);
	TF1 *fBDTb2 = new TF1("fBDTb2",fBkgBDT2,xbdtmin,xbdtmax,nparbdt2); fBDTb2->SetLineColor(kBlue);
	fBDTb1->SetParameters(tf1Bkg1->GetParameters());
	fBDTb2->SetParameters(tf1Bkg2->GetParameters());
	
	
	
	BDTb->GetListOfFunctions()->Add(fBDTb1);
	BDTb->GetListOfFunctions()->Add(fBDTb2);
	// float iBDT1_0 = fBDTb1->Integral(xbdtmin,xbdtopt)/BDTb->GetBinWidth(1);
	float iBDT1_0 = fBDTb1->Integral(xbdtmin,xbdtmax)/BDTb->GetBinWidth(1);
	float iBDT1_1 = fBDTb1->Integral(xbdtopt,xbdtmax)/BDTb->GetBinWidth(1);
	float transferFactor_1 = iBDT1_1/iBDT1_0;
	// float iBDT2_0 = fBDTb2->Integral(xbdtmin,xbdtopt)/BDTb->GetBinWidth(1);
	float iBDT2_0 = fBDTb2->Integral(xbdtmin,xbdtmax)/BDTb->GetBinWidth(1);
	float iBDT2_1 = fBDTb2->Integral(xbdtopt,xbdtmax)/BDTb->GetBinWidth(1);
	float transferFactor_2 = iBDT2_1/iBDT2_0;
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	BDTb->Draw("p0 e0");
	tf1BDTub->Draw("same");
	cnv->SaveAs("figures/bdtfit.ROOT.BDT.png");
	cnv->SaveAs("figures/bdtfit.ROOT.BDT.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	
	
	tD->SetEventList(0);
	tD->SetEventList(elistD);
	tD->Draw("m3body>>mD0");
	for(int b=1 ; b<=mD0->GetNbinsX() ; ++b) { if(mD0->GetBinContent(b)<1) mD0->SetBinError(b,1); }
	for(int b=1 ; b<=mD0->GetNbinsX() ; ++b) { if(mD0->GetBinContent(b)<1 && mD0->GetXaxis()->GetBinLowEdge(b)>=xbmin && mD0->GetXaxis()->GetBinUpEdge(b)<=xbmax) mD0->SetBinContent(b,-10); }
	mD0->SetMinimum(0.);
	if(mD0->GetMaximum()<3.) mD0->SetMaximum(3.);

	//// first SB0 function (nominal)
	TF1* fbkg0 = new TF1("Background0",fBkg0,m3bodyMin,m3bodyMax, npar0);
	fbkg0->SetParName(0,"Nconst");  /*fbkg0->SetParameter(0,...);*/   fbkg0->SetParLimits(0,0,1000);
	fbkg0->SetParName(1,"Nexpo");   /*fbkg0->SetParameter(1,...);*/   fbkg0->SetParLimits(1,0,1000);
	fbkg0->SetParName(2,"Gexpo");   /*fbkg0->SetParameter(2,...);*/   fbkg0->SetParLimits(2,-10,0);
	// fbkg0->SetParName(3,"Ngaus");   /*fbkg0->SetParameter(3,10);*/    fbkg0->SetParLimits(3,0,50);
	// fbkg0->SetParName(4,"Mean");    /*fbkg0->SetParameter(4,1600);*/  fbkg0->SetParLimits(4,1520,1680);
	// fbkg0->SetParName(5,"Sigma");   /*fbkg0->SetParameter(5,20.);*/   fbkg0->SetParLimits(5,1,100);
	fbkg0->SetLineColor(kBlack);
	reject = kTRUE;
	mD0->Fit("Background0","L0 EMR");
	reject = kFALSE;
	chi2   = fbkg0->GetChisquare();
	ndf    = fbkg0->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR0 = fbkg0->Integral(ximin,ximax)/mD0->GetBinWidth(1);
	float iBkgSB0 = (fbkg0->Integral(m3bodyMin,xbmin)+fbkg0->Integral(xbmax,m3bodyMax))/mD0->GetBinWidth(1);
	cout << "Sidebands0   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fleft0 = new TF1("fleft0",fBkg0,m3bodyMin,xbmin,npar0); fleft0->SetLineColor(kBlack);
	fleft0->SetParameters(fbkg0->GetParameters());
	// mD0->GetListOfFunctions()->Add(fleft0);
	// gROOT->GetListOfFunctions()->Remove(fleft0);
	TF1 *fright0 = new TF1("fright0",fBkg0,xbmax,m3bodyMax,npar0); fright0->SetLineColor(kBlack);
	fright0->SetParameters(fbkg0->GetParameters());
	// mD0->GetListOfFunctions()->Add(fright0);
	// gROOT->GetListOfFunctions()->Remove(fright0);


	//// first SB0 function
	TF1* fbgCR0_1 = new TF1("BackgroundCR0_1",fBkgMassCR0_1,m3bodyMin,m3bodyMax,nparmass1);
	fbgCR0_1->SetLineColor(kRed);
	fbgCR0_1->SetLineStyle(2);
	reject = kTRUE;
	mD0->Fit("BackgroundCR0_1","L0 EMR");
	reject = kFALSE;
	chi2   = fbgCR0_1->GetChisquare();
	ndf    = fbgCR0_1->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR0_1 = fbgCR0_1->Integral(ximin,ximax)/mD0->GetBinWidth(1);
	float iBkgSB0_1 = (fbgCR0_1->Integral(m3bodyMin,xbmin)+fbgCR0_1->Integral(xbmax,m3bodyMax))/mD0->GetBinWidth(1);
	cout << "Sidebands0   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fleftCR0_1 = new TF1("fleftCR0_1",fBkgMassCR0_1,m3bodyMin,xbmin,nparmass1); fleftCR0_1->SetLineColor(kRed);
	fleftCR0_1->SetParameters(fbgCR0_1->GetParameters());
	// mD0->GetListOfFunctions()->Add(fleftCR0_1);
	// gROOT->GetListOfFunctions()->Remove(fleftCR0_1);
	TF1 *frightCR0_1 = new TF1("frightCR0_1",fBkgMassCR0_1,xbmax,m3bodyMax,nparmass1); frightCR0_1->SetLineColor(kRed);
	frightCR0_1->SetParameters(fbgCR0_1->GetParameters());
	// mD0->GetListOfFunctions()->Add(frightCR0_1);
	// gROOT->GetListOfFunctions()->Remove(frightCR0_1);
	
	//// second SB0 function
	TF1* fbgCR0_2 = new TF1("BackgroundCR0_2",fBkgMassCR0_2,m3bodyMin,m3bodyMax,nparmass2);
	fbgCR0_2->SetLineColor(kBlue);
	fbgCR0_2->SetLineStyle(2);
	reject = kTRUE;
	mD0->Fit("BackgroundCR0_2","L0 EMR");
	reject = kFALSE;
	chi2   = fbgCR0_2->GetChisquare();
	ndf    = fbgCR0_2->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR0_2 = fbgCR0_2->Integral(ximin,ximax)/mD0->GetBinWidth(1);
	float iBkgSB0_2 = (fbgCR0_2->Integral(m3bodyMin,xbmin)+fbgCR0_2->Integral(xbmax,m3bodyMax))/mD0->GetBinWidth(1);
	cout << "Sidebands0   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fleftCR0_2 = new TF1("fleftCR0_2",fBkgMassCR0_2,m3bodyMin,xbmin,nparmass2); fleftCR0_2->SetLineColor(kBlue);
	fleftCR0_2->SetParameters(fbgCR0_2->GetParameters());
	// mD0->GetListOfFunctions()->Add(fleftCR0_2);
	// gROOT->GetListOfFunctions()->Remove(fleftCR0_2);
	TF1 *frightCR0_2 = new TF1("frightCR0_2",fBkgMassCR0_2,xbmax,m3bodyMax,nparmass2); frightCR0_2->SetLineColor(kBlue);
	frightCR0_2->SetParameters(fbgCR0_2->GetParameters());
	// mD0->GetListOfFunctions()->Add(frightCR0_2);
	// gROOT->GetListOfFunctions()->Remove(frightCR0_2);
	
	//// third SB0 function
	TF1* fbgCR0_3 = new TF1("BackgroundCR0_3",fBkgMassCR0_3,m3bodyMin,m3bodyMax,nparmass3);
	fbgCR0_3->SetLineColor(kGreen);
	fbgCR0_3->SetLineStyle(2);
	reject = kTRUE;
	mD0->Fit("BackgroundCR0_3","L0 EMR");
	reject = kFALSE;
	chi2   = fbgCR0_3->GetChisquare();
	ndf    = fbgCR0_3->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR0_3 = fbgCR0_3->Integral(ximin,ximax)/mD0->GetBinWidth(1);
	float iBkgSB0_3 = (fbgCR0_3->Integral(m3bodyMin,xbmin)+fbgCR0_3->Integral(xbmax,m3bodyMax))/mD0->GetBinWidth(1);
	cout << "Sidebands0   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fleftCR0_3 = new TF1("fleftCR0_3",fBkgMassCR0_3,m3bodyMin,xbmin,nparmass3); fleftCR0_3->SetLineColor(kGreen);
	fleftCR0_3->SetParameters(fbgCR0_3->GetParameters());
	// mD0->GetListOfFunctions()->Add(fleftCR0_3);
	// gROOT->GetListOfFunctions()->Remove(fleftCR0_3);
	TF1 *frightCR0_3 = new TF1("frightCR0_3",fBkgMassCR0_3,xbmax,m3bodyMax,nparmass3); frightCR0_3->SetLineColor(kGreen);
	frightCR0_3->SetParameters(fbgCR0_3->GetParameters());
	// mD0->GetListOfFunctions()->Add(frightCR0_3);
	// gROOT->GetListOfFunctions()->Remove(frightCR0_3);
	
	//// fourth SB0 function
	TF1* fbgCR0_4 = new TF1("BackgroundCR0_4",fBkgMassCR0_4,m3bodyMin,m3bodyMax,nparmass4);
	fbgCR0_4->SetLineColor(kOrange);
	fbgCR0_4->SetLineStyle(2);
	reject = kTRUE;
	mD0->Fit("BackgroundCR0_4","L0 EMR");
	reject = kFALSE;
	chi2   = fbgCR0_4->GetChisquare();
	ndf    = fbgCR0_4->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR0_4 = fbgCR0_4->Integral(ximin,ximax)/mD0->GetBinWidth(1);
	float iBkgSB0_4 = (fbgCR0_4->Integral(m3bodyMin,xbmin)+fbgCR0_4->Integral(xbmax,m3bodyMax))/mD0->GetBinWidth(1);
	cout << "Sidebands0   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TF1 *fleftCR0_4 = new TF1("fleftCR0_4",fBkgMassCR0_4,m3bodyMin,xbmin,nparmass4); fleftCR0_4->SetLineColor(kOrange);
	fleftCR0_4->SetParameters(fbgCR0_4->GetParameters());
	// mD0->GetListOfFunctions()->Add(fleftCR0_4);
	// gROOT->GetListOfFunctions()->Remove(fleftCR0_4);
	TF1 *frightCR0_4 = new TF1("frightCR0_4",fBkgMassCR0_4,xbmax,m3bodyMax,nparmass4); frightCR0_4->SetLineColor(kOrange);
	frightCR0_4->SetParameters(fbgCR0_4->GetParameters());
	// mD0->GetListOfFunctions()->Add(frightCR0_4);
	// gROOT->GetListOfFunctions()->Remove(frightCR0_4);
	
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	mD0->Draw("p0 e0");
	fbkg0->Draw("same");
	fbgCR0_1->Draw("same");
	fbgCR0_2->Draw("same");
	fbgCR0_3->Draw("same");
	fbgCR0_4->Draw("same");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:EventNumber","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.events.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.events.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:mOS1","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.m3body_vs_mOS1.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.m3body_vs_mOS1.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:mOS2","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.m3body_vs_mOS2.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB0.m3body_vs_mOS2.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	


	
	
	tD->SetEventList(0);
	tD->SetEventList(elistDopt);
	tD->Draw("m3body>>mD1");
	for(int b=1 ; b<=mD1->GetNbinsX() ; ++b) { if(mD1->GetBinContent(b)<1) mD1->SetBinError(b,1); }
	for(int b=1 ; b<=mD1->GetNbinsX() ; ++b) { if(mD1->GetBinContent(b)<1 && mD1->GetXaxis()->GetBinLowEdge(b)>=xbmin && mD1->GetXaxis()->GetBinUpEdge(b)<=xbmax) mD1->SetBinContent(b,-10); }
	mD1->SetMinimum(0.);
	if(mD1->GetMaximum()<3.) mD1->SetMaximum(3.);
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:EventNumber","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.events.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.events.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:mOS1","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.m3body_vs_mOS2.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.m3body_vs_mOS2.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	tD->Draw("m3body:mOS2","","box text0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.m3body_vs_mOS2.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.m3body_vs_mOS2.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf");
	
	
	
	
	// TF1* fbgCR1 = new TF1("BackgroundCR1",fBkgMassCR1,m3bodyMin,m3bodyMax,nparmass);
	TF1* fbkg1 = new TF1("BackgroundCR1",fBkg1,m3bodyMin,m3bodyMax,npar1);
	fbkg1->SetParName(0,"Nconst");  /*fbkg1->SetParameter(0,...);*/   fbkg1->SetParLimits(0,0,1000);
	fbkg1->SetParName(1,"Nexpo");   /*fbkg1->SetParameter(1,...);*/   fbkg1->SetParLimits(1,0,1000);
	fbkg1->SetParName(2,"Gexpo");   /*fbkg1->SetParameter(2,...);*/   fbkg1->SetParLimits(2,-10,0);
	reject = kTRUE;
	mD1->Fit("BackgroundCR1","L0 EMR");
	reject = kFALSE;
	chi2   = fbkg1->GetChisquare();
	ndf    = fbkg1->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	float iBkgSR1 = fbkg1->Integral(ximin,ximax)/mD1->GetBinWidth(1);
	float iBkgSB1 = (fbkg1->Integral(m3bodyMin,xbmin)+fbkg1->Integral(xbmax,m3bodyMax))/mD1->GetBinWidth(1);
	cout << "Sidebands1   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	// TF1 *fleft1 = new TF1("fleft1",fBkgMassCR1,m3bodyMin,xbmin,nparmass);
	TF1 *fleft1 = new TF1("fleft1",fbkg1,m3bodyMin,xbmin,npar1);
	fleft1->SetParameters(fbkg1->GetParameters());
	mD1->GetListOfFunctions()->Add(fleft1);
	gROOT->GetListOfFunctions()->Remove(fleft1);
	// TF1 *fright1 = new TF1("fright1",fBkgMassCR1,xbmax,m3bodyMax,nparmass);
	TF1 *fright1 = new TF1("fright1",fbkg1,xbmax,m3bodyMax,npar1);
	fright1->SetParameters(fbkg1->GetParameters());
	mD1->GetListOfFunctions()->Add(fright1);
	gROOT->GetListOfFunctions()->Remove(fright1);
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	mD1->Draw("p0 e0");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.png");
	cnv->SaveAs("figures/bdtfit.ROOT.SB1.eps");
	cnv->SaveAs("figures/bdtfit.ROOT.pdf)");
	
	
	
	
	
	
	
	
	vector<float> vF01;
	vector<float> vNSB0;
	vector<float> vDNSB0;
	vector<float> vNSR0;
	vector<float> vDNSR0;
	vector<float> vNSR1;
	vector<float> vDNSR1;
	
	float transferFactorNominal = transferFactor_2;
	
	vF01.push_back(transferFactorNominal);
	vF01.push_back(transferFactor_1);
	vF01.push_back(transferFactor_2);
	
	vNSB0.push_back(iBkgSB0);
	vNSB0.push_back(iBkgSB0_1);
	vNSB0.push_back(iBkgSB0_2);
	vNSB0.push_back(iBkgSB0_3);
	vNSB0.push_back(iBkgSB0_4);
	
	vDNSB0.push_back(sqrt(iBkgSB0));
	vDNSB0.push_back(sqrt(iBkgSB0_1));
	vDNSB0.push_back(sqrt(iBkgSB0_2));
	vDNSB0.push_back(sqrt(iBkgSB0_3));
	vDNSB0.push_back(sqrt(iBkgSB0_4));
	
	vNSR0.push_back(iBkgSR0);
	vNSR0.push_back(iBkgSR0_1);
	vNSR0.push_back(iBkgSR0_2);
	vNSR0.push_back(iBkgSR0_3);
	vNSR0.push_back(iBkgSR0_4);
	
	vDNSR0.push_back((iBkgSR0/iBkgSB0)*sqrt(iBkgSB0));
	vDNSR0.push_back((iBkgSR0_1/iBkgSB0_1)*sqrt(iBkgSB0_1));
	vDNSR0.push_back((iBkgSR0_2/iBkgSB0_2)*sqrt(iBkgSB0_2));
	vDNSR0.push_back((iBkgSR0_3/iBkgSB0_3)*sqrt(iBkgSB0_3));
	vDNSR0.push_back((iBkgSR0_4/iBkgSB0_4)*sqrt(iBkgSB0_4));
	
	vNSR1.push_back(iBkgSR0*transferFactorNominal);
	vNSR1.push_back(iBkgSR0_1*transferFactorNominal);
	vNSR1.push_back(iBkgSR0_2*transferFactorNominal);
	vNSR1.push_back(iBkgSR0_3*transferFactorNominal);
	vNSR1.push_back(iBkgSR0_4*transferFactorNominal);
	
	vDNSR1.push_back((iBkgSR0/iBkgSB0)*sqrt(iBkgSB0)*transferFactorNominal);
	vDNSR1.push_back((iBkgSR0_1/iBkgSB0_1)*sqrt(iBkgSB0_1)*transferFactorNominal);
	vDNSR1.push_back((iBkgSR0_2/iBkgSB0_2)*sqrt(iBkgSB0_2)*transferFactorNominal);
	vDNSR1.push_back((iBkgSR0_3/iBkgSB0_3)*sqrt(iBkgSB0_3)*transferFactorNominal);
	vDNSR1.push_back((iBkgSR0_4/iBkgSB0_4)*sqrt(iBkgSB0_4)*transferFactorNominal);
	
	
	////////////
	//// summary
	cout << "========================================================================================" << endl;
	cout << "SR                       " << ximin << "->" << ximax << " MeV" << endl;
	cout << "Minimal BDT cut (x0)     " << xbdtmin << endl;
	cout << "Optimal BDT cut (x1)     " << xbdtopt << endl;
	cout << "Maximal BDT cut (x2)     " << xbdtmax << endl;
	cout << "TF1 normalization        " << ncandidates/dummyIntegral << endl;
	cout << "Transfer factor (f01)    "; for(unsigned int i=0 ; i<vF01.size() ; ++i) cout << tstr(vF01[i],5) << " | "; cout << endl;
	cout << "N events in SB0 (tree)   " << tstr(ncandidates,0) << endl;
	cout << "N events in SB0 (BDTfit) " << iBDT1_0 << " | " << iBDT2_0 << endl;
	cout << "N events in SB0          "; for(unsigned int i=0 ; i<vNSB0.size() ; ++i)  cout << tstr(vNSB0[i],2)  << " | "; cout << endl;
	cout << "N events in SR0          "; for(unsigned int i=0 ; i<vNSR0.size() ; ++i)  cout << tstr(vNSR0[i],2)  << " | "; cout << endl; // << tstr(iBkgSR0_1,2) << endl;
	cout << "Stat error in SR0        "; for(unsigned int i=0 ; i<vDNSR0.size() ; ++i) cout << tstr(vDNSR0[i],2) << " | "; cout << endl; // << tstr((iBkgSR0_1/iBkgSB0_1)*sqrt(iBkgSB0_1),2) << endl;
	cout << "N events in SB1 (cut)    " << npassedDopt << endl;
	cout << "N events in SB1 (fit)    " << tstr(iBkgSB1,3) << endl;
	cout << "N events in SR1 (cut)    " << tstr(iBkgSR1,3) << endl;
	cout << "N events in SR1 (fit)    "; for(unsigned int i=0 ; i<vNSR1.size() ; ++i) cout << tstr(vNSR1[i],2) << " | "; cout << endl; // << tstr(iBkgSR0_1*transferFactor_1,3) << endl;
	cout << "Stat error in SR1 (cut)  " << tstr((iBkgSR1/iBkgSB1)*sqrt(iBkgSB1),3) << " (dSR1/SR1="<< tstr((iBkgSR1/iBkgSB1)*sqrt(iBkgSB1)/iBkgSR1*100,0) <<"\%)" << endl;
	cout << "Stat error in SR1 (fit)  "; for(unsigned int i=0 ; i<vDNSR1.size() ; ++i) cout << tstr(vDNSR1[i],2) << " | "; cout << endl; //  << tstr((iBkgSR0_1/iBkgSB0_1)*sqrt(iBkgSB0_1)*transferFactor_1,3) << " (i.e. dSR0*f01)" << " (dSR1/SR1="<< tstr((iBkgSR0_1/iBkgSB0_1)*sqrt(iBkgSB0_1)*transferFactor_1/(iBkgSR0_1*transferFactor_1)*100,0) <<"\%)" << endl;
	cout << "Total Acc*Eff in SR1     " << tstr(npassedSopt/ninitS*100,2) << "\%" << endl;
	cout << "========================================================================================" << endl;
}
