/////////////////////////////////////////////////////////
//// root -b -l -q m3muFit.C++\(\"muons\",\"cuts\"\) ////
/////////////////////////////////////////////////////////

#include "std.h"

Bool_t reject;
Double_t xmin;
Double_t xmax;
Double_t xbmin;
Double_t xbmax;
Int_t npar;
TString hname;

Double_t fBg(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	// return (hname.Contains("lin_zoom")) ? (par[0]+x*par[1])/exp(par[2]*x) : (par[0]+x*par[1]);
	return par[0];
}

void fit3body(TString master, TString method)
{
	cout << "master = " << master << endl;
	cout << "method = " << method << endl;
	
	TString figures = "figures/"+method;
	hname = "_m3body_sigregion";
	TString hnamesigregion = "_m3body_sigregion";
	
	Double_t ximin = -1.;
	Double_t ximax = -1.;
	xmin  = (hname.Contains("lin_zoom")) ? 500.  : 1300.;
	xmax  = (hname.Contains("lin_zoom")) ? 3500. : 2300.;
	xbmin = (hname.Contains("lin_zoom")) ? 1600. : 1650.;
	xbmax = (hname.Contains("lin_zoom")) ? 1940. : 1920.;
	// npar  = (hname.Contains("lin_zoom")) ? 3     : 2;
	npar  = 1;
	
	// TFile* f = new TFile(figures+".m3bodyFit.root","READ");
	TFile* f = new TFile("histos."+master+"."+method+".root","READ");
	TH1D* hData      = (TH1D*)f->Get("Data"+hname);
	TH1D* hSignal    = (TH1D*)f->Get("bbTotau10_3mu"+hname);
	TH1D* hSignal_bb = (TH1D*)f->Get("bbTotau10_3mu"+hname);
	TH1D* hSignal_cc = (TH1D*)f->Get("ccTotau10_3mu"+hname);
	TH1D* hSignal_W  = (TH1D*)f->Get("Wtaunu_3mu"+hname);
	TH1D* hSignal_W_sigregion  = (TH1D*)f->Get("Wtaunu_3mu"+hnamesigregion);
	hSignal_W_sigregion->Fit("gaus","0");
	TF1* fGauss = hSignal_W_sigregion->GetFunction("gaus");
	Double_t chi2 = fGauss->GetChisquare();
	Int_t    ndf = fGauss->GetNDF();
	Double_t pvalue = TMath::Prob(chi2,ndf);
	Double_t N  = fGauss->GetParameter(0);
	Double_t dN = fGauss->GetParError(0);
	Double_t avg  = fGauss->GetParameter(1);
	Double_t davg = fGauss->GetParError(1);
	Double_t sigma  = fGauss->GetParameter(2);
	Double_t dsigma = fGauss->GetParError(2);
	cout << "Gauss chi2/ndf: " << chi2/ndf  << endl;
	cout << "Gauss pvalue:   " << pvalue  << endl;
	cout << "Gauss Average:  " << avg   << " +- " << davg << endl;
	cout << "Gauss Sigma:    " << sigma << " +- " << dsigma << endl;
	
	ximin = avg-2.*sigma;
	ximax = avg+2.*sigma;
	TLine* lmin = new TLine(ximin,hSignal_W_sigregion->GetMaximum()/2.,ximin,0.);
	TLine* lmax = new TLine(ximax,hSignal_W_sigregion->GetMaximum()/2.,ximax,0.);
		
	
	hSignal->Reset();
	hSignal->Add(hSignal_bb);
	hSignal->Add(hSignal_cc);
	hSignal->Add(hSignal_W);
	////////////////////////////////////////////////////////////////
	// (2.1e-8 * 5e2) * (1/5) = (1e-8 * 1e2) = 2.1e-6 = Br*1e2 ////
	hSignal->Scale(1./5.); /////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	TH1F* hSignalBg = (TH1F*)hSignal->Clone();
	hSignal->SetLineColor(kBlue);
	hSignal->SetFillColor(kBlue);
	hSignalBg->SetFillStyle(0);
	hSignalBg->SetFillColor(0);
	
	TF1 *fbg = new TF1("Background", fBg, xmin,xmax, npar);
	reject = kTRUE;
	hData->Fit("Background","L0");
	reject = kFALSE;
	//////////////////////////
	hSignalBg->Add(fbg); /////
	//////////////////////////
	chi2   = fbg->GetChisquare();
	ndf    = fbg->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	Double_t iBkg     = fbg->Integral(ximin,ximax)/hData->GetBinWidth(1);
	Double_t iSigNorm = fGauss->Integral(ximin,ximax)/hData->GetBinWidth(1);
	Double_t iSig     = iSigNorm/(2.1e-8 * 1.e4);
	cout << "Sidebands   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	cout << "Backgrounds integral(" << ximin << "," << ximax << "): " << iBkg << endl;
	cout << "Signals     integral(" << ximin << "," << ximax << "): " << iSig << " (for Br=1)" << endl;	
	
	
	// store 2 separate functions for visualization	
	TF1 *fleft = new TF1("fleft",fBg,xmin,xbmin,npar);
	fleft->SetParameters(fbg->GetParameters());
	hData->GetListOfFunctions()->Add(fleft);
	gROOT->GetListOfFunctions()->Remove(fleft);
	TF1 *fright = new TF1("fright",fBg,xbmax,xmax,npar);
	fright->SetParameters(fbg->GetParameters());
	hData->GetListOfFunctions()->Add(fright);
	gROOT->GetListOfFunctions()->Remove(fright);
	hData->SetLineWidth(2);
	hData->SetMarkerStyle(20);
	hData->SetMarkerSize(1);
	
	
	TLegend* leg = new TLegend(0.57,0.62,0.85,0.85,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(hSignal,"#SigmaSignals @ Br#times10^{2}","f");
	leg->AddEntry(hData,"Data @ 8.51 fb^{-1}","p");
	if(hname.Contains("lin_zoom")) leg->AddEntry(fbg,"Data fit (a_{0}+a_{1}x)/e^{a_{2}x}","f");
	else                           leg->AddEntry(fbg,"Data fit (flat)","f");
	leg->AddEntry(hSignalBg,"Fit+#SigmaSignals","f");
	
	
	
	TCanvas* cnv = NULL;
	
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	hSignal_W_sigregion->Draw("p");
	fGauss->Draw("same");
	lmin->Draw("same");
	lmax->Draw("same");
	cnv->Update();
	cnv->SaveAs(figures+".m3bodyFit.pdf(");
	cnv->SaveAs(figures+".m3bodyFit.sigGausFit.eps");
	
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	Double_t max = (hSignalBg->GetMaximum()>hData->GetMaximum()) ? hSignalBg->GetMaximum() : hData->GetMaximum();
	hSignalBg->SetMaximum(max*1.7);
	hSignalBg->SetMinimum(0.);
	hSignalBg->Draw("hist");
	hSignal->Draw("same hist");
	hData->Draw("same p");
	leg->Draw("same");
	cnv->Update();
	cnv->SaveAs(figures+".m3bodyFit.pdf)");
	cnv->SaveAs(figures+".m3bodyFit.sidBndFit.eps");
}
