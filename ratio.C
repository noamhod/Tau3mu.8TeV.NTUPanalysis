TCanvas* hratio(TString name, TH1* hNumerator, TH1* hDenominator,
				TLegend* leg=NULL, TLegend* legR=NULL, TString ratioLabel="Ratio",
				Double_t rmin=0., Double_t rmax=3.,
				Bool_t logx=false, Bool_t logy=false, TString drawopt="hist same",
				TGraphAsymmErrors* g=NULL, TGraphAsymmErrors* gr=NULL, TString errortype="")
{	
	TCanvas* cnv = new TCanvas(name,name,600,550);
	cnv->Divide(1,2);
	TVirtualPad* tvp_hists = cnv->cd(1);
	TVirtualPad* tvp_ratio = cnv->cd(2);
	
	if(logx) tvp_ratio->SetLogx();
	if(logx) tvp_hists->SetLogx();
	if(logy) tvp_hists->SetLogy();
	
	tvp_hists->SetPad(0.00, 0.35, 1.00, 1.00);
	tvp_ratio->SetPad(0.00, 0.00, 1.00, 0.355);

	tvp_hists->SetBottomMargin(0.012);
	tvp_ratio->SetBottomMargin(0.20);
	tvp_ratio->SetTopMargin(0.012);
	
	tvp_hists->SetTicks(1,1);
	tvp_ratio->SetTicks(1,1);
	
	TString cloneName_n = hNumerator->GetName();
	TString cloneName_d = hDenominator->GetName();
	TH1D* th1n_tmp = (TH1D*)hNumerator->Clone(cloneName_n+"_th1n_tmp");
	TH1D* th1d_tmp = (TH1D*)hDenominator->Clone(cloneName_d+"_th1d_tmp");
	

	TH1D* hr = (TH1D*)hNumerator->Clone();
	TString sXtitle = (TString)hNumerator->GetXaxis()->GetTitle();
	TString sTitle = ";"+sXtitle+";"+ratioLabel;
	hr->SetTitle(sTitle);
	
	if(errortype=="") hr->Divide(th1n_tmp,th1d_tmp,1.,1.);
	else
	{
		for(Int_t b=1 ; b<=hr->GetNbinsX() ; b++)
		{
			float v1  = th1n_tmp->GetBinContent(b);
			float v2  = th1d_tmp->GetBinContent(b);
			float d1  = th1n_tmp->GetBinError(b);
			float d2  = th1d_tmp->GetBinError(b);
			
			////////////////////////
			if(v2==0.) continue; ///
			////////////////////////
			
			// set the bin value
			hr->SetBinContent(b,(v1/v2));
			
			// set the bin error
			if     (errortype=="combined")    hr->SetBinError(b,(v1/v2)*TMath::Sqrt((d1/v1)*(d1/v1) + (d2/v2)*(d2/v2))); // combined error 
			else if(errortype=="numerator")   hr->SetBinError(b,(1./v2)*d1);      // numerator error only
			else if(errortype=="denominator") hr->SetBinError(b,(v1/(v2*v2))*d2); // denominator error only
			else { cout << "error..." << endl; return; }
		}
	}
	
	hr->SetMarkerStyle(20);
	hr->SetMarkerSize(0.8);
	hr->SetMarkerColor(kBlack);
	hr->SetLineColor(kBlack);
	hr->SetLineStyle(1);
	hr->SetLineWidth(2);
	Double_t xLabelSize = hNumerator->GetXaxis()->GetLabelSize()*1.5;
	Double_t yLabelSize = hNumerator->GetYaxis()->GetLabelSize()*1.5;
	Double_t xTitleSize = hNumerator->GetXaxis()->GetTitleSize()*1.5;
	Double_t yTitleSize = hNumerator->GetYaxis()->GetTitleSize()*1.5;
	Double_t titleSize = hNumerator->GetTitleSize()*1.5;
	hr->GetXaxis()->SetLabelSize(xLabelSize);
	hr->GetYaxis()->SetLabelSize(yLabelSize);
	hr->GetXaxis()->SetTitleSize(xTitleSize);
	hr->GetYaxis()->SetTitleSize(yTitleSize);
	hr->SetTitleSize(titleSize);
	hr->GetYaxis()->SetTitleOffset(0.5);
	hr->SetMinimum(rmin);
	hr->SetMaximum(rmax);
	
	if(logx) setlogx(hr);
	
	TLine* line = new TLine(hr->GetXaxis()->GetXmin(),1.,hr->GetXaxis()->GetXmax(),1.);

	tvp_hists->cd();
	setMinMax(th1d_tmp, th1n_tmp, true);
	th1d_tmp->Draw("hist");
	th1d_tmp->Draw("e1x0 SAMES");
	if(g!=NULL) g->Draw("3 SAMES");
	th1n_tmp->Draw(drawopt);
	// th1n_tmp->Draw("e1x0 SAMES");
	if(leg!=NULL) leg->Draw("SAMES");
	tvp_hists->Update();
	tvp_hists->RedrawAxis();

	tvp_ratio->cd();
	tvp_ratio->SetGridy();
	hr->Draw("epx0");
	if(gr!=NULL) gr->Draw("3 SAMES");
	line->Draw("SAMES");
	if(legR!=NULL) legR->Draw("SAMES");
	tvp_ratio->Update();
	tvp_ratio->RedrawAxis();
	
	cnv->Update();
	
	return cnv;
}
