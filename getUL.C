///////////////////////////
//// root -l getUL.C++ ////
///////////////////////////
#include "std.h"

void getUL()
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
	gStyle->SetTitleX(0.25); //title X location 
	gStyle->SetTitleY(0.96); //title Y location 
	gStyle->SetTitleW(0.5); //title width 
	gStyle->SetTitleH(0.05); //title height
	gStyle->SetTitleBorderSize(0);
	
	
	
	// Example macro of using the TFeldmanCousins class in root.
	// get a FeldmanCousins calculation object with the default limits
	// of calculating a 90% CL with the minimum signal value scanned
	// = 0.0 and the maximum signal value scanned of 50.0
	//Author : Adrian John Bevan <bevan@SLAC.Stanford.EDU>

	// if (!gROOT->GetClass("TFeldmanCousins")) gSystem->Load("libPhysics");
	// TFeldmanCousins fc;

	// // calculate either the upper or lower limit for 10 observerd
	// // events with an estimated background of 3.  The calculation of
	// // either upper or lower limit will return that limit and fill
	// // data members with both the upper and lower limit for you.
	// Double_t Nobserved1   = 10.0;
	// Double_t Nbackground1 = 3.0;
	// Double_t ul = fc.CalculateUpperLimit(Nobserved1, Nbackground1);
	// Double_t ll = fc.GetLowerLimit();
	// cout << "For " <<  Nobserved1 << " data observed with and estimated background"<<endl;
	// cout << "of " << Nbackground1 << " candidates, the Feldman-Cousins method of "<<endl;
	// cout << "calculating confidence limits gives:"<<endl;
	// cout << "\tUpper Limit = " <<  ul << endl;
	// cout << "\tLower Limit = " <<  ll << endl;
	// cout << "at the 90% CL"<< endl;
 
	


	TFile* fin  = new TFile("TMVA.root","READ");
	TFile* fout = new TFile("UL.root","RECREATE");
	TH1D* hEffS = (TH1D*)fin->Get("Method_BDT/BDTG/MVA_BDTG_trainingEffS");
	TH1D* hEffB = (TH1D*)fin->Get("Method_BDT/BDTG/MVA_BDTG_trainingEffB");
	TH1D* hUL   = (TH1D*)hEffS->Clone("UL optimization curve");
	hUL->Reset();
	
	hEffS->SetLineColor(kBlue);
	hEffB->SetLineColor(kRed);
	hEffS->SetLineStyle(2);
	hEffB->SetLineStyle(2);
	hEffS->SetTitle("");
	hEffB->SetTitle("");
	
	Double_t Nsidebands = 6000; // the number of sideband events befor MVA cut
	Double_t SRfraction = 0.82; // the fraction of events in the SR

	Double_t Ymin = +1.e20;
	Double_t Xmin = -999;
	for(Int_t i=1 ; i<=hUL->GetNbinsX() ; i++)
	{
		Double_t effS = hEffS->GetBinContent(i);
		Double_t effB = hEffB->GetBinContent(i); 

		Double_t n0   = effB;
		Double_t nUL  = TMath::Sqrt(n0);
		Double_t UL   = nUL/effS;
		
		// Double_t nSR = effB; // floor(effB*Nsidebands*SRfraction);
		// Double_t Nobserved   = nSR;
		// Double_t Nbackground = nSR;
		// Double_t nUL = fc.CalculateUpperLimit(Nobserved,Nbackground);
		// Double_t nLL = fc.GetLowerLimit();
		// // cout << "For nSR=" << nSR << ", UL=" << nUL << endl;
		// Double_t UL = nUL/effS;
		
		if(effS!=0.) hUL->SetBinContent(i,UL);

		if(UL>0. && UL<Ymin)
		{
			Xmin = hUL->GetBinCenter(i);
			Ymin = UL;
		}
	}
	fout->cd();
	hEffS->Write();
	hEffB->Write();
	hUL->Write();
	fout->Write();
	// fout->Close();
	fin->cd();
	
	
	TLegend* leg = new TLegend(0.55,0.5,0.83,0.73);
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(hEffS,"#it{#epsilon}_{S}","f");
	leg->AddEntry(hEffB,"#it{#epsilon}_{B}","f");
	leg->AddEntry(hUL,"#sqrt{#it{#epsilon}_{B}}/#it{#epsilon}_{S}","f");
	
	TCanvas* c1 = new TCanvas("cnv","",1200,1000);
	c1->cd();
	c1->Draw();
	
	TPad *p1 = new TPad("pad1","",0,0,1,1);
	TPad *p2 = new TPad("pad2","",0,0,1,1);
	p2->SetFillStyle(4000); //will be transparent
	p1->SetTicks(1,0);
	p2->SetTicks(1,0);
	// p1->SetGridx();
	
	p1->Draw();
	p1->cd();
	hUL->SetTitle("Optimization on Br(#it{#tau#rightarrow3#mu}) #propto #sqrt{#it{#epsilon}_{B}}/#it{#epsilon}_{S}");
	hUL->GetXaxis()->SetTitle("BDTG score cut");
	hUL->GetYaxis()->SetTitle("#sqrt{#it{#epsilon}_{B}}/#it{#epsilon}_{S} [arbitrary units]");
	hUL->DrawNormalized("hist");
	// leg->Draw("same");
	p1->Modified();
	c1->cd();
	
	//compute the pad range with suitable margins
	Double_t ymin = 0;
	Double_t ymax = 1;
	Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
	Double_t xmin = hUL->GetXaxis()->GetXmin();
	Double_t xmax = hUL->GetXaxis()->GetXmax();
	Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
	p2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
	p2->Draw();
	p2->cd();
	leg->Draw("same");
	hEffS->Draw("][hist sames");
	hEffB->Draw("][hist sames");

	// draw axis on the right side of the pad
	TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
	axis->SetTitle("Efficiency");
	axis->SetLabelColor(kBlack);
	axis->SetLabelFont(font);
	axis->SetLabelSize(tsize);
	axis->SetTextFont(font);
	axis->SetTextSize(tsize);
	axis->Draw();
	
	c1->cd();
	// c1->Update();
	c1->SaveAs("ToyULonBR.pdf");
	c1->SaveAs("ToyULonBR.eps");
	c1->SaveAs("ToyULonBR.png");
	c1->SaveAs("ToyULonBR.root");
	cout << "Optimal BDT score is " << Xmin << " for Sqrt(effB)/effS=" << Ymin << endl;
	cout << "and efficiency of " << hEffS->GetBinContent(hEffS->FindBin(Xmin)) << " and rejection of " << 1-hEffB->GetBinContent(hEffB->FindBin(Xmin)) << endl;
}
