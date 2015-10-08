//////////////////////////////////////////////////////////
//// root -b -l -q syst.C++\(\"Data\"\) or Wtaunu_3mu ////
//////////////////////////////////////////////////////////
#include "std.h"
#include "type.h"
#include "const.h"

vector<TFile*> files;
vector<TString> modes;
TMapTSP2TH1 histos1;
TMapTSP2TH2 histos2;

TPaveText* ptxt;
void makeAtlasLabel()
{	
	ptxt = new TPaveText(0.50,0.30,0.75,0.52,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

void addHistos(TString hname)
{
	cout << "adding " << hname << endl;
	for(unsigned int i=0 ; i<modes.size() ; ++i)
	{
		TString fullname = hname+"_"+modes[i];
		histos1.insert(make_pair(fullname,(TH1F*)files[i]->Get(hname)));
		if(modes[i]=="uncalib")
		{
			histos1[fullname]->SetMarkerColor(kGray+2);
			histos1[fullname]->SetMarkerStyle(22);
			histos1[fullname]->SetMarkerSize(0.8);
			histos1[fullname]->SetLineColor(kGray+2);
		}
		if(modes[i]=="calib")
		{
	   		histos1[fullname]->SetMarkerColor(kBlack);
			histos1[fullname]->SetMarkerStyle(20);
			histos1[fullname]->SetMarkerSize(0.8);
			histos1[fullname]->SetLineColor(kBlack);
		}
	}
}

void plot(TString hname, TString legpos, TString pdfname, TString pdfmode="", Int_t rebinFactor=2)
{
	TGraphAsymmErrors* gJES = new TGraphAsymmErrors(); gJES->GetXaxis()->SetTitle(histos1[hname+"_calib"]->GetXaxis()->GetTitle()); gJES->GetYaxis()->SetTitle("Normalized"); gJES->SetTitle(histos1[hname+"_calib"]->GetTitle()); 
	TGraphAsymmErrors* gJER = new TGraphAsymmErrors(); gJER->GetXaxis()->SetTitle(histos1[hname+"_calib"]->GetXaxis()->GetTitle()); gJER->GetYaxis()->SetTitle("Normalized"); gJER->SetTitle(histos1[hname+"_calib"]->GetTitle()); 
	TGraphAsymmErrors* gTOT = new TGraphAsymmErrors(); gTOT->GetXaxis()->SetTitle(histos1[hname+"_calib"]->GetXaxis()->GetTitle()); gTOT->GetYaxis()->SetTitle("Normalized"); gTOT->SetTitle(histos1[hname+"_calib"]->GetTitle()); 

	Double_t entries_calib = 0;
	for(Int_t i=1 ; i<=histos1[hname+"_calib"]->GetNbinsX() ; i++) entries_calib+=histos1[hname+"_calib"]->GetBinContent(i);
	Double_t entries_uncalib = 0;
	for(Int_t i=1 ; i<=histos1[hname+"_uncalib"]->GetNbinsX() ; i++) entries_uncalib+=histos1[hname+"_uncalib"]->GetBinContent(i);

	histos1[hname+"_uncalib"]->Rebin(rebinFactor);
	histos1[hname+"_calib"]->Rebin(rebinFactor);
	histos1[hname+"_JESUP"]->Rebin(rebinFactor);
	histos1[hname+"_JESDWN"]->Rebin(rebinFactor);
	histos1[hname+"_JERUP"]->Rebin(rebinFactor);
	histos1[hname+"_JERDWN"]->Rebin(rebinFactor);

	for(Int_t i=1 ; i<=histos1[hname+"_calib"]->GetNbinsX() ; i++)
	{
		float x = histos1[hname+"_calib"]->GetBinCenter(i);
		float y = histos1[hname+"_calib"]->GetBinContent(i)/entries_calib;
		float dx = histos1[hname+"_calib"]->GetBinWidth(i)/2;

		gJES->SetPoint(i-1, x, y);
		gJER->SetPoint(i-1, x, y);
		gTOT->SetPoint(i-1, x, y);
		
		float dyJES1 = fabs(y-histos1[hname+"_JESUP"]->GetBinContent(i)/entries_calib);
		float dyJES2 = fabs(y-histos1[hname+"_JESDWN"]->GetBinContent(i)/entries_calib);
		float dyJES  = (dyJES1>dyJES2) ? dyJES1 : dyJES2;
		float yErrDwnJES = dyJES;
		float yErrUpJES  = dyJES;

		if((y-yErrDwnJES)<0.) yErrDwnJES = 0.;
		gJES->SetPointError(i-1, dx,dx, yErrDwnJES,yErrUpJES);
		
		float dyJER1 = fabs(y-histos1[hname+"_JERUP"]->GetBinContent(i)/entries_calib);
		float dyJER2 = fabs(y-histos1[hname+"_JERDWN"]->GetBinContent(i)/entries_calib);
		float dyJER  = (dyJER1>dyJER2) ? dyJER1 : dyJER2;
		float yErrDwnJER = dyJER;
		float yErrUpJER  = dyJER;

		if((y-yErrDwnJER)<0.) yErrDwnJER = 0.;
		gJER->SetPointError(i-1, dx,dx, yErrDwnJER,yErrUpJER);
		
		
		float yErrDwnTOT = sqrt(yErrDwnJES*yErrDwnJES + yErrDwnJER*yErrDwnJER);
		float yErrUpTOT  = sqrt(yErrUpJES*yErrUpJES + yErrUpJER*yErrUpJER);
		if((y-yErrDwnTOT)<0.) yErrDwnTOT = 0.;
		gTOT->SetPointError(i-1, dx,dx, yErrDwnTOT,yErrUpTOT);
	}
	
	gTOT->SetFillColor(kRed); gTOT->SetLineColor(kRed);
	gTOT->SetFillColor(kRed); gTOT->SetLineColor(kRed);
	
	gJER->SetFillColor(kYellow); gJER->SetLineColor(kYellow);
	gJER->SetFillColor(kYellow); gJER->SetLineColor(kYellow);

	gJES->SetFillColor(kGreen); gJES->SetLineColor(kGreen);
	gJES->SetFillColor(kGreen); gJES->SetLineColor(kGreen);
	
	TMultiGraph *mgBands = new TMultiGraph();
	mgBands->Add(gTOT);
	mgBands->Add(gJER);
	mgBands->Add(gJES);
	
	
	TLegend* legLeft = new TLegend(0.12,0.65,0.42,0.85,NULL,"brNDC");
	legLeft->SetFillStyle(4000); //will be transparent
	legLeft->SetFillColor(0);
	legLeft->SetTextFont(42);
	legLeft->SetBorderSize(0);
	legLeft->AddEntry(histos1[hname+"_calib"], "Calibrated jets","ple");
	legLeft->AddEntry(histos1[hname+"_uncalib"], "Non calibrated jets","ple");
	legLeft->AddEntry(gJES, "JES uncertainty","f");
	legLeft->AddEntry(gJER, "JER uncertainty","f");
	legLeft->AddEntry(gTOT, "Total uncertainty","f");
	
	
	TLegend* legRight = new TLegend(0.5,0.65,0.8,0.85,NULL,"brNDC");
	legRight->SetFillStyle(4000); //will be transparent
	legRight->SetFillColor(0);
	legRight->SetTextFont(42);
	legRight->SetBorderSize(0);
	legRight->AddEntry(histos1[hname+"_calib"], "Calibrated jets","ple");
	legRight->AddEntry(histos1[hname+"_uncalib"], "Non calibrated jets","ple");
	legRight->AddEntry(gJES, "JES uncertainty","f");
	legRight->AddEntry(gJER, "JER uncertainty","f");
	legRight->AddEntry(gTOT, "Total uncertainty","f");
	
	histos1[hname+"_calib"]->Scale(1./entries_calib);
	histos1[hname+"_uncalib"]->Scale(1./entries_uncalib);
	
	TCanvas* cnv = new TCanvas("c","c",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	cnv->cd();
	histos1[hname+"_calib"]->Draw("p");
	mgBands->Draw("a2");
	mgBands->GetXaxis()->SetTitle(histos1[hname+"_calib"]->GetXaxis()->GetTitle());
	mgBands->GetYaxis()->SetTitle("Normalized");
	mgBands->SetTitle(histos1[hname+"_calib"]->GetTitle());
	histos1[hname+"_uncalib"]->GetYaxis()->SetTitle("Normalized");
	histos1[hname+"_uncalib"]->Draw("p same");
	histos1[hname+"_calib"]->GetYaxis()->SetTitle("Normalized");
	histos1[hname+"_calib"]->Draw("p same");
	if(legpos=="left")  legLeft->Draw("same");
	if(legpos=="right") legRight->Draw("same");
	ptxt->Draw("same");
	
	cnv->RedrawAxis();
	cnv->Update();
	TString pngname = "figures/syst."+hname+".png"; pngname.ReplaceAll("calib_","");
	TString epsname = "figures/syst."+hname+".eps"; epsname.ReplaceAll("calib_","");
	cnv->SaveAs(pngname);
	cnv->SaveAs(epsname);
	cnv->SaveAs(pdfname+pdfmode);
	
	delete gJES;
	delete gJER;
	delete gTOT;
	delete mgBands;
	delete legLeft;
	delete legRight;
	delete cnv;
}


void syst(TString name)
{
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);

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

	modes.push_back("uncalib");  files.push_back(new TFile(name+".histos.muons.mva.uncalib.root", "READ"));
	modes.push_back("calib");    files.push_back(new TFile(name+".histos.muons.mva.calib.root", "READ"));
	modes.push_back("JESUP");    files.push_back(new TFile(name+".histos.muons.mva.JESUP.root",   "READ"));
	modes.push_back("JESDWN");   files.push_back(new TFile(name+".histos.muons.mva.JESDWN.root",  "READ"));
	modes.push_back("JERUP");    files.push_back(new TFile(name+".histos.muons.mva.JERUP.root",   "READ"));
	modes.push_back("JERDWN");   files.push_back(new TFile(name+".histos.muons.mva.JERDWN.root",  "READ"));
	
	TString hname, hnewname;
	TString varname;
	TString pdfname = "figures/"+name+".syst.pdf";
	
	varname = "jet_pt_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname,"(");
	
	varname = "dPhi3bodyJ1_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"left",pdfname);
	
	varname = "met_et_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname);
	
	varname = "met_mt_et3mu_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"left",pdfname);
	
	varname = "met_dphi3mu_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname);
	
	varname = "mettrk_dphimet_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname);
	
	varname = "ht_pt_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname);
	
	varname = "ht_mT_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"left",pdfname);
	
	varname = "ht_dphimet_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname);
	
	varname = "ht_dR3body_calib_all";
	hname = name+"_"+varname;
	addHistos(hname);
	plot(hname,"right",pdfname+")");
}