////////////////////////////////////////
//// root -b -l -q sidebandsDef.C++ ////
////////////////////////////////////////
#include "std.h"
#include "const.h"


TPaveText* ptxt;
void makeAtlasLabel()
{	
	ptxt = new TPaveText(0.55,0.60,0.90,0.82,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

void sidebandsDef()
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
	
	// TFile* f = new TFile("histos.muons.cuts.root","READ");
	// TFile* f = new TFile("hM3bodyLoose.root","READ");

	// TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.loose.root","READ");
	TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.veryloose.root","READ");
	TTree* tD = (TTree*)f->Get("flatout_Data");
	TCut    isLoose    = "pass_loose==1";
	TCut    oneCand    = "@vtx_mass->size()==1";
	TCut    masses2    = "vtx_mOS1>220 && vtx_mOS2>220 && vtx_mSS>220";
	TCut    distances  = "geo_lxySig>-10 && geo_lxySig<50 && geo_a0xySig<25";
	TCut    pvalues    = "trks_fitprob>1e-9";
	TCut    kinematics = "vtx_pt>10000 && (met_muons_et>10000 && met_muons_et<250000) && (met_track_et>10000 && met_track_et<250000)  &&  met_muons_mT>20000 && met_track_mT>20000";
	TCut    iso        = "vtx_isolation020<0.3 && vtx_isolation030<1";
	TCut preTraining   = isLoose+oneCand+masses2+iso+distances+pvalues+kinematics;
	TH1F* hM3body  = new TH1F("hM3body",";#it{m}_{3body} [MeV];Events", 100,700,2700); hM3body->SetMarkerStyle(20); hM3body->SetBinErrorOption(TH1::kPoisson); 
	tD->Draw("vtx_mass>>hM3body",preTraining);


	makeAtlasLabel();

	TCanvas* cnv = new TCanvas("","",600,400);
	cnv->Draw();

	// TH1F* hM3body = (TH1F*)f->Get("Data_m3body_lin_zoom_after_triplet");
	// TH1F* hM3body = (TH1F*)f->Get("hM3body");
	hM3body->Rebin(5);
	hM3body->SetMinimum(0);
	hM3body->SetMaximum(2.*hM3body->GetMaximum());
	hM3body->SetLineColor(kBlack);
	hM3body->SetLineWidth(1);
	hM3body->SetLineStyle(1);
	hM3body->SetMarkerStyle(20);
	hM3body->SetMarkerSize(0.6);
	hM3body->SetMarkerColor(kBlack);
	hM3body->Draw("p0 e0");
	ptxt->Draw("same");
	
	cnv->Update();
	Double_t ymax = cnv->GetUymax();
	Double_t y = -1;
	
	y = ymax/1.7;
	TLine* lTrainingMin = new TLine(750,0,750,y);
	lTrainingMin->SetLineColor(kAzure+10);
	lTrainingMin->SetLineWidth(2);
	lTrainingMin->SetLineStyle(1);
	lTrainingMin->Draw("L same");
	TLine* lTrainingMax = new TLine(2500,0,2500,y);
	lTrainingMax->SetLineColor(kAzure+10);
	lTrainingMax->SetLineWidth(2);
	lTrainingMax->SetLineStyle(1);
	lTrainingMax->Draw("L same");
	TLine* lTraining = new TLine(750,y,2500,y);
	lTraining->SetLineColor(kAzure+10);
	lTraining->SetLineWidth(2);
	lTraining->SetLineStyle(1);
	lTraining->Draw("L same");
	
	y = ymax/1.8;
	TLine* lSidebandsMin = new TLine(1450,0,1450,y);
	lSidebandsMin->SetLineColor(kGreen+1);
	lSidebandsMin->SetLineWidth(2);
	lSidebandsMin->SetLineStyle(2);
	lSidebandsMin->Draw("L same");
	TLine* lSidebandsMax = new TLine(2110,0,2110,y);
	lSidebandsMax->SetLineColor(kGreen+1);
	lSidebandsMax->SetLineWidth(2);
	lSidebandsMax->SetLineStyle(2);
	lSidebandsMax->Draw("L same");
	TLine* lSidebands = new TLine(1450,y,2110,y);
	lSidebands->SetLineColor(kGreen+1);
	lSidebands->SetLineWidth(2);
	lSidebands->SetLineStyle(2);
	lSidebands->Draw("L same");
	
	y = ymax/1.9;
	TLine* lBlindedMin = new TLine(1690,0,1690,y);
	lBlindedMin->SetLineColor(kBlack);
	lBlindedMin->SetLineWidth(2);
	lBlindedMin->SetLineStyle(3);
	lBlindedMin->Draw("L same");
	TLine* lBlindedMax = new TLine(1870,0,1870,y);
	lBlindedMax->SetLineColor(kBlack);
	lBlindedMax->SetLineWidth(2);
	lBlindedMax->SetLineStyle(3);
	lBlindedMax->Draw("L same");
	TLine* lBlinded = new TLine(1690,y,1870,y);
	lBlinded->SetLineColor(kBlack);
	lBlinded->SetLineWidth(2);
	lBlinded->SetLineStyle(3);
	lBlinded->Draw("L same");
	
	y = ymax/2.0;
	TLine* lSignalRegionMin = new TLine(1713,0,1713,y);
	lSignalRegionMin->SetLineColor(kRed);
	lSignalRegionMin->SetLineWidth(2);
	lSignalRegionMin->SetLineStyle(5);
	lSignalRegionMin->Draw("L same");
	TLine* lSignalRegionMax = new TLine(1841,0,1841,y);
	lSignalRegionMax->SetLineColor(kRed);
	lSignalRegionMax->SetLineWidth(2);
	lSignalRegionMax->SetLineStyle(5);
	lSignalRegionMax->Draw("L same");
	TLine* lSignalRegion = new TLine(1713,y,1841,y);
	lSignalRegion->SetLineColor(kRed);
	lSignalRegion->SetLineWidth(2);
	lSignalRegion->SetLineStyle(5);
	lSignalRegion->Draw("L same");
	
	
	hM3body->Draw("p0 e0 same");
	
	TLegend* leg = new TLegend(0.12,0.60,0.65,0.89);
	leg->SetFillStyle(4000); //will be transparent
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(hM3body,"Data (loose)","ple");
	leg->AddEntry(lTraining,"Training","l");
	leg->AddEntry(lSidebands,"Sidebands","l");
	leg->AddEntry(lBlinded,"Blinded","l");
	leg->AddEntry(lSignalRegion,"Signal","l");
	leg->Draw("same");
	
	
	cnv->SaveAs("figures/sidebandsDef.png");
	cnv->SaveAs("figures/sidebandsDef.eps");
	cnv->SaveAs("figures/sidebandsDef.pdf");
}