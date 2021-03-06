//////////////////////////////////////////////////////
//// root -b -l -q allFitSystAuto.C++\(63,0.933\) ////
//////////////////////////////////////////////////////

#include "postBDTcuts.h"
#include "allFitSystAuto.h"

TPaveText* ptxt;
TPaveText* ptxtR;
void makeAtlasLabel()
{
	// ptxt = new TPaveText(0.62,0.70,0.87,0.87,"NDC");
	ptxt = new TPaveText(0.15,0.70,0.40,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
	
	ptxtR = new TPaveText(0.60,0.70,0.85,0.87,"NDC");
	ptxtR->SetFillStyle(4000); //will be transparent
	ptxtR->SetFillColor(0);
	ptxtR->SetTextFont(42);
	ptxtR->SetBorderSize(0);
	ptxtR->AddText("#bf{#it{ATLAS}} internal");
	ptxtR->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxtR->AddText("#sqrt{s}=8 TeV");
}

TLegend* leg;
TLegend* legR;
void makeLegend()
{
	// leg = new TLegend(0.15,0.55,0.45,0.88,NULL,"brNDC");
	leg = new TLegend(0.15,0.55,0.43,0.87,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	
	legR = new TLegend(0.4,0.55,0.88,0.86,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
}



void allFitSystAuto(unsigned int iSR, float currentBDTcut=optBDTcut)
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
	
	
	bool print = false;
	bool write = true;

	
	makeAtlasLabel();
	makeLegend();

	setFiles();
	
	//// get the exact values for the optimal cut
	Result* res = getExtrapolationResult(currentBDTcut,iSR,print,write);
	TString SRname = "iSR="+tstr(iSR,0)+": ["+tstr(res->xSRmin,0)+","+tstr(res->xSRmax,0)+"] #rightarrow "+tstr(res->xSRmax-res->xSRmin,0)+" MeV";

	////////////
	return; ////
	////////////

	TString smMinSBleft    = tstr(mSideBandLeftLowerMeVGlob,0);
	TString smMaxSBleft    = tstr(mBlindMinGlob,0);
	TString smMinSBright   = tstr(mBlindMaxGlob,0);
	TString smMaxSBright   = tstr(mSideBandRightUpperMeVGlob,0);
	TString sxSRmin        = tstr(res->xSRmin,0);
	TString sxSRmax        = tstr(res->xSRmax,0);
	TString scurrentBDTcut = tstr(currentBDTcut,4);
	TString sminBDTcut     = tstr(minBDTcut,4);

	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	// TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_200k_3mu");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");
	
	
	
	float dx      = 0.001; // 0.005;
	float xmin    = 0.5; // minBDTcut;
	float xmax    = +1.0; // maxBDTcut;
	Float_t hXmin = xmin-dx/2.;
	Float_t hXmax = xmax+dx/2.;
	Int_t hNbins  = (hXmax-hXmin)/dx;
	Int_t hbin    = 1;
	
	TH1F* hSBshape  = new TH1F("SBshape",  ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hSBshape->SetLineWidth(2);  hSBshape->SetLineColor(kRed+1);  hSBshape->SetLineStyle(2); 
	TH1F* hSBrange  = new TH1F("SBrange",  ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hSBrange->SetLineWidth(2);  hSBrange->SetLineColor(kRed+1);  hSBrange->SetLineStyle(3); 
	TH1F* hSBcutoff = new TH1F("SBcutoff", ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hSBcutoff->SetLineWidth(2); hSBcutoff->SetLineColor(kRed+1); hSBcutoff->SetLineStyle(5); 
	
	TH1F* hBDTshape  = new TH1F("BDTshape",  ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hBDTshape->SetLineWidth(2);  hBDTshape->SetLineColor(kAzure+1);  hBDTshape->SetLineStyle(2);  
	TH1F* hBDTrange  = new TH1F("BDTrange",  ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hBDTrange->SetLineWidth(2);  hBDTrange->SetLineColor(kAzure+1);  hBDTrange->SetLineStyle(3);  
	TH1F* hBDTcutoff = new TH1F("BDTcutoff", ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hBDTcutoff->SetLineWidth(2); hBDTcutoff->SetLineColor(kAzure+1); hBDTcutoff->SetLineStyle(5); 

	TH1F* hStat = new TH1F("Stat", ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hStat->SetLineWidth(2); hStat->SetLineColor(kGreen+1); hStat->SetLineStyle(1); 
	TH1F* hQuad = new TH1F("Quad", ";BDT score cut;Relative uncertainty [\%]",hNbins,hXmin,hXmax); hQuad->SetLineWidth(2); hQuad->SetLineColor(kBlack); hQuad->SetLineStyle(1); 
	
	TH1F* hCorrection = new TH1F("Correction", ";BDT score cut;Correction factor",hNbins,hXmin,hXmax); hCorrection->SetLineWidth(2); hCorrection->SetLineColor(kBlack); hCorrection->SetLineStyle(1); 
	
	TH1F* hEff     = new TH1F("AccEff",              ";BDT score cut;#it{A}_{s}#times#it{#epsilon}_{s}",hNbins,hXmin,hXmax);
	TH1F* hBkg     = new TH1F("Nbkg",                ";BDT score cut;#it{N}_{bkg}",hNbins,hXmin,hXmax);
	TH1F* hBkgSyst = new TH1F("Nbkg_with_syst_err",  ";BDT score cut;#it{N}_{bkg}",hNbins,hXmin,hXmax);
	TH1F* hBkgStat = new TH1F("Nbkg_with_stat_err",  ";BDT score cut;#it{N}_{bkg}",hNbins,hXmin,hXmax);
	
	leg->Clear();
	leg->AddEntry(hQuad,     "Systematic (total)","l");
	leg->AddEntry(hStat,     "Statistic","l");
	leg->AddEntry(hBDTshape, "BDT shape #Delta#it{R}","l");
	leg->AddEntry(hSBshape,  "SB shape #DeltaN_{SR0}","l");
	leg->AddEntry(hBDTrange, "SB range #Delta#it{R}","l");
	leg->AddEntry(hSBrange,  "SB range #DeltaN_{SR0}","l");
	leg->AddEntry(hBDTcutoff,"BDT cutoff #Delta#it{R}","l");
	leg->AddEntry(hSBcutoff, "BDT cutoff #DeltaN_{SR0}","l");
	
	ofstream fbkg;
	fbkg.open("fbkg.txt");
	
	for(float x=xmin ; x<xmax ; x+=dx, ++hbin)
	{
		Result* R = getExtrapolationResult(x,iSR,print);
		
		hSBshape->SetBinContent(hbin,  R->VARdnSR0shape/R->VARnSR0*100);
		hSBrange->SetBinContent(hbin,  R->VARdnSR0range/R->VARnSR0*100);
		hSBcutoff->SetBinContent(hbin, R->VARdnSR0cutoff/R->VARnSR0*100);
		// hSBshape->SetBinContent(hbin,  R->dnSR0shape/R->nSR0*100);
		// hSBrange->SetBinContent(hbin,  R->dnSR0range/R->nSR0*100);
		// hSBcutoff->SetBinContent(hbin, R->dnSR0cutoff/R->nSR0*100);
		
		hBDTshape->SetBinContent(hbin,  R->df01shape/R->f01*100);	
		hBDTrange->SetBinContent(hbin,  R->df01range/R->f01*100);	
		hBDTcutoff->SetBinContent(hbin, R->df01cutoff/R->f01*100);
		
		// hQuad->SetBinContent(hbin, R->dnSR1quad/R->nSR1*100);
		hQuad->SetBinContent(hbin, R->VARdnSR1quad/R->VARnSR1*100);
	
		// hStat->SetBinContent(hbin, R->dnSR1stat/R->nSR1*100);
		hStat->SetBinContent(hbin, R->VARdnSR1stat/R->VARnSR1*100);
		
		hCorrection->SetBinContent(hbin, R->correction);

		TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, sminBDTcut,tstr(x,4), "sig", sxSRmin,sxSRmax, false,true,"postTraining");
		float npassedS   = tS->GetEntries(cuts_sig);
		float ninitS     = (((TString)tS->GetName()).Contains("Wtaunu_200k_3mu")) ? 200000 : 99900;
		float AccEffSig  = npassedS/ninitS;
		hEff->SetBinContent(hbin,AccEffSig);
		hBkg->SetBinContent(hbin,R->VARnSR1);
		hBkgSyst->SetBinContent(hbin,R->VARnSR1); hBkgSyst->SetBinError(hbin,R->VARdnSR1quad);
		hBkgStat->SetBinContent(hbin,R->VARnSR1); hBkgStat->SetBinError(hbin,R->VARdnSR1stat);
		
		fbkg << "x1=" << tstr(x,3) << " Acc*Eff=" << tstr(AccEffSig,4) << " Nbkg=" << tstr(R->VARnSR1,4) << " +- " << tstr(R->VARdnSR1quad,4) << "(syst) +- " << tstr(R->VARdnSR1stat,4) << "(stat)" << endl;
		
		delete R;
	}
	
	TFile* fBkg = new TFile("Background.root","RECREATE");
	fBkg->cd();
	hEff->Write();
	hBkg->Write();
	hBkgSyst->Write();
	hBkgStat->Write();
	fBkg->Write();
	fBkg->Close();
	fbkg.close();
	
	
	
	
	
	TCanvas* cnv = NULL;
	
	if(cnv) delete cnv;
	cnv = new TCanvas("","",800,600);
	cnv->Draw();
	cnv->SetTicks(1,1);
	// cnv->SetGridy();
	// cnv->SetLogy();
	hQuad->SetMinimum(0);
	hQuad->SetMaximum(100);
	hQuad->Draw();
	hBDTshape->Draw("same");
	hBDTrange->Draw("same");
	hBDTcutoff->Draw("same");
	hSBshape->Draw("same");
	hSBrange->Draw("same");
	hSBcutoff->Draw("same");
	hStat->Draw("same");
	ptxtR->Draw("same");
	leg->Draw("same");
	cnv->Update();
	cnv->SaveAs("figures/nSR1fitSyst.png");
	cnv->SaveAs("figures/nSR1fitSyst.eps");
	cnv->SaveAs("figures/nSR1fitSyst.pdf(");
	
	if(cnv) delete cnv;
	cnv = new TCanvas("","",800,600);
	cnv->Draw();
	cnv->SetTicks(1,1);
	hCorrection->SetMinimum(0.75);
	hCorrection->SetMaximum(1.0);
	hCorrection->Draw();
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs("figures/CorrectionFactor.png");
	cnv->SaveAs("figures/CorrectionFactor.eps");
	cnv->SaveAs("figures/nSR1fitSyst.pdf");
	
	
	//// calculate the signal efficiency
	TString loose_test_cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, "-0.5","-0.5", "bkg", "",     "",      true,true,"preTraining");
	cout << "loose_test_cuts_bkg=" << loose_test_cuts_bkg << endl;
	TString loose_cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, sminBDTcut,sminBDTcut, "bkg", "",     "",      true,true,"preTraining");
	cout << "loose_cuts_bkg=" << loose_cuts_bkg << endl;
	TString tight_test_cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, "-0.5","-0.5", "bkg", "",     "",      true,true,"postTraining");
	cout << "tight_test_cuts_bkg=" << tight_test_cuts_bkg << endl;
	
	TString tight_cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, sminBDTcut,sminBDTcut, "bkg", "",     "",      true,true,"postTraining");
	cout << "tight_cuts_bkg=" << tight_cuts_bkg << endl;
	Double_t m3bodyMin = mSideBandLeftLowerMeVGlob;
	Double_t m3bodyMax = mSideBandRightUpperMeVGlob;
	Double_t m3bodyBinSize = mBinSize;
	Int_t nm3bodybins = (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize);
	TH1* hBDTloose  = new TH1F("hBDTloose",";BDTloose score;Events",nBDTbins,minBDTcut,+1); hBDTloose->Sumw2(); hBDTloose->SetLineColor(kRed); hBDTloose->SetLineWidth(1); hBDTloose->SetMarkerColor(kRed); hBDTloose->SetMarkerStyle(21); hBDTloose->SetMarkerSize(1.2); hBDTloose->SetBinErrorOption(TH1::kPoisson);
	TH1* hBDTloosezoom  = new TH1F("hBDTloosezoom",";BDT score;Events",nBDTbins,minBDTcut,+1); hBDTloosezoom->Sumw2(); hBDTloosezoom->SetLineColor(kRed); hBDTloosezoom->SetLineWidth(1); hBDTloosezoom->SetMarkerColor(kRed); hBDTloosezoom->SetMarkerStyle(21); hBDTloosezoom->SetMarkerSize(1.2); hBDTloosezoom->SetBinErrorOption(TH1::kPoisson);
	TH1* hBDT  = new TH1F("hBDT",";BDT score;Events",nBDTbins,minBDTcut,+1); hBDT->Sumw2(); hBDT->SetLineColor(kBlack); hBDT->SetLineWidth(1); hBDT->SetMarkerColor(kBlack); hBDT->SetMarkerStyle(20); hBDT->SetMarkerSize(1.2); hBDT->SetBinErrorOption(TH1::kPoisson);
	TH1* hBDTzoom  = new TH1F("hBDTzoom",";BDT score;Events",nBDTbins,minBDTcut,+1); hBDT->Sumw2(); hBDTzoom->SetLineColor(kBlack); hBDTzoom->SetLineWidth(1); hBDTzoom->SetMarkerColor(kBlack); hBDTzoom->SetMarkerStyle(20); hBDTzoom->SetMarkerSize(1.2); hBDTzoom->SetBinErrorOption(TH1::kPoisson);
	TH1* hdBDT = new TH1F("hdBDT",";BDT score;Events",nBDTbins,minBDTcut,+1); hdBDT->Sumw2(); hdBDT->SetLineColor(kBlack); hdBDT->SetFillColor(kBlack); hdBDT->SetFillStyle(3344);
	TH1* hdBDTzoom = new TH1F("hdBDTzoom",";BDT score;Events",nBDTbins,minBDTcut,+1); hdBDTzoom->Sumw2(); hdBDTzoom->SetLineColor(kBlack); hdBDTzoom->SetFillColor(kBlack); hdBDTzoom->SetFillStyle(3344);
	TH1* hSB   = new TH1F("hSB",";#it{m}_{3body} [MeV];Events",nm3bodybins,m3bodyMin,m3bodyMax); hSB->Sumw2(); hSB->SetLineColor(kBlack); hSB->SetLineWidth(1); hSB->SetMarkerColor(kBlack); hSB->SetMarkerStyle(20); hSB->SetMarkerSize(1.2); hSB->SetBinErrorOption(TH1::kPoisson);
	TH1* hdSB  = new TH1F("hdSB",";#it{m}_{3body} [MeV];Events",nm3bodybins,m3bodyMin,m3bodyMax); hdSB->Sumw2(); hdSB->SetLineColor(kBlack); hdSB->SetFillColor(kBlack); hdSB->SetFillStyle(3344);
	
	
	Double_t NlooseTest = tD->GetEntries(loose_test_cuts_bkg);
	Double_t NtightTest = tD->GetEntries(tight_test_cuts_bkg);
	Double_t Nloose = tD->Draw("score>>hBDTloose",loose_cuts_bkg);
	Double_t Ntight = tD->Draw("score>>hBDT",tight_cuts_bkg);
	// tD->Draw("score>>hBDTloosezoom",loose_cuts_bkg);
	// tD->Draw("score>>hBDTzoom",tight_cuts_bkg);
	tD->Draw("m3body>>hSB",tight_cuts_bkg);
	getEnvelopes(hdBDT,hdSB,hBDT,hSB,false);
	// hBDTloose->Scale(Ntight/Nloose);
	hBDTloose->Scale(NtightTest/NlooseTest);
	
	cout << "NtightTest/NlooseTest = " << NtightTest/NlooseTest << endl;
	cout << "Ntight/Nloose = " << Ntight/Nloose << endl;
	
	TH1* hlBDT = (TH1*)hdBDT->Clone("lBDT"); hlBDT->SetLineColor(kBlack); hlBDT->SetLineWidth(2); hlBDT->SetFillColor(0); hlBDT->Smooth();
	TH1* hlSB = (TH1*)hdSB->Clone("lSB");    hlSB->SetLineColor(kBlack);  hlSB->SetLineWidth(2);  hlSB->SetFillColor(0);  hlSB->Smooth();
	
	
	// for(Int_t b=1 ; b<=hBDT->GetNbinsX() ; ++b) { if(hBDT->GetBinContent(b)<1) { hBDTzro->SetBinContent(b,0); hBDTzro->SetBinError(b,0); } else { hBDTzro->SetBinContent(b,-1); hBDTzro->SetBinError(b,0); } }
	// for(Int_t b=1 ; b<=hSB->GetNbinsX() ; ++b)  { if(hSB->GetBinContent(b)<1)  { hSBzro->SetBinContent(b,0);  hSBzro->SetBinError(b,0);  } else { hSBzro->SetBinContent(b,-1);  hSBzro->SetBinError(b,0);  } }
	
	legR->Clear();
	legR->AddEntry(hBDT,"SB data (tight+x>x_{0} )","lpe");
	legR->AddEntry(hBDTloose,"SB data (loose+x>x_{0})","lpe");
	legR->AddEntry((TObject*)0, "Normalized to tight in x>-0.5", "");
	legR->AddEntry(hlBDT,"Fit to the data","l");
	legR->AddEntry(hdBDT,"Fit uncertainty","f");
	if(cnv) delete cnv;
	cnv = new TCanvas("","",800,600);
	cnv->Draw();
	cnv->SetTicks(1,1);
	hdBDT->SetMinimum(0);
	hdBDT->SetMaximum(12);
	hdBDT->Draw("e2");
	hlBDT->Draw("hist same");
	hBDT->Draw("p0 e same");
	hBDTloose->Draw("p0 e same");
	ptxt->Draw("same");
	legR->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/BDTenvelope.pdf");
	cnv->SaveAs("figures/BDTenvelope.png");
	cnv->SaveAs("figures/BDTenvelope.eps");
	cnv->SaveAs("figures/nSR1fitSyst.pdf");
	
	legR->Clear();
	legR->AddEntry(hSB,"SB data (tight+x>x_{0})","lpe");
	legR->AddEntry(hlSB,"Nominal fit","l");
	legR->AddEntry(hdSB,"Fit uncertainty","f");
	if(cnv) delete cnv;
	cnv = new TCanvas("","",800,600);
	cnv->Draw();
	cnv->SetTicks(1,1);
	hdSB->SetMinimum(0);
	hdSB->SetMaximum(12);
	hdSB->Draw("e2");
	hlSB->Draw("hist same");
	hSB->Draw("p0 e same");
	ptxt->Draw("same");
	legR->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/SBenvelope.pdf");
	cnv->SaveAs("figures/SBenvelope.png");
	cnv->SaveAs("figures/SBenvelope.eps");
	cnv->SaveAs("figures/nSR1fitSyst.pdf");
	
	
	if(cnv) delete cnv;
	cnv = new TCanvas("","",800,600);
	cnv->Draw();
	cnv->SaveAs("figures/nSR1fitSyst.pdf)");
	
	
	TFile* fUncert = new TFile("Uncert.root","RECREATE");
	fUncert->cd();
	hdBDT->Write();
	hlBDT->Write();
	hBDT->Write();
	hdSB->Write();
	hlSB->Write();
	hSB->Write();
	fUncert->Write();
	fUncert->Close();
	
	
	TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, sminBDTcut,scurrentBDTcut, "sig", sxSRmin,sxSRmax, false,true,"postTraining");
	TString cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright, sminBDTcut,scurrentBDTcut, "bkg", "",     "",      true,true,"postTraining");
	
	Float_t npassedS = tS->GetEntries(cuts_sig);
	Float_t ninitS   = (((TString)tS->GetName()).Contains("Wtaunu_200k_3mu")) ? 200000 : 99900;
	float AccEffSig  = npassedS/ninitS*100;
	cout << "cuts_sig=" << cuts_sig << endl;
	
	// tD->SetEventList(0);
	// tD->Draw(">>elistD",cuts_bkg);
	// TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
	// Float_t npassedD   = elistD->GetN();  // number of events to pass cuts
	// cout << "cuts_bkg=" << cuts_bkg << endl;
	
	Float_t npassedD   = tD->GetEntries(cuts_bkg);
	cout << "cuts_bkg=" << cuts_bkg << endl;
	
	
	cout << "npassedS=" << npassedS << endl;
	cout << "npassedD=" << npassedD << endl;
	
	
	cout << "=================================" << endl;
	cout << "currentBDTcut      : " << currentBDTcut << endl;
	cout << "iSR                : " << iSR << endl;
	cout << "xSRmin             : " << res->xSRmin << " MeV" << endl;
	cout << "xSRmax             : " << res->xSRmax << " MeV" << endl;
	cout << "SR width           : " << res->xSRmax-res->xSRmin << " MeV" << endl;
	cout << "Signal A*e         : " << AccEffSig << "\%" << endl;
	cout << "=================================" << endl;
	cout << "nSB0               : " << tstr(res->nSB0,3) << endl;
	// cout << "nSR0               : " << tstr(res->nSR0,3) << endl;
	cout << "nSR0               : " << tstr(res->VARnSR0,3) << endl;
	cout << "f01                : " << tstr(res->f01,7) << endl;
	cout << "Correction         : " << tstr(res->correction,3) << endl;
	cout << "nSB1(cut)          : " << tstr(npassedD,3) << endl;
	cout << "nSB1(fit)          : " << tstr(res->nSB1,3) << endl;
	// cout << "nSR1               : " << tstr(res->nSR1,3) << endl;
	cout << "nSR1               : " << tstr(res->VARnSR1,3) << "+-" << res->VARdnSR1stat << "(stat)" << endl;
	cout << "=================================" << endl;
	// cout << "dnSR0 (SB shape)   : " << tstr(res->dnSR0shape/res->nSR0*100 ,1) << "\%" << endl;
	// cout << "dnSR0 (SB range)   : " << tstr(res->dnSR0range/res->nSR0*100 ,1) << "\%" << endl;
	// cout << "dnSR0 (BDT cutoff) : " << tstr(res->dnSR0cutoff/res->nSR0*100,1) << "\%" << endl;
	cout << "dnSR0 (SB shape)   : " << tstr(res->VARdnSR0shape/res->VARnSR0*100 ,1) << "\%" << endl;
	cout << "dnSR0 (SB range)   : " << tstr(res->VARdnSR0range/res->VARnSR0*100 ,1) << "\%" << endl;
	cout << "dnSR0 (BDT cutoff) : " << tstr(res->VARdnSR0cutoff/res->VARnSR0*100,1) << "\%" << endl;
	cout << "df01 (BDT shape)   : " << tstr(res->df01shape/res->f01*100   ,1) << "\%" << endl;	
	cout << "df01 (SB range)    : " << tstr(res->df01range/res->f01*100   ,1) << "\%" << endl;	
	cout << "df01 (BDT cutoff)  : " << tstr(res->df01cutoff/res->f01*100  ,1) << "\%" << endl;
	// cout << "dnSR1 (quad.)      : " << tstr(res->dnSR1quad/res->nSR1*100  ,1) << "\%" << endl;
	// cout << "dnSR1 (stat.)      : " << res->dnSR1stat/res->nSR1*100   << "\%" << endl;
	cout << "dnSR1 (quad.)      : " << tstr(res->VARdnSR1quad/res->VARnSR1*100  ,1) << "\%" << endl;
	cout << "dnSR1 (stat.)      : " << res->VARdnSR1stat/res->VARnSR1*100   << "\%" << endl;
	cout << "=================================" << endl;
	delete res;
}
