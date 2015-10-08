/////////////////////////////////////
//// root -b -l -q phitest.C++ //////
/////////////////////////////////////
#include "std.h"
#include "const.h"
#include <iterator>
#include <cstdlib>

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}

Int_t npar = 5;
Double_t func(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	return par[0]*exp(x*par[1]) + par[2]*exp(-pow((x-par[3])/par[4],2));
}

TPaveText* ptxt;
void makeAtlasLabel()
{	
	ptxt = new TPaveText(0.35,0.65,0.60,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

void phitest()
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


	// TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.loose.root","READ"); TString cuts = "pass_loose==1 && @vtx_mass->size()==1";
	TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.root","READ");    TString cuts = "TMath::Abs(vtx_charge)==1";
	TTree* tD = (TTree*)f->Get("flatout_Data");
	TTree* tS = (TTree*)f->Get("flatout_Wtaunu_3mu");

	makeAtlasLabel();

	tS->SetEventList(0);
	tS->Draw(">>elistS",cuts);
	TEventList* elistS = (TEventList*)gDirectory->Get("elistS");
	Float_t npassedS = elistS->GetN();  // number of events to pass cuts

	tD->SetEventList(0);
	tD->Draw(">>elistD",cuts);
	TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
	Float_t npassedD = elistD->GetN();  // number of events to pass cuts

	tD->SetEventList(elistD);
	tS->SetEventList(elistS);

	TH1F* mPhiD = new TH1F("mPhiD",";m_{OS} [MeV];Events",80,700,1200);
	TH1F* mPhiD1 = new TH1F("mPhiD1",";m_{OS1} [MeV];Events",80,700,1200);
	TH1F* mPhiD2 = new TH1F("mPhiD2",";m_{OS2} [MeV];Events",80,700,1200);
	TH1F* mPhiD3 = new TH1F("mPhiD3",";m_{SS} [MeV];Events",80,700,1200);
	TH1F* mPhiS = new TH1F("mPhiS",";m_{OS} [MeV];Events",80,700,1200);
	TH1F* mPhiS1 = new TH1F("mPhiS1",";m_{OS1} [MeV];Events",80,700,1200);
	TH1F* mPhiS2 = new TH1F("mPhiS2",";m_{OS2} [MeV];Events",80,700,1200);
	TH1F* mPhiS3 = new TH1F("mPhiS3",";m_{SS} [MeV];Events",80,700,1200);
	
	tD->Draw("vtx_mOS1>>mPhiD1");
	tD->Draw("vtx_mOS2>>mPhiD2");
	tD->Draw("vtx_mSS>>mPhiD3");
	tD->Draw("vtx_mOS1>>mPhiD"); mPhiD->Add(mPhiD2);
	tS->Draw("vtx_mOS1>>mPhiS1");
	tS->Draw("vtx_mOS2>>mPhiS2");
	tS->Draw("vtx_mSS>>mPhiS3");
	tS->Draw("vtx_mOS1>>mPhiS"); mPhiS->Add(mPhiS2);
	
	
	

	TH2F* mSS_vs_pcnHits1D  = new TH2F("mSS_vs_pcnHits1D", ";m_{SS} [MeV];Precision hits trk1;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_pcnHits1S  = new TH2F("mSS_vs_pcnHits1S", ";m_{SS} [MeV];Precision hits trk1;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits1D = new TH2F("mOS1_vs_pcnHits1D",";m_{OS1} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits1S = new TH2F("mOS1_vs_pcnHits1S",";m_{OS1} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits1D = new TH2F("mOS2_vs_pcnHits1D",";m_{OS2} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits1S = new TH2F("mOS2_vs_pcnHits1S",";m_{OS2} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);

	TH2F* mSS_vs_pcnHits2D  = new TH2F("mSS_vs_pcnHits2D", ";m_{SS} [MeV];Precision hits trk2;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_pcnHits2S  = new TH2F("mSS_vs_pcnHits2S", ";m_{SS} [MeV];Precision hits trk2;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits2D = new TH2F("mOS1_vs_pcnHits2D",";m_{OS1} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits2S = new TH2F("mOS1_vs_pcnHits2S",";m_{OS1} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits2D = new TH2F("mOS2_vs_pcnHits2D",";m_{OS2} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits2S = new TH2F("mOS2_vs_pcnHits2S",";m_{OS2} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);

	TH2F* mSS_vs_pcnHits3D  = new TH2F("mSS_vs_pcnHits3D", ";m_{SS} [MeV];Precision hits trk3;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_pcnHits3S  = new TH2F("mSS_vs_pcnHits3S", ";m_{SS} [MeV];Precision hits trk3;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits3D = new TH2F("mOS1_vs_pcnHits3D",";m_{OS1} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_pcnHits3S = new TH2F("mOS1_vs_pcnHits3S",";m_{OS1} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits3D = new TH2F("mOS2_vs_pcnHits3D",";m_{OS2} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_pcnHits3S = new TH2F("mOS2_vs_pcnHits3S",";m_{OS2} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	
	tD->Draw("mu_nPIXhits1:vtx_mSS>>mSS_vs_pcnHits1D");   tS->Draw("mu_nPIXhits1:vtx_mSS>>mSS_vs_pcnHits1S");
	tD->Draw("mu_nPIXhits1:vtx_mOS1>>mOS1_vs_pcnHits1D"); tS->Draw("mu_nPIXhits1:vtx_mOS1>>mOS1_vs_pcnHits1S");
	tD->Draw("mu_nPIXhits1:vtx_mOS2>>mOS2_vs_pcnHits1D"); tS->Draw("mu_nPIXhits1:vtx_mOS2>>mOS2_vs_pcnHits1S");
                                                                    
	tD->Draw("mu_nPIXhits2:vtx_mSS>>mSS_vs_pcnHits2D");   tS->Draw("mu_nPIXhits2:vtx_mSS>>mSS_vs_pcnHits2S");
	tD->Draw("mu_nPIXhits2:vtx_mOS1>>mOS1_vs_pcnHits2D"); tS->Draw("mu_nPIXhits2:vtx_mOS1>>mOS1_vs_pcnHits2S");
	tD->Draw("mu_nPIXhits2:vtx_mOS2>>mOS2_vs_pcnHits2D"); tS->Draw("mu_nPIXhits2:vtx_mOS2>>mOS2_vs_pcnHits2S");
                                                                    
	tD->Draw("mu_nPIXhits3:vtx_mSS>>mSS_vs_pcnHits3D");   tS->Draw("mu_nPIXhits3:vtx_mSS>>mSS_vs_pcnHits3S");
	tD->Draw("mu_nPIXhits3:vtx_mOS1>>mOS1_vs_pcnHits3D"); tS->Draw("mu_nPIXhits3:vtx_mOS1>>mOS1_vs_pcnHits3S");
	tD->Draw("mu_nPIXhits3:vtx_mOS2>>mOS2_vs_pcnHits3D"); tS->Draw("mu_nPIXhits3:vtx_mOS2>>mOS2_vs_pcnHits3S");
	
	
	
	TH2F* mSS_vs_etaphiLayers1D  = new TH2F("mSS_vs_etaphiLayers1D", ";m_{SS} [MeV];Precision hits trk1;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_etaphiLayers1S  = new TH2F("mSS_vs_etaphiLayers1S", ";m_{SS} [MeV];Precision hits trk1;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers1D = new TH2F("mOS1_vs_etaphiLayers1D",";m_{OS1} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers1S = new TH2F("mOS1_vs_etaphiLayers1S",";m_{OS1} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers1D = new TH2F("mOS2_vs_etaphiLayers1D",";m_{OS2} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers1S = new TH2F("mOS2_vs_etaphiLayers1S",";m_{OS2} [MeV];Precision hits trk1;Normalized",100,250,1250, 10,0,10);
	
	TH2F* mSS_vs_etaphiLayers2D  = new TH2F("mSS_vs_etaphiLayers2D", ";m_{SS} [MeV];Precision hits trk2;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_etaphiLayers2S  = new TH2F("mSS_vs_etaphiLayers2S", ";m_{SS} [MeV];Precision hits trk2;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers2D = new TH2F("mOS1_vs_etaphiLayers2D",";m_{OS1} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers2S = new TH2F("mOS1_vs_etaphiLayers2S",";m_{OS1} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers2D = new TH2F("mOS2_vs_etaphiLayers2D",";m_{OS2} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers2S = new TH2F("mOS2_vs_etaphiLayers2S",";m_{OS2} [MeV];Precision hits trk2;Normalized",100,250,1250, 10,0,10);
	
	TH2F* mSS_vs_etaphiLayers3D  = new TH2F("mSS_vs_etaphiLayers3D", ";m_{SS} [MeV];Precision hits trk3;Normalized", 100,250,1250, 10,0,10);
	TH2F* mSS_vs_etaphiLayers3S  = new TH2F("mSS_vs_etaphiLayers3S", ";m_{SS} [MeV];Precision hits trk3;Normalized", 100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers3D = new TH2F("mOS1_vs_etaphiLayers3D",";m_{OS1} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS1_vs_etaphiLayers3S = new TH2F("mOS1_vs_etaphiLayers3S",";m_{OS1} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers3D = new TH2F("mOS2_vs_etaphiLayers3D",";m_{OS2} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);
	TH2F* mOS2_vs_etaphiLayers3S = new TH2F("mOS2_vs_etaphiLayers3S",";m_{OS2} [MeV];Precision hits trk3;Normalized",100,250,1250, 10,0,10);

	tD->Draw("mu_nEtaPhiLayers1:vtx_mSS>>mSS_vs_etaphiLayers1D");   tS->Draw("mu_nEtaPhiLayers1:vtx_mSS>>mSS_vs_etaphiLayers1S");
	tD->Draw("mu_nEtaPhiLayers1:vtx_mOS1>>mOS1_vs_etaphiLayers1D"); tS->Draw("mu_nEtaPhiLayers1:vtx_mOS1>>mOS1_vs_etaphiLayers1S");
	tD->Draw("mu_nEtaPhiLayers1:vtx_mOS2>>mOS2_vs_etaphiLayers1D"); tS->Draw("mu_nEtaPhiLayers1:vtx_mOS2>>mOS2_vs_etaphiLayers1S");
	
	tD->Draw("mu_nEtaPhiLayers2:vtx_mSS>>mSS_vs_etaphiLayers2D");   tS->Draw("mu_nEtaPhiLayers2:vtx_mSS>>mSS_vs_etaphiLayers2S");
	tD->Draw("mu_nEtaPhiLayers2:vtx_mOS1>>mOS1_vs_etaphiLayers2D"); tS->Draw("mu_nEtaPhiLayers2:vtx_mOS1>>mOS1_vs_etaphiLayers2S");
	tD->Draw("mu_nEtaPhiLayers2:vtx_mOS2>>mOS2_vs_etaphiLayers2D"); tS->Draw("mu_nEtaPhiLayers2:vtx_mOS2>>mOS2_vs_etaphiLayers2S");
	
	tD->Draw("mu_nEtaPhiLayers3:vtx_mSS>>mSS_vs_etaphiLayers3D");   tS->Draw("mu_nEtaPhiLayers2:vtx_mSS>>mSS_vs_etaphiLayers3S");
	tD->Draw("mu_nEtaPhiLayers3:vtx_mOS1>>mOS1_vs_etaphiLayers3D"); tS->Draw("mu_nEtaPhiLayers2:vtx_mOS1>>mOS1_vs_etaphiLayers3S");
	tD->Draw("mu_nEtaPhiLayers3:vtx_mOS2>>mOS2_vs_etaphiLayers3D"); tS->Draw("mu_nEtaPhiLayers2:vtx_mOS2>>mOS2_vs_etaphiLayers3S");



	TF1* fcnRhoOmega = new TF1("fcnRhoOmega",func,700,890,npar);
	fcnRhoOmega->SetParName(0,"Nexpo_rhoOmega");   /*fcnRhoOmega->SetParameter(0,...);*/   fcnRhoOmega->SetParLimits(0,0,1e5);
	fcnRhoOmega->SetParName(1,"Gexpo_rhoOmega");   /*fcnRhoOmega->SetParameter(1,...);*/   fcnRhoOmega->SetParLimits(1,-10,0);
	fcnRhoOmega->SetParName(2,"Ngaus_rhoOmega");   /*fcnRhoOmega->SetParameter(2,10);*/    fcnRhoOmega->SetParLimits(2,0,1e4);
	fcnRhoOmega->SetParName(3,"Mean_rhoOmega");    fcnRhoOmega->SetParameter(3,782);       fcnRhoOmega->SetParLimits(3,750,800);
	fcnRhoOmega->SetParName(4,"Sigma_rhoOmega");   fcnRhoOmega->SetParameter(4,30.);       fcnRhoOmega->SetParLimits(4,0,100);
	fcnRhoOmega->SetLineColor(kRed);
	mPhiD->Fit("fcnRhoOmega","L0 EMR");
	float chi2   = fcnRhoOmega->GetChisquare();
	float ndf    = fcnRhoOmega->GetNDF();
	float pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptRhoOmega = new TPaveText(0.12,0.12,0.37,0.37,"NDC");
	ptRhoOmega->SetFillStyle(4000); //will be transparent
	ptRhoOmega->SetFillColor(0);
	ptRhoOmega->SetTextFont(42);
	ptRhoOmega->SetBorderSize(0);
	ptRhoOmega->AddText("#scale[1.5]{m_{OS1} & m_{OS2}, #rho/#omega:}");
	ptRhoOmega->AddText("#sigma: "+tstr(fcnRhoOmega->GetParameter(4),2)+"#pm"+tstr(fcnRhoOmega->GetParError(4),2)+" MeV");
	ptRhoOmega->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptRhoOmega->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnRhoOmega1 = new TF1("fcnRhoOmega1",func,700,890,npar);
	fcnRhoOmega1->SetParName(0,"Nexpo_rhoOmega1");   /*fcnRhoOmega1->SetParameter(0,...);*/   fcnRhoOmega1->SetParLimits(0,0,1e5);
	fcnRhoOmega1->SetParName(1,"Gexpo_rhoOmega1");   /*fcnRhoOmega1->SetParameter(1,...);*/   fcnRhoOmega1->SetParLimits(1,-10,0);
	fcnRhoOmega1->SetParName(2,"Ngaus_rhoOmega1");   /*fcnRhoOmega1->SetParameter(2,10);*/    fcnRhoOmega1->SetParLimits(2,0,1e4);
	fcnRhoOmega1->SetParName(3,"Mean_rhoOmega1");    fcnRhoOmega1->SetParameter(3,782);       fcnRhoOmega1->SetParLimits(3,750,800);
	fcnRhoOmega1->SetParName(4,"Sigma_rhoOmega1");   fcnRhoOmega1->SetParameter(4,30.);       fcnRhoOmega1->SetParLimits(4,0,100);
	fcnRhoOmega1->SetLineColor(kRed);
	mPhiD1->Fit("fcnRhoOmega1","L0 EMR");
	chi2   = fcnRhoOmega1->GetChisquare();
	ndf    = fcnRhoOmega1->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptRhoOmega1 = new TPaveText(0.12,0.12,0.37,0.37,"NDC");
	ptRhoOmega1->SetFillStyle(4000); //will be transparent
	ptRhoOmega1->SetFillColor(0);
	ptRhoOmega1->SetTextFont(42);
	ptRhoOmega1->SetBorderSize(0);
	ptRhoOmega1->AddText("#scale[1.5]{m_{OS1}, #rho/#omega:}");
	ptRhoOmega1->AddText("#sigma: "+tstr(fcnRhoOmega1->GetParameter(4),2)+"#pm"+tstr(fcnRhoOmega1->GetParError(4),2)+" MeV");
	ptRhoOmega1->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptRhoOmega1->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnRhoOmega2 = new TF1("fcnRhoOmega2",func,700,890,npar);
	fcnRhoOmega2->SetParName(0,"Nexpo_rhoOmega2");   /*fcnRhoOmega2->SetParameter(0,...);*/   fcnRhoOmega2->SetParLimits(0,0,1e5);
	fcnRhoOmega2->SetParName(1,"Gexpo_rhoOmega2");   /*fcnRhoOmega2->SetParameter(1,...);*/   fcnRhoOmega2->SetParLimits(1,-10,0);
	fcnRhoOmega2->SetParName(2,"Ngaus_rhoOmega2");   /*fcnRhoOmega2->SetParameter(2,10);*/    fcnRhoOmega2->SetParLimits(2,0,1e4);
	fcnRhoOmega2->SetParName(3,"Mean_rhoOmega2");    fcnRhoOmega2->SetParameter(3,782);       fcnRhoOmega2->SetParLimits(3,750,800);
	fcnRhoOmega2->SetParName(4,"Sigma_rhoOmega2");   fcnRhoOmega2->SetParameter(4,30.);       fcnRhoOmega2->SetParLimits(4,0,100);
	fcnRhoOmega2->SetLineColor(kRed);
	mPhiD2->Fit("fcnRhoOmega2","L0 EMR");
	chi2   = fcnRhoOmega2->GetChisquare();
	ndf    = fcnRhoOmega2->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptRhoOmega2 = new TPaveText(0.12,0.12,0.37,0.37,"NDC");
	ptRhoOmega2->SetFillStyle(4000); //will be transparent
	ptRhoOmega2->SetFillColor(0);
	ptRhoOmega2->SetTextFont(42);
	ptRhoOmega2->SetBorderSize(0);
	ptRhoOmega2->AddText("#scale[1.5]{m_{OS2}, #rho/#omega:}");
	ptRhoOmega2->AddText("#sigma: "+tstr(fcnRhoOmega2->GetParameter(4),2)+"#pm"+tstr(fcnRhoOmega2->GetParError(4),2)+" MeV");
	ptRhoOmega2->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptRhoOmega2->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnRhoOmega3 = new TF1("fcnRhoOmega3",func,700,890,npar);
	fcnRhoOmega3->SetParName(0,"Nexpo_rhoOmega3");   /*fcnRhoOmega3->SetParameter(0,...);*/   fcnRhoOmega3->SetParLimits(0,0,1e5);
	fcnRhoOmega3->SetParName(1,"Gexpo_rhoOmega3");   /*fcnRhoOmega3->SetParameter(1,...);*/   fcnRhoOmega3->SetParLimits(1,-10,0);
	fcnRhoOmega3->SetParName(2,"Ngaus_rhoOmega3");   /*fcnRhoOmega3->SetParameter(2,10);*/    fcnRhoOmega3->SetParLimits(2,0,1e4);
	fcnRhoOmega3->SetParName(3,"Mean_rhoOmega3");    fcnRhoOmega3->SetParameter(3,782);       fcnRhoOmega3->SetParLimits(3,750,800);
	fcnRhoOmega3->SetParName(4,"Sigma_rhoOmega3");   fcnRhoOmega3->SetParameter(4,30.);       fcnRhoOmega3->SetParLimits(4,0,100);
	fcnRhoOmega3->SetLineColor(kRed);
	mPhiD3->Fit("fcnRhoOmega3","L0 EMR");
	chi2   = fcnRhoOmega3->GetChisquare();
	ndf    = fcnRhoOmega3->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptRhoOmega3 = new TPaveText(0.12,0.12,0.37,0.37,"NDC");
	ptRhoOmega3->SetFillStyle(4000); //will be transparent
	ptRhoOmega3->SetFillColor(0);
	ptRhoOmega3->SetTextFont(42);
	ptRhoOmega3->SetBorderSize(0);
	ptRhoOmega3->AddText("#scale[1.5]{m_{OS2}, #rho/#omega:}");
	ptRhoOmega3->AddText("#sigma: "+tstr(fcnRhoOmega3->GetParameter(4),2)+"#pm"+tstr(fcnRhoOmega3->GetParError(4),2)+" MeV");
	ptRhoOmega3->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptRhoOmega3->AddText("#it{p}-value: "+tstr(pvalue,5));



	TF1* fcnPhi = new TF1("fcnPhi",func,910,1200,npar);
	fcnPhi->SetParName(0,"Nexpo_phi");   /*fcnPhi->SetParameter(0,...);*/   fcnPhi->SetParLimits(0,0,1e5);
	fcnPhi->SetParName(1,"Gexpo_phi");   /*fcnPhi->SetParameter(1,...);*/   fcnPhi->SetParLimits(1,-10,0);
	fcnPhi->SetParName(2,"Ngaus_phi");   /*fcnPhi->SetParameter(2,10);*/    fcnPhi->SetParLimits(2,0,1e4);
	fcnPhi->SetParName(3,"Mean_phi");    fcnPhi->SetParameter(3,1020);      /*fcnPhi->SetParLimits(3,980,1060);*/
	fcnPhi->SetParName(4,"Sigma_phi");   fcnPhi->SetParameter(4,30.);       fcnPhi->SetParLimits(4,1,100);
	fcnPhi->SetLineColor(kRed);
	mPhiD->Fit("fcnPhi","L0 EMR");
	chi2   = fcnPhi->GetChisquare();
	ndf    = fcnPhi->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptPhi = new TPaveText(0.6,0.5,0.85,0.75,"NDC");
	ptPhi->SetFillStyle(4000); //will be transparent
	ptPhi->SetFillColor(0);
	ptPhi->SetTextFont(42);
	ptPhi->SetBorderSize(0);
	ptPhi->AddText("#scale[1.5]{m_{OS1} & m_{OS2}, #phi:}");
	ptPhi->AddText("#sigma: "+tstr(fcnPhi->GetParameter(4),2)+"#pm"+tstr(fcnPhi->GetParError(4),2)+" MeV");
	ptPhi->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptPhi->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnPhi1 = new TF1("fcnPhi1",func,910,1200,npar);
	fcnPhi1->SetParName(0,"Nexpo_phi1");   /*fcnPhi1->SetParameter(0,...);*/   fcnPhi1->SetParLimits(0,0,1e5);
	fcnPhi1->SetParName(1,"Gexpo_phi1");   /*fcnPhi1->SetParameter(1,...);*/   fcnPhi1->SetParLimits(1,-10,0);
	fcnPhi1->SetParName(2,"Ngaus_phi1");   /*fcnPhi1->SetParameter(2,10);*/    fcnPhi1->SetParLimits(2,0,1e4);
	fcnPhi1->SetParName(3,"Mean_phi1");    fcnPhi1->SetParameter(3,1020);      /*fcnPhi1->SetParLimits(3,980,1060);*/
	fcnPhi1->SetParName(4,"Sigma_phi1");   fcnPhi1->SetParameter(4,30.);       fcnPhi1->SetParLimits(4,1,100);
	fcnPhi1->SetLineColor(kRed);
	mPhiD1->Fit("fcnPhi1","L0 EMR");
	chi2   = fcnPhi1->GetChisquare();
	ndf    = fcnPhi1->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptPhi1 = new TPaveText(0.6,0.5,0.85,0.75,"NDC");
	ptPhi1->SetFillStyle(4000); //will be transparent
	ptPhi1->SetFillColor(0);
	ptPhi1->SetTextFont(42);
	ptPhi1->SetBorderSize(0);
	ptPhi1->AddText("#scale[1.5]{m_{OS1}, #phi1:}");
	ptPhi1->AddText("#sigma: "+tstr(fcnPhi1->GetParameter(4),2)+"#pm"+tstr(fcnPhi1->GetParError(4),2)+" MeV");
	ptPhi1->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptPhi1->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnPhi2 = new TF1("fcnPhi2",func,910,1200,npar);
	fcnPhi2->SetParName(0,"Nexpo_phi2");   /*fcnPhi2->SetParameter(0,...);*/   fcnPhi2->SetParLimits(0,0,1e5);
	fcnPhi2->SetParName(1,"Gexpo_phi2");   /*fcnPhi2->SetParameter(1,...);*/   fcnPhi2->SetParLimits(1,-10,0);
	fcnPhi2->SetParName(2,"Ngaus_phi2");   /*fcnPhi2->SetParameter(2,10);*/    fcnPhi2->SetParLimits(2,0,1e4);
	fcnPhi2->SetParName(3,"Mean_phi2");    fcnPhi2->SetParameter(3,1020);      /*fcnPhi2->SetParLimits(3,980,1060);*/
	fcnPhi2->SetParName(4,"Sigma_phi2");   fcnPhi2->SetParameter(4,30.);       fcnPhi2->SetParLimits(4,1,100);
	fcnPhi2->SetLineColor(kRed);
	mPhiD2->Fit("fcnPhi2","L0 EMR");
	chi2   = fcnPhi2->GetChisquare();
	ndf    = fcnPhi2->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptPhi2 = new TPaveText(0.6,0.5,0.85,0.75,"NDC");
	ptPhi2->SetFillStyle(4000); //will be transparent
	ptPhi2->SetFillColor(0);
	ptPhi2->SetTextFont(42);
	ptPhi2->SetBorderSize(0);
	ptPhi2->AddText("#scale[1.5]{m_{OS2}, #phi2:}");
	ptPhi2->AddText("#sigma: "+tstr(fcnPhi2->GetParameter(4),2)+"#pm"+tstr(fcnPhi2->GetParError(4),2)+" MeV");
	ptPhi2->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptPhi2->AddText("#it{p}-value: "+tstr(pvalue,5));

	TF1* fcnPhi3 = new TF1("fcnPhi3",func,910,1200,npar);
	fcnPhi3->SetParName(0,"Nexpo_phi2");   /*fcnPhi3->SetParameter(0,...);*/   fcnPhi3->SetParLimits(0,0,1e5);
	fcnPhi3->SetParName(1,"Gexpo_phi2");   /*fcnPhi3->SetParameter(1,...);*/   fcnPhi3->SetParLimits(1,-10,0);
	fcnPhi3->SetParName(2,"Ngaus_phi2");   /*fcnPhi3->SetParameter(2,10);*/    fcnPhi3->SetParLimits(2,0,1e4);
	fcnPhi3->SetParName(3,"Mean_phi2");    fcnPhi3->SetParameter(3,1020);      /*fcnPhi3->SetParLimits(3,980,1060);*/
	fcnPhi3->SetParName(4,"Sigma_phi2");   fcnPhi3->SetParameter(4,30.);       fcnPhi3->SetParLimits(4,1,100);
	fcnPhi3->SetLineColor(kRed);
	mPhiD3->Fit("fcnPhi3","L0 EMR");
	chi2   = fcnPhi3->GetChisquare();
	ndf    = fcnPhi3->GetNDF();
	pvalue = TMath::Prob(chi2,ndf);
	cout << "fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
	TPaveText *ptPhi3 = new TPaveText(0.6,0.5,0.85,0.75,"NDC");
	ptPhi3->SetFillStyle(4000); //will be transparent
	ptPhi3->SetFillColor(0);
	ptPhi3->SetTextFont(42);
	ptPhi3->SetBorderSize(0);
	ptPhi3->AddText("#scale[1.5]{m_{OS2}, #phi2:}");
	ptPhi3->AddText("#sigma: "+tstr(fcnPhi3->GetParameter(4),2)+"#pm"+tstr(fcnPhi3->GetParError(4),2)+" MeV");
	ptPhi3->AddText("#chi^{2}_{DOF}: "+tstr(chi2/ndf,3));
	ptPhi3->AddText("#it{p}-value: "+tstr(pvalue,5));



	TCanvas* cnv = new TCanvas("cnv","",600,400);
	gPad->SetTicks(1,1);
	mPhiD->SetMarkerStyle(20); mPhiD->SetMarkerSize(0.8); mPhiD->SetMarkerColor(kBlack); mPhiD->Draw("p0 x0");
	fcnRhoOmega->Draw("same");
	fcnPhi->Draw("same");
	ptRhoOmega->Draw("same");
	ptPhi->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFit.png");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFit.eps");
	cnv->SaveAs("figures/phitest.pdf(");

	delete cnv;
	cnv = new TCanvas("cnv","",600,400);
	gPad->SetTicks(1,1);
	mPhiD1->SetMarkerStyle(20); mPhiD1->SetMarkerSize(0.8); mPhiD1->SetMarkerColor(kBlack); mPhiD1->Draw("p0 x0");
	fcnRhoOmega1->Draw("same");
	fcnPhi1->Draw("same");
	ptRhoOmega1->Draw("same");
	ptPhi1->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitOS1.png");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitOS1.eps");
	cnv->SaveAs("figures/phitest.pdf");

	delete cnv;
	cnv = new TCanvas("cnv","",600,400);
	gPad->SetTicks(1,1);
	mPhiD2->SetMarkerStyle(20); mPhiD2->SetMarkerSize(0.8); mPhiD2->SetMarkerColor(kBlack); mPhiD2->Draw("p0 x0");
	fcnRhoOmega2->Draw("same");
	fcnPhi2->Draw("same");
	ptRhoOmega2->Draw("same");
	ptPhi2->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitOS2.png");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitOS2.eps");
	cnv->SaveAs("figures/phitest.pdf");

	delete cnv;
	cnv = new TCanvas("cnv","",600,400);
	gPad->SetTicks(1,1);
	mPhiD3->SetMarkerStyle(20); mPhiD3->SetMarkerSize(0.8); mPhiD3->SetMarkerColor(kBlack); mPhiD3->Draw("p0 x0");
	fcnRhoOmega3->Draw("same");
	fcnPhi3->Draw("same");
	ptRhoOmega3->Draw("same");
	ptPhi3->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitSS.png");
	cnv->SaveAs("figures/RhoOmegaPhiShapeFitSS.eps");
	cnv->SaveAs("figures/phitest.pdf");






	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mSS_vs_pcnHits1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mSS_vs_pcnHits1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mSS_vs_pcnHits2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mSS_vs_pcnHits2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mSS_vs_pcnHits3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mSS_vs_pcnHits3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");
	
	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mOS1_vs_pcnHits1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mOS1_vs_pcnHits1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mOS1_vs_pcnHits2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mOS1_vs_pcnHits2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mOS1_vs_pcnHits3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mOS1_vs_pcnHits3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");
	
	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mOS2_vs_pcnHits1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mOS2_vs_pcnHits1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mOS2_vs_pcnHits2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mOS2_vs_pcnHits2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mOS2_vs_pcnHits3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mOS2_vs_pcnHits3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");






	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mSS_vs_etaphiLayers1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mSS_vs_etaphiLayers1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mSS_vs_etaphiLayers2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mSS_vs_etaphiLayers2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mSS_vs_etaphiLayers3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mSS_vs_etaphiLayers3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");
	
	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mOS1_vs_etaphiLayers3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");
	
	delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(2,3);
	cnv->cd(1); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers1D->Draw("col");
	cnv->cd(2); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers1S->Draw("col");
	cnv->cd(3); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers2D->Draw("col");
	cnv->cd(4); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers2S->Draw("col");
	cnv->cd(5); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers3D->Draw("col");
	cnv->cd(6); gPad->SetTicks(1,1); mOS2_vs_etaphiLayers3S->Draw("col");
	cnv->SaveAs("figures/phitest.pdf");
	
	
	
	
	
	
	


	delete cnv;
	cnv = new TCanvas("cnv","",600,400);
	cnv->SaveAs("figures/phitest.pdf)");
}
