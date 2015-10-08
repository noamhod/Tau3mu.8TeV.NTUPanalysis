/////////////////////////////////////////
//// root -b -l -q sidebandsBias.C++ ////
/////////////////////////////////////////
#include "std.h"
#include "const.h"
#include "type.h"
#include "postBDTcuts.h"
#include <iterator>
#include <cstdlib>

bool blind = false;
bool doBkg = false;

Double_t m3bodyMin = mTrainingMinGlob;
Double_t m3bodyMax = mTrainingMaxGlob;
Double_t xbmin  = mBlindMinGlob;
Double_t xbmax  = mBlindMaxGlob;
Double_t xbdtmin  = minBDTcut;
Double_t xbdtmax  = +1.0;
Double_t xbdtopt  = optBDTcut;
Int_t    nbdtbins = 30;

const int     nmbins = 8;
const Double_t mbins[nmbins+1] = {m3bodyMin,980,1220,1450,1690,1870,2110,2320,m3bodyMax};

const int     nmbinsSB = 3;
const Double_t mbinsSB[nmbinsSB+1] = {m3bodyMin,1450,2110,m3bodyMax};

const int     nmbinsSR = 3;
const Double_t mbinsSR[nmbinsSR+1] = {m3bodyMin,1713,1841,m3bodyMax};

const int     nmbinsBL = 3;
const Double_t mbinsBL[nmbinsBL+1] = {m3bodyMin,1690,1870,m3bodyMax};

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

void reBL(TH1* h)
{	
	Int_t bMin = h->FindBin(xbmin);
	Int_t bMax = h->FindBin(xbmax);
	for(Int_t b=bMin ; b<=bMax && b>0 ; b++)
	{	
		if(h->GetBinContent(b)>0)
		{
			h->SetBinContent(b,0);
			h->SetBinError(b,0);
		}
	}
}

float Sum(TH1* h, bool addUunderFlow=false, bool addOverFlow=false)
{
	float I=0.;
	for(int i=1 ; i<=h->GetNbinsX() ; i++)// without under- and over-flow
	{
		I += h->GetBinContent(i);
	}
	if(addUunderFlow) I+=h->GetBinContent(0);
	if(addOverFlow)   I+=h->GetBinContent(h->GetNbinsX()+1);
	return I;
}

void NormToEntries(TH1* h)
{
	Double_t entries = Sum(h);
	Double_t scale = (entries==0.) ? 1. : 1./entries;
	for(int i=1 ; i<=h->GetNbinsX() ; i++)// without under- and over-flow
	{
		h->SetBinContent(i,h->GetBinContent(i)*scale);
		h->SetBinError(i,h->GetBinError(i)*scale);
	}
	// h->Scale(scale);
	// h->GetYaxis()->SetTitle("Normalized");
}

void Residuals(TH1* h1, TH1* h2, TH1* hR)
{
	for(Int_t b=1 ; b<=h1->GetNbinsX() ; b++)
	{
		float x = h1->GetBinContent(b); float dx = h1->GetBinError(b);
		float y = h2->GetBinContent(b); float dy = h2->GetBinError(b);
		if(y<=0 || x<=0) { /*hR->SetBinContent(b,0);*/ hR->SetBinError(b,0); continue; }
		float residual = fabs(x-y)/y;
		float error    = (x/y)*sqrt((dx/x)*(dx/x) + (dy/y)*(dy/y));
		hR->SetBinContent(b,residual);
		hR->SetBinError(b,error);
	}
}

void Ratio(TH1* h1, TH1* h2, TH1* hR)
{
	for(Int_t b=1 ; b<=h1->GetNbinsX() ; b++)
	{
		float x = h1->GetBinContent(b); float dx = h1->GetBinError(b);
		float y = h2->GetBinContent(b); float dy = h2->GetBinError(b);
		if(y<=0 || x<=0) { /*hR->SetBinContent(b,0);*/ hR->SetBinError(b,0); continue; }
		float rat   = x/y;
		float error = (x/y)*sqrt((dx/x)*(dx/x) + (dy/y)*(dy/y));
		hR->SetBinContent(b,rat);
		hR->SetBinError(b,error);
	}
}

float getDistErr(float da, float db)
{
	return sqrt(da*da+db*db);
}
float getRatErr(float a, float b, float da)
{
	// return (a/b)*sqrt((da/a)*(da/a)+(db/b)*(db/b));
	return (b!=0) ? da/b : 0.;
}

TMapTSf gaus, dgaus, slopes, dslopes, distances, ddistances, ratios, dratios;
void plot(TString name, TString title,
		  TH1* hSBres, TH1* hBLres, TH1* hSRres,
		  TH1* hSBrat, TH1* hBLrat, TH1* hSRrat,
		  TH1* mDres, TH1* mDrat,
		  TH1* mBres, TH1* mBrat,
		  TString pdfstatus="")
{
	
	TCanvas* cnv = new TCanvas("","",1200,400); cnv->Draw(); cnv->cd(); cnv->Divide(2,1);
	
	float xmin = mDres->GetXaxis()->GetXmin();
	float xmax = mDres->GetXaxis()->GetXmax();
	
	cnv->cd(1); gPad->SetTicks(1,1);
	hSBres->SetTitle("Residuals "+title);
	hSBres->Draw();
	hBLres->Draw("same");
	hSRres->Draw("same");
	TF1* fDres = new TF1("fDres"+name,"[0]+[1]*x+[2]*TMath::Gaus(x,1777,40)",xmin,xmax); TFitResultPtr fitptrDres = mDres->Fit(fDres,"EM0","",xmin,xmax); fDres->SetLineColor(mDres->GetLineColor()); fDres->SetLineWidth(2);
	TF1* fDres1 = new TF1("fDres1"+name,"[0]+[1]*x",xmin,xmax); fDres1->SetParameter(0,fDres->GetParameter(0)); fDres1->SetParameter(1,fDres->GetParameter(1)); fDres1->SetLineColor(kGreen); fDres1->SetLineWidth(2); fDres1->SetLineStyle(2);
	TF1* fBres = new TF1("fBres"+name,"[0]+[1]*x+[2]*TMath::Gaus(x,1777,40)",xmin,xmax); TFitResultPtr fitptrBres = mBres->Fit(fBres,"EM0","",xmin,xmax); fBres->SetLineColor(mBres->GetLineColor()); fBres->SetLineWidth(2); 
	TF1* fBres1 = new TF1("fBres1"+name,"[0]+[1]*x",xmin,xmax); fBres1->SetParameter(0,fBres->GetParameter(0)); fBres1->SetParameter(1,fBres->GetParameter(1)); fBres1->SetLineColor(kGreen); fBres1->SetLineWidth(2); fBres1->SetLineStyle(2);
	
	mDres->Draw("p0 e0 same"); fDres1->Draw("same"); fDres->Draw("same");
	if(doBkg) { mBres->Draw("p0 e0 same"); fBres1->Draw("same"); fBres->Draw("same"); }
	ptxt->Draw("smae");
	gPad->Update();
	gPad->RedrawAxis();
	
	cnv->cd(2); gPad->SetTicks(1,1);
	hSBrat->SetTitle("Ratio "+title);
	hSBrat->Draw();
	hBLrat->Draw("same");
	hSRrat->Draw("same");
	TF1* fDrat = new TF1("fDrat"+name,"[0]+[1]*x+[2]*TMath::Gaus(x,1777,40)",xmin,xmax); TFitResultPtr fitptrDrat = mDrat->Fit(fDrat,"EM0","",xmin,xmax); fDrat->SetLineColor(mDrat->GetLineColor()); fDrat->SetLineWidth(2);
	TF1* fDrat1 = new TF1("fDrat1"+name,"[0]+[1]*x",xmin,xmax); fDrat1->SetParameter(0,fDrat->GetParameter(0)); fDrat1->SetParameter(1,fDrat->GetParameter(1)); fDrat1->SetLineColor(kGreen); fDrat1->SetLineWidth(2); fDrat1->SetLineStyle(2);
	TF1* fBrat = new TF1("fBrat"+name,"[0]+[1]*x+[2]*TMath::Gaus(x,1777,40)",xmin,xmax); TFitResultPtr fitptrBrat = mBrat->Fit(fBrat,"EM0","",xmin,xmax); fBrat->SetLineColor(mBrat->GetLineColor()); fBrat->SetLineWidth(2);
	TF1* fBrat1 = new TF1("fBrat1"+name,"[0]+[1]*x",xmin,xmax); fBrat1->SetParameter(0,fBrat->GetParameter(0)); fBrat1->SetParameter(1,fBrat->GetParameter(1)); fBrat1->SetLineColor(kGreen); fBrat1->SetLineWidth(2); fBrat1->SetLineStyle(2);
	mDrat->Draw("p0 e0 same");  fDrat1->Draw("same"); fDrat->Draw("same");
	if(doBkg) { mBrat->Draw("p0 e0 same");  fBrat1->Draw("same"); fBrat->Draw("same"); }
	ptxt->Draw("smae");
	gPad->Update();
	gPad->RedrawAxis();

	gaus.insert(make_pair(mDres->GetName(),fDres->GetParameter(2))); dgaus.insert(make_pair(mDres->GetName(),fDres->GetParError(2)));  
	gaus.insert(make_pair(mBres->GetName(),fBres->GetParameter(2))); dgaus.insert(make_pair(mBres->GetName(),fBres->GetParError(2)));
	gaus.insert(make_pair(mDrat->GetName(),fDrat->GetParameter(2))); dgaus.insert(make_pair(mDrat->GetName(),fDrat->GetParError(2)));  
	gaus.insert(make_pair(mBrat->GetName(),fBrat->GetParameter(2))); dgaus.insert(make_pair(mBrat->GetName(),fBrat->GetParError(2)));
	
	slopes.insert(make_pair(mDres->GetName(),fDres->GetParameter(1))); dslopes.insert(make_pair(mDres->GetName(),fDres->GetParError(1)));  
	slopes.insert(make_pair(mBres->GetName(),fBres->GetParameter(1))); dslopes.insert(make_pair(mBres->GetName(),fBres->GetParError(1)));
	slopes.insert(make_pair(mDrat->GetName(),fDrat->GetParameter(1))); dslopes.insert(make_pair(mDrat->GetName(),fDrat->GetParError(1)));  
	slopes.insert(make_pair(mBrat->GetName(),fBrat->GetParameter(1))); dslopes.insert(make_pair(mBrat->GetName(),fBrat->GetParError(1)));
	
	distances.insert(make_pair(mDres->GetName(),mDres->GetBinContent(4)-fDres->Eval(mDres->GetBinCenter(4))));  ddistances.insert(make_pair(mDres->GetName(),getDistErr(mDres->GetBinContent(4),fDres->Eval(mDres->GetBinCenter(4)))));
	distances.insert(make_pair(mBres->GetName(),mBres->GetBinContent(4)-fBres->Eval(mBres->GetBinCenter(4))));  ddistances.insert(make_pair(mBres->GetName(),getDistErr(mBres->GetBinContent(4),fBres->Eval(mBres->GetBinCenter(4)))));
	distances.insert(make_pair(mDrat->GetName(),mDrat->GetBinContent(4)-fDrat->Eval(mDrat->GetBinCenter(4))));  ddistances.insert(make_pair(mDrat->GetName(),getDistErr(mDrat->GetBinContent(4),fDrat->Eval(mDrat->GetBinCenter(4)))));
	distances.insert(make_pair(mBrat->GetName(),mBrat->GetBinContent(4)-fBrat->Eval(mBrat->GetBinCenter(4))));  ddistances.insert(make_pair(mBrat->GetName(),getDistErr(mBrat->GetBinContent(4),fBrat->Eval(mBrat->GetBinCenter(4)))));
	
	ratios.insert(make_pair(mDres->GetName(),mDres->GetBinContent(4)/fDres->Eval(mDres->GetBinCenter(4))));  ratios.insert(make_pair(mDres->GetName(), getRatErr(mDres->GetBinContent(4),fDres->Eval(mDres->GetBinCenter(4)),mDres->GetBinError(4)) ));
	ratios.insert(make_pair(mBres->GetName(),mBres->GetBinContent(4)/fBres->Eval(mBres->GetBinCenter(4))));  ratios.insert(make_pair(mBres->GetName(), getRatErr(mBres->GetBinContent(4),fBres->Eval(mBres->GetBinCenter(4)),mBres->GetBinError(4)) ));
	ratios.insert(make_pair(mDrat->GetName(),mDrat->GetBinContent(4)/fDrat->Eval(mDrat->GetBinCenter(4))));  ratios.insert(make_pair(mDrat->GetName(), getRatErr(mDrat->GetBinContent(4),fDrat->Eval(mDrat->GetBinCenter(4)),mDrat->GetBinError(4)) )); 
	ratios.insert(make_pair(mBrat->GetName(),mBrat->GetBinContent(4)/fBrat->Eval(mBrat->GetBinCenter(4))));  ratios.insert(make_pair(mBrat->GetName(), getRatErr(mBrat->GetBinContent(4),fBrat->Eval(mBrat->GetBinCenter(4)),mBrat->GetBinError(4)) ));
	
	cnv->SaveAs("figures/sidebandsBias."+name+".png"+pdfstatus);
	cnv->SaveAs("figures/sidebandsBias."+name+".pdf"+pdfstatus);
	cnv->SaveAs("figures/sidebandsBias."+name+".eps"+pdfstatus);
	cnv->SaveAs("figures/sidebandsBias.pdf"+pdfstatus);
	
	delete cnv; cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hSBres->Draw();
	hBLres->Draw("same");
	hSRres->Draw("same");
	mDres->Draw("p0 e0 same"); fDres1->Draw("same"); fDres->Draw("same");
	if(doBkg) { mBres->Draw("p0 e0 same"); fBres1->Draw("same"); fBres->Draw("same"); }
	ptxt->Draw("smae");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.residuals."+name+".png");
	cnv->SaveAs("figures/sidebandsBias.residuals."+name+".pdf");
	cnv->SaveAs("figures/sidebandsBias.residuals."+name+".eps");
	
	delete cnv; cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hSBrat->Draw();
	hBLrat->Draw("same");
	hSRrat->Draw("same");
	mDrat->Draw("p0 e0 same"); fDrat1->Draw("same"); fDrat->Draw("same");
	if(doBkg) { mBrat->Draw("p0 e0 same"); fBrat1->Draw("same"); fBrat->Draw("same"); }
	ptxt->Draw("smae");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figurat/sidebandsBias.ratio."+name+".png");
	cnv->SaveAs("figurat/sidebandsBias.ratio."+name+".pdf");
	cnv->SaveAs("figurat/sidebandsBias.ratio."+name+".eps");
}

void addHist(TMapTSP2TH1& histos, TString name, TString titles)
{
	histos.insert(make_pair(name, new TH1F(name,titles,nmbins,mbins)));
	histos[name]->Sumw2();
	if(name.Contains("mD"))
	{
		histos[name]->SetLineColor(kBlack);
		histos[name]->SetLineWidth(1);
		histos[name]->SetLineStyle(1);
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(0.8);
		histos[name]->SetMarkerColor(kBlack);
		
		histos.insert(make_pair(name+"rat", (TH1F*)histos[name]->Clone(name+"rat")));
		histos.insert(make_pair(name+"res", (TH1F*)histos[name]->Clone(name+"res")));
		histos[name+"rat"]->Reset();
		histos[name+"res"]->Reset();
	}
	else
	{
		histos[name]->SetLineColor(kRed);
		histos[name]->SetLineWidth(1);
		histos[name]->SetLineStyle(1);
		histos[name]->SetMarkerStyle(21);
		histos[name]->SetMarkerSize(0.8);
		histos[name]->SetMarkerColor(kRed);
		
		histos.insert(make_pair(name+"rat", (TH1F*)histos[name]->Clone(name+"rat")));
		histos.insert(make_pair(name+"res", (TH1F*)histos[name]->Clone(name+"res")));
		histos[name+"rat"]->Reset();
		histos[name+"res"]->Reset();
	}	
}


void sidebandsBias()
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
	
	makeAtlasLabel();
	
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");
	TTree* tB = (TTree*)fmvaout->Get("fltmva_bb_mu4mu4");
	
	TString sm3bodyMin = tstr(m3bodyMin,0);
	TString sm3bodyMax = tstr(m3bodyMax,0);
	TString smBLMin = tstr(xbmin,0);
	TString smBLMax = tstr(xbmax,0);	
	
	Double_t chi2 = 0;
	Double_t pvalue = 0;
	Int_t ndf  = 0;
	
	TH1F* hSBrat = new TH1F("hSBrat", ";m_{3body} [MeV];Ratio", nmbinsSB,mbinsSB); hSBrat->SetFillColor(kGray+2);  hSBrat->SetLineColor(kGray+2);
	TH1F* hBLrat = new TH1F("hBLrat", ";m_{3body} [MeV];Ratio", nmbinsBL,mbinsBL); hBLrat->SetFillColor(kGray+1);  hBLrat->SetLineColor(kGray+1);
	TH1F* hSRrat = new TH1F("hSRrat", ";m_{3body} [MeV];Ratio", nmbinsSR,mbinsSR); hSRrat->SetFillColor(kGray+0);  hSRrat->SetLineColor(kGray+0);

	TH1F* hSBres = new TH1F("hSBres", ";m_{3body} [MeV];Residuals", nmbinsSB,mbinsSB); hSBres->SetFillColor(kGray+2); hSBres->SetLineColor(kGray+2);
	TH1F* hBLres = new TH1F("hBLres", ";m_{3body} [MeV];Residuals", nmbinsBL,mbinsBL); hBLres->SetFillColor(kGray+1); hBLres->SetLineColor(kGray+1);
	TH1F* hSRres = new TH1F("hSRres", ";m_{3body} [MeV];Residuals", nmbinsSR,mbinsSR); hSRres->SetFillColor(kGray+0); hSRres->SetLineColor(kGray+0);
	
	TMapTSP2TH1 histos1;
	TMapTSTS    cuts;
	
	TMapfTS xnames;
	xnames.insert(make_pair(-1.00, "neg100"));
	xnames.insert(make_pair(-0.90, "neg090")); 
	// xnames.insert(make_pair(-0.80, "neg080")); 
	xnames.insert(make_pair(-0.70, "neg070")); 
	// xnames.insert(make_pair(-0.60, "neg060")); 
	xnames.insert(make_pair(-0.50, "neg050")); 
	// xnames.insert(make_pair(-0.40, "neg040")); 
	xnames.insert(make_pair(-0.30, "neg030")); 
	// xnames.insert(make_pair(-0.20, "neg020")); 
	xnames.insert(make_pair(-0.10, "neg010")); 
	// xnames.insert(make_pair( 0.00, "zro000")); 
	xnames.insert(make_pair(+0.10, "pos010")); 
	// xnames.insert(make_pair(+0.20, "pos020")); 
	xnames.insert(make_pair(+0.30, "pos030")); 
	// xnames.insert(make_pair(+0.40, "pos040")); 
	xnames.insert(make_pair(+0.50, "pos050")); 
	// xnames.insert(make_pair(+0.60, "pos060"));
	xnames.insert(make_pair(+0.70, "pos070")); 
	// xnames.insert(make_pair(+0.80, "pos080")); 
	xnames.insert(make_pair(+0.90, "pos090")); 
	xnames.insert(make_pair(+0.95, "pos095")); 
	
	const int     nbdtbins = 11;
	const Double_t bdtbins[nbdtbins+1] = {-1,-0.9,-0.7,-0.5,-0.3,-0.1,+0.1,+0.3,+0.5,+0.7,+0.9,+0.95};
	
	
	for(TMapfTS::iterator it=xnames.begin() ; it!=xnames.end() ; ++it)
	{
		TString xname = it->second;
		TString xtitle = tstr(it->first,2);
		
		addHist(histos1,"mD"+xname,"x_{1}>"+xtitle+";m_{3body} [MeV];Events");
		addHist(histos1,"mB"+xname,"x_{1}>"+xtitle+";m_{3body} [MeV];Events");
		
		cout << "Adding: hname=" << "mD"+xname << ", htitles=" << "x_{1}>"+xtitle+";m_{3body} [MeV];Events" << endl;
		
		TString cutD = postBDTcut(sm3bodyMin,smBLMin,smBLMax,sm3bodyMax,xtitle,xtitle,"bkg","-1","-1",blind,true,"postTraining");
		// cutD.ReplaceAll("(@m3body->size()==1) && (m3body>750. && m3body<2500.) && (score>-1.00 && score>-1.00) && ((iso020<0.3 && iso030<1) && (trkspval>1.e-9) && (lxysig>-10 && lxysig<50 && a0xysig<25) && ((mettrk>10000 && mettrk<250000 && mttrk>20000) && (metcal>10000 && metcal<250000 && mtcal>20000)) && (mOS1>220. && mOS2>220. && mSS>220.)) && (lxysig>1 && trkspval>8.e-9 && pval>0.2) && (mttrk>60000 && mtcal>60000 && dphihttrk>2 && dphihtcal>2) && !( ((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50) || (TMath::Abs(mOS1-782)<50 || TMath::Abs(mOS2-782)<50)) && (mettrk<35000 || metcal<35000 || pt3body<35000) ) && (mOS1>300. && mOS2>300. && mSS>300.) && !((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50) && TMath::Abs(m3body-1968)<100)","");
		// if(!blind) cutD.ReplaceAll(" && (m3body<1690. || m3body>1870.)","");
		cout << "cutD: " << cutD << endl;
		// TString cutB = "(score>"+xtitle+") && (pval>0.2 && lxysig>1 && trkspval>8.e-9)";
		TString cutB = "(score>"+xtitle+")";
		cout << "cutB: " << cutB << "\n" << endl;
		
		cuts.insert(make_pair("mD"+xname,cutD));
		cuts.insert(make_pair("mB"+xname,cutB));
	}
	

	float ymin = 0;
	float ymax = +2;
	float ymid = 1;
	
	hSBrat->SetBinContent(2,ymax*0.9); hSBrat->SetMinimum(ymin); hSBrat->SetMaximum(ymax);
	hBLrat->SetBinContent(2,ymax*0.7); hBLrat->SetMinimum(ymin); hBLrat->SetMaximum(ymax);
	hSRrat->SetBinContent(2,ymax*0.5); hSRrat->SetMinimum(ymin); hSRrat->SetMaximum(ymax);
	
	hSBres->SetBinContent(2,ymax*0.9); hSBres->SetMinimum(ymin); hSBres->SetMaximum(ymax);
	hBLres->SetBinContent(2,ymax*0.7); hBLres->SetMinimum(ymin); hBLres->SetMaximum(ymax);
	hSRres->SetBinContent(2,ymax*0.5); hSRres->SetMinimum(ymin); hSRres->SetMaximum(ymax);
	
	TCanvas* cnv = new TCanvas("","",1200,400); cnv->Draw(); cnv->cd();
	cnv->SaveAs("figures/sidebandsBias.pdf(");
	
	TString prevxname  = "";
	TString prevxtitle = "";
	for(TMapfTS::iterator it=xnames.begin() ; it!=xnames.end() ; ++it)
	{
		TString xname = it->second;
		TString xtitle = tstr(it->first,2);
		
		tD->Draw("m3body>>mD"+xname, cuts["mD"+xname]);
		tB->Draw("m3body>>mB"+xname, cuts["mB"+xname]);
		
		if(it==xnames.begin())
		{
			prevxname  = xname;
			prevxtitle = xtitle;
			continue;
		}
		Ratio(histos1["mD"+xname],histos1["mD"+prevxname],histos1["mD"+xname+"rat"]);
		Ratio(histos1["mB"+xname],histos1["mB"+prevxname],histos1["mB"+xname+"rat"]);
		Residuals(histos1["mD"+xname],histos1["mD"+prevxname],histos1["mD"+xname+"res"]);
		Residuals(histos1["mB"+xname],histos1["mB"+prevxname],histos1["mB"+xname+"res"]);
		
		histos1["mD"+xname+"rat"]->SetMinimum(ymin);
		histos1["mD"+xname+"res"]->SetMinimum(ymin);
		histos1["mB"+xname+"rat"]->SetMinimum(ymin);
		histos1["mB"+xname+"res"]->SetMinimum(ymin);
		
		histos1["mD"+xname+"rat"]->SetMaximum(ymax);
		histos1["mD"+xname+"res"]->SetMaximum(ymax);
		histos1["mB"+xname+"rat"]->SetMaximum(ymax);
		histos1["mB"+xname+"res"]->SetMaximum(ymax);
		
		cout << "xname=" << xname << ", xtitle=" << xtitle << endl;
		
		plot(xname+"_vs_"+prevxname,"x_{1}>"+xtitle+" vs x_{1}>"+prevxtitle, hSBres,hBLres,hSRres,hSBrat,hBLrat,hSRrat, histos1["mD"+xname+"res"],histos1["mD"+xname+"rat"], histos1["mB"+xname+"res"],histos1["mB"+xname+"rat"]);
		
		prevxname  = xname;
		prevxtitle = xtitle;
	}
	
	
	
	// Int_t nbins = xnames.size();
	TH1F* hResDgaus     = new TH1F("hResDgaus",     "Residuals #muGauss;x_{1} cut;Residuals #muGauss",nbdtbins,bdtbins);
	TH1F* hResBgaus     = new TH1F("hResBgaus",     "Residuals #muGauss;x_{1} cut;Residuals #muGauss",nbdtbins,bdtbins);
	TH1F* hRatDgaus     = new TH1F("hRatDgaus",     "Ratios #muGauss;x_{1} cut;Ratios #muGauss",nbdtbins,bdtbins);
	TH1F* hRatBgaus     = new TH1F("hRatBgaus",     "Ratios #muGauss;x_{1} cut;Ratios #muGauss",nbdtbins,bdtbins);
	TH1F* hResDslopes   = new TH1F("hResDslopes",   "Residuals slopes;x_{1} cut;Residuals slope",nbdtbins,bdtbins);
	TH1F* hResBslopes   = new TH1F("hResBslopes",   "Residuals slopes;x_{1} cut;Residuals slope",nbdtbins,bdtbins);
	TH1F* hRatDslopes   = new TH1F("hRatDslopes",   "Ratios slopes;x_{1} cut;Ratios slope",nbdtbins,bdtbins);
	TH1F* hRatBslopes   = new TH1F("hRatBslopes",   "Ratios slopes;x_{1} cut;Ratios slope",nbdtbins,bdtbins);
	TH1F* hResDdistance = new TH1F("hResDdistance", "Residuals distances;x_{1} cut;Residuals distance",nbdtbins,bdtbins);
	TH1F* hResBdistance = new TH1F("hResBdistance", "Residuals distances;x_{1} cut;Residuals distance",nbdtbins,bdtbins);
	TH1F* hRatDdistance = new TH1F("hRatDdistance", "Ratios distances;x_{1} cut;Ratios distance",nbdtbins,bdtbins);
	TH1F* hRatBdistance = new TH1F("hRatBdistance", "Ratios distances;x_{1} cut;Ratios distance",nbdtbins,bdtbins);
	TH1F* hResDratio    = new TH1F("hResDratio",    "Residuals ratio;x_{1} cut;Residuals ratio",nbdtbins,bdtbins);
	TH1F* hResBratio    = new TH1F("hResBratio",    "Residuals ratio;x_{1} cut;Residuals ratio",nbdtbins,bdtbins);	
	TH1F* hRatDratio    = new TH1F("hRatDratio",    "Ratios ratio;x_{1} cut;Ratios ratio",nbdtbins,bdtbins);
	TH1F* hRatBratio    = new TH1F("hRatBratio",    "Ratios ratio;x_{1} cut;Ratios ratio",nbdtbins,bdtbins);
	
	Int_t b = 1;
	for(TMapfTS::iterator it=xnames.begin() ; it!=xnames.end() ; ++it)
	{
		TString xname = it->second;
		TString xtitle = tstr(it->first,2);
		
		if(it==xnames.begin()) continue;
		
		TString name = "";
		
		// name = "mD"+xname+"res"; cout << name << ": slope=" << slopes[name]  << " +- " << dslopes[name] << "   distance=" << distances[name] << "   ratios=" << ratios[name] << endl;
		// name = "mB"+xname+"res"; cout << name << ": slope=" << slopes[name]  << " +- " << dslopes[name] << "   distance=" << distances[name] << "   ratios=" << ratios[name] << endl;
		// name = "mD"+xname+"rat"; cout << name << ": slope=" << slopes[name]  << " +- " << dslopes[name] << "   distance=" << distances[name] << "   ratios=" << ratios[name] << endl;
		// name = "mB"+xname+"rat"; cout << name << ": slope=" << slopes[name]  << " +- " << dslopes[name] << "   distance=" << distances[name] << "   ratios=" << ratios[name] << endl;
		
		name = "mD"+xname+"res"; hResDgaus->SetBinContent(b,gaus[name]);                                 hResDgaus->SetBinError(b,dgaus[name]);
		name = "mB"+xname+"res"; hResBgaus->SetBinContent(b,gaus[name]);                                 hResBgaus->SetBinError(b,dgaus[name]);
		name = "mD"+xname+"rat"; hRatDgaus->SetBinContent(b,gaus[name]);                                 hRatDgaus->SetBinError(b,dgaus[name]);
		name = "mB"+xname+"rat"; hRatBgaus->SetBinContent(b,gaus[name]);                                 hRatBgaus->SetBinError(b,dgaus[name]);
		name = "mD"+xname+"res"; hResDslopes->SetBinContent(b,slopes[name]);                             hResDslopes->SetBinError(b,dslopes[name]);
		name = "mB"+xname+"res"; hResBslopes->SetBinContent(b,slopes[name]);                             hResBslopes->SetBinError(b,dslopes[name]);
		name = "mD"+xname+"rat"; hRatDslopes->SetBinContent(b,slopes[name]);                             hRatDslopes->SetBinError(b,dslopes[name]);
		name = "mB"+xname+"rat"; hRatBslopes->SetBinContent(b,slopes[name]);                             hRatBslopes->SetBinError(b,dslopes[name]);
		name = "mD"+xname+"res"; hResDdistance->SetBinContent(b,distances[name]);                        hResDdistance->SetBinError(b,ddistances[name]);
		name = "mB"+xname+"res"; hResBdistance->SetBinContent(b,distances[name]);                        hResBdistance->SetBinError(b,ddistances[name]);
		name = "mD"+xname+"rat"; hRatDdistance->SetBinContent(b,distances[name]);                        hRatDdistance->SetBinError(b,ddistances[name]);
		name = "mB"+xname+"rat"; hRatBdistance->SetBinContent(b,distances[name]);                        hRatBdistance->SetBinError(b,ddistances[name]);
		name = "mD"+xname+"res"; hResDratio->SetBinContent(b, (isnan(ratios[name])) ? 0 : ratios[name]); hResDratio->SetBinError(b,dratios[name]);
		name = "mB"+xname+"res"; hResBratio->SetBinContent(b, (isnan(ratios[name])) ? 0 : ratios[name]); hResBratio->SetBinError(b,dratios[name]);
		name = "mD"+xname+"rat"; hRatDratio->SetBinContent(b, (isnan(ratios[name])) ? 0 : ratios[name]); hRatDratio->SetBinError(b,dratios[name]);
		name = "mB"+xname+"rat"; hRatBratio->SetBinContent(b, (isnan(ratios[name])) ? 0 : ratios[name]); hRatBratio->SetBinError(b,dratios[name]);
		
		
		cout << "bincenter["<<b<<"] = " << hResDslopes->GetBinCenter(b) << " --> xname=" << xname << endl;
		
		b++; 
	}
	
	
	
	
	
	Double_t xmin = hResDslopes->GetXaxis()->GetXmin();
	Double_t xmax = hResDslopes->GetXaxis()->GetXmax();
	
	TLine* line = new TLine(xmin,0.,xmax,0.);
	line->SetLineColor(kGreen); line->SetLineStyle(2); line->SetLineWidth(2);
	
	TString drawopt = "p0 e0"; // "hist"
	
	
	cnv = new TCanvas("","",1200,400); cnv->Draw(); cnv->cd(); cnv->Divide(2,1);
	cnv->cd(1); gPad->SetTicks(1,1);
	if(doBkg) { hResDslopes->SetMinimum( (hResDslopes->GetMinimum()<=hResBslopes->GetMinimum()) ? hResDslopes->GetMinimum()*2 : hResBslopes->GetMinimum()*2 ); } else { hResDslopes->SetMinimum(-0.002); }
	if(doBkg) { hResDslopes->SetMaximum( (hResDslopes->GetMaximum()>=hResBslopes->GetMaximum()) ? hResDslopes->GetMaximum()*2 : hResBslopes->GetMaximum()*2 ); } else { hResDslopes->SetMaximum(+0.002); }
	TF1* fResDslopes = new TF1("fResDslopes","[0]+[1]*x",xmin,xmax); TFitResultPtr fitptrDresSlopes = hResDslopes->Fit(fResDslopes,"EM0","",xmin,xmax); fResDslopes->SetLineColor(hResDslopes->GetLineColor()); fResDslopes->SetLineWidth(2);
	hResDslopes->SetLineColor(kBlack); hResDslopes->Draw(drawopt); line->Draw("same"); hResDslopes->Draw(drawopt+" same"); fResDslopes->Draw("same");
	if(doBkg) { hResBslopes->SetLineColor(kRed);   hResBslopes->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	gPad->Update();
	gPad->RedrawAxis();
	cnv->cd(2); gPad->SetTicks(1,1);
	line->Draw();
	if(doBkg) { hRatDslopes->SetMinimum( (hRatDslopes->GetMinimum()<=hRatBslopes->GetMinimum()) ? hRatDslopes->GetMinimum()*2 : hRatBslopes->GetMinimum()*2 ); } else { hRatDslopes->SetMinimum(-0.002); }
	if(doBkg) { hRatDslopes->SetMaximum( (hRatDslopes->GetMaximum()>=hRatBslopes->GetMaximum()) ? hRatDslopes->GetMaximum()*2 : hRatBslopes->GetMaximum()*2 ); } else { hRatDslopes->SetMaximum(+0.002); }
	TF1* fRatDslopes = new TF1("fRatDslopes","[0]+[1]*x",xmin,xmax); TFitResultPtr fitptrDratSlopes = hRatDslopes->Fit(fRatDslopes,"EM0","",xmin,xmax); fRatDslopes->SetLineColor(hRatDslopes->GetLineColor()); fRatDslopes->SetLineWidth(2);
	hRatDslopes->SetLineColor(kBlack); hRatDslopes->Draw(drawopt); line->Draw("same"); hRatDslopes->Draw(drawopt+" same"); fRatDslopes->Draw("same");
	if(doBkg) { hRatBslopes->SetLineColor(kRed);   hRatBslopes->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	gPad->Update();
	gPad->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.slopes.png");
	cnv->SaveAs("figures/sidebandsBias.slopes.pdf");
	cnv->SaveAs("figures/sidebandsBias.slopes.eps");
	cnv->SaveAs("figures/sidebandsBias.pdf");
	
	cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hResDslopes->Draw(drawopt); line->Draw("same"); hResDslopes->Draw(drawopt+" same"); fResDslopes->Draw("same");
	if(doBkg) { hResBslopes->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.residuals.Fit_slopes.eps");
	cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hRatDslopes->Draw(drawopt); line->Draw("same"); hRatDslopes->Draw(drawopt+" same"); fRatDslopes->Draw("same");
	if(doBkg) { hRatBslopes->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.ratios.Fit_slopes.eps");
	
	
	
	
	cnv = new TCanvas("","",1200,400); cnv->Draw(); cnv->cd(); cnv->Divide(2,1);
	cnv->cd(1); gPad->SetTicks(1,1);
	if(doBkg) { hResDgaus->SetMinimum( (hResDgaus->GetMinimum()<=hResBgaus->GetMinimum()) ? hResDgaus->GetMinimum()*2 : hResBgaus->GetMinimum()*2 ); } else { hResDgaus->SetMinimum(-2); }
	if(doBkg) { hResDgaus->SetMaximum( (hResDgaus->GetMaximum()>=hResBgaus->GetMaximum()) ? hResDgaus->GetMaximum()*2 : hResBgaus->GetMaximum()*2 ); } else { hResDgaus->SetMaximum(+2); }
	TF1* fResDgaus = new TF1("fResDgaus","[0]+[1]*x",xmin,xmax); TFitResultPtr fitptrDresGaus = hResDgaus->Fit(fResDgaus,"EM0","",xmin,xmax); fResDgaus->SetLineColor(hResDgaus->GetLineColor()); fResDgaus->SetLineWidth(2);
	hResDgaus->SetLineColor(kBlack);  hResDgaus->Draw(drawopt); line->Draw("smae"); hResDgaus->Draw(drawopt+" same"); fResDgaus->Draw("same");
	if(doBkg) { hResBgaus->SetLineColor(kRed);   hResBgaus->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	gPad->Update();
	gPad->RedrawAxis();
	cnv->cd(2); gPad->SetTicks(1,1);
	if(doBkg) { hRatDgaus->SetMinimum( (hRatDgaus->GetMinimum()<=hRatBgaus->GetMinimum()) ? hRatDgaus->GetMinimum()*2 : hRatBgaus->GetMinimum()*2 ); } else { hRatDgaus->SetMinimum(-2); }
	if(doBkg) { hRatDgaus->SetMaximum( (hRatDgaus->GetMaximum()>=hRatBgaus->GetMaximum()) ? hRatDgaus->GetMaximum()*2 : hRatBgaus->GetMaximum()*2 ); } else { hRatDgaus->SetMaximum(+2); }
	TF1* fRatDgaus = new TF1("fRatDgaus","[0]+[1]*x",xmin,xmax); TFitResultPtr fitptrDratGaus = hRatDgaus->Fit(fRatDgaus,"EM0","",xmin,xmax); fRatDgaus->SetLineColor(hRatDgaus->GetLineColor()); fRatDgaus->SetLineWidth(2);
	hRatDgaus->SetLineColor(kBlack); hRatDgaus->Draw(drawopt); line->Draw("same"); hRatDgaus->Draw(drawopt+" same"); fRatDgaus->Draw("same");
	if(doBkg) { hRatBgaus->SetLineColor(kRed);   hRatBgaus->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	gPad->Update();
	gPad->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.Gaus.png");
	cnv->SaveAs("figures/sidebandsBias.Gaus.pdf");
	cnv->SaveAs("figures/sidebandsBias.Gaus.eps");
	
	cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hResDgaus->Draw(drawopt); line->Draw("same"); hResDgaus->Draw(drawopt+" same"); fResDgaus->Draw("same");
	if(doBkg) { hResBgaus->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.residuals.Fit_gaus.eps");
	cnv = new TCanvas("","",600,400); cnv->Draw(); cnv->cd();
	hRatDgaus->Draw(drawopt); line->Draw("same"); hRatDgaus->Draw(drawopt+" same"); fRatDgaus->Draw("same");
	if(doBkg) { hRatBgaus->Draw(drawopt+" same"); }
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/sidebandsBias.ratios.Fit_gaus.eps");
	
	
	
	
	cnv = new TCanvas("","",1200,400); cnv->Draw(); cnv->cd();
	cnv->SaveAs("figures/sidebandsBias.pdf)");
}