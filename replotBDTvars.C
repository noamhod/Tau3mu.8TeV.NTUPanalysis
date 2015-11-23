//////////////////////////////////////////////////////////////////////
//// root -b -l -q replotBDTvars.C++\(1450.,1690.,1870.,2110.\) //////
//////////////////////////////////////////////////////////////////////
#include "std.h"
#include "const.h"
#include "type.h"
#include "postBDTcuts.h"
#include <iterator>
#include <cstdlib>

bool isBlinded = true;
bool resonVeto = true;

enum twoBody
{
	INCLUSIVE,
	ONPHI1,   OFFPHI1,
	ONPHI2,   OFFPHI2,
	ONOMEGA1, OFFOMEGA1,
	ONOMEGA2, OFFOMEGA2,
	ONOMEGAPHI1, OFFOMEGAPHI1,
	ONOMEGAPHI2, OFFOMEGAPHI2,
	ONLOWOS1, OFFLOWOS1,
	ONLOWOS2, OFFLOWOS2,
	ONLOWSS, OFFLOWSS,
	MAXRANGE
};

TString getRangeOS(int k)
{
	TString  rangeOS = "";
	switch(k)
	{
		case INCLUSIVE: rangeOS = "";                break;
		case ONPHI1:    rangeOS = "_OS1onPhi";       break;
		case OFFPHI1:   rangeOS = "_OS1offPhi";      break;
		case ONPHI2:    rangeOS = "_OS2onPhi";       break;
		case OFFPHI2:   rangeOS = "_OS2offPhi";      break;
		case ONOMEGA1:  rangeOS = "_OS1onOmega";     break;
		case OFFOMEGA1: rangeOS = "_OS1offOmega";    break;
		case ONOMEGA2:  rangeOS = "_OS2onOmega";     break;
		case OFFOMEGA2: rangeOS = "_OS2offOmega";    break;
		case ONOMEGAPHI1:  rangeOS = "_OS1onOmegaPhi";  break;
		case OFFOMEGAPHI1: rangeOS = "_OS1offOmegaPhi"; break;
		case ONOMEGAPHI2:  rangeOS = "_OS2onOmegaPhi";  break;
		case OFFOMEGAPHI2: rangeOS = "_OS2offOmegaPhi"; break;
		case ONLOWOS1:     rangeOS = "_OS1onLow";       break;
		case OFFLOWOS1:    rangeOS = "_OS1OffLow";      break;
		case ONLOWOS2:     rangeOS = "_OS2onLow";       break;
		case OFFLOWOS2:    rangeOS = "_OS2OffLow";      break;
		case ONLOWSS:      rangeOS = "_SSonLow";        break;
		case OFFLOWSS:     rangeOS = "_SSOffLow";       break;
		default: _FAT("k="<<k<<" is not supported"); break;
	}
	return rangeOS;
}

void getRangeOSdecisions(int k, vector<bool>& rangeOSflags, float mOS1, float mOS2, float mSS)
{
	rangeOSflags.push_back(0);
	
	bool onPhi1    = (fabs(mOS1-1020)<30);
	bool offPhi1   = (fabs(mOS1-1020)>=30 && fabs(mOS1-1020)<120);
	bool onPhi2    = (fabs(mOS2-1020)<30);
	bool offPhi2   = (fabs(mOS2-1020)>=30 && fabs(mOS2-1020)<120);
	bool onOmega1  = (fabs(mOS1-782)<30);
	bool offOmega1 = (fabs(mOS1-782)>=30 && fabs(mOS1-782)<120);
	bool onOmega2  = (fabs(mOS2-782)<30);
	bool offOmega2 = (fabs(mOS2-782)>=30 && fabs(mOS2-782)<120);
	bool onLowOS1  = (mOS1<500);
	bool offLowOS1 = (mOS1>=500 && mOS1<720);
	bool onLowOS2  = (mOS2<500);
	bool offLowOS2 = (mOS2>=500 && mOS2<720);
	bool onLowSS   = (mSS<500);
	bool offLowSS  = (mSS>=500 && mSS<720);

	switch(k)
	{
		case INCLUSIVE:     rangeOSflags[INCLUSIVE]    = 1;                                  break;
		case ONPHI1:        rangeOSflags[ONPHI1   ]    = (onPhi1 && !onPhi2 && onOmega2);    break;
		case OFFPHI1:       rangeOSflags[OFFPHI1  ]    = (offPhi1 && !onPhi2 && onOmega2);   break;
		case ONPHI2:        rangeOSflags[ONPHI2   ]    = (onPhi2 && !onPhi1 && onOmega1);    break;
		case OFFPHI2:       rangeOSflags[OFFPHI2  ]    = (offPhi2 && !onPhi1 && onOmega1);   break;
		case ONOMEGA1:      rangeOSflags[ONOMEGA1 ]    = (onOmega1 && !onPhi2 && onOmega2);  break;
		case OFFOMEGA1:     rangeOSflags[OFFOMEGA1]    = (offOmega1 && !onPhi2 && onOmega2); break;
		case ONOMEGA2:      rangeOSflags[ONOMEGA2 ]    = (onOmega2 && !onPhi1 && onOmega1);  break;
		case OFFOMEGA2:     rangeOSflags[OFFOMEGA2]    = (offOmega2 && !onPhi1 && onOmega1); break;
		case ONOMEGAPHI1:   rangeOSflags[ONOMEGAPHI1 ] = ((onPhi1 || onOmega1) && !(onPhi2 || onOmega2));  break;
		case OFFOMEGAPHI1:  rangeOSflags[OFFOMEGAPHI1] = ((offPhi1 || offOmega1) && !onPhi2 && !onOmega2); break;
		case ONOMEGAPHI2:   rangeOSflags[ONOMEGAPHI2 ] = ((onPhi2 || onOmega2) && !(onPhi1 || onOmega1));  break;
		case OFFOMEGAPHI2:  rangeOSflags[OFFOMEGAPHI2] = ((offPhi2 || offOmega2) && !onPhi1 && !onOmega1); break;
		case ONLOWOS1:      rangeOSflags[ONLOWOS1]     = onLowOS1;  break;
		case OFFLOWOS1:     rangeOSflags[OFFLOWOS1]    = offLowOS1; break;
		case ONLOWOS2:      rangeOSflags[ONLOWOS2]     = onLowOS2;  break;
		case OFFLOWOS2:     rangeOSflags[OFFLOWOS2]    = offLowOS2; break;
		case ONLOWSS:       rangeOSflags[ONLOWSS]      = onLowSS;   break;
		case OFFLOWSS:      rangeOSflags[OFFLOWSS]     = offLowSS;  break;
		default: _FAT("k="<<k<<" is not supported"); break;
	}
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

void rgbPalette(Double_t r, Double_t g, Double_t b, Int_t nb=50)
{
	const UInt_t Number = 3;
	Double_t Red[Number]    = { r,   0.0, 0.0 };
	Double_t Green[Number]  = { g,   0.0, 0.0 };
	Double_t Blue[Number]   = { b,   0.0, 0.0 };
	Double_t Length[Number] = { 0.1, 0.5, 1.0 };
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}
void gryPalette(Int_t nb=50)
{
	const UInt_t Number = 2;
	Double_t Red[Number]    = { 1.0, 0.0 };
	Double_t Green[Number]  = { 1.0, 0.0 };
	Double_t Blue[Number]   = { 1.0, 0.0 };
	Double_t Length[Number] = { 0.0, 1.0,};
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}
TExec* exeRed   = new TExec("exeRed",   "rgbPalette(0.7,0,0,50);");
TExec* exeGreen = new TExec("exeGreen", "rgbPalette(0,0.5,0,50);");
// TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0.2,0.9,50);");
TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0,0.8,10);");
TExec* exeGray  = new TExec("exeGray",  "rgbPalette(1,1,1,50);");

TPaveText* ptxt;
TPaveText* ptxtLumi;
void makeAtlasLabel()
{	
	// ptxt = new TPaveText(0.65,0.75,0.90,0.95,"NDC");
	// ptxt->SetFillStyle(4000); //will be transparent
	// ptxt->SetFillColor(0);
	// ptxt->SetTextFont(42);
	// ptxt->SetBorderSize(0);
	// ptxt->AddText("#bf{#it{ATLAS}} internal");
	// ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	// ptxt->AddText("#sqrt{s}=8 TeV");

	ptxt = new TPaveText(0.65,0.75,0.90,0.95,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} Internal");
	ptxt->AddText("#sqrt{s}=8 TeV, "+slumi);
}

void plotAtlasLabel()
{
	Double_t x = 0.60;
	Double_t y = 0.88;
	
	TLatex t1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
	t1.SetNDC();
	t1.SetTextFont(72);
	t1.SetTextColor(kBlack);
	double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
	// double dely = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
	t1.DrawLatex(x,y,"ATLAS");
	
	TLatex t2; 
	t2.SetNDC();
	t2.SetTextFont(42);
	t2.SetTextColor(kBlack);
	t2.DrawLatex(x+delx,y,"Internal");
	
	TLatex t3; 
	t3.SetNDC();
	t3.SetTextFont(42);
	t3.SetTextColor(kBlack);
	t3.DrawLatex(x,y-0.7*delx,"#sqrt{s}=8 TeV, "+slumi);
}

TLegend* legR;
TLegend* legL;
TLegend* legTL;
TLegend* legTM;
void makeLegend()
{	
	legR = new TLegend(0.58,0.57,0.88,0.72,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
	
	legTM = new TLegend(0.25,0.78,0.55,0.93,NULL,"brNDC");
	legTM->SetFillStyle(4000); //will be transparent
	legTM->SetFillColor(0);
	legTM->SetTextFont(42);
	legTM->SetBorderSize(0);
	
	legTL = new TLegend(0.15,0.78,0.45,0.93,NULL,"brNDC");
	legTL->SetFillStyle(4000); //will be transparent
	legTL->SetFillColor(0);
	legTL->SetTextFont(42);
	legTL->SetBorderSize(0);
	
	legL  = new TLegend(0.15,0.78,0.45,0.93,NULL,"brNDC");
	legL->SetFillStyle(4000); //will be transparent
	legL->SetFillColor(0);
	legL->SetTextFont(42);
	legL->SetBorderSize(0);
}

TPaletteAxis* palette;
void fixPalette(TH2* h)
{
	palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.855); palette->SetX2NDC(0.905); palette->SetY1NDC(0.1); palette->SetY2NDC(0.9);
}



void addHist(TMapTSP2TH1& histos, TMapTSP2TH2& histos2, TMapTSP2TProfile& profiles, TString channel, TString rangeOS, TString name, TString titles, int nbins, float xmin, float xmax)
{
	if(nbins<=0) _FAT("nbins="<<nbins<<" for histogram: "<<name);
	TString hname = channel+"_"+name+rangeOS;
	histos.insert(make_pair(hname, new TH1F(hname,titles,nbins,xmin,xmax)));
	histos.insert(make_pair(hname+"_frame", new TH1F(hname+"_frame",titles,nbins,xmin,xmax)));
	histos[hname]->Sumw2();
	if(channel=="Data")
	{
		histos[hname]->SetMarkerStyle(24);
		histos[hname]->SetMarkerSize(1.2);
		histos[hname]->SetLineWidth(1);
		histos[hname]->SetMarkerColor(kBlack);
		histos[hname]->SetLineColor(kBlack);
	}
	else
	{
		histos[hname]->SetFillColor(kGray);
		histos[hname]->SetLineColor(kGray);
	}
	
	TString titlesBDT2d = titles; titlesBDT2d.ReplaceAll(";Events",";BDT score;Events");
	TString hnameBDT2d = "bdt_vs_"+hname;
	histos2.insert(make_pair(hnameBDT2d, new TH2F(hnameBDT2d,titlesBDT2d,nbins,xmin,xmax,nBDTbins,minBDTcut,+1.0)));
	histos2[hnameBDT2d]->Sumw2();
	
	TString titlesProfBDT = titles; titlesProfBDT.ReplaceAll(";Events",";BDT score;Events");
	TString pnameBDT = "prof_bdt_"+hname;
	profiles.insert(make_pair(pnameBDT, new TProfile(pnameBDT,titlesProfBDT,(int)((float)nbins/2.),xmin,xmax,minBDTcut,+1.0)));
	profiles[pnameBDT]->SetLineColor(kBlack); profiles[pnameBDT]->SetMarkerColor(kBlack); profiles[pnameBDT]->SetMarkerStyle(20);
	
	TString titlesMass2d = titles; titlesMass2d.ReplaceAll(";Events",";#it{m}_{3#mu} [MeV];Events");
	TString hnameMass2d = "m3body_vs_"+hname;
	histos2.insert(make_pair(hnameMass2d, new TH2F(hnameMass2d,titlesMass2d,nbins,xmin,xmax,22,1450,2110)));
	histos2[hnameMass2d]->Sumw2();
	
	TString titlesProfMass = titles; titlesProfMass.ReplaceAll(";Events",";#it{m}_{3#mu} [MeV];Events");
	TString pnameMass = "prof_m3body_"+hname;
	profiles.insert(make_pair(pnameMass, new TProfile(pnameMass,titlesProfMass,(int)((float)nbins/2.),xmin,xmax,1450,2110)));
	profiles[pnameMass]->SetLineColor(kBlack); profiles[pnameMass]->SetMarkerColor(kBlack); profiles[pnameMass]->SetMarkerStyle(20);
}

void addHist(TMapTSP2TH2& histos, TString channel, TString name, TString titles, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax)
{
	TString hname = channel+"_"+name;
	histos.insert(make_pair(hname, new TH2F(hname,titles,nxbins,xmin,xmax,nybins,ymin,ymax)));
	histos.insert(make_pair(hname+"_frame", new TH2F(hname+"_frame",titles,nxbins,xmin,xmax,nybins,ymin,ymax)));
	histos[hname]->Sumw2();
	if(channel=="Data")
	{
		// histos[hname]->SetMarkerStyle(20);
		// histos[hname]->SetMarkerSize(1.2);
		// histos[hname]->SetMarkerColor(kBlack);
		// histos[hname]->SetLineColor(kBlack);
	}
	else
	{
		// histos[hname]->SetFillColor(kGray);
		// histos[hname]->SetLineColor(kGray);
	}
}

void setMax(TH1* h1, TH1* h2, float scale=0.9)
{
	Int_t maxbin1 = h1->GetMaximumBin();
	Int_t maxbin2 = h2->GetMaximumBin();
	float max1 = h1->GetBinContent(maxbin1)+h1->GetBinError(maxbin1);
	float max2 = h2->GetBinContent(maxbin2)+h2->GetBinError(maxbin2);
	float max = (max1>max2) ? max1 : max2;
	h1->Scale(scale/max);
	h2->Scale(scale/max);
}
void setMax3(TH1* h1, TH1* h2, TH1* h3, float scale=0.9)
{
	Int_t maxbin1 = h1->GetMaximumBin();
	Int_t maxbin2 = h2->GetMaximumBin();
	float max1 = h1->GetBinContent(maxbin1)+h1->GetBinError(maxbin1);
	float max2 = h2->GetBinContent(maxbin2)+h2->GetBinError(maxbin2);
	float max = (max1>max2) ? max1 : max2;
	// float max = (h1->GetMaximum()>h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();
	h1->Scale(scale/max);
	h2->Scale(scale/max);
	h3->Scale(scale/max);
}
void setMax4(TH1* h1, TH1* h2, TH1* h3, TH1* h4)
{
	Int_t maxbin1 = h1->GetMaximumBin();
	Int_t maxbin2 = h2->GetMaximumBin();
	float max1 = h1->GetBinContent(maxbin1)+h1->GetBinError(maxbin1);
	float max2 = h2->GetBinContent(maxbin2)+h2->GetBinError(maxbin2);
	float max = (max1>max2) ? max1 : max2;
	
	float scale = max2;
	
	h1->Scale(scale/max);    // Signal loose
	// h2->Scale(scale/max); // Data loose
	h3->Scale(scale/max);    // Signal tight
	// h4->Scale(scale/max); // Data tight
}

void setMax(TH2* h, float scale=0.9)
{
	float max = h->GetMaximum();
	h->Scale(scale/max);
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
float Sum(TH2* h, bool addUunderFlow=false, bool addOverFlow=false)
{
	float I=0.;
	for(int i=1 ; i<=h->GetNbinsX() ; i++)// without under- and over-flow
	{
		for(int j=1 ; j<=h->GetNbinsY() ; j++)
		{
			I += h->GetBinContent(i,j);
		}
	}
	if(addUunderFlow) I+=h->GetBinContent(0);
	if(addOverFlow)   I+=h->GetBinContent(h->GetNbinsX()*h->GetNbinsY()+1);
	return I;
}
void NormToEntries(TH1* h, Double_t factor=-1)
{
	Double_t entries = (factor<0) ? Sum(h) : factor;
	Double_t scale = (entries==0.) ? 1. : 1./entries;
	h->Scale(scale);
	h->GetYaxis()->SetTitle("Events");
}
void NormToEntries(TH2* h)
{
	Double_t entries = Sum(h);
	Double_t scale = (entries==0.) ? 1. : 1./entries;
	h->Scale(scale);
	h->GetZaxis()->SetTitle("Events");
}


void plot(TString prefix, TString var, TMapTSP2TH1& histos1, TMapTSP2TH2& histos2, TMapTSP2TProfile& profiles, TMapTSf& rescales, TLegend* leg, TString rangeOS="")
{
	TString varname = var+rangeOS;
	TCanvas* cnv = new TCanvas("cnv","",800,600);
	gPad->SetTicks(1,1);
	
	float scaleLoose = Sum(histos1["Data_"+varname])/Sum(histos1["Wtaunu_3mu_"+varname]);
	float scaleTight = scaleLoose*rescales["Wtaunu_3mu_tight_"+varname];
	histos1["Wtaunu_3mu_"+varname]->Scale(scaleLoose);
	histos1["Wtaunu_3mu_tight_"+varname]->Scale(scaleTight);
	
	histos1["Data_"+varname]->SetBinErrorOption(TH1::kPoisson);
	histos1["Data_tight_"+varname]->SetBinErrorOption(TH1::kPoisson);
	
	Int_t maxbinSig = histos1["Wtaunu_3mu_"+varname]->GetMaximumBin();
	Int_t maxbinDat = histos1["Data_"+varname]->GetMaximumBin();
	float maxSig = histos1["Wtaunu_3mu_"+varname]->GetBinContent(maxbinSig)+histos1["Wtaunu_3mu_"+varname]->GetBinError(maxbinSig);
	float maxDat = histos1["Data_"+varname]->GetBinContent(maxbinDat)+histos1["Data_"+varname]->GetBinError(maxbinDat);
	float max = (maxSig>maxDat) ? maxSig : maxDat;
	histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMaximum(max*1.1);
	
	histos1["Wtaunu_3mu_"+varname+"_frame"]->Draw();
	histos1["Wtaunu_3mu_"+varname]->Draw("hist same");
	histos1["Wtaunu_3mu_tight_"+varname]->Draw("hist same");
	histos1["Data_"+varname]->Draw("p1x1 same");
	histos1["Data_tight_"+varname]->Draw("p1x1 same");	
	
	plotAtlasLabel();
	leg->Clear();
	leg->AddEntry(histos1["Data_"+varname],"SB data (loose)","ple");
	leg->AddEntry(histos1["Data_tight_"+varname],"SB data (tight+#it{x}>#it{x}_{0})","ple");
	leg->AddEntry(histos1["Wtaunu_3mu_"+varname],"Signal (loose)","f");
	leg->AddEntry(histos1["Wtaunu_3mu_tight_"+varname],"Signal (tight+#it{x}>#it{x}_{0})","f");
	leg->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	if(rangeOS!="") cnv->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".pdf");
	else
	{
		cnv->SaveAs("figures/"+prefix+".pdf");
		cnv->SaveAs("figures/"+prefix+"."+varname+".png");
		cnv->SaveAs("figures/"+prefix+"."+varname+".eps");
		cnv->SaveAs("figures/"+prefix+"."+varname+".pdf");
	}
	
	///////////////////// Plot in log Y axis
	histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMaximum(max*20);

	histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMinimum(5e-1);
	histos1["Wtaunu_3mu_"+varname]->SetMinimum(5e-1);
	histos1["Wtaunu_3mu_tight_"+varname]->SetMinimum(5e-1);
	histos1["Data_"+varname]->SetMinimum(5e-1);
	histos1["Data_tight_"+varname]->SetMinimum(5e-1);
	
	TCanvas* cnvlog = new TCanvas("cnvlog","",800,600);
	gPad->SetTicks(1,1);
	gPad->SetLogy();
	histos1["Wtaunu_3mu_"+varname+"_frame"]->Draw();
	histos1["Wtaunu_3mu_"+varname]->Draw("hist same");
	histos1["Wtaunu_3mu_tight_"+varname]->Draw("hist same");
	histos1["Data_"+varname]->Draw("p1x1 same");
	histos1["Data_tight_"+varname]->Draw("p1x1 same");
	plotAtlasLabel();
	leg->Draw("same");
	cnvlog->Update();
	cnvlog->RedrawAxis();
	if(rangeOS!="") cnvlog->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".logy.pdf");
	else
	{
		cnvlog->SaveAs("figures/"+prefix+".pdf");
		cnvlog->SaveAs("figures/"+prefix+"."+varname+".logy.png");
		cnvlog->SaveAs("figures/"+prefix+"."+varname+".logy.eps");
		cnvlog->SaveAs("figures/"+prefix+"."+varname+".logy.pdf");
	}
	
	//// set it back to the default
	histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMaximum(max*1.1);

	histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMinimum(0);
	histos1["Wtaunu_3mu_"+varname]->SetMinimum(0);
	histos1["Wtaunu_3mu_tight_"+varname]->SetMinimum(0);
	histos1["Data_"+varname]->SetMinimum(0);
	histos1["Data_tight_"+varname]->SetMinimum(0);
	/////////////////////
	
	TCanvas* cnv2 = new TCanvas("cnv2","",1200,600);
	cnv2->Divide(2,1);
	cnv2->cd(1); histos2["bdt_vs_Data_"+varname]->Draw("col");       profiles["prof_bdt_Data_"+varname]->Draw("same");
	cnv2->cd(2); histos2["bdt_vs_Wtaunu_3mu_"+varname]->Draw("col"); profiles["prof_bdt_Wtaunu_3mu_"+varname]->Draw("same");
	if(rangeOS!="") cnv2->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".pdf");
	else
	{
		cnv2->SaveAs("figures/"+prefix+".pdf");
		cnv2->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_loose.pdf");
		cnv2->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_loose.eps");
		cnv2->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_loose.png");
	}
	
	TCanvas* cnv3 = new TCanvas("cnv3","",1200,600);
	cnv3->Divide(2,1);
	cnv3->cd(1); histos2["bdt_vs_Data_tight_"+varname]->Draw("col");       profiles["prof_bdt_Data_tight_"+varname]->Draw("same");
	cnv3->cd(2); histos2["bdt_vs_Wtaunu_3mu_tight_"+varname]->Draw("col"); profiles["prof_bdt_Wtaunu_3mu_tight_"+varname]->Draw("same");
	if(rangeOS!="") cnv3->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".pdf");
	else
	{
		cnv3->SaveAs("figures/"+prefix+".pdf");
		cnv3->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_tight.pdf");
		cnv3->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_tight.eps");
		cnv3->SaveAs("figures/"+prefix+".bdt_vs_"+varname+"_tight.png");
	}
	
	TCanvas* cnv4 = new TCanvas("cnv4","",1200,600);
	cnv4->Divide(2,1);
	cnv4->cd(1); histos2["m3body_vs_Data_"+varname]->Draw("col");       profiles["prof_m3body_Data_"+varname]->Draw("same");
	cnv4->cd(2); histos2["m3body_vs_Wtaunu_3mu_"+varname]->Draw("col"); profiles["prof_m3body_Wtaunu_3mu_"+varname]->Draw("same");
	if(rangeOS!="") cnv4->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".pdf");
	else
	{
		cnv4->SaveAs("figures/"+prefix+".pdf");
		cnv4->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_loose.pdf");
		cnv4->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_loose.eps");
		cnv4->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_loose.png");
	}
	
	TCanvas* cnv5 = new TCanvas("cnv5","",1200,600);
	cnv5->Divide(2,1);
	cnv5->cd(1); histos2["m3body_vs_Data_tight_"+varname]->Draw("col");       profiles["prof_m3body_Data_tight_"+varname]->Draw("same");
	cnv5->cd(2); histos2["m3body_vs_Wtaunu_3mu_tight_"+varname]->Draw("col"); profiles["prof_m3body_Wtaunu_3mu_tight_"+varname]->Draw("same");
	if(rangeOS!="") cnv5->SaveAs("figures/"+prefix+"_resonances"+rangeOS+".pdf");
	else
	{
		cnv5->SaveAs("figures/"+prefix+".pdf");
		cnv5->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_tight.pdf");
		cnv5->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_tight.eps");
		cnv5->SaveAs("figures/"+prefix+".m3body_vs_"+varname+"_tight.png");
	}
	
	delete cnv;
	delete cnv2;
	delete cnv3;
	delete cnv4;
	delete cnv5;
}


void plot9(TString prefix, TString name, vector<TString>& vars, vector<TLegend*>& legs, TMapTSP2TH1& histos1, bool doLogy=false)
{
	TCanvas* cnv = new TCanvas("cnv","",1200,900);
	cnv->Divide(3,3);
	
	for(unsigned int pad=0 ; pad<vars.size() ; ++pad)
	{
		TString varname = vars[pad];
		cnv->cd(pad+1);
		gPad->SetTicks(1,1);
		if(doLogy) gPad->SetLogy();
	
		Int_t maxbinSig = histos1["Wtaunu_3mu_"+varname]->GetMaximumBin();
		Int_t maxbinDat = histos1["Data_"+varname]->GetMaximumBin();
		float maxSig = histos1["Wtaunu_3mu_"+varname]->GetBinContent(maxbinSig)+histos1["Wtaunu_3mu_"+varname]->GetBinError(maxbinSig);
		float maxDat = histos1["Data_"+varname]->GetBinContent(maxbinDat)+histos1["Data_"+varname]->GetBinError(maxbinDat);
		float max = (maxSig>maxDat) ? maxSig : maxDat;
		histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMaximum(max*20);
	
	
		histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMinimum(5e-1);
		histos1["Wtaunu_3mu_"+varname]->SetMinimum(5e-1);
		histos1["Wtaunu_3mu_tight_"+varname]->SetMinimum(5e-1);
		histos1["Data_"+varname]->SetMinimum(5e-1);
		histos1["Data_tight_"+varname]->SetMinimum(5e-1);
	
		histos1["Wtaunu_3mu_"+varname+"_frame"]->Draw();
		histos1["Wtaunu_3mu_"+varname]->Draw("hist same");
		histos1["Wtaunu_3mu_tight_"+varname]->Draw("hist same");
		histos1["Data_"+varname]->Draw("p1x1 same");
		histos1["Data_tight_"+varname]->Draw("p1x1 same");	
		
		plotAtlasLabel();
		
		legs[pad]->Clear();
		legs[pad]->AddEntry(histos1["Data_"+varname],"SB data (loose)","ple");
		legs[pad]->AddEntry(histos1["Data_tight_"+varname],"SB data (tight+#it{x}>#it{x}_{0})","ple");
		legs[pad]->AddEntry(histos1["Wtaunu_3mu_"+varname],"Signal (loose)","f");
		legs[pad]->AddEntry(histos1["Wtaunu_3mu_tight_"+varname],"Signal (tight+#it{x}>#it{x}_{0})","f");
		legs[pad]->Draw("same");
		
		gPad->Update();
		gPad->RedrawAxis();
	}
	
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/"+prefix+"."+name+".png");
	cnv->SaveAs("figures/"+prefix+"."+name+".eps");
	cnv->SaveAs("figures/"+prefix+"."+name+".pdf");
	cnv->SaveAs("figures/"+prefix+".pdf");
	delete cnv;
}

void plot2(TString prefix, TString name, vector<TString>& vars, vector<TLegend*>& legs, TMapTSP2TH1& histos1, bool doLogy=false)
{
	TCanvas* cnv = new TCanvas("cnv","",800,600);
	cnv->Divide(2,1);
	
	for(unsigned int pad=0 ; pad<vars.size() ; ++pad)
	{
		TString varname = vars[pad];
		cnv->cd(pad+1);
		gPad->SetTicks(1,1);
		if(doLogy) gPad->SetLogy();
	
		Int_t maxbinSig = histos1["Wtaunu_3mu_"+varname]->GetMaximumBin();
		Int_t maxbinDat = histos1["Data_"+varname]->GetMaximumBin();
		float maxSig = histos1["Wtaunu_3mu_"+varname]->GetBinContent(maxbinSig)+histos1["Wtaunu_3mu_"+varname]->GetBinError(maxbinSig);
		float maxDat = histos1["Data_"+varname]->GetBinContent(maxbinDat)+histos1["Data_"+varname]->GetBinError(maxbinDat);
		float max = (maxSig>maxDat) ? maxSig : maxDat;
		histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMaximum(max*20);
	
	
		histos1["Wtaunu_3mu_"+varname+"_frame"]->SetMinimum(5e-1);
		histos1["Wtaunu_3mu_"+varname]->SetMinimum(5e-1);
		histos1["Wtaunu_3mu_tight_"+varname]->SetMinimum(5e-1);
		histos1["Data_"+varname]->SetMinimum(5e-1);
		histos1["Data_tight_"+varname]->SetMinimum(5e-1);
	
		histos1["Wtaunu_3mu_"+varname+"_frame"]->Draw();
		histos1["Wtaunu_3mu_"+varname]->Draw("hist same");
		histos1["Wtaunu_3mu_tight_"+varname]->Draw("hist same");
		histos1["Data_"+varname]->Draw("p1x1 same");
		histos1["Data_tight_"+varname]->Draw("p1x1 same");	
		
		plotAtlasLabel();
		
		legs[pad]->Clear();
		legs[pad]->AddEntry(histos1["Data_"+varname],"SB data (loose)","ple");
		legs[pad]->AddEntry(histos1["Data_tight_"+varname],"SB data (tight+#it{x}>#it{x}_{0})","ple");
		legs[pad]->AddEntry(histos1["Wtaunu_3mu_"+varname],"Signal (loose)","f");
		legs[pad]->AddEntry(histos1["Wtaunu_3mu_tight_"+varname],"Signal (tight+#it{x}>#it{x}_{0})","f");
		legs[pad]->Draw("same");
		
		gPad->Update();
		gPad->RedrawAxis();
	}
	
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/"+prefix+"."+name+".png");
	cnv->SaveAs("figures/"+prefix+"."+name+".eps");
	cnv->SaveAs("figures/"+prefix+"."+name+".pdf");
	cnv->SaveAs("figures/"+prefix+".pdf");
	delete cnv;
}


int readData(TTree* t, TMapTSP2TH1& histos1, TMapTSP2TH2& histos2, TMapTSP2TProfile& profiles, TString channel, float xFullMin,float xBlindMin,float xBlindMax,float xFullMax, bool isTight=false)
{
	vector<int>*   loose           = 0;
	vector<int>*   tight           = 0;
	vector<int>*   pvntrk          = 0;
	vector<int>*   njets           = 0;
	vector<int>*   muonauthor      = 0;
	vector<float>* drmax           = 0;
	vector<float>* mass            = 0;
	vector<float>* pt              = 0;
	vector<float>* mSS             = 0;
	vector<float>* mOS1            = 0;
	vector<float>* mOS2            = 0;
	vector<float>* trkprob         = 0;
	vector<float>* iso003          = 0;
	vector<float>* iso010          = 0;
	vector<float>* iso014          = 0;
	vector<float>* iso020          = 0;
	vector<float>* iso030          = 0;
	vector<float>* maxpbalsig      = 0;
	vector<float>* pvalue          = 0;
	vector<float>* lxy             = 0;
	vector<float>* dlxy            = 0;
	vector<float>* a0xy            = 0;
	vector<float>* da0xy           = 0;
	float          calo_met        = 0;
	float          trk_met         = 0;
	vector<float>* calo_mt         = 0;
	vector<float>* trk_mt          = 0;
	vector<float>* calo_dphimet    = 0;
	vector<float>* trk_dphimet     = 0;
	vector<float>* dphimets        = 0;
	vector<float>* ht              = 0;
	vector<float>* calo_mht        = 0;
	vector<float>* trk_mht         = 0;
	vector<float>* calo_dphimet_ht = 0;
	vector<float>* trk_dphimet_ht  = 0;
	vector<float>* dr_ht           = 0;
	vector<float>* score           = 0;
	vector<float>* lxySig          = 0;
	vector<float>* a0xySig         = 0;
	vector<float>* dptreltrk       = 0;
	vector<float>* dptrelcal       = 0;
	
	t->SetBranchAddress("pass_loose",       &loose);
	t->SetBranchAddress("pass_tight",       &tight);
	t->SetBranchAddress("vtx_pvNtrk",       &pvntrk);
	t->SetBranchAddress("vtx_dRmax",        &drmax);
	t->SetBranchAddress("vtx_pt",           &pt);
	t->SetBranchAddress("vtx_mass",         &mass);
	t->SetBranchAddress("vtx_mSS",          &mSS);
	t->SetBranchAddress("vtx_mOS1",         &mOS1);
	t->SetBranchAddress("vtx_mOS2",         &mOS2);
	t->SetBranchAddress("trks_fitprob",     &trkprob);
	t->SetBranchAddress("vtx_isolation003", &iso003);
	t->SetBranchAddress("vtx_isolation010", &iso010);
	t->SetBranchAddress("vtx_isolation014", &iso014);
	t->SetBranchAddress("vtx_isolation020", &iso020);
	t->SetBranchAddress("vtx_isolation030", &iso030);
	t->SetBranchAddress("muons_maxpbalsig", &maxpbalsig);
	t->SetBranchAddress("vtx_pval",         &pvalue);
	t->SetBranchAddress("vtx_lxy",          &lxy);
	t->SetBranchAddress("vtx_lxyErr",       &dlxy);
	t->SetBranchAddress("geo_lxySig",       &lxySig);
	t->SetBranchAddress("vtx_a0xy",         &a0xy);
	t->SetBranchAddress("vtx_a0xyErr",      &da0xy);
	t->SetBranchAddress("geo_a0xySig",      &a0xySig);
	t->SetBranchAddress("met_muons_et",     &calo_met);
	t->SetBranchAddress("met_track_et",     &trk_met);
	t->SetBranchAddress("met_muons_mT",     &calo_mt);
	t->SetBranchAddress("met_track_mT",     &trk_mt);
	t->SetBranchAddress("met_muons_dPhi3mu",&calo_dphimet);
	t->SetBranchAddress("met_track_dPhi3mu",&trk_dphimet);
	t->SetBranchAddress("mets_dphi",        &dphimets);
	t->SetBranchAddress("mets_dptreltrk",   &dptreltrk);
	t->SetBranchAddress("mets_dptrelcal",   &dptrelcal);
	t->SetBranchAddress("ht_pt",            &ht);
	t->SetBranchAddress("ht_mT",            &calo_mht);
	t->SetBranchAddress("ht_mT_mettrk",     &trk_mht);
	t->SetBranchAddress("ht_dphimet_muons", &calo_dphimet_ht);
	t->SetBranchAddress("ht_dphimet_track", &trk_dphimet_ht);
	t->SetBranchAddress("ht_dr3body",       &dr_ht);
	t->SetBranchAddress("jets_n",           &njets);
	t->SetBranchAddress("vtx_code",         &muonauthor);
	t->SetBranchAddress("mva_score",        &score);

	vector<bool> rangeOSflags;
	int ncandidates = 0;
	for(Int_t entry=1 ; entry<=t->GetEntries() ; entry++)
	{
		t->GetEntry(entry);
		
		if(mass->size()!=1) continue;
		
		bool isFilled = false;
		
		TString type = (channel.Contains("_3mu")) ? "sig" : "bkg";
		TString mode = (isTight) ? "postTraining" : "preTraining";
		bool doBDT   = (isTight);
		vector<float>* metcal = new vector<float>; metcal->push_back(calo_met);
		vector<float>* mettrk = new vector<float>; mettrk->push_back(trk_met);
		TMapTSP2vi vi;
		TMapTSP2vf vf;
		vi.insert(make_pair("pass_loose", loose));
		vi.insert(make_pair("njets"     , njets));
		vf.insert(make_pair("pt3body"   , pt));
		vf.insert(make_pair("m3body"    , mass));
		vf.insert(make_pair("mSS"       , mSS));
		vf.insert(make_pair("mOS1"      , mOS1));
		vf.insert(make_pair("mOS2"      , mOS2));
		vf.insert(make_pair("score"     , score));
		vf.insert(make_pair("pval"      , pvalue));
		vf.insert(make_pair("iso003"    , iso003));
		vf.insert(make_pair("iso020"    , iso020));
		vf.insert(make_pair("iso030"    , iso030));
		vf.insert(make_pair("trkspval"  , trkprob));
		vf.insert(make_pair("lxysig"    , lxySig));
		vf.insert(make_pair("a0xysig"   , a0xySig));
		vf.insert(make_pair("mettrk"    , mettrk));
		vf.insert(make_pair("metcal"    , metcal));
		vf.insert(make_pair("mttrk"     , trk_mt));
		vf.insert(make_pair("mtcal"     , calo_mt));
		vf.insert(make_pair("dphihttrk" , trk_dphimet_ht));
		vf.insert(make_pair("dphihtcal" , calo_dphimet_ht));
		bool pass = passPostBDTcut(0,vf,vi,xFullMin,xBlindMin,xBlindMax,xFullMax, minBDTcut,minBDTcut, type, -1,-1,isBlinded,doBDT,mode,resonVeto);
		delete metcal;
		delete mettrk;
		if(!pass) continue;
		
		
		rangeOSflags.clear();
		for(int k=0 ; k<MAXRANGE ; ++k)
		{
			TString rangeOS = getRangeOS(k);				
			getRangeOSdecisions(k,rangeOSflags,mOS1->at(0),mOS2->at(0),mSS->at(0));
		
			// if(k==INCLUSIVE && !rangeOSflags[INCLUSIVE])    continue;
			if(k==ONPHI1    && !rangeOSflags[ONPHI1   ])       continue;
			if(k==OFFPHI1   && !rangeOSflags[OFFPHI1  ])       continue;
			if(k==ONPHI2    && !rangeOSflags[ONPHI2   ])       continue;
			if(k==OFFPHI2   && !rangeOSflags[OFFPHI2  ])       continue;
			if(k==ONOMEGA1  && !rangeOSflags[ONOMEGA1 ])       continue;
			if(k==OFFOMEGA1 && !rangeOSflags[OFFOMEGA1])       continue;
			if(k==ONOMEGA2  && !rangeOSflags[ONOMEGA2 ])       continue;
			if(k==OFFOMEGA2 && !rangeOSflags[OFFOMEGA2])       continue;
			if(k==ONOMEGAPHI1  && !rangeOSflags[ONOMEGAPHI1 ]) continue;
			if(k==OFFOMEGAPHI1 && !rangeOSflags[OFFOMEGAPHI1]) continue;
			if(k==ONOMEGAPHI2  && !rangeOSflags[ONOMEGAPHI2 ]) continue;
			if(k==OFFOMEGAPHI2 && !rangeOSflags[OFFOMEGAPHI2]) continue;
			if(k==ONLOWOS1  && !rangeOSflags[ONLOWOS1 ])       continue;
			if(k==OFFLOWOS1 && !rangeOSflags[OFFLOWOS1])       continue;
			if(k==ONLOWOS2  && !rangeOSflags[ONLOWOS2 ])       continue;
			if(k==OFFLOWOS2 && !rangeOSflags[OFFLOWOS2])       continue;
			if(k==ONLOWSS   && !rangeOSflags[ONLOWSS  ])       continue;
			if(k==OFFLOWSS  && !rangeOSflags[OFFLOWSS ])       continue;
			
		
			TString hname = "";
			
			hname = channel+"_score"+rangeOS;        histos1[hname]->Fill(score->at(0));      histos2["bdt_vs_"+hname]->Fill(score->at(0),score->at(0));      histos2["m3body_vs_"+hname]->Fill(score->at(0),mass->at(0));       profiles["prof_bdt_"+hname]->Fill(score->at(0),score->at(0));      profiles["prof_m3body_"+hname]->Fill(score->at(0),mass->at(0));
			hname = channel+"_m3body"+rangeOS;       histos1[hname]->Fill(mass->at(0));       histos2["bdt_vs_"+hname]->Fill(mass->at(0),score->at(0));       histos2["m3body_vs_"+hname]->Fill(mass->at(0),mass->at(0));        profiles["prof_bdt_"+hname]->Fill(mass->at(0),score->at(0));       profiles["prof_m3body_"+hname]->Fill(mass->at(0),mass->at(0));
			hname = channel+"_PVNtrk"+rangeOS;       histos1[hname]->Fill(pvntrk->at(0));     histos2["bdt_vs_"+hname]->Fill(pvntrk->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(pvntrk->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(pvntrk->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(pvntrk->at(0),mass->at(0));
			hname = channel+"_dRmax"+rangeOS;        histos1[hname]->Fill(drmax->at(0));      histos2["bdt_vs_"+hname]->Fill(drmax->at(0),score->at(0));      histos2["m3body_vs_"+hname]->Fill(drmax->at(0),mass->at(0));       profiles["prof_bdt_"+hname]->Fill(drmax->at(0),score->at(0));      profiles["prof_m3body_"+hname]->Fill(drmax->at(0),mass->at(0));
			hname = channel+"_pT3body"+rangeOS;      histos1[hname]->Fill(pt->at(0)*MeV2GeV);         histos2["bdt_vs_"+hname]->Fill(pt->at(0)*MeV2GeV,score->at(0));         histos2["m3body_vs_"+hname]->Fill(pt->at(0)*MeV2GeV,mass->at(0));          profiles["prof_bdt_"+hname]->Fill(pt->at(0)*MeV2GeV,score->at(0));         profiles["prof_m3body_"+hname]->Fill(pt->at(0)*MeV2GeV,mass->at(0));
			hname = channel+"_mSS"+rangeOS;          histos1[hname]->Fill(mSS->at(0));        histos2["bdt_vs_"+hname]->Fill(mSS->at(0),score->at(0));        histos2["m3body_vs_"+hname]->Fill(mSS->at(0),mass->at(0));         profiles["prof_bdt_"+hname]->Fill(mSS->at(0),score->at(0));        profiles["prof_m3body_"+hname]->Fill(mSS->at(0),mass->at(0));
			hname = channel+"_mOS2"+rangeOS;         histos1[hname]->Fill(mOS2->at(0));       histos2["bdt_vs_"+hname]->Fill(mOS2->at(0),score->at(0));       histos2["m3body_vs_"+hname]->Fill(mOS2->at(0),mass->at(0));        profiles["prof_bdt_"+hname]->Fill(mOS2->at(0),score->at(0));       profiles["prof_m3body_"+hname]->Fill(mOS2->at(0),mass->at(0));
			hname = channel+"_mOS1"+rangeOS;         histos1[hname]->Fill(mOS1->at(0));       histos2["bdt_vs_"+hname]->Fill(mOS1->at(0),score->at(0));       histos2["m3body_vs_"+hname]->Fill(mOS1->at(0),mass->at(0));        profiles["prof_bdt_"+hname]->Fill(mOS1->at(0),score->at(0));       profiles["prof_m3body_"+hname]->Fill(mOS1->at(0),mass->at(0));
			hname = channel+"_isolation003"+rangeOS; histos1[hname]->Fill(iso003->at(0));     histos2["bdt_vs_"+hname]->Fill(iso003->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(iso003->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(iso003->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(iso003->at(0),mass->at(0));
			hname = channel+"_isolation010"+rangeOS; histos1[hname]->Fill(iso010->at(0));     histos2["bdt_vs_"+hname]->Fill(iso010->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(iso010->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(iso010->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(iso010->at(0),mass->at(0));
			hname = channel+"_isolation020"+rangeOS; histos1[hname]->Fill(iso020->at(0));     histos2["bdt_vs_"+hname]->Fill(iso020->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(iso020->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(iso020->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(iso020->at(0),mass->at(0));
			hname = channel+"_isolation030"+rangeOS; histos1[hname]->Fill(iso030->at(0));     histos2["bdt_vs_"+hname]->Fill(iso030->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(iso030->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(iso030->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(iso030->at(0),mass->at(0));
			hname = channel+"_trksfitprob"+rangeOS;  histos1[hname]->Fill(trkprob->at(0));    histos2["bdt_vs_"+hname]->Fill(trkprob->at(0),score->at(0));    histos2["m3body_vs_"+hname]->Fill(trkprob->at(0),mass->at(0));     profiles["prof_bdt_"+hname]->Fill(trkprob->at(0),score->at(0));    profiles["prof_m3body_"+hname]->Fill(trkprob->at(0),mass->at(0));
			hname = channel+"_maxpbalsig"+rangeOS;   histos1[hname]->Fill(maxpbalsig->at(0)); histos2["bdt_vs_"+hname]->Fill(maxpbalsig->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(maxpbalsig->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(maxpbalsig->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(maxpbalsig->at(0),mass->at(0));
			hname = channel+"_pvalue"+rangeOS;       histos1[hname]->Fill(pvalue->at(0));     histos2["bdt_vs_"+hname]->Fill(pvalue->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(pvalue->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(pvalue->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(pvalue->at(0),mass->at(0));
			hname = channel+"_pvalue_zoom"+rangeOS;  histos1[hname]->Fill(pvalue->at(0));     histos2["bdt_vs_"+hname]->Fill(pvalue->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(pvalue->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(pvalue->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(pvalue->at(0),mass->at(0));
	
			hname = channel+"_Lxy"+rangeOS;        histos1[hname]->Fill(lxy->at(0));      histos2["bdt_vs_"+hname]->Fill(lxy->at(0),score->at(0));      histos2["m3body_vs_"+hname]->Fill(lxy->at(0),mass->at(0));       profiles["prof_bdt_"+hname]->Fill(lxy->at(0),score->at(0));      profiles["prof_m3body_"+hname]->Fill(lxy->at(0),mass->at(0));
			hname = channel+"_dLxy"+rangeOS;       histos1[hname]->Fill(dlxy->at(0));     histos2["bdt_vs_"+hname]->Fill(dlxy->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(dlxy->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(dlxy->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(dlxy->at(0),mass->at(0));
			hname = channel+"_SLxy"+rangeOS;       histos1[hname]->Fill(lxySig->at(0));   histos2["bdt_vs_"+hname]->Fill(lxySig->at(0),score->at(0));   histos2["m3body_vs_"+hname]->Fill(lxySig->at(0),mass->at(0));    profiles["prof_bdt_"+hname]->Fill(lxySig->at(0),score->at(0));   profiles["prof_m3body_"+hname]->Fill(lxySig->at(0),mass->at(0));
			hname = channel+"_SLxy_zoom"+rangeOS;  histos1[hname]->Fill(lxySig->at(0));   histos2["bdt_vs_"+hname]->Fill(lxySig->at(0),score->at(0));   histos2["m3body_vs_"+hname]->Fill(lxySig->at(0),mass->at(0));    profiles["prof_bdt_"+hname]->Fill(lxySig->at(0),score->at(0));   profiles["prof_m3body_"+hname]->Fill(lxySig->at(0),mass->at(0));
			hname = channel+"_a0xy"+rangeOS;       histos1[hname]->Fill(a0xy->at(0));     histos2["bdt_vs_"+hname]->Fill(a0xy->at(0),score->at(0));     histos2["m3body_vs_"+hname]->Fill(a0xy->at(0),mass->at(0));      profiles["prof_bdt_"+hname]->Fill(a0xy->at(0),score->at(0));     profiles["prof_m3body_"+hname]->Fill(a0xy->at(0),mass->at(0));
			hname = channel+"_da0xy"+rangeOS;      histos1[hname]->Fill(da0xy->at(0));    histos2["bdt_vs_"+hname]->Fill(da0xy->at(0),score->at(0));    histos2["m3body_vs_"+hname]->Fill(da0xy->at(0),mass->at(0));     profiles["prof_bdt_"+hname]->Fill(da0xy->at(0),score->at(0));    profiles["prof_m3body_"+hname]->Fill(da0xy->at(0),mass->at(0));
			hname = channel+"_Sa0xy"+rangeOS;      histos1[hname]->Fill(a0xySig->at(0));  histos2["bdt_vs_"+hname]->Fill(a0xySig->at(0),score->at(0));  histos2["m3body_vs_"+hname]->Fill(a0xySig->at(0),mass->at(0));   profiles["prof_bdt_"+hname]->Fill(a0xySig->at(0),score->at(0));  profiles["prof_m3body_"+hname]->Fill(a0xySig->at(0),mass->at(0));
			hname = channel+"_Sa0xy_zoom"+rangeOS; histos1[hname]->Fill(a0xySig->at(0));  histos2["bdt_vs_"+hname]->Fill(a0xySig->at(0),score->at(0));  histos2["m3body_vs_"+hname]->Fill(a0xySig->at(0),mass->at(0));   profiles["prof_bdt_"+hname]->Fill(a0xySig->at(0),score->at(0));  profiles["prof_m3body_"+hname]->Fill(a0xySig->at(0),mass->at(0));

			if(!isFilled) { hname = channel+"_calo_met"+rangeOS; histos1[hname]->Fill(calo_met*MeV2GeV);             histos2["bdt_vs_"+hname]->Fill(calo_met*MeV2GeV,score->at(0));            histos2["m3body_vs_"+hname]->Fill(calo_met*MeV2GeV,mass->at(0));             profiles["prof_bdt_"+hname]->Fill(calo_met*MeV2GeV,score->at(0));            profiles["prof_m3body_"+hname]->Fill(calo_met*MeV2GeV,mass->at(0)); }
			hname = channel+"_calo_mt"+rangeOS;                  histos1[hname]->Fill(calo_mt->at(0)*MeV2GeV);       histos2["bdt_vs_"+hname]->Fill(calo_mt->at(0)*MeV2GeV,score->at(0));      histos2["m3body_vs_"+hname]->Fill(calo_mt->at(0)*MeV2GeV,mass->at(0));       profiles["prof_bdt_"+hname]->Fill(calo_mt->at(0)*MeV2GeV,score->at(0));      profiles["prof_m3body_"+hname]->Fill(calo_mt->at(0)*MeV2GeV,mass->at(0));
			hname = channel+"_calo_dphi3mu"+rangeOS;             histos1[hname]->Fill(calo_dphimet->at(0));  histos2["bdt_vs_"+hname]->Fill(calo_dphimet->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(calo_dphimet->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(calo_dphimet->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(calo_dphimet->at(0),mass->at(0));

			if(!isFilled) { hname = channel+"_trk_met"+rangeOS; histos1[hname]->Fill(trk_met*MeV2GeV); histos2["bdt_vs_"+hname]->Fill(trk_met*MeV2GeV,score->at(0));            histos2["m3body_vs_"+hname]->Fill(trk_met*MeV2GeV,mass->at(0));             profiles["prof_bdt_"+hname]->Fill(trk_met*MeV2GeV,score->at(0));            profiles["prof_m3body_"+hname]->Fill(trk_met*MeV2GeV,mass->at(0));             }
			hname = channel+"_trk_mt"+rangeOS;      histos1[hname]->Fill(trk_mt->at(0)*MeV2GeV);       histos2["bdt_vs_"+hname]->Fill(trk_mt->at(0)*MeV2GeV,score->at(0));      histos2["m3body_vs_"+hname]->Fill(trk_mt->at(0)*MeV2GeV,mass->at(0));       profiles["prof_bdt_"+hname]->Fill(trk_mt->at(0)*MeV2GeV,score->at(0));      profiles["prof_m3body_"+hname]->Fill(trk_mt->at(0)*MeV2GeV,mass->at(0));
			hname = channel+"_trk_dphi3mu"+rangeOS; histos1[hname]->Fill(trk_dphimet->at(0));  histos2["bdt_vs_"+hname]->Fill(trk_dphimet->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(trk_dphimet->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(trk_dphimet->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(trk_dphimet->at(0),mass->at(0));

			if(!isFilled) { hname = channel+"_calo_trk_dphi"+rangeOS; histos1[hname]->Fill(dphimets->at(0));  histos2["bdt_vs_"+hname]->Fill(dphimets->at(0),score->at(0));  histos2["m3body_vs_"+hname]->Fill(dphimets->at(0),mass->at(0));   profiles["prof_bdt_"+hname]->Fill(dphimets->at(0),score->at(0));  profiles["prof_m3body_"+hname]->Fill(dphimets->at(0),mass->at(0));  }
			hname = channel+"_dptreltrk"+rangeOS;                     histos1[hname]->Fill(dptreltrk->at(0)); histos2["bdt_vs_"+hname]->Fill(dptreltrk->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(dptreltrk->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(dptreltrk->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(dptreltrk->at(0),mass->at(0));
			hname = channel+"_dptrelcal"+rangeOS;                     histos1[hname]->Fill(dptrelcal->at(0)); histos2["bdt_vs_"+hname]->Fill(dptrelcal->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(dptrelcal->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(dptrelcal->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(dptrelcal->at(0),mass->at(0));

			hname = channel+"_ht"+rangeOS;              histos1[hname]->Fill(ht->at(0)*MeV2GeV);               histos2["bdt_vs_"+hname]->Fill(ht->at(0)*MeV2GeV,score->at(0));              histos2["m3body_vs_"+hname]->Fill(ht->at(0)*MeV2GeV,mass->at(0));               profiles["prof_bdt_"+hname]->Fill(ht->at(0)*MeV2GeV,score->at(0));              profiles["prof_m3body_"+hname]->Fill(ht->at(0)*MeV2GeV,mass->at(0));              
			hname = channel+"_ht_dphimet_calo"+rangeOS; histos1[hname]->Fill(calo_dphimet_ht->at(0));  histos2["bdt_vs_"+hname]->Fill(calo_dphimet_ht->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(calo_dphimet_ht->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(calo_dphimet_ht->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(calo_dphimet_ht->at(0),mass->at(0));
			hname = channel+"_ht_dphimet_trk"+rangeOS;  histos1[hname]->Fill(trk_dphimet_ht->at(0));   histos2["bdt_vs_"+hname]->Fill(trk_dphimet_ht->at(0),score->at(0));  histos2["m3body_vs_"+hname]->Fill(trk_dphimet_ht->at(0),mass->at(0));   profiles["prof_bdt_"+hname]->Fill(trk_dphimet_ht->at(0),score->at(0));  profiles["prof_m3body_"+hname]->Fill(trk_dphimet_ht->at(0),mass->at(0));
			hname = channel+"_calo_mht"+rangeOS;        histos1[hname]->Fill(calo_mht->at(0)*MeV2GeV);         histos2["bdt_vs_"+hname]->Fill(calo_mht->at(0)*MeV2GeV,score->at(0));        histos2["m3body_vs_"+hname]->Fill(calo_mht->at(0)*MeV2GeV,mass->at(0));         profiles["prof_bdt_"+hname]->Fill(calo_mht->at(0)*MeV2GeV,score->at(0));        profiles["prof_m3body_"+hname]->Fill(calo_mht->at(0)*MeV2GeV,mass->at(0));
			hname = channel+"_trk_mht"+rangeOS;         histos1[hname]->Fill(trk_mht->at(0)*MeV2GeV);          histos2["bdt_vs_"+hname]->Fill(trk_mht->at(0)*MeV2GeV,score->at(0));         histos2["m3body_vs_"+hname]->Fill(trk_mht->at(0)*MeV2GeV,mass->at(0));          profiles["prof_bdt_"+hname]->Fill(trk_mht->at(0)*MeV2GeV,score->at(0));         profiles["prof_m3body_"+hname]->Fill(trk_mht->at(0)*MeV2GeV,mass->at(0));
			
			hname = channel+"_dR3body_ht"+rangeOS; histos1[hname]->Fill(dr_ht->at(0));      histos2["bdt_vs_"+hname]->Fill(dr_ht->at(0),score->at(0));      histos2["m3body_vs_"+hname]->Fill(dr_ht->at(0),mass->at(0));       profiles["prof_bdt_"+hname]->Fill(dr_ht->at(0),score->at(0));      profiles["prof_m3body_"+hname]->Fill(dr_ht->at(0),mass->at(0));      
			hname = channel+"_njets"+rangeOS;      histos1[hname]->Fill(njets->at(0));      histos2["bdt_vs_"+hname]->Fill(njets->at(0),score->at(0));      histos2["m3body_vs_"+hname]->Fill(njets->at(0),mass->at(0));       profiles["prof_bdt_"+hname]->Fill(njets->at(0),score->at(0));      profiles["prof_m3body_"+hname]->Fill(njets->at(0),mass->at(0));
			hname = channel+"_muonauthor"+rangeOS; histos1[hname]->Fill(muonauthor->at(0)); histos2["bdt_vs_"+hname]->Fill(muonauthor->at(0),score->at(0)); histos2["m3body_vs_"+hname]->Fill(muonauthor->at(0),mass->at(0));  profiles["prof_bdt_"+hname]->Fill(muonauthor->at(0),score->at(0)); profiles["prof_m3body_"+hname]->Fill(muonauthor->at(0),mass->at(0));
			
			if(k>0) continue;
			hname = channel+"_Dalitz0"+rangeOS; histos2[hname]->Fill(mOS1->at(0)*MeV2GeV*mOS1->at(0)*MeV2GeV,mOS2->at(0)*MeV2GeV*mOS2->at(0)*MeV2GeV);
			hname = channel+"_Dalitz1"+rangeOS; histos2[hname]->Fill(mOS1->at(0)*MeV2GeV*mOS1->at(0)*MeV2GeV,mSS->at(0)*MeV2GeV*mSS->at(0)*MeV2GeV);
			hname = channel+"_Dalitz2"+rangeOS; histos2[hname]->Fill(mOS2->at(0)*MeV2GeV*mOS2->at(0)*MeV2GeV,mSS->at(0)*MeV2GeV*mSS->at(0)*MeV2GeV);
		}
		
		isFilled = true; // fill floats only once.
		
		ncandidates++;
	}
	
	cout << t->GetName() << ": ncandidates is " << ncandidates << " for x>" << minBDTcut << endl;
	
	return ncandidates;
}


void replotBDTvars(float mMinSBleft, float mMaxSBleft, float mMinSBright, float mMaxSBright)
{	
	// use plain black on white colors
	Int_t icol=0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects
	
	// set the paper & margin sizes
	gStyle->SetPaperSize(20,26);
	
	// set margin sizes
	// gStyle->SetPadTopMargin(0.03);
	// gStyle->SetPadRightMargin(0.07);
	// gStyle->SetPadBottomMargin(0.125);
	// gStyle->SetPadLeftMargin(0.115);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.12);
	
	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.05);
	gStyle->SetTitleYOffset(0.95);
	
	// use large fonts
	//Int_t font=72; // Helvetica italics
	Int_t font=42; // Helvetica
	Double_t tsize=0.05;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	
	gStyle->SetLabelSize(tsize*0.85,"x");
	gStyle->SetTitleSize(tsize*1.10,"x");
	gStyle->SetLabelSize(tsize*0.85,"y");
	gStyle->SetTitleSize(tsize*1.10,"y");
	gStyle->SetLabelSize(tsize*0.85,"z");
	gStyle->SetTitleSize(tsize*1.10,"z");
	
	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2.);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
	// get rid of X error bars 
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);
	
	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
	
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	
	
	
	// TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.loose.root","READ");
	// TTree* tD = (TTree*)f->Get("flatout_Data");
	// TTree* tS = (TTree*)f->Get("flatout_Wtaunu_3mu");
	// 
	// TFile* fTight = new TFile("flatout.periodall+MC.muons.mva.analysis.n0.j0.loose.root","READ");
	// TTree* tDtight = (TTree*)fTight->Get("flatout_Data");
	// TTree* tStight = (TTree*)fTight->Get("flatout_Wtaunu_3mu");
	
	TFile* f = new TFile("flatout.periodall+MC.muons.mva.analysis.n0.j0.loose.root","READ");
	TTree* tD = (TTree*)f->Get("flatout_Data");
	TTree* tS = (TTree*)f->Get("flatout_Wtaunu_200k_3mu");

	makeAtlasLabel();
	makeLegend();
	

	// TString mSBleft  = "(vtx_mass>"+mMinSBleft+" && vtx_mass<"+mMaxSBleft+")";
	// TString mSBright = "(vtx_mass>"+mMinSBright+" && vtx_mass<"+mMaxSBright+")";
	// TString mSB      = "("+mSBleft+" || "+mSBright+")";
	// TString isLoose  = "(pass_loose==1)";
	// TString oneCand  = "(@vtx_mass->size()==1)"; // for siganl only
	// // TString basecuts = isLoose+" && "+oneCand;
	// TString basecuts = isLoose;
	// TString mycuts = "("+basecuts+")";
	// TString mycutb = "("+basecuts+" && "+mSB+")";
	
	// TString mSBleft  = "(vtx_mass>"+mMinSBleft+" && vtx_mass<"+mMaxSBleft+")";
	// TString mSBright = "(vtx_mass>"+mMinSBright+" && vtx_mass<"+mMaxSBright+")";
	// TString mSB      = mSBleft+" || "+mSBright;
	// TCut    sidebandsOnly(mSB.Data());
	// TCut    isLoose  = "pass_loose==1";
	// TCut    oneCand  = "@vtx_mass->size()==1"; // for siganl only
	// TCut mycuts = isLoose;
	// TCut mycutb = isLoose+sidebandsOnly;
	// // TCut mycuts = isLoose+oneCand;
	// // TCut mycutb = isLoose+sidebandsOnly;
	// cout << "mycuts=" << mycuts << endl;
	// cout << "mycutb=" << mycutb << endl;
	// 
	// tS->SetEventList(0);
	// tS->Draw(">>elistS",mycuts);
	// TEventList* elistS = (TEventList*)gDirectory->Get("elistS");
	// Float_t npassedS = elistS->GetN();  // number of events to pass cuts
	// 
	// tD->SetEventList(0);
	// tD->Draw(">>elistD",mycutb);
	// TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
	// Float_t npassedD = elistD->GetN();  // number of events to pass cuts
	// 
	// tD->SetEventList(elistD);
	// tS->SetEventList(elistS);

	
	TMapTSP2TH1 histos1;
	TMapTSP2TH2 histos2;
	TMapTSP2TProfile profiles;
	
	vector<TString> channels;
	channels.push_back("Data");
	channels.push_back("Wtaunu_3mu");
	
	for(unsigned int i=0 ; i<channels.size() ; ++i)
	{
		for(unsigned int j=0 ; j<2 ; ++j)
		{
			TString channel = (j==1) ? channels[i] : channels[i]+"_tight";
			
			for(int k=0 ; k<MAXRANGE ; ++k)
			{
				TString rangeOS = getRangeOS(k);
				
				addHist(histos1,histos2,profiles,channel,rangeOS,"score",        ";BDT score;Events",40,-1,+1);
				addHist(histos1,histos2,profiles,channel,rangeOS,"m3body",       ";#it{m}_{3#mu} [MeV];Events",22,1450,2110);
				
				addHist(histos1,histos2,profiles,channel,rangeOS,"PVNtrk",       ";#it{N}_{trk}^{PV};Events",50,0,200);
				
				addHist(histos1,histos2,profiles,channel,rangeOS,"dRmax",        ";3body #Delta#it{R}_{max};Events", 30,0.,0.3);
				addHist(histos1,histos2,profiles,channel,rangeOS,"pT3body",      ";#it{p}_{T}^{3#mu} [GeV];Events", 50,0.,100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"mSS",          ";#it{m}_{SS} [MeV];Events", 50,0.,2100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"mOS2",         ";#it{m}_{OS2} [MeV];Events", 50,0.,2100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"mOS1",         ";#it{m}_{OS1} [MeV];Events", 50,0.,2100.);
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"isolation003", ";#Sigma#it{p}_{T}^{trk}(cone #Delta#it{R}_{max}+0.03)/#it{p}_{T}^{3#mu};Events", 60,0.,0.3);
				addHist(histos1,histos2,profiles,channel,rangeOS,"isolation010", ";#Sigma#it{p}_{T}^{trk}(cone #Delta#it{R}_{max}+0.10)/#it{p}_{T}^{3#mu};Events", 60,0.,0.3);
				addHist(histos1,histos2,profiles,channel,rangeOS,"isolation020", ";#Sigma#it{p}_{T}^{trk}(cone #Delta#it{R}_{max}+0.20)/#it{p}_{T}^{3#mu};Events", 60,0.,0.3);
				addHist(histos1,histos2,profiles,channel,rangeOS,"isolation030", ";#Sigma#it{p}_{T}^{trk}(cone #Delta#it{R}_{max}+0.30)/#it{p}_{T}^{3#mu};Events", 50,0.,1.0);
				addHist(histos1,histos2,profiles,channel,rangeOS,"trksfitprob",  ";#it{P}_{trks};Events",50,0.,1.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"maxpbalsig",   ";#it{#sigma}_{#it{p}-balance}^{max};Events", 70,-3.,+4.);
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"pvalue",       ";#it{p}-value (3#mu vertex);Events", 50,0.,1.);		
				addHist(histos1,histos2,profiles,channel,rangeOS,"pvalue_zoom",  ";#it{p}-value (3#mu vertex);Events", 100,0.,0.5);		
				addHist(histos1,histos2,profiles,channel,rangeOS,"Lxy",          ";#it{L}_{xy} [#mum];Events", 52,-1.,+12.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"a0xy",         ";#it{a}_{0}^{xy} [#mum];Events", 50,0.,0.1);
				addHist(histos1,histos2,profiles,channel,rangeOS,"dLxy",         ";#Delta#it{L}_{xy} [#mum];Events", 50,0.,+1.5);
				addHist(histos1,histos2,profiles,channel,rangeOS,"SLxy",         ";#it{S}(#it{L}_{xy});Events", 60,-10.,+50.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"SLxy_zoom",    ";#it{S}(#it{L}_{xy});Events", 60,-10.,+20.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"da0xy",        ";#Delta#it{a}_{0}^{xy} [#mum];Events", 50,0.,+0.05);
				addHist(histos1,histos2,profiles,channel,rangeOS,"Sa0xy",        ";#it{S}(#it{a}_{0}^{xy});Events", 50,0.,25.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"Sa0xy_zoom",   ";#it{S}(#it{a}_{0}^{xy});Events", 50,0.,3.);
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"calo_met",      ";#it{E}_{T,cal}^{miss} [GeV];Events",45,10.,100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"calo_mt",       ";#it{m}_{T}^{cal} [GeV];Events",65,20.,150.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"calo_dphi3mu",  ";#Delta#it{#phi}_{3#mu}^{cal};Events",32,0.,TMath::Pi());
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"trk_met",      ";#it{E}_{T,trk}^{miss} [GeV];Events",45,10.,100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"trk_mt",       ";#it{m}_{T}^{trk} [GeV];Events",60,30.,150.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"trk_dphi3mu",  ";#Delta#it{#phi}_{3#mu}^{trk};Events",32,0.,TMath::Pi());
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"calo_trk_dphi", ";#Delta#it{#phi}_{trk}^{cal};Events",32,0.,TMath::Pi());
				addHist(histos1,histos2,profiles,channel,rangeOS,"dptreltrk",     ";p_{T}^{3#mu}/#it{E}_{T,trk}^{miss}-1;Events",60,-1.,+5.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"dptrelcal",     ";p_{T}^{3#mu}/#it{E}_{T,cal}^{miss}-1;Events",60,-1.,+5.);
				                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"ht",              ";#it{#Sigma}_{T} [GeV];Events",100,0.,100.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"ht_dphimet_calo", ";#Delta#it{#phi}_{#it{#Sigma}_{T}}^{cal};Events",32,0.,TMath::Pi());
				addHist(histos1,histos2,profiles,channel,rangeOS,"ht_dphimet_trk",  ";#Delta#it{#phi}_{#it{#Sigma}_{T}}^{trk};Events",32,0.,TMath::Pi());
				addHist(histos1,histos2,profiles,channel,rangeOS,"calo_mht",        ";#it{m}_{#it{#Sigma}_{T}}^{cal} [GeV];Events",65,20.,150.);
				addHist(histos1,histos2,profiles,channel,rangeOS,"trk_mht",         ";#it{m}_{#it{#Sigma}_{T}}^{trk} [GeV];Events",60,30.,150.);
        		                       
				addHist(histos1,histos2,profiles,channel,rangeOS,"dR3body_ht",      ";#Delta#it{R}(#it{p}_{3#mu},#it{#Sigma});Events",32,0.,2.*TMath::Pi());
				addHist(histos1,histos2,profiles,channel,rangeOS,"njets",           ";#it{N}_{jets};Events",5,-0.5,4.5);
				addHist(histos1,histos2,profiles,channel,rangeOS,"muonauthor",      ";Muon author;Events",2,0.,2.); 
            	
				if(k>0) continue;
				addHist(histos2,channel,"Dalitz0",          ";m_{OS1}^{2} [GeV^{2}];m_{OS2}^{2} [GeV^{2}];Events",50,0.,1.9*1.9, 50,0.,1.9*1.9); 
				addHist(histos2,channel,"Dalitz1",          ";m_{OS1}^{2} [GeV^{2}];m_{SS}^{2} [GeV^{2}];Events",50,0.,1.9*1.9, 50,0.,1.9*1.9); 
				addHist(histos2,channel,"Dalitz2",          ";m_{OS2}^{2} [GeV^{2}];m_{SS}^{2} [GeV^{2}];Events",50,0.,1.9*1.9, 50,0.,1.9*1.9);
			}
		}
	}
		
	// tD->Draw("vtx_pvNtrk>>Data_PVNtrk");
	// tS->Draw("vtx_pvNtrk>>Wtaunu_3mu_PVNtrk");
	// 
	// tD->Draw("vtx_mOS2>>Data_mOS2");
	// tS->Draw("vtx_mOS2>>Wtaunu_3mu_mOS2");
	// 
	// tD->Draw("vtx_mOS1>>Data_mOS1");
	// tS->Draw("vtx_mOS1>>Wtaunu_3mu_mOS1");
	// 
	// tD->Draw("vtx_isolation003>>Data_isolation003");
	// tS->Draw("vtx_isolation003>>Wtaunu_3mu_isolation003");
	// 
	// tD->Draw("trks_fitprob>>Data_trksfitprob");
	// tS->Draw("trks_fitprob>>Wtaunu_3mu_trksfitprob");
	// 
	// tD->Draw("vtx_pval>>Data_pvalue");
	// tS->Draw("vtx_pval>>Wtaunu_3mu_pvalue");
	
	readData(tD, histos1,histos2,profiles,"Data",             mMinSBleft,mMaxSBleft,mMinSBright,mMaxSBright);
	readData(tS, histos1,histos2,profiles,"Wtaunu_3mu",       mMinSBleft,mMaxSBleft,mMinSBright,mMaxSBright);
	readData(tD, histos1,histos2,profiles,"Data_tight",       mMinSBleft,mMaxSBleft,mMinSBright,mMaxSBright,true);
	readData(tS, histos1,histos2,profiles,"Wtaunu_3mu_tight", mMinSBleft,mMaxSBleft,mMinSBright,mMaxSBright,true);
	
	
	TMapTSf rescales;
	for(TMapTSP2TH1::iterator hit=histos1.begin() ; hit!=histos1.end() ; ++hit)
	{
		if(!hit->first.Contains("tight")) continue;
		TString name = hit->first; name.ReplaceAll("_tight","");
		float rescale = Sum(hit->second)/Sum(histos1[name]);
		rescales.insert(make_pair(hit->first,rescale));
	}
	
	
	for(TMapTSP2TH1::iterator hit=histos1.begin() ; hit!=histos1.end() ; ++hit) { if(hit->first.Contains("Data")) continue; NormToEntries(hit->second); }
	// for(TMapTSP2TH2::iterator hit=histos2.begin() ; hit!=histos2.end() ; ++hit) { if(hit->first.Contains("Data")) continue; NormToEntries(hit->second); }

	for(TMapTSP2TH1::iterator hit=histos1.begin() ; hit!=histos1.end() ; ++hit)
	{
		if(!hit->first.Contains("tight")) continue;
		if(hit->first.Contains("Data"))
		{
			hit->second->SetLineColor(kBlack);
			hit->second->SetMarkerColor(kBlack);
			hit->second->SetMarkerStyle(20);
			hit->second->SetLineWidth(1);
			hit->second->SetMarkerSize(1.2);
		}
		else
		{
			// hit->second->SetFillStyle(3005);
			hit->second->SetFillColor(kGray+1);
			hit->second->SetLineColor(kGray+1);
			hit->second->SetLineWidth(1);
		}
	}
	
	
	float max = 0;
	
	TCanvas* cnv = NULL;
	
	TString pdffilename = "paperplots";
	if(!isBlinded) pdffilename += ".unblinded";
	if(!resonVeto) pdffilename += ".noresveto";
	
	cnv = new TCanvas("cnv","",800,600);
	cnv->SaveAs("figures/"+pdffilename+".pdf(");
	
	for(int k=0 ; k<MAXRANGE ; ++k)
	{
		TString rangeOS = getRangeOS(k);
		
		if(k==INCLUSIVE) continue;
		
		delete cnv;
		cnv = new TCanvas("cnv","",800,600);
		cnv->SaveAs("figures/"+pdffilename+"_resonances"+rangeOS+".pdf(");
	}
	
	_INF(1,"");	
	
	for(int k=0 ; k<MAXRANGE ; ++k)
	{
		TString rangeOS = getRangeOS(k);
		
		plot(pdffilename,"score", histos1, histos2,profiles, rescales, legTM,rangeOS);
		plot(pdffilename,"m3body", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"PVNtrk", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"dRmax", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"pT3body", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"mSS", histos1, histos2,profiles, rescales, legTM,rangeOS);
		plot(pdffilename,"mOS1", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"mOS2", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"isolation003", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"isolation010", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"isolation020", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"isolation030", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"trksfitprob", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"maxpbalsig", histos1, histos2,profiles, rescales, legTL,rangeOS);
		plot(pdffilename,"pvalue", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"pvalue_zoom", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"Lxy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"dLxy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"SLxy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"SLxy_zoom", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"a0xy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"da0xy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"Sa0xy", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"Sa0xy_zoom", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"calo_met", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"trk_met", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"calo_mt", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"trk_mt", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"calo_dphi3mu", histos1, histos2,profiles, rescales, legL,rangeOS);
		plot(pdffilename,"trk_dphi3mu", histos1, histos2,profiles, rescales, legL,rangeOS);
		plot(pdffilename,"calo_trk_dphi", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"dptreltrk", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"dptrelcal", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"ht", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"ht_dphimet_calo", histos1, histos2,profiles, rescales, legL,rangeOS);
		plot(pdffilename,"ht_dphimet_trk", histos1, histos2,profiles, rescales, legL,rangeOS);
		plot(pdffilename,"calo_mht", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"trk_mht", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"dR3body_ht", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"njets", histos1, histos2,profiles, rescales, legR,rangeOS);
		plot(pdffilename,"muonauthor", histos1, histos2,profiles, rescales, legR,rangeOS);
	}
	
	vector<TString>  vars;
	vector<TLegend*> legs;
	
	vars.clear(); legs.clear();
	vars.push_back("calo_mt");  vars.push_back("trk_met");   vars.push_back("isolation020");
	vars.push_back("ht");       vars.push_back("trk_mt");    vars.push_back("calo_trk_dphi");
	vars.push_back("calo_met"); vars.push_back("dptreltrk"); vars.push_back("calo_dphi3mu");
	legs.push_back(legR); legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR); legs.push_back(legL);
	plot9(pdffilename,"BDTinputs01to09",vars,legs,histos1);
	plot9(pdffilename,"BDTinputs01to09.logy",vars,legs,histos1,true);
	
	vars.clear(); legs.clear();
	vars.push_back("pvalue");    vars.push_back("Sa0xy");   vars.push_back("trksfitprob");
	vars.push_back("pT3body");   vars.push_back("PVNtrk");  vars.push_back("SLxy");
	vars.push_back("dptrelcal");
	legs.push_back(legR); legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR); legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR);
	plot9(pdffilename,"BDTinputs10to16",vars,legs,histos1);
	plot9(pdffilename,"BDTinputs10to16.logy",vars,legs,histos1,true);
	
	vars.clear(); legs.clear();
	vars.push_back("mSS");          vars.push_back("mOS1");            vars.push_back("mOS2");
	vars.push_back("isolation030"); vars.push_back("ht_dphimet_calo"); vars.push_back("ht_dphimet_trk");
	vars.push_back("dRmax");        vars.push_back("Lxy");             vars.push_back("a0xy");
	legs.push_back(legTM); legs.push_back(legR); legs.push_back(legR);
	legs.push_back(legR);  legs.push_back(legL); legs.push_back(legL);
	legs.push_back(legR);  legs.push_back(legR); legs.push_back(legR);
	plot9(pdffilename,"Others17to25",vars,legs,histos1);
	plot9(pdffilename,"Others17to25.logy",vars,legs,histos1,true);
	
	vars.clear(); legs.clear();
	vars.push_back("m3body");   vars.push_back("score");
	legs.push_back(legR); legs.push_back(legTM);
	plot2(pdffilename,"Others26to27",vars,legs,histos1);
	plot2(pdffilename,"Others26to27.logy",vars,legs,histos1,true);
	
	
	
	
	// delete cnv;
	cnv = new TCanvas("cnv","",800,600);
	cnv->SaveAs("figures/"+pdffilename+".pdf)");
	
	for(int k=0 ; k<MAXRANGE ; ++k)
	{
		TString rangeOS = getRangeOS(k);
		
		if(k==INCLUSIVE) continue;
		
		delete cnv;
		cnv = new TCanvas("cnv","",800,600);
		cnv->SaveAs("figures/"+pdffilename+"_resonances"+rangeOS+".pdf)");
	}
}
