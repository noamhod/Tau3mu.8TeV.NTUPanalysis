#ifndef HIST_H
#define HIST_H

#include "std.h"
#include "type.h"
#include "const.h"
#include "enums.h"
#include "count.h"

TMapTSP2TH1 histos1;
TMapTSP2TH2 histos2;

TMapuiTS channels;
TMapTSTS labels;
TMapTSTS legoptions;
TMapTSTC colors;
TMapTSI  patterns;
TMapTSI  drawchannels;
TMapTSi  categories;
TMapiTS  ordered_categories;

vector<TString> triggers;

TDirectory* olddir = gDirectory;


float mT2(float MET, float phiMET, float pTX, float phiX)
{
	return 2*pTX*MET*(1-TMath::Cos(phiX-phiMET));
}
float mT(float MET, float phiMET, float pTX, float phiX)
{
	return TMath::Sqrt(mT2(MET,phiMET,pTX,phiX));
}
// float deltaR(float eta1, float phi1, float eta2, float phi2)
// {
// 	return TMath::Sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
// }



void setStyle()
{
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.12);
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
	gStyle->SetTitleY(0.94); //title Y location
	gStyle->SetTitleW(0.5); //title width
	gStyle->SetTitleH(0.05); //title height
	gStyle->SetTitleBorderSize(0);
}

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm << fixed << x;
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
TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0.2,0.9,50);");
TExec* exeGray  = new TExec("exeGray",  "rgbPalette(1,1,1,50);");

void makeTriggers()
{
	//////////////////
	//// Trigger names
	triggers.push_back("EF_3mu4T");
	triggers.push_back("EF_3mu6");
	triggers.push_back("EF_3mu6_MSonly");
	triggers.push_back("EF_2mu13");
	triggers.push_back("EF_mu18_tight_mu8_EFFS");
	triggers.push_back("EF_mu18_tight_2mu4_EFFS");
	triggers.push_back("EF_2mu8_EFxe30_tclcw");
	// triggers.push_back("EF_mu24_tight_EFxe40");
	// triggers.push_back("EF_mu24i_tight");
	// triggers.push_back("EF_mu36_tight");
}

void makeCategories()
{
	int     order = -1;
	TString name  = "";
	
	name = "3muons";              order = MUONS3;          categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muons1tpmuA";        order = MUONS2TPA1;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muons1tpmuB";        order = MUONS2TPB1;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muons1calomu";       order = MUONS2CALO1;     categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muons2tpmuA";        order = MUONS1TPA2;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muons2tpmuB";        order = MUONS1TPB2;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muons1tpmuA1tpmuB";  order = MUONS1TPA1TPB1;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muons3tpmuA";        order = MUONS0TPA3;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muons3tpmuB";        order = MUONS0TPB3;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muons2tpmuA1tmpmuB"; order = MUONS0TPA2TPB1;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muons2tpmuA1calomu"; order = MUONS0TPA2CALO1; categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muons1tpmuA2tmpmuB"; order = MUONS0TPA1TPB2;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	                                                                                                                              
	name = "3muid";              order = MUID3;          categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muid1tpmuA";        order = MUID2TPA1;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muid1tpmuB";        order = MUID2TPB1;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muid2tpmuA";        order = MUID2CALO1;     categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muid2tpmuB";        order = MUID1TPA2;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "1muid1tpmuA1tpmuB";  order = MUID1TPB2;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "2muid1calomu";       order = MUID1TPA1TPB1;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muid3tpmuA";        order = MUID0TPA3;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muid3tpmuB";        order = MUID0TPB3;      categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muid2tpmuA1tmpmuB"; order = MUID0TPA2TPB1;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muid1tpmuA2tmpmuB"; order = MUID0TPA2CALO1; categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
	name = "0muid2tpmuA1calomu"; order = MUID0TPA1TPB2;  categories.insert(make_pair(name,order)); ordered_categories.insert(make_pair(order,name));
}

void properties(TString channel, TString label, Color_t color, Int_t pattern, TString legoption, Int_t isdraw = 1)
{
	labels.insert(make_pair(channel,     label));
	colors.insert(make_pair(channel,     color));
	patterns.insert(make_pair(channel,   pattern));
	legoptions.insert(make_pair(channel, legoption));
	drawchannels.insert(make_pair(channel, isdraw));
}
bool isWsignal(TString name)
{
	if(name=="Wtaunu_3mu") return true;
	return false;
}
bool isSignal(TString name)
{
	if(name.Contains("3mu")) return true;
	return false;
}
bool isData(TString name)
{
	if(name.Contains("Data") || name.Contains("period")) return true;
	return false;
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
void NormToEntries(TH1* h)
{
	Double_t entries = Sum(h);
	Double_t scale = (entries==0.) ? 1. : 1./entries;
	h->Scale(scale);
	h->GetYaxis()->SetTitle("Normalized");
}
void NormToEntries(TH2* h)
{
	Double_t entries = Sum(h);
	Double_t scale = (entries==0.) ? 1. : 1./entries;
	h->Scale(scale);
	h->GetZaxis()->SetTitle("Normalized");
}
void NormTo1stBin(TH1* h)
{
	Double_t frstbin = h->GetBinContent(1);
	Double_t scale = (frstbin==0.) ? 1. : 1./frstbin;
	h->Scale(scale);
	h->GetYaxis()->SetTitle("Normalized to 1st bin");
}
void NormToTriggers(TH2* h2, TH1* h1)
{
	for(Int_t A=1 ; A<=h2->GetNbinsX() ; ++A)
	{
		for(Int_t B=1 ; B<=h2->GetNbinsY() ; ++B)
		{
			Double_t nAandB = h2->GetBinContent(A,B);
			Double_t nAorB = h1->GetBinContent(A)+h1->GetBinContent(B);
			Double_t correlation = (nAorB>0) ? nAandB/nAorB : 0;
			h2->SetBinContent(A,B,correlation);
		}
	}
}

TCanvas* separatorPage(TString title, Int_t nlines, TString* lines)
{
	TCanvas *cs = new TCanvas("txt","txt",0,0,800,800);
	cs->Range(0,0,1,1);
	float x   = 0.05;
	float y   = 0.9;
	float dy  = (nlines>8) ? 0.08 : 0.1;
	
	TLatex* txt0 = new TLatex(x,y,"#scale[1.5]{#font[22]{"+title+"}}");
	txt0->SetTextAlign(11);
	txt0->SetTextAngle(0);
	txt0->SetTextColor(kRed+2);
	txt0->Draw();
	
	vector<TLatex*> vtxt;
	for(Int_t i=0 ; i<nlines; i++)
	{
		if(lines[i]!="")
		{
			vtxt.push_back(new TLatex(x,y-(i+1)*dy,"#font[132]{"+lines[i]+"}"));
			vtxt[i]->SetTextAlign(11);
			vtxt[i]->SetTextAngle(0);
			vtxt[i]->SetTextColor(kRed+2);
			vtxt[i]->Draw();
		}
	}
	
	return cs;
}
void setLines(TString* lines,
			  TString l0="",
			  TString l1="",
			  TString l2="",
			  TString l3="",
			  TString l4="",
			  TString l5="",
			  TString l6="",
			  TString l7="",
			  TString l8="",
			  TString l9="")
{
	lines[0] = l0;
	lines[1] = l1;
	lines[2] = l2;
	lines[3] = l3;
	lines[4] = l4;
	lines[5] = l5;
	lines[6] = l6;
	lines[7] = l7;
	lines[8] = l8;
	lines[9] = l9;
}
void resetLines(Int_t nlines, TString* lines)
{
	for(Int_t i=0 ; i<nlines; i++) lines[i] = "";
}
void setPage(TString pdffilename, TString pgtitle,
			 Int_t nlines, TString* lines,
			 TString l0="",
			 TString l1="",
			 TString l2="",
			 TString l3="",
			 TString l4="",
			 TString l5="",
			 TString l6="",
			 TString l7="",
			 TString l8="",
			 TString l9="")
{
	setLines(lines,l0,l1,l2,l3,l4,l5,l6,l7,l8,l9);
	separatorPage(pgtitle,nlines,lines)->SaveAs(pdffilename);
	resetLines(nlines,lines);
}

void setLegendDefaults(TLegend* l)
{
	l->SetFillStyle(4000); //will be transparent
	l->SetFillColor(0);
	l->SetTextFont(42);
	l->SetBorderSize(0);
}
double getYmin(TH1* h)
{
	if(h==NULL)
	{
		_ERR(1,"Histogram is null, getYmin(TH1* h) returning 0.001");
		return 0.01;
	}
	double min = 1.e20;
	double binval = 0.;
	for(int b=2 ; b<h->GetNbinsX() ; b++) // don't count the first and last bins
	{	
		binval = h->GetBinContent(b);
		min = (binval<min  &&  binval>0.) ? binval : min;
	}
	return min;
}

double getYmax(TH1* h)
{
	if(h==NULL)
	{
		_ERR(1,"Histogram is null, getYmax(TH1* h) returning -1");
		return -1.;
	}
	double max    = 0.;
	double binval = 0.;
	for(int b=1 ; b<=h->GetNbinsX() ; b++)
	{
		binval = h->GetBinContent(b);
		max = (binval>max) ? binval : max;
	}
	return max;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


void addHist(TMapTSP2TH1& histos, TString channel, TString name, TString titles,
			 int nbins, float xmin, float xmax, bool addchannelName=true)
{
	if(nbins<=0) _FAT("nbins="<<nbins<<" for histogram: "<<name);
	TString hname = channel+"_"+name;
	histos.insert(make_pair(hname, new TH1F(channel+"_"+name,titles,nbins,xmin,xmax)));
	TString title = histos[hname]->GetTitle();
	title = (title!="" && addchannelName) ? labels[channel]+" "+titles : titles;
	histos[hname]->SetTitle(title);
	histos[hname]->Sumw2();
}
void addHist(TMapTSP2TH2& histos, TString channel, TString name, TString titles,
			 int nbinsx, float xmin, float xmax,
			 int nbinsy, float ymin, float ymax, bool addchannelName=true)
{
	if(nbinsx<=0) _FAT("nbins="<<nbinsx<<" for histogram: "<<name);
	if(nbinsy<=0) _FAT("nbins="<<nbinsy<<" for histogram: "<<name);
	TString hname = channel+"_"+name;
	histos.insert(make_pair(hname, new TH2F(channel+"_"+name,titles,nbinsx,xmin,xmax,nbinsy,ymin,ymax)));
	TString title = histos[hname]->GetTitle();
	title = (title!="" && addchannelName) ? labels[channel]+" "+titles : titles;
	histos[hname]->SetTitle(title);
	histos[hname]->Sumw2();
}
void makeCutflowHisto(TString channel, bool addchannelName=true)
{
	vector<TString> cutnames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it)
	{
		if(it->first>=0) cutnames.push_back(it->second);
	}
	Int_t ncuts = cutnames.size();
	
	TString label = "after 2nd skim";
	
	addHist(histos1,channel,"cutflow_normalized","Normalized cutflow "+label+";;Normalized",ncuts,0,ncuts,false);
	addHist(histos1,channel,"cutflow_absolute","Absolute cutflow "+label+";;Events",ncuts,0,ncuts,addchannelName);
	addHist(histos1,channel,"cutflow_weighted","Weighted cutflow "+label+";;Events",ncuts,0,ncuts,addchannelName);
	
	for(Int_t b=1 ; b<=ncuts ; b++)
	{
		TString cutname = cutnames[b-1];
		cutname.ReplaceAll("nPassing_","");
		histos1[channel+"_cutflow_normalized"]->GetXaxis()->SetBinLabel(b,cutname);
		histos1[channel+"_cutflow_absolute"]->GetXaxis()->SetBinLabel(b,cutname);
		histos1[channel+"_cutflow_weighted"]->GetXaxis()->SetBinLabel(b,cutname);
	}
}
void makeCountersHisto(TString channel, bool addchannelName=true)
{
	vector<TString> counternames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it)
	{
		if(it->first<0) counternames.push_back(it->second);
	}
	Int_t ncounters = counternames.size();
	
	TString label = "after 2nd skim";
	
	addHist(histos1,channel,"counters_normalized","Normalized counters "+label+";;Normalized",ncounters,0,ncounters,false);
	addHist(histos1,channel,"counters_absolute","Absolute counters "+label+";;Events",ncounters,0,ncounters,addchannelName);
	addHist(histos1,channel,"counters_weighted","Weighted counters "+label+";;Events",ncounters,0,ncounters,addchannelName);
	
	for(Int_t b=1 ; b<=ncounters ; b++)
	{
		TString countername = counternames[b-1];
		countername.ReplaceAll("nPassing_","");
		histos1[channel+"_counters_normalized"]->GetXaxis()->SetBinLabel(b,countername);
		histos1[channel+"_counters_absolute"]->GetXaxis()->SetBinLabel(b,countername);
		histos1[channel+"_counters_weighted"]->GetXaxis()->SetBinLabel(b,countername);
	}
}
void makeCategoriesHisto(TString channel, bool addchannelName=false)
{
	vector<TString> catnames;
	catnames.push_back("3mu");
	catnames.push_back("2mu+1tpA");
	catnames.push_back("2mu+1tpB | 0tpA");
	catnames.push_back("2mu+1calo(|#eta|<0.1) | 0tpA");
	catnames.push_back("1mu+2tpA");
	catnames.push_back("1mu+2tpB | 0tpA");
	catnames.push_back("1mu+1tpA+1tpB");
	catnames.push_back("3tpA | 0mu");
	catnames.push_back("3tpB | 0mu | 0tpA");
	catnames.push_back("2tpA+1tpB | 0mu");
	catnames.push_back("2tpA+1calo | 0mu");
	catnames.push_back("1tpA+1tpB | 0mu");
	unsigned int ncat = catnames.size();
	TString hname = "";
	TString fullname = "";
	
	
	// hname = "tripletCategories_noVertexing";
	// fullname = channel+"_"+hname;
	// addHist(histos1,channel,hname, "3body categories no vertexing;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	// for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	// 
	// hname = "tripletCategories_norm_noVertexing";
	// fullname = channel+"_"+hname;
	// addHist(histos1,channel,hname, "3body categories no vertexing;;Normalized", ncat,0,ncat, addchannelName);
	// for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_afterVertexing";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after vertexing;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_afterVertexing";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after vertexing;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);

	
	hname = "tripletCategories_after_muons";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after muon cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_muons";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after muon cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_after_triplet";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after triplet cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_triplet";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after triplet cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_after_hadclean";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after hadronic cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_hadclean";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after hadronic cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_after_met";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #it{E}_{T}^{miss} cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_met";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #it{E}_{T}^{miss} cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);


	hname = "tripletCategories_after_ht";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #it{H}_{T} cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);

	hname = "tripletCategories_norm_after_ht";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #it{H}_{T} cuts;;Normalized", ncat,0,ncat, addchannelName);
	
	
	hname = "tripletCategories_after_hadcleanTight";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after tight hadronic cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_hadcleanTight";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after tight hadronic cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_after_tripletTight";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after tight triplet cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_tripletTight";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after tight triplet cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);


	hname = "tripletCategories_after_isolation";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after isolation cuts;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_isolation";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after isolation cuts;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	
	hname = "tripletCategories_after_lowmassres";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #rho/#omega,#phi rejection;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_after_lowmassres";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories after #rho/#omega,#phi rejection;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
		
	
	hname = "tripletCategories_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories end of selection;;Normalized to 1st bin", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
	
	hname = "tripletCategories_norm_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "3body categories end of selection;;Normalized", ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos1[fullname]->GetXaxis()->SetBinLabel(i,catnames[i-1]);
}
void makeTriggerHisto(TString channel)
{
	unsigned int ntrig = triggers.size();
	TString hname = "";
	TString fullname = "";
	
	hname = "triggers_absoluteEvents_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Trigger absolute contribution after vertexing;;Events", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_absolute_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Fractional trigger absolute contribution after vertexing;;Normalized", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_unique_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Fractional trigger unique contribution after vertexing;;Normalized", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}

	hname = "triggers_absoluteEvents_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Trigger absolute contribution end of selection;;Events", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_absolute_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Fractional trigger absolute contribution end of selection;;Normalized", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_unique_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos1,channel,hname, "Fractional trigger unique contribution end of selection;;Normalized", ntrig,0,ntrig,false);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos1[fullname]->GetXaxis()->SetBinLabel(i,label);
	}
	
	
	
	hname = "triggers_correlation_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation after vertexing;;;Events", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_correlation_norm_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation after vertexing;;;Normalized", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_correlation_normAB_after_vertexing";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation after vertexing;;;Correlation", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_",""); label.ReplaceAll("EF",""); label.ReplaceAll("_tclcw",""); label.ReplaceAll("tight_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
	
	hname = "triggers_correlation_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation end of selection;;;Events", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_correlation_norm_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation end of selection;;;Normalized", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
	hname = "triggers_correlation_normAB_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos2,channel,hname, "Trigger correlation end of selection;;;Correlation", ntrig,0,ntrig, ntrig,0,ntrig);
	for(unsigned int i=1 ; i<=ntrig ; ++i)
	{
		TString label = triggers[i-1]; label.ReplaceAll("EF_",""); label.ReplaceAll("EF",""); label.ReplaceAll("_tclcw",""); label.ReplaceAll("tight_","");
		histos2[fullname]->GetXaxis()->SetBinLabel(i,label);
		histos2[fullname]->GetYaxis()->SetBinLabel(i,label);
	}
}


void fillCutFlowHisto(TString channel, TMapTSP2TH1& histos)
{	
	vector<TString> cutnames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it) { if(it->first>=0) cutnames.push_back(it->second); }
	Int_t ncuts = cutnames.size();
	for(Int_t b=1 ; b<=ncuts ; b++)
	{
		unsigned int counter  = getCounter(cutnames[b-1]);
		unsigned int wcounter = getCounter(cutnames[b-1],true);
		counter  = (std::isinf(counter)  || std::isnan(counter)  || counter<1)  ? 0/*1.e-3*/ : counter;
		wcounter = (std::isinf(wcounter) || std::isnan(wcounter) || wcounter<1) ? 0/*1.e-3*/ : wcounter;
		Double_t bincontent = histos[channel+"_cutflow_absolute"]->GetBinContent(b); // that is nonzero only for channel="Data"
		histos[channel+"_cutflow_normalized"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_cutflow_absolute"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_cutflow_weighted"]->SetBinContent(b,wcounter+bincontent);
	}
}
void fillCountersHisto(TString channel, TMapTSP2TH1& histos)
{	
	vector<TString> counternames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it) { if(it->first<0) counternames.push_back(it->second); }
	Int_t ncounters = counternames.size();
	for(Int_t b=1 ; b<=ncounters ; b++)
	{
		unsigned int counter  = getCounter(counternames[b-1]);
		unsigned int wcounter = getCounter(counternames[b-1],true);
		counter  = (std::isinf(counter)  || std::isnan(counter)  || counter<1)  ? 0/*1.e-3*/ : counter;
		wcounter = (std::isinf(wcounter) || std::isnan(wcounter) || wcounter<1) ? 0/*1.e-3*/ : wcounter;
		Double_t bincontent = histos[channel+"_counters_absolute"]->GetBinContent(b); // that is nonzero only for channel="Data"
		histos[channel+"_counters_normalized"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_counters_absolute"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_counters_weighted"]->SetBinContent(b,wcounter+bincontent);
	}
}

void fillCategories(int category, TString channel, TString hname)
{
	TString histname = channel+"_"+hname;
	int bin = (category<MUID3) ? category+1 : category+1-MUONS0TPA1TPB2;
	histos1[histname]->AddBinContent(bin);
}

void fillTriggerHistos(TString channel, TString suffix, TMapTSP2vi& vi, float wgt)
{
	if(0) cout << "vi.size()=" << vi.size() << endl;
	if(0) cout << "wgt="       << wgt       << endl;
	
	unsigned int ntriggers = triggers.size();
	unsigned int binx = 0;
	unsigned int biny = 0;
	TString histname  = "";
	TString histname2 = "";
	TString histname3 = "";
	TString histname4 = "";
	
	histname = channel+"_triggers_absolute_"+suffix;
	for(unsigned int t=0 ; t<ntriggers ; ++t)
	{
		if(ints[triggers[t]]==1) histos1[histname]->AddBinContent(t+1);
	}
	
	histname = channel+"_triggers_absoluteEvents_"+suffix;
	for(unsigned int t=0 ; t<ntriggers ; ++t)
	{
		if(ints[triggers[t]]==1) histos1[histname]->AddBinContent(t+1);
	}
	
	//triggers_unique_after_vertexing
	histname = channel+"_triggers_unique_"+suffix;
	histname2 = channel+"_triggers_correlation_"+suffix;
	histname3 = channel+"_triggers_correlation_norm_"+suffix;
	histname4 = channel+"_triggers_correlation_normAB_"+suffix;
	for(unsigned int t=0 ; t<ntriggers ; ++t)
	{
		if(ints[triggers[t]]==1)
		{
			bool unique  = true;
			for(unsigned int tt=0 ; tt<ntriggers ; ++tt)
			{
				// if(tt==t) continue; // dont compare with itself
				if(ints[triggers[tt]]==1)
				{
					if(tt!=t) unique = false; // triggers[t] is not unique
					Int_t binxy = histos2[histname2]->GetBin(t+1,tt+1);
					Float_t z = histos2[histname2]->GetBinContent(binxy);
					histos2[histname2]->SetBinContent(binxy,z+1);
					histos2[histname3]->SetBinContent(binxy,z+1);
					histos2[histname4]->SetBinContent(binxy,z+1);
				}
			}
			if(unique) histos1[histname]->AddBinContent(t+1);
		}
	}
}

void fillHitHistos(unsigned int vtx, TString name, TMapTSP2vi& vi, float wgt)
{	
	float htTRTfrac1 = (vi["mu_nTRThits1"]->at(vtx)>0) ? (float)vi["mu_htTRThits1"]->at(vtx)/(float)vi["mu_nTRThits1"]->at(vtx) : 0.;
	float htTRTfrac2 = (vi["mu_nTRThits2"]->at(vtx)>0) ? (float)vi["mu_htTRThits2"]->at(vtx)/(float)vi["mu_nTRThits2"]->at(vtx) : 0.;
	float htTRTfrac3 = (vi["mu_nTRThits3"]->at(vtx)>0) ? (float)vi["mu_htTRThits3"]->at(vtx)/(float)vi["mu_nTRThits3"]->at(vtx) : 0.;
	int ltTRThits1 = vi["mu_nTRThits1"]->at(vtx) - vi["mu_htTRThits1"]->at(vtx);
	int ltTRThits2 = vi["mu_nTRThits2"]->at(vtx) - vi["mu_htTRThits2"]->at(vtx);
	int ltTRThits3 = vi["mu_nTRThits3"]->at(vtx) - vi["mu_htTRThits3"]->at(vtx);
	float htTRTrat1 = (ltTRThits1>0) ? (float)vi["mu_htTRThits1"]->at(vtx)/(float)ltTRThits1 : 0.;
	float htTRTrat2 = (ltTRThits2>0) ? (float)vi["mu_htTRThits2"]->at(vtx)/(float)ltTRThits2 : 0.;
	float htTRTrat3 = (ltTRThits3>0) ? (float)vi["mu_htTRThits3"]->at(vtx)/(float)ltTRThits3 : 0.;
	
	if     (vi["mu_type1"]->at(vtx)==MUON) { histos1[name+"_mu_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits1"]->at(vtx),wgt);   histos1[name+"_mu_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac1,wgt);    histos1[name+"_mu_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat1,wgt);     }
	else if(vi["mu_type1"]->at(vtx)==TPA)  { histos1[name+"_TPa_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits1"]->at(vtx),wgt);  histos1[name+"_TPa_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac1,wgt);   histos1[name+"_TPa_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat1,wgt);    }
	else if(vi["mu_type1"]->at(vtx)==TPB)  { histos1[name+"_TPb_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits1"]->at(vtx),wgt);  histos1[name+"_TPb_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac1,wgt);   histos1[name+"_TPb_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat1,wgt);    }
	if     (vi["mu_type2"]->at(vtx)==MUON) { histos1[name+"_mu_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits2"]->at(vtx),wgt);   histos1[name+"_mu_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac2,wgt);    histos1[name+"_mu_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat2,wgt);     }
	else if(vi["mu_type2"]->at(vtx)==TPA)  { histos1[name+"_TPa_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits2"]->at(vtx),wgt);  histos1[name+"_TPa_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac2,wgt);   histos1[name+"_TPa_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat2,wgt);    }
	else if(vi["mu_type2"]->at(vtx)==TPB)  { histos1[name+"_TPb_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits2"]->at(vtx),wgt);  histos1[name+"_TPb_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac2,wgt);   histos1[name+"_TPb_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat2,wgt);    }
	if     (vi["mu_type3"]->at(vtx)==MUON) { histos1[name+"_mu_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits3"]->at(vtx),wgt);   histos1[name+"_mu_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac3,wgt);    histos1[name+"_mu_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat3,wgt);     }
	else if(vi["mu_type3"]->at(vtx)==TPA)  { histos1[name+"_TPa_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits3"]->at(vtx),wgt);  histos1[name+"_TPa_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac3,wgt);   histos1[name+"_TPa_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat3,wgt);    }
	else if(vi["mu_type3"]->at(vtx)==TPB)  { histos1[name+"_TPb_nHighThresholdTRTHits"]->Fill(vi["mu_htTRThits3"]->at(vtx),wgt);  histos1[name+"_TPb_nHighThresholdTRTHitsFraction"]->Fill(htTRTfrac3,wgt);   histos1[name+"_TPb_nHighThresholdTRTHitsRatio"]->Fill(htTRTrat3,wgt);    }
	
	
	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_outliersOnTrack1"]->Fill(vi["mu_nOutliersOnTrack1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_outliersOnTrack2"]->Fill(vi["mu_nOutliersOnTrack2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_outliersOnTrack3"]->Fill(vi["mu_nOutliersOnTrack3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_outliersOnTrack1"]->Fill(vi["mu_nOutliersOnTrack1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_outliersOnTrack2"]->Fill(vi["mu_nOutliersOnTrack2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_outliersOnTrack3"]->Fill(vi["mu_nOutliersOnTrack3"]->at(vtx),wgt);
	
	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_MDTHits1"]->Fill(vi["mu_nMDThits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_MDTHits2"]->Fill(vi["mu_nMDThits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_MDTHits3"]->Fill(vi["mu_nMDThits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_MDTHits1"]->Fill(vi["mu_nMDThits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_MDTHits2"]->Fill(vi["mu_nMDThits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_MDTHits3"]->Fill(vi["mu_nMDThits3"]->at(vtx),wgt);
	
	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_RPCHits1"]->Fill(vi["mu_nRPCPhiHits1"]->at(vtx)+vi["mu_nRPCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_RPCHits2"]->Fill(vi["mu_nRPCPhiHits2"]->at(vtx)+vi["mu_nRPCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_RPCHits3"]->Fill(vi["mu_nRPCPhiHits3"]->at(vtx)+vi["mu_nRPCEtaHits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_RPCHits1"]->Fill(vi["mu_nRPCPhiHits1"]->at(vtx)+vi["mu_nRPCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_RPCHits2"]->Fill(vi["mu_nRPCPhiHits2"]->at(vtx)+vi["mu_nRPCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_RPCHits3"]->Fill(vi["mu_nRPCPhiHits3"]->at(vtx)+vi["mu_nRPCEtaHits3"]->at(vtx),wgt);

	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_TGCHits1"]->Fill(vi["mu_nTGCPhiHits1"]->at(vtx)+vi["mu_nTGCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_TGCHits2"]->Fill(vi["mu_nTGCPhiHits2"]->at(vtx)+vi["mu_nTGCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_TGCHits3"]->Fill(vi["mu_nTGCPhiHits3"]->at(vtx)+vi["mu_nTGCEtaHits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_TGCHits1"]->Fill(vi["mu_nTGCPhiHits1"]->at(vtx)+vi["mu_nTGCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_TGCHits2"]->Fill(vi["mu_nTGCPhiHits2"]->at(vtx)+vi["mu_nTGCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_TGCHits3"]->Fill(vi["mu_nTGCPhiHits3"]->at(vtx)+vi["mu_nTGCEtaHits3"]->at(vtx),wgt);

	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_precisionHits1"]->Fill(vi["mu_nPrecisionHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_precisionHits2"]->Fill(vi["mu_nPrecisionHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_precisionHits3"]->Fill(vi["mu_nPrecisionHits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHits1"]->Fill(vi["mu_nPrecisionHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHits2"]->Fill(vi["mu_nPrecisionHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHits3"]->Fill(vi["mu_nPrecisionHits3"]->at(vtx),wgt);

	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_precisionHoles1"]->Fill(vi["mu_nPrecisionHoles1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_precisionHoles2"]->Fill(vi["mu_nPrecisionHoles2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_precisionHoles3"]->Fill(vi["mu_nPrecisionHoles3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHoles1"]->Fill(vi["mu_nPrecisionHoles1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHoles2"]->Fill(vi["mu_nPrecisionHoles2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_precisionHoles3"]->Fill(vi["mu_nPrecisionHoles3"]->at(vtx),wgt);
    
	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_precisionOutliers1"]->Fill(vi["mu_nPrecisionOutliers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_precisionOutliers2"]->Fill(vi["mu_nPrecisionOutliers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_precisionOutliers3"]->Fill(vi["mu_nPrecisionOutliers3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_precisionOutliers1"]->Fill(vi["mu_nPrecisionOutliers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_precisionOutliers2"]->Fill(vi["mu_nPrecisionOutliers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_precisionOutliers3"]->Fill(vi["mu_nPrecisionOutliers3"]->at(vtx),wgt);
    
	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_phiLayers1"]->Fill(vi["mu_nPhiLayers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_phiLayers2"]->Fill(vi["mu_nPhiLayers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_phiLayers3"]->Fill(vi["mu_nPhiLayers3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_phiLayers1"]->Fill(vi["mu_nPhiLayers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_phiLayers2"]->Fill(vi["mu_nPhiLayers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_phiLayers3"]->Fill(vi["mu_nPhiLayers3"]->at(vtx),wgt);

	if(vi["mu_type1"]->at(vtx)==MUON) histos1[name+"_mu_etaphiLayers1"]->Fill(vi["mu_nEtaPhiLayers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos1[name+"_mu_etaphiLayers2"]->Fill(vi["mu_nEtaPhiLayers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos1[name+"_mu_etaphiLayers3"]->Fill(vi["mu_nEtaPhiLayers3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos1[name+"_TPa_etaphiLayers1"]->Fill(vi["mu_nEtaPhiLayers1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos1[name+"_TPa_etaphiLayers2"]->Fill(vi["mu_nEtaPhiLayers2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos1[name+"_TPa_etaphiLayers3"]->Fill(vi["mu_nEtaPhiLayers3"]->at(vtx),wgt);
	
	if(vi["mu_type1"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_phiHits1"]->Fill(vi["mu_nTGCPhiHits1"]->at(vtx),vi["mu_nRPCPhiHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_phiHits2"]->Fill(vi["mu_nTGCPhiHits2"]->at(vtx),vi["mu_nRPCPhiHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_phiHits3"]->Fill(vi["mu_nTGCPhiHits3"]->at(vtx),vi["mu_nRPCPhiHits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_phiHits1"]->Fill(vi["mu_nTGCPhiHits1"]->at(vtx),vi["mu_nRPCPhiHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_phiHits2"]->Fill(vi["mu_nTGCPhiHits2"]->at(vtx),vi["mu_nRPCPhiHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_phiHits3"]->Fill(vi["mu_nTGCPhiHits3"]->at(vtx),vi["mu_nRPCPhiHits3"]->at(vtx),wgt);
    
	if(vi["mu_type1"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_etaHits1"]->Fill(vi["mu_nTGCEtaHits1"]->at(vtx),vi["mu_nRPCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_etaHits2"]->Fill(vi["mu_nTGCEtaHits2"]->at(vtx),vi["mu_nRPCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_etaHits3"]->Fill(vi["mu_nTGCEtaHits3"]->at(vtx),vi["mu_nRPCEtaHits3"]->at(vtx),wgt);
	if(vi["mu_type1"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_etaHits1"]->Fill(vi["mu_nTGCEtaHits1"]->at(vtx),vi["mu_nRPCEtaHits1"]->at(vtx),wgt);
	if(vi["mu_type2"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_etaHits2"]->Fill(vi["mu_nTGCEtaHits2"]->at(vtx),vi["mu_nRPCEtaHits2"]->at(vtx),wgt);
	if(vi["mu_type3"]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_etaHits3"]->Fill(vi["mu_nTGCEtaHits3"]->at(vtx),vi["mu_nRPCEtaHits3"]->at(vtx),wgt);
}

void fillMETHistos(unsigned int vtx, TString name, TString suffix, TMapTSf& f, TMapTSP2vf& vf, TString metType, float wgt)
{	
	histos1[name+"_met_et"+suffix]->Fill(f["met_"+metType+"_et"],wgt);
	histos1[name+"_met_mt_et3mu"+suffix]->Fill(vf["met_"+metType+"_mT"]->at(vtx),wgt);
	histos1[name+"_met_dphi3mu"+suffix]->Fill(vf["met_"+metType+"_dPhi3mu"]->at(vtx),wgt);
	
	histos1[name+"_mettrk_et"+suffix]->Fill(f["met_track_et"],wgt);
	histos1[name+"_mettrk_mt_et3mu"+suffix]->Fill(vf["met_track_mT"]->at(vtx),wgt);
	histos1[name+"_mettrk_dphi3mu"+suffix]->Fill(vf["met_track_dPhi3mu"]->at(vtx),wgt);
	
	histos1[name+"_mettrk_dphimet"+suffix]->Fill(fabs(TVector2::Phi_mpi_pi(f["met_"+metType+"_phi"]-f["met_track_phi"])),wgt);
}

void fill3BodyHistos(unsigned int vtx, TString name, TString suffix, TMapTSP2vf& vf, float wgt)
{	
	histos1[name+"_m3body"+suffix]->Fill(vf["vtx_mass"]->at(vtx),wgt);
	histos1[name+"_m3body_lin"+suffix]->Fill(vf["vtx_mass"]->at(vtx),wgt);
	histos1[name+"_m3body_lin_zoom"+suffix]->Fill(vf["vtx_mass"]->at(vtx),wgt);
	histos1[name+"_m3body_sigregion"+suffix]->Fill(vf["vtx_mass"]->at(vtx),wgt);
	histos1[name+"_pT3body"+suffix]->Fill(vf["vtx_pt"]->at(vtx),wgt);
	if(fabs(vf["vtx_charge"]->at(vtx))==1)
	{
		histos1[name+"_mOS"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),wgt/2.);
		histos1[name+"_mOS"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx),wgt/2.);
		histos1[name+"_mOS1"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),wgt);
		histos1[name+"_mOS2"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx),wgt);
		histos1[name+"_mSS"+suffix]->Fill(vf["vtx_mSS"]->at(vtx),wgt);

		histos2[name+"_mOS2_vs_mOS1"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),vf["vtx_mOS2"]->at(vtx),wgt);
		histos2[name+"_mSS_vs_mOS"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),vf["vtx_mSS"]->at(vtx),wgt/2.);
		histos2[name+"_mSS_vs_mOS"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx),vf["vtx_mSS"]->at(vtx),wgt/2.);		
		histos2[name+"_mSSsq_vs_mOS1sq"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx)*vf["vtx_mOS1"]->at(vtx),vf["vtx_mSS"]->at(vtx)*vf["vtx_mSS"]->at(vtx),wgt);
		histos2[name+"_mSSsq_vs_mOS2sq"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx)*vf["vtx_mOS2"]->at(vtx),vf["vtx_mSS"]->at(vtx)*vf["vtx_mSS"]->at(vtx),wgt);
		histos2[name+"_mSSsq_vs_mOSsq"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx)*vf["vtx_mOS1"]->at(vtx),vf["vtx_mSS"]->at(vtx)*vf["vtx_mSS"]->at(vtx),wgt/2.);
		histos2[name+"_mSSsq_vs_mOSsq"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx)*vf["vtx_mOS2"]->at(vtx),vf["vtx_mSS"]->at(vtx)*vf["vtx_mSS"]->at(vtx),wgt/2.);
		histos2[name+"_m3body_vs_mOS"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),vf["vtx_mass"]->at(vtx),wgt/2.);
		histos2[name+"_m3body_vs_mOS"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx),vf["vtx_mass"]->at(vtx),wgt/2.);
		histos2[name+"_pT3body_vs_mOS"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),vf["vtx_pt"]->at(vtx),wgt/2.);
		histos2[name+"_pT3body_vs_mOS"+suffix]->Fill(vf["vtx_mOS2"]->at(vtx),vf["vtx_pt"]->at(vtx),wgt/2.);
		
		TLorentzVector p1, p2, p3, p, p123, p12, p13, p23;
		p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
		p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
		p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
		p123 = p1+p2+p3;
		p12 = p1+p2;
		p13 = p1+p3;
		p23 = p2+p3;
		float q1 = (vf["mu_trkqoverp1"]->at(vtx)<0) ? -1. : +1;
		float q2 = (vf["mu_trkqoverp2"]->at(vtx)<0) ? -1. : +1;
		float q3 = (vf["mu_trkqoverp3"]->at(vtx)<0) ? -1. : +1;
		float q12 = fabs(q1+q2);
		float q13 = fabs(q1+q3);
		float q23 = fabs(q2+q3);
		float dr12 = p1.DeltaR(p2); //deltaR(vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),f["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx));
		float dr13 = p1.DeltaR(p3); //deltaR(vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),f["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx));
		float dr23 = p2.DeltaR(p3); //deltaR(vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),f["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx));
		if(q12==0.)
		{
			float mOS = p12.M();
			p.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
			TLorentzVector pBoosted = p;
			TLorentzVector pOSBoosted = p12;
			pBoosted.Boost(-1.*p123.BoostVector());
			pOSBoosted.Boost(-1.*p123.BoostVector());
			float angle = pBoosted.Vect().Angle(pOSBoosted.Vect());
			histos2[name+"_3rdTrk_pt_vs_mOS"+suffix]->Fill(mOS,vf["mu_pt3"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_chi2dof_vs_mOS"+suffix]->Fill(mOS,vf["mu_chi2ndftrkfit3"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_pval_vs_mOS"+suffix]->Fill(mOS,vf["mu_pvaltrkfit3"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_angle_vs_mOS"+suffix]->Fill(mOS,angle,wgt);
			histos2[name+"_3rdTrk_PtOverPtOS_vs_mOS"+suffix]->Fill(mOS,p.Pt()/p12.Pt(),wgt);
			
			histos1[name+"_dROS1"+suffix]->Fill(dr12,wgt);
			histos2[name+"_dROS1_vs_mOS1"+suffix]->Fill(vf["vtx_mOS1"]->at(vtx),dr12,wgt);
			if(q13==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr13,wgt); histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p13.M(),dr13,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr23,wgt);  histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p23.M(),dr23,wgt);
			}
			if(q23==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr23,wgt); histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p23.M(),dr23,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr13,wgt);  histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p13.M(),dr13,wgt);
			}
		}
		if(q13==0.)
		{
			float mOS = p13.M();
			p.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
			TLorentzVector pBoosted = p;
			TLorentzVector pOSBoosted = p13;
			pBoosted.Boost(-1.*p123.BoostVector());
			pOSBoosted.Boost(-1.*p123.BoostVector());
			float angle = pBoosted.Vect().Angle(pOSBoosted.Vect());
			histos2[name+"_3rdTrk_pt_vs_mOS"+suffix]->Fill(mOS,vf["mu_pt2"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_chi2dof_vs_mOS"+suffix]->Fill(mOS,vf["mu_chi2ndftrkfit2"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_pval_vs_mOS"+suffix]->Fill(mOS,vf["mu_pvaltrkfit2"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_angle_vs_mOS"+suffix]->Fill(mOS,angle,wgt);
			histos2[name+"_3rdTrk_PtOverPtOS_vs_mOS"+suffix]->Fill(mOS,p.Pt()/p13.Pt(),wgt);
			
			histos1[name+"_dROS1"+suffix]->Fill(dr13,wgt);
			if(q12==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr12,wgt);  histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p12.M(),dr12,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr23,wgt);   histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p23.M(),dr23,wgt);
			}
			if(q23==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr23,wgt); histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p23.M(),dr23,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr12,wgt);  histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p12.M(),dr12,wgt);
			}
		}
		if(q23==0.)
		{
			float mOS = p23.M();
			p.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
			TLorentzVector pBoosted = p;
			TLorentzVector pOSBoosted = p23;
			pBoosted.Boost(-1.*p123.BoostVector());
			pOSBoosted.Boost(-1.*p123.BoostVector());
			float angle = pBoosted.Vect().Angle(pOSBoosted.Vect());
			histos2[name+"_3rdTrk_pt_vs_mOS"+suffix]->Fill(mOS,vf["mu_pt1"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_chi2dof_vs_mOS"+suffix]->Fill(mOS,vf["mu_chi2ndftrkfit1"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_pval_vs_mOS"+suffix]->Fill(mOS,vf["mu_pvaltrkfit1"]->at(vtx),wgt);
			histos2[name+"_3rdTrk_angle_vs_mOS"+suffix]->Fill(mOS,angle,wgt);
			histos2[name+"_3rdTrk_PtOverPtOS_vs_mOS"+suffix]->Fill(mOS,p.Pt()/p23.Pt(),wgt);
			
			histos1[name+"_dROS1"+suffix]->Fill(dr23,wgt);
			if(q12==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr12,wgt);  histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p12.M(),dr12,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr13,wgt);   histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p13.M(),dr13,wgt);
			}
			if(q13==0.)
			{
				histos1[name+"_dROS2"+suffix]->Fill(dr13,wgt);  histos2[name+"_dROS2_vs_mOS2"+suffix]->Fill(p13.M(),dr13,wgt);
				histos1[name+"_dRSS"+suffix]->Fill(dr12,wgt);   histos2[name+"_dRSS_vs_mSS"+suffix]->Fill(p12.M(),dr12,wgt);
			}
		}
	}
}

void fillIsolationHistos(unsigned int vtx, TString name, TMapTSP2vf& vf, float wgt)
{
	histos1[name+"_dRmin"]->Fill(vf["vtx_dRmin"]->at(vtx),wgt);
	histos1[name+"_dRmax"]->Fill(vf["vtx_dRmax"]->at(vtx),wgt);
	
	histos1[name+"_isolation000"]->Fill(vf["vtx_isolation000"]->at(vtx),wgt);
	histos1[name+"_isolation001"]->Fill(vf["vtx_isolation001"]->at(vtx),wgt);
	histos1[name+"_isolation002"]->Fill(vf["vtx_isolation002"]->at(vtx),wgt);
	histos1[name+"_isolation003"]->Fill(vf["vtx_isolation003"]->at(vtx),wgt);
	histos1[name+"_isolation004"]->Fill(vf["vtx_isolation004"]->at(vtx),wgt);
	histos1[name+"_isolation005"]->Fill(vf["vtx_isolation005"]->at(vtx),wgt);
	histos1[name+"_isolation006"]->Fill(vf["vtx_isolation006"]->at(vtx),wgt);
	histos1[name+"_isolation007"]->Fill(vf["vtx_isolation007"]->at(vtx),wgt);
	histos1[name+"_isolation008"]->Fill(vf["vtx_isolation008"]->at(vtx),wgt);
	histos1[name+"_isolation009"]->Fill(vf["vtx_isolation009"]->at(vtx),wgt);
	histos1[name+"_isolation010"]->Fill(vf["vtx_isolation010"]->at(vtx),wgt);
	histos1[name+"_isolation012"]->Fill(vf["vtx_isolation012"]->at(vtx),wgt);
	histos1[name+"_isolation014"]->Fill(vf["vtx_isolation014"]->at(vtx),wgt);
	histos1[name+"_isolation016"]->Fill(vf["vtx_isolation016"]->at(vtx),wgt);
	histos1[name+"_isolation018"]->Fill(vf["vtx_isolation018"]->at(vtx),wgt);
	histos1[name+"_isolation020"]->Fill(vf["vtx_isolation020"]->at(vtx),wgt);
	histos1[name+"_isolation022"]->Fill(vf["vtx_isolation022"]->at(vtx),wgt);
	histos1[name+"_isolation024"]->Fill(vf["vtx_isolation024"]->at(vtx),wgt);
	histos1[name+"_isolation026"]->Fill(vf["vtx_isolation026"]->at(vtx),wgt);
	histos1[name+"_isolation028"]->Fill(vf["vtx_isolation028"]->at(vtx),wgt);
	histos1[name+"_isolation030"]->Fill(vf["vtx_isolation030"]->at(vtx),wgt);
}

void fillVertexHistos(unsigned int vtx, TString name, TMapTSP2vf& vf, float wgt)
{
	histos1[name+"_chi2ndf"]->Fill(vf["vtx_chi2ndf"]->at(vtx),wgt);
	histos1[name+"_pvalue"]->Fill(vf["vtx_pval"]->at(vtx),wgt);
	histos1[name+"_Lxy"]->Fill(vf["vtx_lxy"]->at(vtx),wgt);
	histos1[name+"_dLxy"]->Fill(vf["vtx_lxyErr"]->at(vtx),wgt);
	histos1[name+"_a0"]->Fill(vf["vtx_a0"]->at(vtx),wgt);
	histos1[name+"_a0xy"]->Fill(vf["vtx_a0xy"]->at(vtx),wgt);
	histos1[name+"_da0xy"]->Fill(vf["vtx_a0xyErr"]->at(vtx),wgt);
	histos1[name+"_cosT"]->Fill(vf["vtx_cosT"]->at(vtx),wgt);
	histos1[name+"_cosTxy"]->Fill(vf["vtx_cosTxy"]->at(vtx),wgt);
}
void fillPVHistos(unsigned int vtx, TString name, TString suffix, TMapTSP2vi& vi, float wgt)
{
	histos1[name+"_PVNtrk"+suffix]->Fill(vi["vtx_pvNtrk"]->at(vtx),wgt);
}

void fillJetCalibrationHistos(unsigned int vtx, TString name, TString suffix, TMapTSf& f, TMapTSP2vf& vf, TString metType, float wgt)
{	
	int njets = getNjetsAll(vtx,vf);
	histos1[name+"_jet_n_calib"+suffix]->Fill(njets,wgt);
	
	int nbjets = getNbjetsAll(vtx,vf,0.3511);
	histos1[name+"_jet_b_calib"+suffix]->Fill(nbjets,wgt);
	
	if(vf["jet_pt1"]->at(vtx)>0.) histos1[name+"_jet_pt_calib"+suffix]->Fill(vf["jet_pt1"]->at(vtx),wgt);
	if(vf["jet_pt2"]->at(vtx)>0.) histos1[name+"_jet_pt_calib"+suffix]->Fill(vf["jet_pt2"]->at(vtx),wgt);
	if(vf["jet_pt3"]->at(vtx)>0.) histos1[name+"_jet_pt_calib"+suffix]->Fill(vf["jet_pt3"]->at(vtx),wgt);
	if(vf["jet_pt4"]->at(vtx)>0.) histos1[name+"_jet_pt_calib"+suffix]->Fill(vf["jet_pt4"]->at(vtx),wgt);
	
	histos1[name+"_jet_pt1_calib"+suffix]->Fill(vf["jet_pt1"]->at(vtx),wgt);
	histos1[name+"_jet_pt2_calib"+suffix]->Fill(vf["jet_pt2"]->at(vtx),wgt);
	histos1[name+"_jet_ptJ1J2_calib"+suffix]->Fill(vf["jet_pt1"]->at(vtx)+vf["jet_pt2"]->at(vtx),wgt);
	
	if(vf["jet_E1"]->at(vtx)>0.) histos1[name+"_jet_E_calib"+suffix]->Fill(vf["jet_E1"]->at(vtx),wgt);
	if(vf["jet_E2"]->at(vtx)>0.) histos1[name+"_jet_E_calib"+suffix]->Fill(vf["jet_E2"]->at(vtx),wgt);
	if(vf["jet_E3"]->at(vtx)>0.) histos1[name+"_jet_E_calib"+suffix]->Fill(vf["jet_E3"]->at(vtx),wgt);
	if(vf["jet_E4"]->at(vtx)>0.) histos1[name+"_jet_E_calib"+suffix]->Fill(vf["jet_E4"]->at(vtx),wgt);
	
	histos1[name+"_jet_E1_calib"+suffix]->Fill(vf["jet_E1"]->at(vtx),wgt);
	histos1[name+"_jet_E2_calib"+suffix]->Fill(vf["jet_E2"]->at(vtx),wgt);
	histos1[name+"_jet_EJ1J2_calib"+suffix]->Fill(vf["jet_E1"]->at(vtx)+vf["jet_E2"]->at(vtx),wgt);
	
	TLorentzVector p3body = getP3body(vtx,vf);
	TLorentzVector pSum = getPsum(vtx,vf);
	float dphiHTMET = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-f["met_"+metType+"_phi"]));
	float dphiHTMETtrk = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-f["met_track_phi"]));
	float dRHT3body = fabs(pSum.DeltaR(p3body));
	float mtHT      = sqrt(2*pSum.Pt()*f["met_"+metType+"_et"]*(1-TMath::Cos(pSum.Phi()-f["met_"+metType+"_phi"])));
	float mtHTtrk   = sqrt(2*pSum.Pt()*f["met_track_et"]*(1-TMath::Cos(pSum.Phi()-f["met_track_phi"])));
	float METsig    = f["met_"+metType+"_et"]/sqrt(pSum.Pt());
	histos1[name+"_ht_pt_calib"+suffix]->Fill(pSum.Pt(),wgt);
	histos1[name+"_ht_mT_calib"+suffix]->Fill(mtHT,wgt);
	histos1[name+"_ht_mT_track_calib"+suffix]->Fill(mtHTtrk,wgt);
	histos1[name+"_ht_metsig_calib"+suffix]->Fill(METsig,wgt);
	histos1[name+"_ht_dphimet_calib"+suffix]->Fill(dphiHTMET,wgt);
	histos1[name+"_ht_dphimet_track_calib"+suffix]->Fill(dphiHTMETtrk,wgt);
	histos1[name+"_ht_dR3body_calib"+suffix]->Fill(dRHT3body,wgt);
	
	histos1[name+"_dPhi3bodyJ1_calib"+suffix]->Fill(vf["jet_dphi3muJ1"]->at(vtx),wgt);
	histos1[name+"_dPhiJ1J2_calib"+suffix]->Fill(vf["jet_dphiJ1J2"]->at(vtx),wgt);
	
	histos1[name+"_met_et_calib"+suffix]->Fill(f["met_"+metType+"_et"],wgt);
	histos1[name+"_met_mt_et3mu_calib"+suffix]->Fill(vf["met_"+metType+"_mT"]->at(vtx),wgt);
	histos1[name+"_met_dphi3mu_calib"+suffix]->Fill(vf["met_"+metType+"_dPhi3mu"]->at(vtx),wgt);
	histos1[name+"_mettrk_dphimet_calib"+suffix]->Fill(fabs(TVector2::Phi_mpi_pi(f["met_"+metType+"_phi"]-f["met_track_phi"])),wgt);
}


void fillJetHistos(unsigned int vtx, TString name, TString suffix, TMapTSP2vf& vf, TMapTSf& f, float minptjet, TString metType, float wgt)
{
	float ht3bodyJets = vf["vtx_pt"]->at(vtx);
	float htJets = 0;
	if(vf["jet_pt1"]->at(vtx)>0. && vf["jet_pt1"]->at(vtx)>minptjet) { htJets += vf["jet_pt1"]->at(vtx); ht3bodyJets += vf["jet_pt1"]->at(vtx); }
	if(vf["jet_pt2"]->at(vtx)>0. && vf["jet_pt2"]->at(vtx)>minptjet) { htJets += vf["jet_pt2"]->at(vtx); ht3bodyJets += vf["jet_pt2"]->at(vtx); }
	if(vf["jet_pt3"]->at(vtx)>0. && vf["jet_pt3"]->at(vtx)>minptjet) { htJets += vf["jet_pt3"]->at(vtx); ht3bodyJets += vf["jet_pt3"]->at(vtx); }
	if(vf["jet_pt4"]->at(vtx)>0. && vf["jet_pt4"]->at(vtx)>minptjet) { htJets += vf["jet_pt4"]->at(vtx); ht3bodyJets += vf["jet_pt4"]->at(vtx); }
	histos1[name+"_jet_3body_ht"+suffix]->Fill(ht3bodyJets,wgt);
	histos1[name+"_jet_ht"+suffix]->Fill(htJets,wgt);
	
	int njets = getNjetsAll(vtx,vf);
	histos1[name+"_jet_n"+suffix]->Fill(njets,wgt);
	
	int nbjets = getNbjetsAll(vtx,vf,0.3511);
	histos1[name+"_jet_b"+suffix]->Fill(nbjets,wgt);
	
	if(vf["jet_pt1"]->at(vtx)>minptjet)
	{	
		histos1[name+"_dPhi3bodyJ1"+suffix]->Fill(vf["jet_dphi3muJ1"]->at(vtx),wgt);
		histos2[name+"_dPhi3bodyJet1_vs_pTjet1"+suffix]->Fill(vf["jet_pt1"]->at(vtx),vf["jet_dphi3muJ1"]->at(vtx),wgt);
		
		float dPhiJet1MET = fabs(TVector2::Phi_mpi_pi(vf["jet_pt1"]->at(vtx)-f["met_"+metType+"_phi"]));
		histos1[name+"_dPhiMETJ1"+suffix]->Fill(dPhiJet1MET,wgt);
		
		if(vf["jet_pt2"]->at(vtx)>minptjet)
		{	
			histos1[name+"_dPhiJ1J2"+suffix]->Fill(vf["jet_dphiJ1J2"]->at(vtx),wgt);
			histos2[name+"_dPhiJet1Jet2_vs_sumpTjet12"+suffix]->Fill(vf["jet_sumpt12"]->at(vtx),vf["jet_dphiJ1J2"]->at(vtx),wgt);
			
			TLorentzVector Jet1, Jet2;
			Jet1.SetPtEtaPhiE(vf["jet_pt1"]->at(vtx),vf["jet_eta1"]->at(vtx),vf["jet_phi1"]->at(vtx),vf["jet_E1"]->at(vtx));
			Jet1.SetPtEtaPhiE(vf["jet_pt2"]->at(vtx),vf["jet_eta2"]->at(vtx),vf["jet_phi2"]->at(vtx),vf["jet_E2"]->at(vtx));
			histos1[name+"_mJ1J2"+suffix]->Fill((Jet1+Jet2).M(),wgt);
		}	
	}
	
	// b-jets
	float scale = 0;
	for(int i=1 ; i<=4 ; ++i)
	{
		TString si = tstr(i,0);
		if(vf["jet_pt"+si]->at(vtx)>minptjet) scale++;
	}
	for(int i=1 ; i<=4 ; ++i)
	{
		TString si = tstr(i,0);
		if(vf["jet_pt"+si]->at(vtx)>minptjet)
		{
			histos2[name+"_MV1_vs_mjet"+suffix]->Fill(vf["jet_m"+si]->at(vtx),vf["jet_MV1w"+si]->at(vtx),wgt/scale);
			histos2[name+"_MV1_vs_pTjet"+suffix]->Fill(vf["jet_pt"+si]->at(vtx),vf["jet_MV1w"+si]->at(vtx),wgt/scale);

			histos1[name+"_JVFall"+suffix]->Fill(vf["jet_vtxf"+si]->at(vtx),wgt/scale);
			if(i==1) histos1[name+"_JVF1"+suffix]->Fill(vf["jet_vtxf"+si]->at(vtx),wgt);
			if(i==2) histos1[name+"_JVF2"+suffix]->Fill(vf["jet_vtxf"+si]->at(vtx),wgt);
			if(i==3) histos1[name+"_JVF3"+suffix]->Fill(vf["jet_vtxf"+si]->at(vtx),wgt);
			if(i==4) histos1[name+"_JVF4"+suffix]->Fill(vf["jet_vtxf"+si]->at(vtx),wgt);
			
			if(i==1) histos2[name+"_JVF1_vs_pTjet1"+suffix]->Fill(vf["jet_pt"+si]->at(vtx),vf["jet_vtxf"+si]->at(vtx),wgt);
			if(i==2) histos2[name+"_JVF2_vs_pTjet2"+suffix]->Fill(vf["jet_pt"+si]->at(vtx),vf["jet_vtxf"+si]->at(vtx),wgt);
		}
	}
	
	TLorentzVector p3body = getP3body(vtx,vf);
	TLorentzVector pSum = getPsum(vtx,vf);
	float dphiHTMET = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-f["met_"+metType+"_phi"]));
	float dphiHTMETtrk = fabs(TVector2::Phi_mpi_pi(pSum.Phi()-f["met_track_phi"]));
	float dRHT3body = fabs(pSum.DeltaR(p3body));
	float mtHT      = sqrt(2*pSum.Pt()*f["met_"+metType+"_et"]*(1-TMath::Cos(pSum.Phi()-f["met_"+metType+"_phi"])));
	float mtHTtrk   = sqrt(2*pSum.Pt()*f["met_track_et"]*(1-TMath::Cos(pSum.Phi()-f["met_track_phi"])));
	float METsig    = f["met_"+metType+"_et"]/sqrt(pSum.Pt());
	histos1[name+"_ht_pt"+suffix]->Fill(pSum.Pt(),wgt);
	histos1[name+"_ht_mT"+suffix]->Fill(mtHT,wgt);
	histos1[name+"_ht_mT_track"+suffix]->Fill(mtHTtrk,wgt);
	histos1[name+"_ht_metsig"+suffix]->Fill(METsig,wgt);
	histos1[name+"_ht_dphimet"+suffix]->Fill(dphiHTMET,wgt);
	histos1[name+"_ht_dphimet_track"+suffix]->Fill(dphiHTMETtrk,wgt);
	histos1[name+"_ht_dR3body"+suffix]->Fill(dRHT3body,wgt);
}

void fillAllObjectHistos(unsigned int vtx, TString name, TMapTSP2vi& vi, TMapTSP2vf& vf, float wgt)
{
	if(vi["mu_order1"]->at(vtx)!=1 ||
	   vi["mu_order2"]->at(vtx)!=2 ||
	   vi["mu_order3"]->at(vtx)!=3 ) _ERR(1,"objects are not ordered correctly in the ntuple");
	
	histos1[name+"_trk_pt1_before_muons"]->Fill(vf["mu_pt1"]->at(vtx),wgt);
	histos1[name+"_trk_pt2_before_muons"]->Fill(vf["mu_pt2"]->at(vtx),wgt);
	histos1[name+"_trk_pt3_before_muons"]->Fill(vf["mu_pt3"]->at(vtx),wgt);
	
	histos1[name+"_trk_eta1_before_muons"]->Fill(vf["mu_eta1"]->at(vtx),wgt);
	histos1[name+"_trk_eta2_before_muons"]->Fill(vf["mu_eta2"]->at(vtx),wgt);
	histos1[name+"_trk_eta3_before_muons"]->Fill(vf["mu_eta3"]->at(vtx),wgt);
	
	histos1[name+"_trk_phi1_before_muons"]->Fill(vf["mu_phi1"]->at(vtx),wgt);
	histos1[name+"_trk_phi2_before_muons"]->Fill(vf["mu_phi2"]->at(vtx),wgt);
	histos1[name+"_trk_phi3_before_muons"]->Fill(vf["mu_phi3"]->at(vtx),wgt);
	
	if(vi["mu_type1"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt_before_muons"]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_mu_eta_before_muons"]->Fill(vf["mu_eta1"]->at(vtx),wgt);
		histos1[name+"_mu_isMedium1"]->Fill(vi["mu_isMedium1"]->at(vtx),wgt);	
		histos1[name+"_mu_isMedium"]->Fill(vi["mu_isMedium1"]->at(vtx),wgt);
		histos1[name+"_mu_sctangsig"]->Fill(vf["mu_sctangsig1"]->at(vtx),wgt);
		histos1[name+"_mu_sctngbsig"]->Fill(vf["mu_sctngbsig1"]->at(vtx),wgt);
		histos1[name+"_mu_pbalsig"]->Fill(vf["mu_pbalsig1"]->at(vtx),wgt);
		histos1[name+"_mu_qoverp"]->Fill(vf["mu_trkqoverp1"]->at(vtx),wgt);
		histos1[name+"_mu_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit1"]->at(vtx),wgt);
		histos1[name+"_mu_pvalue"]->Fill(vf["mu_pvaltrkfit1"]->at(vtx),wgt);	
		// histos1[name+"_mu_matchchi2ndf"]->Fill(vf["mu_trkqoverp1"]->at(vtx),wgt);
	}
	else if(vi["mu_type1"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt_before_muons"]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_TPa_eta_before_muons"]->Fill(vf["mu_eta1"]->at(vtx),wgt);
		histos1[name+"_TPa_qoverp"]->Fill(vf["mu_trkqoverp1"]->at(vtx),wgt);
		histos1[name+"_TPa_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit1"]->at(vtx),wgt);
		histos1[name+"_TPa_pvalue"]->Fill(vf["mu_pvaltrkfit1"]->at(vtx),wgt);
	}
	else if(vi["mu_type1"]->at(vtx)==TPB)
	{	
		histos1[name+"_TPb_pt_before_muons"]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_TPb_eta_before_muons"]->Fill(vf["mu_eta1"]->at(vtx),wgt);
		histos1[name+"_TPb_qoverp"]->Fill(vf["mu_trkqoverp1"]->at(vtx),wgt);
		histos1[name+"_TPb_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit1"]->at(vtx),wgt);
		histos1[name+"_TPb_pvalue"]->Fill(vf["mu_pvaltrkfit1"]->at(vtx),wgt);
	}
	
	if(vi["mu_type2"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt_before_muons"]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_mu_eta_before_muons"]->Fill(vf["mu_eta2"]->at(vtx),wgt);
		histos1[name+"_mu_isMedium2"]->Fill(vi["mu_isMedium2"]->at(vtx),wgt);
		histos1[name+"_mu_isMedium"]->Fill(vi["mu_isMedium2"]->at(vtx),wgt);
		histos1[name+"_mu_sctangsig"]->Fill(vf["mu_sctangsig2"]->at(vtx),wgt);
		histos1[name+"_mu_sctngbsig"]->Fill(vf["mu_sctngbsig2"]->at(vtx),wgt);
		histos1[name+"_mu_pbalsig"]->Fill(vf["mu_pbalsig2"]->at(vtx),wgt);
		histos1[name+"_mu_qoverp"]->Fill(vf["mu_trkqoverp2"]->at(vtx),wgt);
		histos1[name+"_mu_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit2"]->at(vtx),wgt);
		histos1[name+"_mu_pvalue"]->Fill(vf["mu_pvaltrkfit2"]->at(vtx),wgt);
		// histos1[name+"_mu_matchchi2ndf"]->Fill(vf["mu_trkqoverp2"]->at(vtx),wgt);
	}
	else if(vi["mu_type2"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt_before_muons"]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_TPa_eta_before_muons"]->Fill(vf["mu_eta2"]->at(vtx),wgt);
		histos1[name+"_TPa_qoverp"]->Fill(vf["mu_trkqoverp2"]->at(vtx),wgt);
		histos1[name+"_TPa_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit2"]->at(vtx),wgt);
		histos1[name+"_TPa_pvalue"]->Fill(vf["mu_pvaltrkfit2"]->at(vtx),wgt);
	}
	else if(vi["mu_type2"]->at(vtx)==TPB)
	{
		histos1[name+"_TPb_pt_before_muons"]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_TPb_eta_before_muons"]->Fill(vf["mu_eta2"]->at(vtx),wgt);
		histos1[name+"_TPb_qoverp"]->Fill(vf["mu_trkqoverp2"]->at(vtx),wgt);
		histos1[name+"_TPb_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit2"]->at(vtx),wgt);
		histos1[name+"_TPb_pvalue"]->Fill(vf["mu_pvaltrkfit2"]->at(vtx),wgt);
	}
	
	if(vi["mu_type3"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt_before_muons"]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_mu_eta_before_muons"]->Fill(vf["mu_eta3"]->at(vtx),wgt);
		histos1[name+"_mu_isMedium3"]->Fill(vi["mu_isMedium3"]->at(vtx),wgt);
		histos1[name+"_mu_isMedium"]->Fill(vi["mu_isMedium3"]->at(vtx),wgt);
		histos1[name+"_mu_sctangsig"]->Fill(vf["mu_sctangsig3"]->at(vtx),wgt);
		histos1[name+"_mu_sctngbsig"]->Fill(vf["mu_sctngbsig3"]->at(vtx),wgt);
		histos1[name+"_mu_pbalsig"]->Fill(vf["mu_pbalsig3"]->at(vtx),wgt);
		histos1[name+"_mu_qoverp"]->Fill(vf["mu_trkqoverp3"]->at(vtx),wgt);
		histos1[name+"_mu_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit3"]->at(vtx),wgt);
		histos1[name+"_mu_pvalue"]->Fill(vf["mu_pvaltrkfit3"]->at(vtx),wgt);
		// histos1[name+"_mu_matchchi2ndf"]->Fill(vf["mu_trkqoverp3"]->at(vtx),wgt);
	}
	else if(vi["mu_type3"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt_before_muons"]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_TPa_eta_before_muons"]->Fill(vf["mu_eta3"]->at(vtx),wgt);
		histos1[name+"_TPa_qoverp"]->Fill(vf["mu_trkqoverp3"]->at(vtx),wgt);
		histos1[name+"_TPa_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit3"]->at(vtx),wgt);
		histos1[name+"_TPa_pvalue"]->Fill(vf["mu_pvaltrkfit3"]->at(vtx),wgt);
	}
	else if(vi["mu_type3"]->at(vtx)==TPB)
	{
		histos1[name+"_TPb_pt_before_muons"]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_TPb_eta_before_muons"]->Fill(vf["mu_eta3"]->at(vtx),wgt);
		histos1[name+"_TPb_qoverp"]->Fill(vf["mu_trkqoverp3"]->at(vtx),wgt);
		histos1[name+"_TPb_chi2ndf"]->Fill(vf["mu_chi2ndftrkfit3"]->at(vtx),wgt);
		histos1[name+"_TPb_pvalue"]->Fill(vf["mu_pvaltrkfit3"]->at(vtx),wgt);
	}
}



void fillObjectValidationHistos(unsigned int vtx, TString name, TString suffix, TMapTSP2vi& vi, TMapTSP2vf& vf, float wgt)
{	
	bool onJpsi  = suffix.Contains("onJpsi");
	bool offJpsi = suffix.Contains("offJpsi");
	if(!onJpsi && !offJpsi) _FAT("suffix must be either on or off J/psi");

	TLorentzVector p1, p2, p3;
	p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
	p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
	p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
	TLorentzVector p12 = p1+p2;
	TLorentzVector p13 = p1+p3;
	TLorentzVector p23 = p2+p3;
	float m12 = p12.M();
	float m13 = p13.M();
	float m23 = p23.M();
	float q1 = (vf["mu_trkqoverp1"]->at(vtx)<0) ? -1. : +1;
	float q2 = (vf["mu_trkqoverp2"]->at(vtx)<0) ? -1. : +1;
	float q3 = (vf["mu_trkqoverp3"]->at(vtx)<0) ? -1. : +1;
	float q12 = fabs(q1+q2);
	float q13 = fabs(q1+q3);
	float q23 = fabs(q2+q3);
	
	//// safety
	if(m12<=0. || m13<=0. || m23<=0.) return;
	
	bool onJpsi12 = (m12>2850. && m12<3350. && q12==0.);
	bool onJpsi13 = (m13>2850. && m13<3350. && q13==0.);
	bool onJpsi23 = (m23>2850. && m23<3350. && q23==0.);
	bool offJpsi12 = (m12>0. && m12<750. && q12==0.);
	bool offJpsi13 = (m13>0. && m13<750. && q13==0.);
	bool offJpsi23 = (m23>0. && m23<750. && q23==0.);

	// allow only one J/psi doublet per triplet
	if(onJpsi12+onJpsi13+onJpsi23>1) return;
	
	TString i1, i2;
	float m2body = -1.;
	bool isOnJpsi = (onJpsi12 || onJpsi13 || onJpsi23);
	bool isOffJpsi = (offJpsi12 || offJpsi13 || offJpsi23);
	if(isOnJpsi)
	{
		if(!onJpsi) return; /// obviously...
		if     (onJpsi12)  { i1="1"; i2="2"; m2body = m12; }
		else if(onJpsi13)  { i1="1"; i2="3"; m2body = m13; }
		else if(onJpsi23)  { i1="2"; i2="3"; m2body = m23; }
		else _ERR(1,"bla");
	}
	else
	{
		if(isOffJpsi)
		{
			if(!offJpsi) return;
			if     (offJpsi12) { i1="1"; i2="2"; m2body = m12; }
			else if(offJpsi13) { i1="1"; i2="3"; m2body = m13; }
			else if(offJpsi23) { i1="2"; i2="3"; m2body = m23; }
			else _ERR(1,"bla");
		}
		else return;
	}
	
	// the doublet invariant mass histo, just to be sure we're lookin on the right regions
	histos1[name+"_m2body"+suffix]->Fill(m2body,wgt);
	
	// the fit probability product and the sum of dE/dx
	float trksprob = vf["mu_pvaltrkfit"+i1]->at(vtx)*vf["mu_pvaltrkfit"+i2]->at(vtx);
	float trksdEdx = vf["mu_pixeldEdx"+i1]->at(vtx)+vf["mu_pixeldEdx"+i2]->at(vtx);
	histos1[name+"_trks_hitstudy_fitprob"+suffix]->Fill(trksprob,wgt);
	histos1[name+"_trks_hitstudy_pixdEdx"+suffix]->Fill(trksdEdx,wgt);
	
	histos1[name+"_trk_hitstudy_pt1"+suffix]->Fill(vf["mu_pt"+i1]->at(vtx),wgt);
	histos1[name+"_trk_hitstudy_pt2"+suffix]->Fill(vf["mu_pt"+i2]->at(vtx),wgt);
	histos1[name+"_trk_hitstudy_eta1"+suffix]->Fill(vf["mu_eta"+i1]->at(vtx),wgt);
	histos1[name+"_trk_hitstudy_eta2"+suffix]->Fill(vf["mu_eta"+i2]->at(vtx),wgt);
	histos1[name+"_trk_hitstudy_phi1"+suffix]->Fill(vf["mu_phi"+i1]->at(vtx),wgt);
	histos1[name+"_trk_hitstudy_phi2"+suffix]->Fill(vf["mu_phi"+i2]->at(vtx),wgt);
	
	if(vi["mu_type"+i1]->at(vtx)==MUON)
	{
		histos1[name+"_mu_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i1]->at(vtx),wgt);
		histos1[name+"_mu_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i1]->at(vtx),wgt);	
		histos1[name+"_mu_isMedium"+suffix]->Fill(vi["mu_isMedium"+i1]->at(vtx),wgt);
		histos1[name+"_mu_sctangsig"+suffix]->Fill(vf["mu_sctangsig"+i1]->at(vtx),wgt);
		histos1[name+"_mu_sctngbsig"+suffix]->Fill(vf["mu_sctngbsig"+i1]->at(vtx),wgt);
		histos1[name+"_mu_pbalsig"+suffix]->Fill(vf["mu_pbalsig"+i1]->at(vtx),wgt);
		histos1[name+"_mu_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i1]->at(vtx),wgt);
		histos1[name+"_mu_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i1]->at(vtx),wgt);
		histos1[name+"_mu_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i1]->at(vtx),wgt);	
		// histos1[name+"_mu_matchchi2ndf"+suffix]->Fill(vf["mu_trkqoverp"+i1]->at(vtx),wgt);
	}
	else if(vi["mu_type"+i1]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i1]->at(vtx),wgt);
		histos1[name+"_TPa_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i1]->at(vtx),wgt);
		histos1[name+"_TPa_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i1]->at(vtx),wgt);
		histos1[name+"_TPa_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i1]->at(vtx),wgt);
		histos1[name+"_TPa_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i1]->at(vtx),wgt);
	}
	else if(vi["mu_type"+i1]->at(vtx)==TPB)
	{	
		histos1[name+"_TPb_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i1]->at(vtx),wgt);
		histos1[name+"_TPb_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i1]->at(vtx),wgt);
		histos1[name+"_TPb_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i1]->at(vtx),wgt);
		histos1[name+"_TPb_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i1]->at(vtx),wgt);
		histos1[name+"_TPb_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i1]->at(vtx),wgt);
	}
	
	if(vi["mu_type"+i2]->at(vtx)==MUON)
	{
		histos1[name+"_mu_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i2]->at(vtx),wgt);
		histos1[name+"_mu_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i2]->at(vtx),wgt);
		histos1[name+"_mu_isMedium"+suffix]->Fill(vi["mu_isMedium"+i2]->at(vtx),wgt);
		histos1[name+"_mu_sctangsig"+suffix]->Fill(vf["mu_sctangsig"+i2]->at(vtx),wgt);
		histos1[name+"_mu_sctngbsig"+suffix]->Fill(vf["mu_sctngbsig"+i2]->at(vtx),wgt);
		histos1[name+"_mu_pbalsig"+suffix]->Fill(vf["mu_pbalsig"+i2]->at(vtx),wgt);
		histos1[name+"_mu_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i2]->at(vtx),wgt);
		histos1[name+"_mu_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i2]->at(vtx),wgt);
		histos1[name+"_mu_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i2]->at(vtx),wgt);	
		// histos1[name+"_mu_matchchi2ndf"+suffix]->Fill(vf["mu_trkqoverp"+i2]->at(vtx),wgt);
	}
	else if(vi["mu_type"+i2]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i2]->at(vtx),wgt);
		histos1[name+"_TPa_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i2]->at(vtx),wgt);
		histos1[name+"_TPa_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i2]->at(vtx),wgt);
		histos1[name+"_TPa_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i2]->at(vtx),wgt);
		histos1[name+"_TPa_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i2]->at(vtx),wgt);
	}
	else if(vi["mu_type"+i2]->at(vtx)==TPB)
	{	
		histos1[name+"_TPb_hitstudy_pt"+suffix]->Fill(vf["mu_pt"+i2]->at(vtx),wgt);
		histos1[name+"_TPb_hitstudy_eta"+suffix]->Fill(vf["mu_eta"+i2]->at(vtx),wgt);
		histos1[name+"_TPb_qoverp"+suffix]->Fill(vf["mu_trkqoverp"+i2]->at(vtx),wgt);
		histos1[name+"_TPb_chi2ndf"+suffix]->Fill(vf["mu_chi2ndftrkfit"+i2]->at(vtx),wgt);
		histos1[name+"_TPb_pvalue"+suffix]->Fill(vf["mu_pvaltrkfit"+i2]->at(vtx),wgt);
	}
	
	


	float htTRTfrac1 = (vi["mu_nTRThits"+i1]->at(vtx)>0) ? (float)vi["mu_htTRThits"+i1]->at(vtx)/(float)vi["mu_nTRThits"+i1]->at(vtx) : 0.;
	float htTRTfrac2 = (vi["mu_nTRThits"+i2]->at(vtx)>0) ? (float)vi["mu_htTRThits"+i2]->at(vtx)/(float)vi["mu_nTRThits"+i2]->at(vtx) : 0.;
	int ltTRThits1 = vi["mu_nTRThits"+i1]->at(vtx) - vi["mu_htTRThits"+i1]->at(vtx);
	int ltTRThits2 = vi["mu_nTRThits"+i2]->at(vtx) - vi["mu_htTRThits"+i2]->at(vtx);
	float htTRTrat1 = (ltTRThits1>0) ? (float)vi["mu_htTRThits"+i1]->at(vtx)/(float)ltTRThits1 : 0.;
	float htTRTrat2 = (ltTRThits2>0) ? (float)vi["mu_htTRThits"+i2]->at(vtx)/(float)ltTRThits2 : 0.;
	
	if     (vi["mu_type"+i1]->at(vtx)==MUON) { histos1[name+"_mu_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i1]->at(vtx),wgt);   histos1[name+"_mu_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac1,wgt);    histos1[name+"_mu_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat1,wgt);     }
	else if(vi["mu_type"+i1]->at(vtx)==TPA)  { histos1[name+"_TPa_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i1]->at(vtx),wgt);  histos1[name+"_TPa_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac1,wgt);   histos1[name+"_TPa_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat1,wgt);    }
	else if(vi["mu_type"+i1]->at(vtx)==TPB)  { histos1[name+"_TPb_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i1]->at(vtx),wgt);  histos1[name+"_TPb_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac1,wgt);   histos1[name+"_TPb_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat1,wgt);    }
	if     (vi["mu_type"+i2]->at(vtx)==MUON) { histos1[name+"_mu_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i2]->at(vtx),wgt);   histos1[name+"_mu_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac2,wgt);    histos1[name+"_mu_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat2,wgt);     }
	else if(vi["mu_type"+i2]->at(vtx)==TPA)  { histos1[name+"_TPa_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i2]->at(vtx),wgt);  histos1[name+"_TPa_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac2,wgt);   histos1[name+"_TPa_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat2,wgt);    }
	else if(vi["mu_type"+i2]->at(vtx)==TPB)  { histos1[name+"_TPb_nHighThresholdTRTHits"+suffix]->Fill(vi["mu_htTRThits"+i2]->at(vtx),wgt);  histos1[name+"_TPb_nHighThresholdTRTHitsFraction"+suffix]->Fill(htTRTfrac2,wgt);   histos1[name+"_TPb_nHighThresholdTRTHitsRatio"+suffix]->Fill(htTRTrat2,wgt);    }
		
	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_outliersOnTrack1"+suffix]->Fill(vi["mu_nOutliersOnTrack"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_outliersOnTrack2"+suffix]->Fill(vi["mu_nOutliersOnTrack"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_outliersOnTrack1"+suffix]->Fill(vi["mu_nOutliersOnTrack"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_outliersOnTrack2"+suffix]->Fill(vi["mu_nOutliersOnTrack"+i2]->at(vtx),wgt);
	
	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_MDTHits1"+suffix]->Fill(vi["mu_nMDThits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_MDTHits2"+suffix]->Fill(vi["mu_nMDThits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_MDTHits1"+suffix]->Fill(vi["mu_nMDThits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_MDTHits2"+suffix]->Fill(vi["mu_nMDThits"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_RPCHits1"+suffix]->Fill(vi["mu_nRPCPhiHits"+i1]->at(vtx)+vi["mu_nRPCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_RPCHits2"+suffix]->Fill(vi["mu_nRPCPhiHits"+i2]->at(vtx)+vi["mu_nRPCEtaHits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_RPCHits1"+suffix]->Fill(vi["mu_nRPCPhiHits"+i1]->at(vtx)+vi["mu_nRPCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_RPCHits2"+suffix]->Fill(vi["mu_nRPCPhiHits"+i2]->at(vtx)+vi["mu_nRPCEtaHits"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_TGCHits1"+suffix]->Fill(vi["mu_nTGCPhiHits"+i1]->at(vtx)+vi["mu_nTGCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_TGCHits2"+suffix]->Fill(vi["mu_nTGCPhiHits"+i2]->at(vtx)+vi["mu_nTGCEtaHits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_TGCHits1"+suffix]->Fill(vi["mu_nTGCPhiHits"+i1]->at(vtx)+vi["mu_nTGCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_TGCHits2"+suffix]->Fill(vi["mu_nTGCPhiHits"+i2]->at(vtx)+vi["mu_nTGCEtaHits"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_precisionHits1"+suffix]->Fill(vi["mu_nPrecisionHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_precisionHits2"+suffix]->Fill(vi["mu_nPrecisionHits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_precisionHits1"+suffix]->Fill(vi["mu_nPrecisionHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_precisionHits2"+suffix]->Fill(vi["mu_nPrecisionHits"+i2]->at(vtx),wgt);
	   
	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_precisionHoles1"+suffix]->Fill(vi["mu_nPrecisionHoles"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_precisionHoles2"+suffix]->Fill(vi["mu_nPrecisionHoles"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_precisionHoles1"+suffix]->Fill(vi["mu_nPrecisionHoles"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_precisionHoles2"+suffix]->Fill(vi["mu_nPrecisionHoles"+i2]->at(vtx),wgt);
	   
	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_precisionOutliers1"+suffix]->Fill(vi["mu_nPrecisionOutliers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_precisionOutliers2"+suffix]->Fill(vi["mu_nPrecisionOutliers"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_precisionOutliers1"+suffix]->Fill(vi["mu_nPrecisionOutliers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_precisionOutliers2"+suffix]->Fill(vi["mu_nPrecisionOutliers"+i2]->at(vtx),wgt);
	
	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_phiLayers1"+suffix]->Fill(vi["mu_nPhiLayers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_phiLayers2"+suffix]->Fill(vi["mu_nPhiLayers"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_phiLayers1"+suffix]->Fill(vi["mu_nPhiLayers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_phiLayers2"+suffix]->Fill(vi["mu_nPhiLayers"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos1[name+"_mu_etaphiLayers1"+suffix]->Fill(vi["mu_nEtaPhiLayers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos1[name+"_mu_etaphiLayers2"+suffix]->Fill(vi["mu_nEtaPhiLayers"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos1[name+"_TPa_etaphiLayers1"+suffix]->Fill(vi["mu_nEtaPhiLayers"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos1[name+"_TPa_etaphiLayers2"+suffix]->Fill(vi["mu_nEtaPhiLayers"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_phiHits1"+suffix]->Fill(vi["mu_nTGCPhiHits"+i1]->at(vtx),vi["mu_nRPCPhiHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_phiHits2"+suffix]->Fill(vi["mu_nTGCPhiHits"+i2]->at(vtx),vi["mu_nRPCPhiHits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_phiHits1"+suffix]->Fill(vi["mu_nTGCPhiHits"+i1]->at(vtx),vi["mu_nRPCPhiHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_phiHits2"+suffix]->Fill(vi["mu_nTGCPhiHits"+i2]->at(vtx),vi["mu_nRPCPhiHits"+i2]->at(vtx),wgt);

	if(vi["mu_type"+i1]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_etaHits1"+suffix]->Fill(vi["mu_nTGCEtaHits"+i1]->at(vtx),vi["mu_nRPCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==MUON) histos2[name+"_mu_RPC_vs_TGC_etaHits2"+suffix]->Fill(vi["mu_nTGCEtaHits"+i2]->at(vtx),vi["mu_nRPCEtaHits"+i2]->at(vtx),wgt);
	if(vi["mu_type"+i1]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_etaHits1"+suffix]->Fill(vi["mu_nTGCEtaHits"+i1]->at(vtx),vi["mu_nRPCEtaHits"+i1]->at(vtx),wgt);
	if(vi["mu_type"+i2]->at(vtx)==TPA) histos2[name+"_TPa_RPC_vs_TGC_etaHits2"+suffix]->Fill(vi["mu_nTGCEtaHits"+i2]->at(vtx),vi["mu_nRPCEtaHits"+i2]->at(vtx),wgt);
}


void fillObjectHistos(unsigned int vtx, TString name, TString suffix, TMapTSP2vi& vi, TMapTSP2vf& vf, float wgt)
{
	histos1[name+"_trk_pt1"+suffix]->Fill(vf["mu_pt1"]->at(vtx),wgt);
	histos1[name+"_trk_pt2"+suffix]->Fill(vf["mu_pt2"]->at(vtx),wgt);
	histos1[name+"_trk_pt3"+suffix]->Fill(vf["mu_pt3"]->at(vtx),wgt);
	
	histos1[name+"_trk_eta1"+suffix]->Fill(vf["mu_eta1"]->at(vtx),wgt);
	histos1[name+"_trk_eta2"+suffix]->Fill(vf["mu_eta2"]->at(vtx),wgt);
	histos1[name+"_trk_eta3"+suffix]->Fill(vf["mu_eta3"]->at(vtx),wgt);
	
	histos1[name+"_trk_phi1"+suffix]->Fill(vf["mu_phi1"]->at(vtx),wgt);
	histos1[name+"_trk_phi2"+suffix]->Fill(vf["mu_phi2"]->at(vtx),wgt);
	histos1[name+"_trk_phi3"+suffix]->Fill(vf["mu_phi3"]->at(vtx),wgt);
	
	// the fit probability product and the sum of dE/dx
	float trksprob = vf["mu_pvaltrkfit1"]->at(vtx)*vf["mu_pvaltrkfit2"]->at(vtx)*vf["mu_pvaltrkfit3"]->at(vtx);
	float trksdEdx = vf["mu_pixeldEdx1"]->at(vtx)+vf["mu_pixeldEdx2"]->at(vtx)+vf["mu_pixeldEdx3"]->at(vtx);
	histos1[name+"_trks_fitprob"+suffix]->Fill(trksprob,wgt);
	histos1[name+"_trks_pixdEdx"+suffix]->Fill(trksdEdx,wgt);
	
	
	if(vi["mu_type1"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt"+suffix]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_mu_eta"+suffix]->Fill(vf["mu_eta1"]->at(vtx),wgt);
	}
	else if(vi["mu_type1"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt"+suffix]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_TPa_eta"+suffix]->Fill(vf["mu_eta1"]->at(vtx),wgt);
	}
	else if(vi["mu_type1"]->at(vtx)==TPB)
	{	
		histos1[name+"_TPb_pt"+suffix]->Fill(vf["mu_pt1"]->at(vtx),wgt);
		histos1[name+"_TPb_eta"+suffix]->Fill(vf["mu_eta1"]->at(vtx),wgt);
	}
	
	if(vi["mu_type2"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt"+suffix]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_mu_eta"+suffix]->Fill(vf["mu_eta2"]->at(vtx),wgt);
	}
	else if(vi["mu_type2"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt"+suffix]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_TPa_eta"+suffix]->Fill(vf["mu_eta2"]->at(vtx),wgt);
	}
	else if(vi["mu_type2"]->at(vtx)==TPB)
	{
		histos1[name+"_TPb_pt"+suffix]->Fill(vf["mu_pt2"]->at(vtx),wgt);
		histos1[name+"_TPb_eta"+suffix]->Fill(vf["mu_eta2"]->at(vtx),wgt);
	}
	
	if(vi["mu_type3"]->at(vtx)==MUON)
	{
		histos1[name+"_mu_pt"+suffix]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_mu_eta"+suffix]->Fill(vf["mu_eta3"]->at(vtx),wgt);
	}
	else if(vi["mu_type3"]->at(vtx)==TPA)
	{
		histos1[name+"_TPa_pt"+suffix]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_TPa_eta"+suffix]->Fill(vf["mu_eta3"]->at(vtx),wgt);
	}
	else if(vi["mu_type3"]->at(vtx)==TPB)
	{
		histos1[name+"_TPb_pt"+suffix]->Fill(vf["mu_pt3"]->at(vtx),wgt);
		histos1[name+"_TPb_eta"+suffix]->Fill(vf["mu_eta3"]->at(vtx),wgt);
	}
}

void fillMVAHistos(unsigned int vtx, TString name, TMapTSP2vf& vf, float score, float wgt)
{	
	histos1[name+"_MVAscore_all"]->Fill(score,wgt);
	if(vf["vtx_mass"]->at(vtx)<=tauMassMeV) histos1[name+"_MVAscore_below_tauMass"]->Fill(score,wgt);
	if(vf["vtx_mass"]->at(vtx)>tauMassMeV)  histos1[name+"_MVAscore_above_tauMass"]->Fill(score,wgt);
	if(vf["vtx_mass"]->at(vtx)<mBlindMinGlob) histos1[name+"_MVAscore_left_sideband"]->Fill(score,wgt);
	if(vf["vtx_mass"]->at(vtx)>mBlindMaxGlob) histos1[name+"_MVAscore_right_sideband"]->Fill(score,wgt);
	
	histos2[name+"_m3body_vs_MVAscore"]->Fill(score,vf["vtx_mass"]->at(vtx),wgt);
	histos2[name+"_pT3body_vs_MVAscore"]->Fill(score,vf["vtx_pt"]->at(vtx),wgt);
	if(fabs(vf["vtx_charge"]->at(vtx))==1)
	{	
		histos2[name+"_mOS_vs_MVAscore"]->Fill(score,vf["vtx_mOS1"]->at(vtx),wgt/2.);
		histos2[name+"_mOS_vs_MVAscore"]->Fill(score,vf["vtx_mOS2"]->at(vtx),wgt/2.);
		histos2[name+"_mOS1_vs_MVAscore"]->Fill(score,vf["vtx_mOS1"]->at(vtx),wgt);
		histos2[name+"_mOS2_vs_MVAscore"]->Fill(score,vf["vtx_mOS2"]->at(vtx),wgt);
		histos2[name+"_mSS_vs_MVAscore"]->Fill(score,vf["vtx_mSS"]->at(vtx),wgt);
	}
}

/*
void fillMETHistos(unsigned int vtx, TString name, TString suffix, TMapTSf& f, TMapTSP2vf& vf, TString metType, float wgt)
{	
	histos1[name+"_met_et"+suffix]->Fill(f["met_"+metType+"_et"],wgt);
	histos1[name+"_met_mt_et3mu"+suffix]->Fill(vf["met_"+metType+"_mT"]->at(vtx),wgt);
	histos1[name+"_met_dphi3mu"+suffix]->Fill(vf["met_"+metType+"_dPhi3mu"]->at(vtx),wgt);
	
	histos1[name+"_mettrk_et"+suffix]->Fill(f["met_track_et"],wgt);
	histos1[name+"_mettrk_mt_et3mu"+suffix]->Fill(vf["met_track_mT"]->at(vtx),wgt);
	histos1[name+"_mettrk_dphi3mu"+suffix]->Fill(vf["met_track_dPhi3mu"]->at(vtx),wgt);
	
	histos1[name+"_mettrk_dphimet"+suffix]->Fill(fabs(TVector2::Phi_mpi_pi(f["met_"+metType+"_phi"]-f["met_track_phi"])),wgt);
}
*/


void fillRhoOmegaPhiHistos(unsigned int vtx, TString name, TString suffix, TMapTSf& f, TMapTSP2vf& vf, TString metType, float wgt)
{	
	float mOS1 = vf["vtx_mOS1"]->at(vtx);
	float mOS2 = vf["vtx_mOS2"]->at(vtx);
	float mSS  = vf["vtx_mSS"]->at(vtx);
	float m3body = vf["vtx_mass"]->at(vtx);
	float mRho   = 770;
	float mOmega = 782;
	float mPhi   = 1020;
	float range  = 50;
	
	bool isRhoOmega = (fabs(mOS1-mRho)<range || fabs(mOS1-mOmega)<range || fabs(mOS2-mRho)<range || fabs(mOS2-mOmega)<range);
	bool isPhi      = (fabs(mOS1-mPhi)<range || fabs(mOS2-mPhi)<range);
	
	if(isRhoOmega)
	{
		histos1[name+"_m3body_lin_zoom_onRhoOmega"+suffix]->Fill(m3body,wgt);
		histos1[name+"_m3body_sigregion_onRhoOmega"+suffix]->Fill(m3body,wgt);
	}
	if(isPhi)
	{
		histos1[name+"_m3body_lin_zoom_onPhi"+suffix]->Fill(m3body,wgt);
		histos1[name+"_m3body_sigregion_onPhi"+suffix]->Fill(m3body,wgt);
		histos1[name+"_m3body_1600"+suffix]->Fill(m3body,wgt);
	}
	
	TLorentzVector p1, p2, p3, p, p12, p13, p23;
	p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
	p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
	p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
	p12 = p1+p2;
	p13 = p1+p3;
	p23 = p2+p3;
	float m12 = p12.M();
	float m13 = p13.M();
	float m23 = p23.M();
	if(m12<=0. || m13<=0. || m23<=0.) return; //// safety
	bool onPhi12 = (fabs(m12-mPhi)<range);
	bool onPhi13 = (fabs(m13-mPhi)<range);
	bool onPhi23 = (fabs(m23-mPhi)<range);
	bool onRhoOmega12 = (fabs(m12-mPhi)<range);
	bool onRhoOmega13 = (fabs(m13-mPhi)<range);
	bool onRhoOmega23 = (fabs(m23-mPhi)<range);
	
	if(onPhi12)
	{
		p.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),kaonMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p1+p2+p).M(),wgt);
		p.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),pionMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p1+p2+p).M(),wgt);
	}
	if(onPhi13)
	{
		p.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),kaonMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p1+p3+p).M(),wgt);
		p.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),pionMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p1+p3+p).M(),wgt);
	}
	if(onPhi23)
	{
		p.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),kaonMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p2+p3+p).M(),wgt);
		p.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),pionMassMeV);
		histos1[name+"_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p2+p3+p).M(),wgt);
	}
	
	if(onPhi12)
	{
		histos1[name+"_mT_mu1_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p1.Pt(),p1.Phi()),wgt);
		histos1[name+"_mT_mu2_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p2.Pt(),p2.Phi()),wgt);
	}
	if(onPhi13)
	{
		histos1[name+"_mT_mu1_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p1.Pt(),p1.Phi()),wgt);
		histos1[name+"_mT_mu3_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p3.Pt(),p3.Phi()),wgt);
	}
	if(onPhi23)
	{
		histos1[name+"_mT_mu2_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p2.Pt(),p2.Phi()),wgt);
		histos1[name+"_mT_mu3_onPhi"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p3.Pt(),p3.Phi()),wgt);
	}
	
	if(onRhoOmega12)
	{
		histos1[name+"_mT_mu1_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p1.Pt(),p1.Phi()),wgt);
		histos1[name+"_mT_mu2_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p2.Pt(),p2.Phi()),wgt);
	}
	if(onRhoOmega13)
	{
		histos1[name+"_mT_mu1_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p1.Pt(),p1.Phi()),wgt);
		histos1[name+"_mT_mu3_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p3.Pt(),p3.Phi()),wgt);
	}
	if(onRhoOmega23)
	{
		histos1[name+"_mT_mu2_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p2.Pt(),p2.Phi()),wgt);
		histos1[name+"_mT_mu3_onRhoOmega"+suffix]->Fill(mT(f["met_"+metType+"_et"]/2.,f["met_"+metType+"_phi"],p3.Pt(),p3.Phi()),wgt);
	}
}

void fillMVAevoHistos(unsigned int vtx, TString name, TString suffix, TMapTSf& f, TMapTSP2vf& vf, TMapTSP2vi& vid, TMapTSP2vf& vfd, TString metType, float score, float wgt)
{
	float mOS1 = vf["vtx_mOS1"]->at(vtx);
	float mOS2 = vf["vtx_mOS2"]->at(vtx);
	float mSS  = vf["vtx_mSS"]->at(vtx);
	float m3body = vf["vtx_mass"]->at(vtx);
	float pt3body = vf["vtx_pt"]->at(vtx);
	
	float met_cal_et = f["met_"+metType+"_et"];
	float met_trk_et = f["met_track_et"];
	float met_cal_mt = vf["met_"+metType+"_mT"]->at(vtx);
	float met_trk_mt = vf["met_track_mT"]->at(vtx);
	float met_cal_dphi3body = vf["met_"+metType+"_dPhi3mu"]->at(vtx);
	float met_trk_dphi3body = vf["met_track_dPhi3mu"]->at(vtx);
	float met_dphi = vfd["mets_dphi"]->at(vtx);
	
	int   njets  = vid["jets_n"]->at(vtx);
	float ht_pt  = vfd["ht_pt"]->at(vtx);
	float ht_mt  = vfd["ht_mT"]->at(vtx);
	float ht_mt_trk  = vfd["ht_mT_mettrk"]->at(vtx);
	float ht_dphimet_cal = vfd["ht_dphimet_"+metType]->at(vtx);
	float ht_dphimet_trk = vfd["ht_dphimet_track"]->at(vtx);
	float ht_dr3body = vfd["ht_dr3body"]->at(vtx);
	
	float mPhi = 1020;
	float range = 50;
	TLorentzVector p1, p2, p3, p;
	p1.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),muonMassMeV);
	p2.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),muonMassMeV);
	p3.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),muonMassMeV);
	TLorentzVector p12 = p1+p2;
	TLorentzVector p13 = p1+p3;
	TLorentzVector p23 = p2+p3;
	float m12 = p12.M();
	float m13 = p13.M();
	float m23 = p23.M();
	// if(m12<=0. || m13<=0. || m23<=0.) return; //// safety
	bool onPhi12 = (fabs(m12-mPhi)<range);
	bool onPhi13 = (fabs(m13-mPhi)<range);
	bool onPhi23 = (fabs(m23-mPhi)<range);
	if(m3body<mBlindMinGlob || m3body>mBlindMaxGlob) // make sure this is blinded
	{
		if(onPhi12)
		{
			p.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),kaonMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p1+p2+p).M(),wgt);
			p.SetPtEtaPhiM(vf["mu_pt3"]->at(vtx),vf["mu_eta3"]->at(vtx),vf["mu_phi3"]->at(vtx),pionMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p1+p2+p).M(),wgt);
		}
		if(onPhi13)
		{
			p.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),kaonMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p1+p3+p).M(),wgt);
			p.SetPtEtaPhiM(vf["mu_pt2"]->at(vtx),vf["mu_eta2"]->at(vtx),vf["mu_phi2"]->at(vtx),pionMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p1+p3+p).M(),wgt);
		}
		if(onPhi23)
		{
			p.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),kaonMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_kaon_hypothesis"+suffix]->Fill((p2+p3+p).M(),wgt);
			p.SetPtEtaPhiM(vf["mu_pt1"]->at(vtx),vf["mu_eta1"]->at(vtx),vf["mu_phi1"]->at(vtx),pionMassMeV);
			histos1[name+"_MVAevo_m3body_sigregion_onPhi_pion_hypothesis"+suffix]->Fill((p2+p3+p).M(),wgt);
		}
	}
	
	histos1[name+"_MVAevo_score"+suffix]->Fill(score,wgt);
	
	histos1[name+"_MVAevo_m3body_sigregion"+suffix]->Fill(m3body,wgt);
	histos1[name+"_MVAevo_pT3body"+suffix]->Fill(pt3body,wgt);
	histos1[name+"_MVAevo_mOS1"+suffix]->Fill(mOS1,wgt);
	histos1[name+"_MVAevo_mOS2"+suffix]->Fill(mOS2,wgt);
	histos1[name+"_MVAevo_mSS"+suffix]->Fill(mSS,wgt);
	
	histos1[name+"_MVAevo_met_cal_et"+suffix]->Fill(met_cal_et,wgt);
	histos1[name+"_MVAevo_met_track_et"+suffix]->Fill(met_trk_et,wgt);
	histos1[name+"_MVAevo_met_cal_mT"+suffix]->Fill(met_cal_mt,wgt);
	histos1[name+"_MVAevo_met_track_mT"+suffix]->Fill(met_trk_mt,wgt);
	histos1[name+"_MVAevo_met_cal_dPhi3mu"+suffix]->Fill(met_cal_dphi3body,wgt);
	histos1[name+"_MVAevo_met_track_dPhi3mu"+suffix]->Fill(met_trk_dphi3body,wgt);
	histos1[name+"_MVAevo_met_dphiCalTrk"+suffix]->Fill(met_dphi,wgt);
	
	histos1[name+"_MVAevo_jet_n"+suffix]->Fill(njets,wgt);
	histos1[name+"_MVAevo_jet_ht_pt"+suffix]->Fill(ht_pt,wgt);
	histos1[name+"_MVAevo_jet_ht_mT_metcal"+suffix]->Fill(ht_mt,wgt);
	histos1[name+"_MVAevo_jet_ht_mT_mettrk"+suffix]->Fill(ht_mt_trk,wgt);
	histos1[name+"_MVAevo_jet_ht_dphimetcal"+suffix]->Fill(ht_dphimet_cal,wgt);
	histos1[name+"_MVAevo_jet_ht_dphimettrk"+suffix]->Fill(ht_dphimet_trk,wgt);
	histos1[name+"_MVAevo_jet_ht_dr3body"+suffix]->Fill(ht_dr3body,wgt);
	
	histos2[name+"_MVAevo_mOS2_vs_mOS1"+suffix]->Fill(mOS2,mSS,wgt);
	histos2[name+"_MVAevo_mSS_vs_mOS1"+suffix]->Fill(mOS1,mSS,wgt);
	histos2[name+"_MVAevo_mSS_vs_mOS2"+suffix]->Fill(mOS2,mSS,wgt);
	histos2[name+"_MVAevo_m3body_vs_mOS1"+suffix]->Fill(mOS1,m3body,wgt);
	histos2[name+"_MVAevo_m3body_vs_mOS2"+suffix]->Fill(mOS2,m3body,wgt);
	histos2[name+"_MVAevo_m3body_vs_mSS"+suffix]->Fill(mSS,m3body,wgt);
}


void bookHistos(TString name)
{
	olddir->cd();

	makeCutflowHisto(name);
	// makeCountersHisto(name);
	makeCategoriesHisto(name);
	makeTriggerHisto(name);
	
	addHist(histos1,name,"chi2ndf",   ";3body vertex fit #chi^{2};Normalized", 50,0.,50.);
	addHist(histos1,name,"pvalue",    ";3body vertex fit #it{p}-value;Normalized", 200,0.,1.);
	addHist(histos1,name,"Lxy",       ";3body vertex #it{L}_{xy} [#mum] (PV type 1);Normalized", 60,-3.,+23.);
	addHist(histos1,name,"dLxy",      ";3body vertex #Delta#it{L}_{xy} [#mum] (PV type 1);Normalized", 50,0.,+1.5);
	addHist(histos1,name,"a0",        ";3body vertex #it{a}_{0} [#mum];Normalized", 100,0.,+150.);
	addHist(histos1,name,"a0xy",      ";3body vertex #it{a}_{0}^{xy} [#mum];Normalized", 100,0.,+5.);
	addHist(histos1,name,"da0xy",     ";3body vertex #Delta#it{a}_{0}^{xy} [#mum];Normalized", 50,0.,+0.05);
	addHist(histos1,name,"cosT",      ";3body vertex cos#theta;Normalized", 100,-1.,+1.);
	addHist(histos1,name,"cosTxy",    ";3body vertex cos#theta_{xy};Normalized", 100,.0,+1.);
	
	addHist(histos1,name,"PVNtrk_after_vertexing",  ";Refitted PV N_{trk} (after vertexing);Normalized", 200,0,200);
	addHist(histos1,name,"PVNtrk_before_isolation", ";Refitted PV N_{trk} (before isolation);Normalized", 200,0,200);
	addHist(histos1,name,"PVNtrk_after_isolation",  ";Refitted PV N_{trk} (after isolation);Normalized", 200,0,200);
	addHist(histos1,name,"PVNtrk_after_ht",         ";Refitted PV N_{trk} (after loose cuts);Normalized", 200,0,200);
	
	addHist(histos1,name,"trk_eta1_before_muons", ";#eta^{trk1};Normalized",   64,-3.0,+3.0,false);
	addHist(histos1,name,"trk_eta2_before_muons", ";#eta^{trk2} ;Normalized", 64,-3.0,+3.0,false);
	addHist(histos1,name,"trk_eta3_before_muons", ";#eta^{trk3};Normalized",   64,-3.0,+3.0,false);
	addHist(histos1,name,"trk_phi1_before_muons", ";#phi^{trk1};Normalized",   64,-TMath::Pi(),+TMath::Pi(),false);
	addHist(histos1,name,"trk_phi2_before_muons", ";#phi^{trk2} ;Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
	addHist(histos1,name,"trk_phi3_before_muons", ";#phi^{trk3};Normalized",   64,-TMath::Pi(),+TMath::Pi(),false);
	addHist(histos1,name,"trk_pt1_before_muons", ";p_{T}^{trk1} [MeV];Normalized",   50,2000.,102000.,false);
	addHist(histos1,name,"trk_pt2_before_muons", ";p_{T}^{trk2} [MeV];Normalized",   50,2000.,42000.,false);
	addHist(histos1,name,"trk_pt3_before_muons", ";p_{T}^{trk3} [MeV];Normalized",   50,2000.,27000.,false);
	
	addHist(histos1,name,"mu_pt_before_muons",  ";#mu p_{T} [MeV];Normalized", 120,2000.,122000.,false);
	addHist(histos1,name,"TPa_pt_before_muons", ";TPa p_{T} [MeV];Normalized", 60,2000.,62000.,false);
	addHist(histos1,name,"TPb_pt_before_muons", ";TPb p_{T} [MeV];Normalized", 60,2000.,32000.,false);
	addHist(histos1,name,"mu_eta_before_muons", ";#mu #eta;Normalized",   64,-3.0,+3.0,false);
	addHist(histos1,name,"TPa_eta_before_muons",";TPa #eta;Normalized",   64,-3.0,+3.0,false);
	addHist(histos1,name,"TPb_eta_before_muons",";TPb #eta;Normalized",   64,-3.0,+3.0,false);
	
	for(int i=0 ; i<3 ; ++i)
	{
		TString suf = "";
		TString txt = "";
		switch(i)
		{
			case 0: suf = "";         txt = "";                  break;
			case 1: suf = "_onJpsi";  txt = " (on the J/#psi)";  break;
			case 2: suf = "_offJpsi"; txt = " (off the J/#psi)"; break;
			default: break;
		}
		
		addHist(histos1,name,"m2body"+suf, ";m_{2body} [MeV];Normalized", 100,0.,4000.);
		
		addHist(histos1,name,"trks_hitstudy_fitprob"+suf,  ";Tracks #prodp-value"+txt+";Normalized",50,0.,1.,false);
		addHist(histos1,name,"trks_hitstudy_pixdEdx"+suf,  ";Tracks #SigmadE/dx (pixel)"+txt+";Normalized",50,0.,8.,false);
			
		addHist(histos1,name,"trk_hitstudy_eta1"+suf, ";#eta^{trk1}"+txt+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_hitstudy_eta2"+suf, ";#eta^{trk2}"+txt+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_hitstudy_eta3"+suf, ";#eta^{trk3}"+txt+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_hitstudy_phi1"+suf, ";#phi^{trk1}"+txt+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_hitstudy_phi2"+suf, ";#phi^{trk2}"+txt+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_hitstudy_phi3"+suf, ";#phi^{trk3}"+txt+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_hitstudy_pt1"+suf,  ";p_{T}^{trk1} [MeV]"+txt+";Normalized", 50,2000.,102000.,false);
		addHist(histos1,name,"trk_hitstudy_pt2"+suf,  ";p_{T}^{trk2} [MeV]"+txt+";Normalized", 50,2000.,42000., false);
		addHist(histos1,name,"trk_hitstudy_pt3"+suf,  ";p_{T}^{trk3} [MeV]"+txt+";Normalized", 50,2000.,27000., false);
		
		addHist(histos1,name,"mu_hitstudy_pt"+suf,  ";#mu p_{T} [MeV]"+txt+";Normalized", 120,2000.,122000.,false);
		addHist(histos1,name,"TPa_hitstudy_pt"+suf, ";TPa p_{T} [MeV]"+txt+";Normalized", 60,2000.,62000.,false);
		addHist(histos1,name,"TPb_hitstudy_pt"+suf, ";TPb p_{T} [MeV]"+txt+";Normalized", 60,2000.,32000.,false);
		addHist(histos1,name,"mu_hitstudy_eta"+suf, ";#mu #eta"+txt+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"TPa_hitstudy_eta"+suf,";TPa #eta"+txt+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"TPb_hitstudy_eta"+suf,";TPb #eta"+txt+";Normalized", 64,-3.0,+3.0,false);
		
		addHist(histos1,name,"mu_isMedium"+suf,  ";#mu isMedium"+txt+";Normalized", 2,0.,2);
		addHist(histos1,name,"mu_isMedium1"+suf, ";#mu1 isMedium"+txt+";Normalized", 2,0.,2);
		addHist(histos1,name,"mu_isMedium2"+suf, ";#mu2 isMedium"+txt+";Normalized", 2,0.,2);
		addHist(histos1,name,"mu_isMedium3"+suf, ";#mu3 isMedium"+txt+";Normalized", 2,0.,2);

		addHist(histos1,name,"mu_sctangsig"+suf,   ";#mu Scattering curvature significance"+txt+";Normalized",  60,-20.,30.);
		addHist(histos1,name,"mu_sctngbsig"+suf,   ";#mu Scattering neighbours significance"+txt+";Normalized", 60,-20.,30.);
		addHist(histos1,name,"mu_pbalsig"+suf,     ";#mu Momentum balance significance"+txt+";Normalized",      50,-10.,20.);

		addHist(histos1,name,"mu_qoverp"+suf,      ";#mu Q/p [1/MeV]"+txt+";Normalized",100,-7.e-4,+7.e-4);
		addHist(histos1,name,"mu_chi2ndf"+suf,     ";#mu track-fit #chi^{2}/N_{DOF}"+txt+";Normalized", 100,0.,10.);
		addHist(histos1,name,"mu_matchchi2ndf"+suf,";#mu Match #chi^{2}_{DOF}"+txt+";Normalized",100,0.,20.);
		addHist(histos1,name,"mu_pvalue"+suf,      ";#mu track-fit p-value"+txt+";Normalized", 100,0.,1.);

		addHist(histos1,name,"TPa_qoverp"+suf,   ";TPa Q/p [1/MeV]"+txt+";Normalized",100,-7.e-4,+7.e-4);
		addHist(histos1,name,"TPa_chi2ndf"+suf,  ";TPa track-fit #chi^{2}/N_{DOF}"+txt+";Normalized", 100,0.,10.);
		addHist(histos1,name,"TPa_pvalue"+suf,   ";TPa track-fit p-value"+txt+";Normalized", 100,0.,1.);

		addHist(histos1,name,"TPb_qoverp"+suf,   ";TPb Q/p [1/MeV]"+txt+";Normalized",100,-7.e-4,+7.e-4);
		addHist(histos1,name,"TPb_chi2ndf"+suf,  ";TPb track-fit #chi^{2}/N_{DOF}"+txt+";Normalized", 100,0.,10.);
		addHist(histos1,name,"TPb_pvalue"+suf,   ";TPb track-fit p-value"+txt+";Normalized", 100,0.,1.);
		
		
		addHist(histos1,name,"mu_outliersOnTrack1"+suf,   ";#mu1 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_outliersOnTrack2"+suf,   ";#mu2 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_outliersOnTrack3"+suf,   ";#mu3 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_outliersOnTrack1"+suf,   ";TPa1 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_outliersOnTrack2"+suf,   ";TPa2 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_outliersOnTrack3"+suf,   ";TPa3 outliersOnTrack hits"+txt+";Normalized", 20,0.,20.);
		
		addHist(histos1,name,"mu_MDTHits1"+suf,   ";#mu1 MDT hits"+txt+";Normalized", 60,0.,60.);
		addHist(histos1,name,"mu_MDTHits2"+suf,   ";#mu2 MDT hits"+txt+";Normalized", 60,0.,60.);
		addHist(histos1,name,"mu_MDTHits3"+suf,   ";#mu3 MDT hits"+txt+";Normalized", 60,0.,60.);
		addHist(histos1,name,"TPa_MDTHits1"+suf,   ";TPa1 MDT hits"+txt+";Normalized", 60,0.,60.);
		addHist(histos1,name,"TPa_MDTHits2"+suf,   ";TPa2 MDT hits"+txt+";Normalized", 60,0.,60.);
		addHist(histos1,name,"TPa_MDTHits3"+suf,   ";TPa3 MDT hits"+txt+";Normalized", 60,0.,60.);
		
		addHist(histos1,name,"mu_TGCHits1"+suf,   ";#mu1 TGC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_TGCHits2"+suf,   ";#mu2 TGC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_TGCHits3"+suf,   ";#mu3 TGC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_TGCHits1"+suf,   ";TPa1 TGC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_TGCHits2"+suf,   ";TPa2 TGC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_TGCHits3"+suf,   ";TPa3 TGC hits"+txt+";Normalized", 20,0.,20.);
		
		addHist(histos1,name,"mu_RPCHits1"+suf,   ";#mu1 RPC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_RPCHits2"+suf,   ";#mu2 RPC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"mu_RPCHits3"+suf,   ";#mu3 RPC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_RPCHits1"+suf,   ";TPa1 RPC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_RPCHits2"+suf,   ";TPa2 RPC hits"+txt+";Normalized", 20,0.,20.);
		addHist(histos1,name,"TPa_RPCHits3"+suf,   ";TPa3 RPC hits"+txt+";Normalized", 20,0.,20.);
		
		addHist(histos2,name,"mu_RPC_vs_TGC_phiHits1"+suf,   "#mu1 RPC vs TGC #phi hits"+txt+";#mu1 TGC #phi hits;#mu1 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"mu_RPC_vs_TGC_phiHits2"+suf,   "#mu2 RPC vs TGC #phi hits"+txt+";#mu2 TGC #phi hits;#mu2 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"mu_RPC_vs_TGC_phiHits3"+suf,   "#mu3 RPC vs TGC #phi hits"+txt+";#mu3 TGC #phi hits;#mu3 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_phiHits1"+suf,   "TPa1 RPC vs TGC #phi hits"+txt+";TPa1 TGC #phi hits;TPa1 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_phiHits2"+suf,   "TPa2 RPC vs TGC #phi hits"+txt+";TPa2 TGC #phi hits;TPa2 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_phiHits3"+suf,   "TPa3 RPC vs TGC #phi hits"+txt+";TPa3 TGC #phi hits;TPa3 RPC #phi hits;Normalized", 20,0.,20., 20,0.,20.);
		
		addHist(histos2,name,"mu_RPC_vs_TGC_etaHits1"+suf,   "#mu1 RPC vs TGC #eta hits"+txt+";#mu1 TGC #eta hits;#mu1 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"mu_RPC_vs_TGC_etaHits2"+suf,   "#mu2 RPC vs TGC #eta hits"+txt+";#mu2 TGC #eta hits;#mu2 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"mu_RPC_vs_TGC_etaHits3"+suf,   "#mu3 RPC vs TGC #eta hits"+txt+";#mu3 TGC #eta hits;#mu3 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_etaHits1"+suf,   "TPa1 RPC vs TGC #eta hits"+txt+";TPa1 TGC #eta hits;TPa1 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_etaHits2"+suf,   "TPa2 RPC vs TGC #eta hits"+txt+";TPa2 TGC #eta hits;TPa2 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		addHist(histos2,name,"TPa_RPC_vs_TGC_etaHits3"+suf,   "TPa3 RPC vs TGC #eta hits"+txt+";TPa3 TGC #eta hits;TPa3 RPC #eta hits;Normalized", 20,0.,20., 20,0.,20.);
		
		addHist(histos1,name,"mu_precisionHits1"+suf,   ";#mu1 precision hits"+txt+";Normalized", 50,0.,50.);
		addHist(histos1,name,"mu_precisionHits2"+suf,   ";#mu2 precision hits"+txt+";Normalized", 50,0.,50.);
		addHist(histos1,name,"mu_precisionHits3"+suf,   ";#mu3 precision hits"+txt+";Normalized", 50,0.,50.);
		addHist(histos1,name,"TPa_precisionHits1"+suf,   ";TPa1 precision hits"+txt+";Normalized", 50,0.,50.);
		addHist(histos1,name,"TPa_precisionHits2"+suf,   ";TPa2 precision hits"+txt+";Normalized", 50,0.,50.);
		addHist(histos1,name,"TPa_precisionHits3"+suf,   ";TPa3 precision hits"+txt+";Normalized", 50,0.,50.);
		
		addHist(histos1,name,"mu_precisionHoles1"+suf,   ";#mu1 precision holes"+txt+";Normalized", 25,0.,25.);
		addHist(histos1,name,"mu_precisionHoles2"+suf,   ";#mu2 precision holes"+txt+";Normalized", 25,0.,25.);
		addHist(histos1,name,"mu_precisionHoles3"+suf,   ";#mu3 precision holes"+txt+";Normalized", 25,0.,25.);
		addHist(histos1,name,"TPa_precisionHoles1"+suf,   ";TPa1 precision holes"+txt+";Normalized", 25,0.,25.);
		addHist(histos1,name,"TPa_precisionHoles2"+suf,   ";TPa2 precision holes"+txt+";Normalized", 25,0.,25.);
		addHist(histos1,name,"TPa_precisionHoles3"+suf,   ";TPa3 precision holes"+txt+";Normalized", 25,0.,25.);
		
		addHist(histos1,name,"mu_precisionOutliers1"+suf,   ";#mu1 precision outliers"+txt+";Normalized", 15,0.,15.);
		addHist(histos1,name,"mu_precisionOutliers2"+suf,   ";#mu2 precision outliers"+txt+";Normalized", 15,0.,15.);
		addHist(histos1,name,"mu_precisionOutliers3"+suf,   ";#mu3 precision outliers"+txt+";Normalized", 15,0.,15.);
		addHist(histos1,name,"TPa_precisionOutliers1"+suf,   ";TPa1 precision outliers"+txt+";Normalized", 15,0.,15.);
		addHist(histos1,name,"TPa_precisionOutliers2"+suf,   ";TPa2 precision outliers"+txt+";Normalized", 15,0.,15.);
		addHist(histos1,name,"TPa_precisionOutliers3"+suf,   ";TPa3 precision outliers"+txt+";Normalized", 15,0.,15.);
		
		addHist(histos1,name,"mu_phiLayers1"+suf,   ";#mu1 #phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"mu_phiLayers2"+suf,   ";#mu2 #phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"mu_phiLayers3"+suf,   ";#mu3 #phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_phiLayers1"+suf,   ";TPa1 #phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_phiLayers2"+suf,   ";TPa2 #phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_phiLayers3"+suf,   ";TPa3 #phi layers"+txt+";Normalized", 10,0.,10.);
		
		addHist(histos1,name,"mu_etaphiLayers1"+suf,   ";#mu1 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"mu_etaphiLayers2"+suf,   ";#mu2 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"mu_etaphiLayers3"+suf,   ";#mu3 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_etaphiLayers1"+suf,   ";TPa1 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_etaphiLayers2"+suf,   ";TPa2 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		addHist(histos1,name,"TPa_etaphiLayers3"+suf,   ";TPa3 #eta-#phi layers"+txt+";Normalized", 10,0.,10.);
		
		addHist(histos1,name,"mu_nHighThresholdTRTHits"+suf,  ";#mu nHighThresholdTRTHits"+txt+";Normalized", 20,0,20);
		addHist(histos1,name,"TPa_nHighThresholdTRTHits"+suf, ";TPa nHighThresholdTRTHits"+txt+";Normalized", 20,0,20);
		addHist(histos1,name,"TPb_nHighThresholdTRTHits"+suf, ";TPb nHighThresholdTRTHits"+txt+";Normalized", 20,0,20);
		addHist(histos1,name,"mu_nHighThresholdTRTHitsFraction"+suf,  ";#mu nHighThresholdTRTHits/nTRTHits"+txt+";Normalized", 50,0.,1.);
		addHist(histos1,name,"TPa_nHighThresholdTRTHitsFraction"+suf, ";TPa nHighThresholdTRTHits/nTRTHits"+txt+";Normalized", 50,0.,1.);
		addHist(histos1,name,"TPb_nHighThresholdTRTHitsFraction"+suf, ";TPb nHighThresholdTRTHits/nTRTHits"+txt+";Normalized", 50,0.,1.);
		addHist(histos1,name,"mu_nHighThresholdTRTHitsRatio"+suf,  ";#mu nHighThresholdTRTHits/nLowThresholdTRTHits"+txt+";Normalized", 50,0.,1.);
		addHist(histos1,name,"TPa_nHighThresholdTRTHitsRatio"+suf, ";TPa nHighThresholdTRTHits/nLowThresholdTRTHits"+txt+";Normalized", 50,0.,1.);
		addHist(histos1,name,"TPb_nHighThresholdTRTHitsRatio"+suf, ";TPb nHighThresholdTRTHits/nLowThresholdTRTHits"+txt+";Normalized", 50,0.,1.);
	}
	
	
	addHist(histos1,name,"jet_n_calib_all",       ";N_{jets} (all);Normalized",5,0.,5.);
	addHist(histos1,name,"jet_b_calib_all",       ";N_{b-jets} (all);Normalized",5,0.,5.);
	addHist(histos1,name,"jet_pt_calib_all",      ";Jets #it{p}_{T} [MeV] (all);Events",40,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_pt1_calib_all",     ";#it{p}_{T}(jet1) [MeV] (all);Events",40,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_pt2_calib_all",     ";#it{p}_{T}(jet2) [MeV] (all);Events",40,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_ptJ1J2_calib_all",  ";#Sigma#it{p}_{T}(jet1+jet2) [MeV] (all);Events",50,0.,150.*GeV2MeV);
	addHist(histos1,name,"jet_E_calib_all",       ";Jets #it{E} [MeV] (all);Events",50,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_E1_calib_all",      ";#it{E}(jet1) [MeV] (all);Events",50,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_E2_calib_all",      ";#it{E}(jet2) [MeV] (all);Events",50,0.,100.*GeV2MeV);
	addHist(histos1,name,"jet_EJ1J2_calib_all",   ";#Sigma#it{E}(jet1+jet2) [MeV] (all);Events",50,0.,150.*GeV2MeV);
	
	addHist(histos1,name,"ht_pt_calib_all",      ";#it{H}_{T} [MeV] (all);Normalized",100,0.,100.*GeV2MeV);
	addHist(histos1,name,"ht_dphimet_calib_all", ";#Delta#phi(#it{H}_{T},#it{E}_{T}^{miss}) (all);Normalized",32,0.,TMath::Pi());
	addHist(histos1,name,"ht_dphimet_track_calib_all", ";#Delta#phi(#it{H}_{T},Track #it{E}_{T}^{miss}) (all);Normalized",32,0.,TMath::Pi());
	addHist(histos1,name,"ht_mT_calib_all",      ";#it{m}_{T}(#it{H}_{T},#it{E}_{T}^{miss}) [MeV] (all);Normalized",60,0,120.*GeV2MeV);
	addHist(histos1,name,"ht_mT_track_calib_all",";#it{m}_{T}(#it{H}_{T},Track #it{E}_{T}^{miss}) [MeV] (all);Normalized",60,0,120.*GeV2MeV);
	addHist(histos1,name,"ht_dR3body_calib_all", ";#DeltaR(#it{H}_{T},3body) (all);Normalized",32,0.,2.*TMath::Pi());
	addHist(histos1,name,"ht_metsig_calib_all",  ";#it{E}_{T}^{miss}/#sqrt{#it{H}_{T}} [#sqrt{MeV}] (all);Normalized",100,0.,600.);
	
	addHist(histos1,name,"met_et_calib_all",         ";#it{E}_{T}^{miss} [MeV] (all);Events",50,0.,100.*GeV2MeV);
	addHist(histos1,name,"met_mt_et3mu_calib_all",   ";m_{T}(3body,#it{E}_{T}^{miss}) [MeV] (all);Events",60,0,120.*GeV2MeV);
	addHist(histos1,name,"met_dphi3mu_calib_all",    ";#Delta#phi(3body,#it{E}_{T}^{miss}) (all);Events",32,0.,TMath::Pi());
	addHist(histos1,name,"mettrk_dphimet_calib_all", ";#Delta#phi(Calo #it{E}_{T}^{miss},Track #it{E}_{T}^{miss}) (all);Normalized",32,0.,TMath::Pi());
	
	addHist(histos1,name,"dPhi3bodyJ1_calib_all", ";#Delta#phi(3body,jet_{1}) (all);Events",32,0.,TMath::Pi());		
	addHist(histos1,name,"dPhiJ1J2_calib_all",    ";#Delta#phi(jet_{1},jet_{2}) (all);Events",32,0.,TMath::Pi());
	
	addHist(histos1,name,"mT_mu1_onPhi_beforeObjectcuts",       "m_{2body} on #phi, before Object cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu1_onRhoOmega_beforeObjectcuts",  "m_{2body} on #rho/#omega, before Object cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onPhi_beforeObjectcuts",       "m_{2body} on #phi, before Object cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onRhoOmega_beforeObjectcuts",  "m_{2body} on #rho/#omega, before Object cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onPhi_beforeObjectcuts",       "m_{2body} on #phi, before Object cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onRhoOmega_beforeObjectcuts",  "m_{2body} on #rho/#omega, before Object cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	
	addHist(histos1,name,"mT_mu1_onPhi_beforeWcuts",       "m_{2body} on #phi, before W cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu1_onRhoOmega_beforeWcuts",  "m_{2body} on #rho/#omega, before W cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onPhi_beforeWcuts",       "m_{2body} on #phi, before W cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onRhoOmega_beforeWcuts",  "m_{2body} on #rho/#omega, before W cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onPhi_beforeWcuts",       "m_{2body} on #phi, before W cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onRhoOmega_beforeWcuts",  "m_{2body} on #rho/#omega, before W cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);

	addHist(histos1,name,"mT_mu1_onPhi_afterLoosecuts",       "m_{2body} on #phi, after loose cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu1_onRhoOmega_afterLoosecuts",  "m_{2body} on #rho/#omega, after loose cuts;m_{T}(#mu_{1},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onPhi_afterLoosecuts",       "m_{2body} on #phi, after loose cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu2_onRhoOmega_afterLoosecuts",  "m_{2body} on #rho/#omega, after loose cuts;m_{T}(#mu_{2},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onPhi_afterLoosecuts",       "m_{2body} on #phi, after loose cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);
	addHist(histos1,name,"mT_mu3_onRhoOmega_afterLoosecuts",  "m_{2body} on #rho/#omega, after loose cuts;m_{T}(#mu_{3},#it{E}_{T}^{miss}/2) [MeV];Events", 150,0.,3000.);

	// addHist(histos1,name,"3rdTrk_onPhi_pt_afterLoosecuts",      "m_{2body} on #phi, after loose cuts;;Events", 150,0.,3000.);
	// addHist(histos1,name,"3rdTrk_onPhi_chi2dof_afterLoosecuts", "m_{2body} on #phi, after loose cuts;;Events", 150,0.,3000.);
	// addHist(histos1,name,"3rdTrk_onPhi_pval_afterLoosecuts",    "m_{2body} on #phi, after loose cuts;;Events", 150,0.,3000.);
	// addHist(histos1,name,"3rdTrk_onPhi_angle_afterLoosecuts",   "m_{2body} on #phi, after loose cuts;;Events", 150,0.,3000.);
	// addHist(histos1,name,"3rdTrk_onPhi_pt_vs_chi2dof",          "m_{2body} on #phi, after loose cuts;MVA score (all range);Normalized",40,-1.,1.);
	// addHist(histos1,name,"3rdTrk_onPhi_pt_vs_pval",             "m_{2body} on #phi, after loose cuts;MVA score (all range);Normalized",40,-1.,1.);



	addHist(histos1,name,"m3body_lin_zoom_onPhi_beforeObjectcuts",       "m_{2body} on #phi, before Object cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onPhi_beforeObjectcuts",      "m_{2body} on #phi, before Object cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_kaon_hypothesis_beforeObjectcuts", "m_{2body} on #phi, before Object cuts, Kaon hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_pion_hypothesis_beforeObjectcuts", "m_{2body} on #phi, before Object cuts, Pion hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_1600_beforeObjectcuts",                 "m_{2body} on #phi, before Object cuts;m_{3body} [MeV];Events", 50,1550.,1700.);
	addHist(histos1,name,"m3body_lin_zoom_onRhoOmega_beforeObjectcuts",  "m_{2body} on #rho/#omega, before Object cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onRhoOmega_beforeObjectcuts", "m_{2body} on #rho/#omega, before Object cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	
	addHist(histos1,name,"m3body_lin_zoom_onPhi_beforeWcuts",       "m_{2body} on #phi, before W cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onPhi_beforeWcuts",      "m_{2body} on #phi, before W cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_kaon_hypothesis_beforeWcuts", "m_{2body} on #phi, before W cuts, Kaon hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_pion_hypothesis_beforeWcuts", "m_{2body} on #phi, before W cuts, Pion hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_1600_beforeWcuts",                 "m_{2body} on #phi, before W cuts;m_{3body} [MeV];Events", 50,1550.,1700.);
	addHist(histos1,name,"m3body_lin_zoom_onRhoOmega_beforeWcuts",  "m_{2body} on #rho/#omega, before W cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onRhoOmega_beforeWcuts", "m_{2body} on #rho/#omega, before W cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);

	addHist(histos1,name,"m3body_lin_zoom_onPhi_afterLoosecuts",       "m_{2body} on #phi, after loose cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onPhi_afterLoosecuts",      "m_{2body} on #phi, after loose cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_kaon_hypothesis_afterLoosecuts", "m_{2body} on #phi, after loose cuts, Kaon hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_sigregion_onPhi_pion_hypothesis_afterLoosecuts", "m_{2body} on #phi, after loose cuts, Pion hypothesis;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos1,name,"m3body_1600_afterLoosecuts",                 "m_{2body} on #phi, after loose cuts;m_{3body} [MeV];Events", 50,1550.,1700.);
	addHist(histos1,name,"m3body_lin_zoom_onRhoOmega_afterLoosecuts",  "m_{2body} on #rho/#omega, after loose cuts;m_{3body} [MeV];Events", 100,800.,2800.);
	addHist(histos1,name,"m3body_sigregion_onRhoOmega_afterLoosecuts", "m_{2body} on #rho/#omega, after loose cuts;m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	
	
	vector<TString> names, xnames;
	names.push_back("muons");     xnames.push_back("(after object cuts)");
	names.push_back("triplet");   xnames.push_back("(after 3body cuts)");
	names.push_back("hadclean");  xnames.push_back("(after hadronic cuts)");
	names.push_back("met");       xnames.push_back("(after #it{E}_{T}^{miss} cuts)");
	names.push_back("ht");        xnames.push_back("(after Jets,#it{E}_{T}^{miss} & isolation cuts)");
	names.push_back("lowmassres");xnames.push_back("(after #rho/#omega,#phi cuts)");
	names.push_back("");          xnames.push_back("");
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="")   ? "_after_"+names[i] : "";
		TString thisxname = (xnames[i]!="") ? " "+xnames[i]      : "";
		
		addHist(histos1,name,"trk_eta1"+thisname, ";#eta^{trk1}"+thisxname+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_eta2"+thisname, ";#eta^{trk2}"+thisxname+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_eta3"+thisname, ";#eta^{trk3}"+thisxname+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"trk_phi1"+thisname, ";#phi^{trk1}"+thisxname+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_phi2"+thisname, ";#phi^{trk2}"+thisxname+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_phi3"+thisname, ";#phi^{trk3}"+thisxname+";Normalized", 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos1,name,"trk_pt1"+thisname,  ";p_{T}^{trk1} [MeV]"+thisxname+";Normalized", 50,2000.,102000.,false);
		addHist(histos1,name,"trk_pt2"+thisname,  ";p_{T}^{trk2} [MeV]"+thisxname+";Normalized", 50,2000.,42000., false);
		addHist(histos1,name,"trk_pt3"+thisname,  ";p_{T}^{trk3} [MeV]"+thisxname+";Normalized", 50,2000.,27000., false);
		
		addHist(histos1,name,"mu_pt"+thisname,  ";#mu p_{T} [MeV]"+thisxname+";Normalized", 120,2000.,122000.,false);
		addHist(histos1,name,"TPa_pt"+thisname, ";TPa p_{T} [MeV]"+thisxname+";Normalized", 60,2000.,62000.,false);
		addHist(histos1,name,"TPb_pt"+thisname, ";TPb p_{T} [MeV]"+thisxname+";Normalized", 60,2000.,32000.,false);
		addHist(histos1,name,"mu_eta"+thisname, ";#mu #eta"+thisxname+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"TPa_eta"+thisname,";TPa #eta"+thisxname+";Normalized", 64,-3.0,+3.0,false);
		addHist(histos1,name,"TPb_eta"+thisname,";TPb #eta"+thisxname+";Normalized", 64,-3.0,+3.0,false);
	
		addHist(histos1,name,"trks_fitprob"+thisname,  ";Tracks #prodp-value"+thisxname+";Normalized",200,0.,1.);
		addHist(histos1,name,"trks_pixdEdx"+thisname,  ";Tracks #SigmadE/dx (pixel)"+thisxname+";Normalized",50,0.,8.);

		addHist(histos1,name,"jet_n_calib"+thisname,  ";N_{jets}"+thisxname+";Normalized",5,0.,5.);
		addHist(histos1,name,"jet_b_calib"+thisname,  ";N_{b-jets}"+thisxname+";Normalized",5,0.,5.);
		addHist(histos1,name,"jet_pt_calib"+thisname,      ";Jets #it{p}_{T} [MeV]"+thisxname+";Events",35,30.*GeV2MeV,100.*GeV2MeV);
		addHist(histos1,name,"jet_pt1_calib"+thisname,     ";#it{p}_{T}(jet1) [MeV]"+thisxname+";Events",35,30.*GeV2MeV,100.*GeV2MeV);
		addHist(histos1,name,"jet_pt2_calib"+thisname,     ";#it{p}_{T}(jet2) [MeV]"+thisxname+";Events",35,30.*GeV2MeV,100.*GeV2MeV);
		addHist(histos1,name,"jet_ptJ1J2_calib"+thisname,  ";#Sigma#it{p}_{T}(jet1+jet2) [MeV]"+thisxname+";Events",50,0.,150.*GeV2MeV);
		addHist(histos1,name,"jet_E_calib"+thisname,       ";Jets #it{E} [MeV]"+thisxname+";Events",50,0.,100.*GeV2MeV);
		addHist(histos1,name,"jet_E1_calib"+thisname,      ";#it{E}(jet1) [MeV]"+thisxname+";Events",50,10.*GeV2MeV,150.*GeV2MeV);
		addHist(histos1,name,"jet_E2_calib"+thisname,      ";#it{E}(jet2) [MeV]"+thisxname+";Events",50,10.*GeV2MeV,150.*GeV2MeV);
		addHist(histos1,name,"jet_EJ1J2_calib"+thisname,   ";#Sigma#it{E}(jet1+jet2) [MeV]"+thisxname+";Events",50,0.,150.*GeV2MeV);
	
		addHist(histos1,name,"met_et_calib"+thisname,         ";#it{E}_{T}^{miss} [MeV]"+thisxname+";Events",50,0.,100.*GeV2MeV);
		addHist(histos1,name,"met_mt_et3mu_calib"+thisname,   ";m_{T}(3body,#it{E}_{T}^{miss}) [MeV]"+thisxname+";Events",60,0,120.*GeV2MeV);
		addHist(histos1,name,"met_dphi3mu_calib"+thisname,    ";#Delta#phi(3body,#it{E}_{T}^{miss})"+thisxname+";Events",32,0.,TMath::Pi());
      	addHist(histos1,name,"mettrk_dphimet_calib"+thisname, ";#Delta#phi(Calo #it{E}_{T}^{miss},Track #it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());

		addHist(histos1,name,"dPhi3bodyJ1_calib"+thisname, ";#Delta#phi(3body,jet_{1})"+thisxname+";Events",32,0.,TMath::Pi());		
		addHist(histos1,name,"dPhiJ1J2_calib"+thisname,    ";#Delta#phi(jet_{1},jet_{2})"+thisxname+";Events",32,0.,TMath::Pi());	


		addHist(histos1,name,"jet_n"+thisname,              ";N_{jets}"+thisxname+";Normalized",5,0.,5.);
		addHist(histos1,name,"jet_b"+thisname,              ";N_{b-jets}"+thisxname+";Normalized",5,0.,5.);
		addHist(histos1,name,"jet_3body_ht"+thisname,       ";#it{H}_{T}(3body,jets) [MeV]"+thisxname+";Normalized",50,0.,200.*GeV2MeV);
		addHist(histos1,name,"jet_ht"+thisname,             ";#it{H}_{T}(jets) [MeV]"+thisxname+";Normalized",50,0.,200.*GeV2MeV);
		
		addHist(histos1,name,"met_et"+thisname,       ";#it{E}_{T}^{miss} [MeV]"+thisxname+";Normalized",50,0.,100.*GeV2MeV);
		addHist(histos1,name,"met_mt_et3mu"+thisname, ";m_{T}(3body,#it{E}_{T}^{miss}) [MeV]"+thisxname+";Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"met_dphi3mu"+thisname,  ";#Delta#phi(3body,#it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());

		addHist(histos1,name,"mettrk_et"+thisname,       ";Track #it{E}_{T}^{miss} [MeV]"+thisxname+";Normalized",50,0.,100.*GeV2MeV);
		addHist(histos1,name,"mettrk_mt_et3mu"+thisname, ";m_{T}(3body,Track #it{E}_{T}^{miss}) [MeV]"+thisxname+";Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"mettrk_dphi3mu"+thisname,  ";#Delta#phi(3body,Track #it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		addHist(histos1,name,"mettrk_dphimet"+thisname,  ";#Delta#phi(Calo #it{E}_{T}^{miss},Track #it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		
		
		addHist(histos1,name,"ht_pt"+thisname,      ";#it{H}_{T}"+thisxname+" [MeV];Normalized",100,0.,100.*GeV2MeV);
		addHist(histos1,name,"ht_dphimet"+thisname, ";#Delta#phi(#it{H}_{T},#it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		addHist(histos1,name,"ht_dphimet_track"+thisname, ";#Delta#phi(#it{H}_{T},Track #it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		addHist(histos1,name,"ht_mT"+thisname,      ";#it{m}_{T}(#it{H}_{T},#it{E}_{T}^{miss})"+thisxname+" [MeV];Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"ht_mT_track"+thisname,";#it{m}_{T}(#it{H}_{T},Track #it{E}_{T}^{miss})"+thisxname+" [MeV];Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"ht_dR3body"+thisname, ";#DeltaR(#it{H}_{T},3body)"+thisxname+";Normalized",32,0.,2.*TMath::Pi());
		addHist(histos1,name,"ht_metsig"+thisname,  ";#it{E}_{T}^{miss}/#sqrt{#it{H}_{T}}"+thisxname+" [#sqrt{MeV}];Normalized",100,0.,600.);

		addHist(histos1,name,"ht_pt_calib"+thisname,      ";#it{H}_{T}"+thisxname+" [MeV];Normalized",100,0.,100.*GeV2MeV);
		addHist(histos1,name,"ht_dphimet_calib"+thisname, ";#Delta#phi(#it{H}_{T},#it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		addHist(histos1,name,"ht_dphimet_track_calib"+thisname, ";#Delta#phi(#it{H}_{T},Track #it{E}_{T}^{miss})"+thisxname+";Normalized",32,0.,TMath::Pi());
		addHist(histos1,name,"ht_mT_calib"+thisname,      ";#it{m}_{T}(#it{H}_{T},#it{E}_{T}^{miss})"+thisxname+" [MeV];Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"ht_mT_track_calib"+thisname,";#it{m}_{T}(#it{H}_{T},Track #it{E}_{T}^{miss})"+thisxname+" [MeV];Normalized",60,0,120.*GeV2MeV);
		addHist(histos1,name,"ht_dR3body_calib"+thisname, ";#DeltaR(#it{H}_{T},3body)"+thisxname+";Normalized",32,0.,2.*TMath::Pi());
		addHist(histos1,name,"ht_metsig_calib"+thisname,  ";#it{E}_{T}^{miss}/#sqrt{#it{H}_{T}}"+thisxname+" [#sqrt{MeV}];Normalized",100,0.,600.);
		
		addHist(histos1,name,"m3body"+thisname,           ";m_{3body} [MeV]"+thisxname+";Events", 50,0.,4000.);
		addHist(histos1,name,"m3body_lin"+thisname,       ";m_{3body} [MeV]"+thisxname+";Events", 50,500.,4000.);
		addHist(histos1,name,"m3body_lin_zoom"+thisname,  ";m_{3body} [MeV]"+thisxname+";Events", 100,800.,2800.);
		addHist(histos1,name,"m3body_sigregion"+thisname, ";m_{3body} [MeV]"+thisxname+";Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		addHist(histos1,name,"pT3body"+thisname,          ";p_{T}^{3body} [MeV]"+thisxname+";Normalized", 80,0.,100000.);
		addHist(histos1,name,"mOS"+thisname,              ";m_{OS} [MeV]"+thisxname+";Events", 100,0.,4000.);
		addHist(histos1,name,"mOS1"+thisname,             ";m_{OS1} [MeV]"+thisxname+";Events", 100,0.,4000.);
		addHist(histos1,name,"mOS2"+thisname,             ";m_{OS2} [MeV]"+thisxname+";Events", 100,0.,4000.);
		addHist(histos1,name,"mSS"+thisname,              ";m_{SS} [MeV]"+thisxname+";Events", 100,0.,4000.);
		
		addHist(histos1,name,"dROS1"+thisname,             ";#DeltaR_{OS1}"+thisxname+";Normalized", 500,0.,0.5);
		addHist(histos1,name,"dROS2"+thisname,             ";#DeltaR_{OS2}"+thisxname+";Normalized", 500,0.,0.5);
		addHist(histos1,name,"dRSS"+thisname,              ";#DeltaR_{SS}"+thisxname+";Normalized", 500,0.,0.5);
		addHist(histos2,name,"dROS1_vs_mOS1"+thisname,     "#DeltaR_{OS1} vs m_{OS1}"+thisxname+";m_{OS1} [MeV];#DeltaR_{OS1};Normalized", 100,0.,4000., 500,0.,0.5);
		addHist(histos2,name,"dROS2_vs_mOS2"+thisname,     "#DeltaR_{OS2} vs m_{OS2}"+thisxname+";m_{OS2} [MeV];#DeltaR_{OS2};Normalized", 100,0.,4000., 500,0.,0.5);
		addHist(histos2,name,"dRSS_vs_mSS"+thisname,       "#DeltaR_{SS} vs m_{SS}"+thisxname+";m_{SS} [MeV];#DeltaR_{SS};Normalized", 100,0.,4000., 500,0.,0.5);
		
		
		addHist(histos2,name,"mOS2_vs_mOS1"+thisname,     "m_{OS2} vs m_{OS1}"+thisxname+";m_{OS1} [MeV];m_{OS2} [MeV];Events", 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS"+thisname,       "m_{SS} vs m_{OS}"+thisxname+";m_{OS} [MeV];m_{SS} [MeV];Events", 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSSsq_vs_mOSsq"+thisname,   "m_{SS}^{2} vs m_{OS}^{2}"+thisxname+";m_{OS}^{2} [MeV^{2}];m_{SS}^{2} [MeV^{2}];Events", 100,0.,2000.*2000., 100,0.,2000.*2000.);
		addHist(histos2,name,"mSSsq_vs_mOS1sq"+thisname,  "m_{SS}^{2} vs m_{OS1}^{2}"+thisxname+";m_{OS1}^{2} [MeV^{2}];m_{SS}^{2} [MeV^{2}];Events", 100,0.,2000.*2000., 100,0.,2000.*2000.);
		addHist(histos2,name,"mSSsq_vs_mOS2sq"+thisname,  "m_{SS}^{2} vs m_{OS2}^{2}"+thisxname+";m_{OS2}^{2} [MeV^{2}];m_{SS}^{2} [MeV^{2}];Events", 100,0.,2000.*2000., 100,0.,2000.*2000.);
		addHist(histos2,name,"m3body_vs_mOS"+thisname,    "m_{3body} vs m_{OS}"+thisxname+";m_{OS} [MeV];m_{3body} [MeV];Events", 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"pT3body_vs_mOS"+thisname,   "p_{T}^{3body} vs m_{OS}"+thisxname+";m_{OS} [MeV];p_{T}^{3body} [MeV];Events", 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"3rdTrk_pt_vs_mOS"+thisname,        "p_{T}(3rd trk) vs m_{OS}"+thisxname+";m_{OS} [MeV];p_{T}(3rd trk) [MeV];Events", 50,0.,2000., 50,0.,50000.);
		addHist(histos2,name,"3rdTrk_chi2dof_vs_mOS"+thisname,   "#chi^{2}_{DOF}(3rd trk) vs m_{OS}"+thisxname+";m_{OS} [MeV];#chi^{2}_{DOF};Events", 50,0.,2000., 50,0.,10.);
		addHist(histos2,name,"3rdTrk_pval_vs_mOS"+thisname,      "p-value(3rd trk) vs m_{OS}"+thisxname+";m_{OS} [MeV];Trk p-value;Events", 50,0.,2000., 50,0.,0.001);
		addHist(histos2,name,"3rdTrk_angle_vs_mOS"+thisname,     "Angle(3rd trk) vs m_{OS}"+thisxname+";m_{OS} [MeV];Trk angle;Events", 50,0.,2000., 50,0.,TMath::Pi());
		addHist(histos2,name,"3rdTrk_PtOverPtOS_vs_mOS"+thisname,"p_{T}(3rd trk)/p_{T}^{OS} vs m_{OS}"+thisxname+";m_{OS} [MeV];p_{T}(3rd trk)/p_{T}^{OS};Events", 50,0.,2000., 50,0.,2.);

		
		addHist(histos2,name,"dPhiJet1Jet2_vs_sumpTjet12"+thisname,  "#Delta#phi(jet_{1},jet_{2}) vs #Sigmap_{T} jet_{1,2}"+thisxname+";#Sigmap_{T} jet_{1,2} [MeV] ;#Delta#phi(jet_{1},jet_{2});Normalized",20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhi3bodyJet1_vs_pTjet1"+thisname,     "#Delta#phi(3body,jet_{1}) vs p_{T} jet_{1}"+thisxname+";p_{T} jet_{1} [MeV];#Delta#phi(3body,jet_{1});Normalized",20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos1,name,"dPhi3bodyJ1"+thisname, ";#Delta#phi(3body,jet_{1})"+thisxname+";Normalized",32,0.,TMath::Pi());		
		addHist(histos1,name,"dPhiJ1J2"+thisname,    ";#Delta#phi(jet_{1},jet_{2})"+thisxname+";Normalized",32,0.,TMath::Pi());	
		addHist(histos1,name,"dPhiMETJ1"+thisname,    ";#Delta#phi(#it{E}_{T}^{miss},jet_{1})"+thisxname+";Normalized",32,0.,TMath::Pi());	
		addHist(histos1,name,"mJ1J2"+thisname,        ";m(jet_{1},jet_{2})"+thisxname+";Normalized",50,0.,150.*GeV2MeV);
		addHist(histos2,name,"MV1_vs_mjet"+thisname, "MV1 flavor weight vs #it{m}_{jet}"+thisxname+";#it{m}_{jet} [MeV];MV1 flavor weight;Normalized",   20,0,20*GeV2MeV, 20,0.,1.);
		addHist(histos2,name,"MV1_vs_pTjet"+thisname,"MV1 flavor weight vs #it{pT}_{jet}"+thisxname+";#it{pT}_{jet} [MeV];MV1 flavor weight;Normalized", 20,0,100*GeV2MeV, 20,0.,1.);			
		
		addHist(histos1,name,"JVFall"+thisname, ";JVFall"+thisxname+";Normalized",50,0.,1.);
		addHist(histos1,name,"JVF1"+thisname, ";JVF1"+thisxname+";Normalized",50,0.,1.);
		addHist(histos1,name,"JVF2"+thisname, ";JVF2"+thisxname+";Normalized",50,0.,1.);
		addHist(histos1,name,"JVF3"+thisname, ";JVF3"+thisxname+";Normalized",50,0.,1.);
		addHist(histos1,name,"JVF4"+thisname, ";JVF4"+thisxname+";Normalized",50,0.,1.);
		
		addHist(histos2,name,"JVF1_vs_pTjet1"+thisname, "JVF1 vs p_{T} jet_{1}"+thisxname+";p_{T} jet_{1} [MeV];JVF1;Normalized",20,0,100*GeV2MeV ,20,0.,1.);
		addHist(histos2,name,"JVF2_vs_pTjet2"+thisname, "JVF2 vs p_{T} jet_{2}"+thisxname+";p_{T} jet_{2} [MeV];JVF2;Normalized",20,0,100*GeV2MeV ,20,0.,1.);
	}
	
	addHist(histos1,name,"isolation000",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.00);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation001",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.01);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation002",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.02);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation003",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.03);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation004",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.04);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation005",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.05);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation006",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.06);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation007",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.07);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation008",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.08);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation009",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.09);Normalized", 100,0.,1.0);
	addHist(histos1,name,"isolation010",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.10);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation012",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.12);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation014",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.14);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation016",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.16);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation018",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.18);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation020",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.20);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation022",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.22);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation024",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.24);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation026",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.26);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation028",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.28);Normalized", 100,0.,2.0);
	addHist(histos1,name,"isolation030",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+0.30);Normalized", 100,0.,2.0);
	
	addHist(histos1,name,"dRmin", ";min[#DeltaR(3body,obj_{i})];Normalized",  100,0.,0.25);
	addHist(histos1,name,"dRmax", ";max[#DeltaR(3body,obj_{i})];Normalized",  100,0.,0.50);	
	
	addHist(histos1,name,"MVAscore_all",            ";MVA score (all range);Normalized",40,-1.,1.);
	addHist(histos1,name,"MVAscore_below_tauMass",  ";MVA score (below the #tau mass);Normalized",40,-1.,1.);
	addHist(histos1,name,"MVAscore_above_tauMass",  ";MVA score (above the #tau mass);Normalized",40,-1.,1.);
	addHist(histos1,name,"MVAscore_left_sideband",  ";MVA score (left sideband);Normalized",40,-1.,1.);
	addHist(histos1,name,"MVAscore_right_sideband", ";MVA score (right sideband);Normalized",40,-1.,1.);
	addHist(histos2,name,"m3body_vs_MVAscore",      "m_{3body} vs MVA score;MVA score;m_{3body} [MeV];Normalized",         30,-1.,+1., (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	addHist(histos2,name,"pT3body_vs_MVAscore",     "p_{T}^{3body} vs MVA score;MVA score;p_{T}^{3body} [MeV];Normalized", 30,-1.,+1., 40,20.*GeV2MeV,100.*GeV2MeV);
	addHist(histos2,name,"mOS_vs_MVAscore",         "m_{OS} vs MVA score;MVA score;m_{OS} [MeV];Normalized",               30,-1.,+1., 100,0.,2500.);
	addHist(histos2,name,"mOS1_vs_MVAscore",        "m_{OS1} vs MVA score;MVA score;m_{OS1} [MeV];Normalized",             30,-1.,+1., 100,0.,2500.);
	addHist(histos2,name,"mOS2_vs_MVAscore",        "m_{OS2} vs MVA score;MVA score;m_{OS2} [MeV];Normalized",             30,-1.,+1., 100,0.,2500.);
	addHist(histos2,name,"mSS_vs_MVAscore",         "m_{SS}  vs MVA score;MVA score;m_{SS} [MeV];Normalized",              30,-1.,+1., 100,0.,2500.);


	vector<TString> mvaevo;
	mvaevo.push_back("neg_veryloose");
	mvaevo.push_back("neg_loose");
	mvaevo.push_back("neg_medium");
	mvaevo.push_back("neg_tight");
	mvaevo.push_back("zero");
	mvaevo.push_back("pos_loose");
	mvaevo.push_back("pos_medium");
	mvaevo.push_back("pos_tight");
	mvaevo.push_back("pos_verytight");	
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
		
		addHist(histos1,name,"MVAevo_score"+suff, mvaevo[i]+";MVA score;Events",40,-1.,1.);
		
		addHist(histos1,name,"MVAevo_m3body_sigregion_onPhi_kaon_hypothesis"+suff, mvaevo[i]+";m_{3body} [MeV] (on #phi, kaon hypo);Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		addHist(histos1,name,"MVAevo_m3body_sigregion_onPhi_pion_hypothesis"+suff, mvaevo[i]+";m_{3body} [MeV] (on #phi, pion hypo);Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		
		addHist(histos1,name,"MVAevo_m3body_sigregion"+suff, mvaevo[i]+";m_{3body} [MeV];Events", (Int_t)((mSideBandRightUpperMeVGlob-mSideBandLeftLowerMeVGlob)/mBinSize),mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		addHist(histos1,name,"MVAevo_pT3body"+suff,          mvaevo[i]+";p_{T}^{3body} [MeV];Events", 100,0.,100000.);
		addHist(histos1,name,"MVAevo_mOS1"+suff,             mvaevo[i]+";m_{OS1} [MeV];Events", 50,0.,2000.);
		addHist(histos1,name,"MVAevo_mOS2"+suff,             mvaevo[i]+";m_{OS2} [MeV];Events", 50,0.,2000.);
		addHist(histos1,name,"MVAevo_mSS"+suff,              mvaevo[i]+";m_{SS} [MeV];Events",  50,0.,2000.);
		
		addHist(histos1,name,"MVAevo_met_cal_et"+suff,         mvaevo[i]+";#it{E}_{T}^{miss} [MeV];Events", 50,0.,100.*GeV2MeV);
		addHist(histos1,name,"MVAevo_met_track_et"+suff,       mvaevo[i]+";Track #it{E}_{T}^{miss} [MeV];Events", 50,0.,100.*GeV2MeV);
		addHist(histos1,name,"MVAevo_met_cal_mT"+suff,         mvaevo[i]+";m_{T}(3body,#it{E}_{T}^{miss}) [MeV];Events", 60,0,120.*GeV2MeV);
		addHist(histos1,name,"MVAevo_met_track_mT"+suff,       mvaevo[i]+";m_{T}(3body,Track #it{E}_{T}^{miss}) [MeV];Events", 60,0,120.*GeV2MeV);
		addHist(histos1,name,"MVAevo_met_cal_dPhi3mu"+suff,    mvaevo[i]+";#Delta#phi(3body,#it{E}_{T}^{miss});Events", 32,0.,TMath::Pi());
		addHist(histos1,name,"MVAevo_met_track_dPhi3mu"+suff,  mvaevo[i]+";#Delta#phi(3body,Track #it{E}_{T}^{miss});Events", 32,0.,TMath::Pi());
		addHist(histos1,name,"MVAevo_met_dphiCalTrk"+suff,     mvaevo[i]+";#Delta#phi(#it{E}_{T}^{miss},Track #it{E}_{T}^{miss});Events", 32,0.,TMath::Pi());

		addHist(histos1,name,"MVAevo_jet_n"+suff,              mvaevo[i]+";N_{jets};Events", 5,0.,5.);
		addHist(histos1,name,"MVAevo_jet_ht_pt"+suff,          mvaevo[i]+";#it{H}_{T} [MeV];Events", 100,0.,100.*GeV2MeV);
		addHist(histos1,name,"MVAevo_jet_ht_mT_metcal"+suff,   mvaevo[i]+";m_{T}(#it{H}_{T},#it{E}_{T}^{miss}) [MeV];Events", 60,0,120.*GeV2MeV);
		addHist(histos1,name,"MVAevo_jet_ht_mT_mettrk"+suff,   mvaevo[i]+";m_{T}(#it{H}_{T},Track #it{E}_{T}^{miss}) [MeV];Events", 60,0,120.*GeV2MeV);
		addHist(histos1,name,"MVAevo_jet_ht_dphimetcal"+suff,  mvaevo[i]+";#Delta#phi(#it{H}_{T},#it{E}_{T}^{miss});Events", 32,0.,TMath::Pi());
		addHist(histos1,name,"MVAevo_jet_ht_dphimettrk"+suff,  mvaevo[i]+";#Delta#phi(#it{H}_{T},Track #it{E}_{T}^{miss});Events", 32,0.,TMath::Pi());
		addHist(histos1,name,"MVAevo_jet_ht_dr3body"+suff,     mvaevo[i]+";#DeltaR(#it{H}_{T},3body);Events", 32,0.,2.*TMath::Pi());
		
		addHist(histos2,name,"MVAevo_mOS2_vs_mOS1"+suff,     "m_{OS2} vs m_{OS1}"+mvaevo[i]+";m_{OS1} [MeV];m_{OS2} [MeV];Events", 50,0.,2000., 50,0.,2000.);
		addHist(histos2,name,"MVAevo_mSS_vs_mOS1"+suff,      "m_{SS} vs m_{OS1}"+mvaevo[i]+";m_{OS1} [MeV];m_{SS} [MeV];Events",   50,0.,2000., 50,0.,2000.);
		addHist(histos2,name,"MVAevo_mSS_vs_mOS2"+suff,      "m_{SS} vs m_{OS2}"+mvaevo[i]+";m_{OS2} [MeV];m_{SS} [MeV];Events",   50,0.,2000., 50,0.,2000.);
		addHist(histos2,name,"MVAevo_m3body_vs_mOS1"+suff,   "m_{3body} vs m_{OS1}"+mvaevo[i]+";m_{OS1} [MeV];m_{3body} [MeV];Events", 50,0.,2000., 60,mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		addHist(histos2,name,"MVAevo_m3body_vs_mOS2"+suff,   "m_{3body} vs m_{OS2}"+mvaevo[i]+";m_{OS2} [MeV];m_{3body} [MeV];Events", 50,0.,2000., 60,mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
		addHist(histos2,name,"MVAevo_m3body_vs_mSS"+suff,    "m_{3body} vs m_{SS}"+mvaevo[i]+";m_{SS} [MeV];m_{3body} [MeV];Events",   50,0.,2000., 60,mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob);
	}
}


void reBlindAllMassHists(TMapTSP2TH1& histos, float mMin, float mMax)
{	
	for(TMapTSP2TH1::iterator it=histos.begin() ; it!=histos.end() ; ++it)
	{	
		TString name = it->first;
		if(!isData(name))            continue;
		if(!name.Contains("m3body")) continue;
		
		if(name.Contains("onPhi"))      continue;
		if(name.Contains("onRhoOmega")) continue;
		
		Int_t bMin = it->second->FindBin(mMin);
		Int_t bMax = it->second->FindBin(mMax);
		for(Int_t b=bMin ; b<=bMax && b>0 ; b++)
		{	
			if(it->second->GetBinContent(b)>0)
			{
				it->second->SetBinContent(b,0);
				it->second->SetBinError(b,0);
			}
		}
	}
}
void reBlindAllMassHists(TMapTSP2TH2& histos, float mMin, float mMax)
{	
	for(TMapTSP2TH2::iterator it=histos.begin() ; it!=histos.end() ; ++it)
	{
		TString name = it->first;
		if(!isData(name))               continue;
		if(!name.Contains("m3body"))    continue;
		if(name.Contains("onPhi"))      continue;
		if(name.Contains("onRhoOmega")) continue;
		
		TString xtitle = it->second->GetXaxis()->GetTitle();
		TString ytitle = it->second->GetYaxis()->GetTitle();
		Bool_t isX = (xtitle.Contains("m_{3body}"));
		if(isX)
		{
			Int_t bMin = it->second->GetXaxis()->FindBin(mMin);
			Int_t bMax = it->second->GetXaxis()->FindBin(mMax);
			Int_t nY = it->second->GetNbinsY();
			for(Int_t b=bMin ; b<=bMax && b>0 ; b++)
			{
				for(Int_t y=1 ; y<=nY ; y++)
				{
					if(it->second->GetBinContent(b,y)>0)
					{
						it->second->SetBinContent(b,y,0);
						it->second->SetBinError(b,y,0);
					}
				}
			}
		}
		else
		{
			Int_t bMin = it->second->GetYaxis()->FindBin(mMin);
			Int_t bMax = it->second->GetYaxis()->FindBin(mMax);
			Int_t nX = it->second->GetNbinsX();
			for(Int_t b=bMin ; b<=bMax && b>0 ; b++)
			{
				for(Int_t x=1 ; x<=nX ; x++)
				{
					if(it->second->GetBinContent(x,b)>0)
					{
						it->second->SetBinContent(x,b,0);
						it->second->SetBinError(x,b,0);
					}
				}
			}
		}
	}
}
bool reBlind(TH1* h, float m, float mMin, float mMax)
{
	TString name = h->GetName();
	if(!isData(name)) return false;
	Int_t bMin = h->FindBin(mMin);
	Int_t bMax = h->FindBin(mMax);
	float xMin = h->GetBinLowEdge(bMin);
	float xMax = h->GetBinLowEdge(bMax)+h->GetBinWidth(bMax);
	if(m>=xMin && m<xMax) return true;
	return false;
}


void finalizeHistos(TString master, bool isBlinded=true, TString channel="all", bool doMVA=false)
{
	olddir->cd();
	
	///////////////////////////////////////
	//// reblind all mass histos1 1d and 2d
	if(isBlinded) reBlindAllMassHists(histos1, mBlindMinGlob, mBlindMaxGlob);
	if(isBlinded) reBlindAllMassHists(histos2, mBlindMinGlob, mBlindMaxGlob);
	
	// histos1 decoration
	for(TMapTSP2TH1::iterator hit=histos1.begin() ; hit!=histos1.end() ; hit++)
	{
		TString hname = hit->first;
		////////////////////////
		// normalize to unity //
		////////////////////////
		TString ytitle = hit->second->GetYaxis()->GetTitle();
		if(ytitle=="Normalized")            NormToEntries(hit->second);
		if(ytitle=="Normalized rate")       NormToEntries(hit->second);
		if(ytitle=="Normalized to 1st bin") NormTo1stBin(hit->second);
		
		// set bin labels:
		// if(hname.EndsWith("mu_type"))
		if(hname.Contains("isTight"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isTight");
			hit->second->GetXaxis()->SetBinLabel(2,"isTight");
		}
		if(hname.Contains("isMedium"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isMedium");
			hit->second->GetXaxis()->SetBinLabel(2,"isMedium");
		}
		if(hname.Contains("isLoose"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isLoose");
			hit->second->GetXaxis()->SetBinLabel(2,"isLoose");
		}
		if(hname.Contains("_type_") || hname.EndsWith("mu_type"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"CB muon");
			hit->second->GetXaxis()->SetBinLabel(2,"TPmuonA");
			hit->second->GetXaxis()->SetBinLabel(3,"TPmuonB");
			hit->second->GetXaxis()->SetBinLabel(4,"CaloMuon");
		}
		if(hname.EndsWith("mu_MCP"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"Fail");
			hit->second->GetXaxis()->SetBinLabel(2,"Pass");
		}
		
		// if(hname.Contains("triggers"))
		// {
		// 	buildTriggerbits();
		// 	for(TMapuiTS::iterator it=triggerorder.begin() ; it!=triggerorder.end() ; ++it)
		// 	{
		// 		unsigned int bin = it->first;
		// 		TString name      = it->second;
		// 		hit->second->GetXaxis()->SetBinLabel(bin,name);
		// 	}
		// }
		
		// styles, colors and min/max
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{	
			TString name  = cit->second;
			Color_t color = colors[name];
			Int_t pattern = patterns[name];
			if(hname.Contains(name))
			{
				// styles, colors
				hit->second->SetLineColor(color);
				hit->second->SetFillColor(color);
				hit->second->SetMarkerColor(color);
				hit->second->SetFillStyle(pattern);
			}
		}
	}
	
	for(TMapTSP2TH2::iterator hit=histos2.begin() ; hit!=histos2.end() ; hit++)
	{
		TString hname = hit->first;
		////////////////////////
		// normalize to unity //
		////////////////////////
		TString ztitle = hit->second->GetZaxis()->GetTitle();
		if(ztitle.Contains("Normalized"))  NormToEntries(hit->second);
		if(ztitle.Contains("Correlation"))
		{
			TString h1name = hname;
			h1name = h1name.ReplaceAll("triggers_correlation_normAB","triggers_absoluteEvents");
			NormToTriggers(hit->second,histos1[h1name]);
		}
		
		// set bin labels:
		if(hname.Contains("_type_") || hname.EndsWith("mu_type"))
		{
			hit->second->GetYaxis()->SetBinLabel(1,"CB muon");
			hit->second->GetYaxis()->SetBinLabel(2,"TPmuonA");
			hit->second->GetYaxis()->SetBinLabel(3,"TPmuonB");
			hit->second->GetYaxis()->SetBinLabel(4,"CaloMuon");
		}

		// styles, colors
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{
			TString name  = cit->second;
			Color_t color = colors[name];
			if(hname.Contains(name))
			{
				// styles, colors
				hit->second->SetLineColor(color);
				hit->second->SetMarkerColor(color);
				hit->second->SetFillColor(color);
			}
		}
	}
	
	olddir->cd();
	TString hfilename = (channel!="all") ? channel+".histos."+master : "histos."+master;
	hfilename += (doMVA) ? ".mva.root" : ".cuts.root";
	TFile* hfile = new TFile(hfilename,"RECREATE");
	for(TMapTSP2TH1::iterator it=histos1.begin() ; it!=histos1.end() ; ++it) { if(it->second) it->second->Write(); else _ERR(1,"Problem with: "<<it->first); }
	for(TMapTSP2TH2::iterator it=histos2.begin() ; it!=histos2.end() ; ++it) { if(it->second) it->second->Write(); else _ERR(1,"Problem with: "<<it->first); }
	hfile->Write();
	hfile->Close();
	delete hfile;
	olddir->cd();
	_INF(1,"Written: "<<hfilename);
}

#endif
