#include <iostream>
#include <stdlib.h>
#include <stdio.h>      // for the sprintf call
#include <string>
#include <sstream>      // for the int to string operation (stringstream call)
#include <cstring>      // for the string functions
#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <time.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TMatrix.h>
#include <TVector.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TApplication.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TThread.h>
#include <TList.h>
#include <TGraphAsymmErrors.h>

using namespace std;

// Declare pointer to data as global (not elegant but TMinuit needs this).
vector<vector<double> *> VVCOSTH;
vector<TH1D*> VHIST_Z0;
vector<TH1D*> VHIST_DATA;
vector<TH1D*> VHIST_ZPRIME;
vector<TVirtualPad*> VPAD;
TCanvas* CNV = new TCanvas("C", "C", 1200, 800);
unsigned int CURRENTBIN;
bool ISZMUMU;
bool ISZPRIME;
bool ISDATA;
string REFNAME = "CosThetaCS";
//string REFNAME = "CosThetaHE";
Float_t IMASS;
Float_t COSTH;
vector<float>* CHARGE;
vector<float>* PX;
vector<float>* PY;
vector<float>* PZ;
vector<float>* E;
TLorentzVector* pa = new TLorentzVector();
TLorentzVector* pb = new TLorentzVector();
TLorentzVector m_pTmp;
float GeV2TeV = 1.e-3;
float MeV2TeV = 1.e-6;
float muonMass = 0.105658367; // GeV


///////////////////////////////////////////////////////////////////////
double error_poisson_up(double data)
{
	double y1 = data + 1.0;
	double d = 1.0 - 1.0/(9.0*y1) + 1.0/(3*TMath::Sqrt(y1));
	return y1*d*d*d-data;
}

double error_poisson_down(double data)
{
	double y = data;
	if (y == 0.0) return 0.0;
	double d = 1.0 - 1.0/(9.0*y) - 1.0/(3.0*TMath::Sqrt(y));
	return data-y*d*d*d;
}

TGraphAsymmErrors* GetPoissonizedGraph(TH1D* histo, bool isXerr)
{
	TGraphAsymmErrors* graph = new TGraphAsymmErrors();
	int j=0;
	for (int i=1;i<=histo->GetNbinsX();++i)
	{
		if (histo->GetBinContent(i)!=0)
		{ 
			graph->SetPoint(j,histo->GetBinCenter(i),histo->GetBinContent(i));
			graph->SetPointError(
								 j,
								 (isXerr) ? histo->GetBinWidth(i)/2. : 0.,
								 (isXerr) ? histo->GetBinWidth(i)/2. : 0.,
								 error_poisson_down(histo->GetBinContent(i)),
								 error_poisson_up(histo->GetBinContent(i))
								);
			++j;
		}
	}
	return graph;
}
///////////////////////////////////////////////////////////////////////



void Scale(TH1D* h, double d)
{ 
	/// scale including over/underflow
	for(int i=0 ; i<=h->GetNbinsX()+1 ; i++)
	{ 
		h->SetBinContent(i,h->GetBinContent(i)*d);
	}
}

void ScaleWerrors(TH1D* h, double d)
{ 
	/// scale including over/underflow
	for ( int i=0 ; i<=h->GetNbinsX()+1 ; i++ )
	{ 
		h->SetBinContent(i,h->GetBinContent(i)*d);
		h->SetBinError(i,h->GetBinError(i)*d);
	}
}

inline float imass( TLorentzVector* pa, TLorentzVector* pb )
{
	m_pTmp = (*pa)+(*pb);
	return m_pTmp.M();
}

inline float QT( TLorentzVector* pa, TLorentzVector* pb )
{
	m_pTmp = (*pa)+(*pb);
	return m_pTmp.Perp();
}

inline float ySystem( TLorentzVector* pa, TLorentzVector* pb )
{
	m_pTmp = (*pa)+(*pb);
	return m_pTmp.Rapidity();
}

inline float cosThetaBoost( TLorentzVector* pa, float ca, TLorentzVector* pb, float cb )
{
	// http://xrootd.slac.stanford.edu/BFROOT/www/doc/workbook_backup_010108/analysis/analysis.html
	// A useful quantity in many analyses is the helicity angle.
	// In the reaction Y -> X -> a + b, the helicity angle of 
	// particle a is the angle measured in the rest frame of the
	//decaying parent particle, X, between the direction of the
	// decay daughter a and the direction of the grandparent particle Y.

	m_pTmp = (*pa)+(*pb); // this is the mumu system (Z) 4vector
	TVector3 ZboostVector = m_pTmp.BoostVector(); // this is the 3vector of the Z
	TLorentzVector p; // this is the muon 4vector
	
	if(ca<0)      p.SetPxPyPzE(pa->Px(),pa->Py(),pa->Pz(),pa->E());
	else if(cb<0) p.SetPxPyPzE(pb->Px(),pb->Py(),pb->Pz(),pb->E());
	p.Boost( -ZboostVector ); // boost p to the dimuon CM (rest) frame
	float cosThetaB = p.Vect()*m_pTmp.Vect()/(p.P()*m_pTmp.P());
	//if (ySystem(pa,pb) < 0) cosThetaB *= -1.; // reclassify ???
	return cosThetaB;
}

inline float cosThetaCollinsSoper( TLorentzVector* pa, float ca, TLorentzVector* pb, float cb )
{
	// this will work only for leptons e, mu and tau
	// by default it is assumed that pa is the lepton
	// if instead pb is the lepton, then the result is
	// reclassified by a (-) sign - see line 4.
	float mass2 = imass(pa,pb)*imass(pa,pb);
	float QT2 = QT(pa,pb)*QT(pa,pb);
	//float cosThetaCS = 2.*( pa->Plus()*pb->Minus() - pa->Minus()*pb->Plus() ) / sqrt( mass2 * (mass2 + QT2) );
	float cosThetaCS = 2.*( pa->Pz()*pb->E() - pa->E()*pb->Pz() ) / sqrt( mass2 * (mass2 + QT2) );
	if (ca>0. && cb<0.)     cosThetaCS *= -1.; // if pb is the lepton
	if (ySystem(pa,pb) < 0) cosThetaCS *= -1.; // reclassify
	return cosThetaCS;
}


// The pdf to be fitted
// First argument needs to be a pointer in order to plot with the TF1 class.
double cosThetaPdf(double* xPtr, double par[])
{
	double x   = *xPtr;
	double pdf = 1.;
	double A4  = par[0];

	if (fabs(x)<=1.)
	{
		// see: http://arxiv.org/PS_cache/arxiv/pdf/1004/1004.1649v1.pdf
		// pdf_CollinsSoper = (1/N)*dN/dcosThetaCS = (3/8)*((1+0.5*A0) + A4*cosThetaCS + (1-3/2*A0)*cosTetaCS^2)
		// pdf_Helicity     = (1/N)*dN/dcosThetaHE = (3/8)*(1+cosTetaHE^2 + Afb*cosThetaHE)
		// if A0=~0 then pdf_CollinsSoper = pdf_Helicity (in the Z/Z'/Z* cases, it is approximately true)
		pdf  = (3./8.)*(1. + x*x + A4*x);
		if (pdf > 0) { return pdf; }
		else { cout << "warning:  pdf<=0, pdf=" << pdf << endl; return 1e-30; }
	}
	else { cout << "warning:  x<-1 || x>+1, x=" << x << endl; return 1e-30; }
}

// fcn passes back f = - 2*ln(L), the function to be minimized.
void fcn(int& npar, double* deriv, double& f, double par[], int flag)
{
	//vector<double> xVec = *VCOSTH; // VCOSTH is global
	vector<double> xVec = *VVCOSTH[CURRENTBIN]; // VVCOSTH is global
	int n = xVec.size();
	double lnL = 0.;
	double x   = 0.;
	double pdf = 0.;
	for (int i=0; i<n; i++){
		x = xVec[i];
		pdf = cosThetaPdf(&x, par);
		if (pdf > 0.) {
			lnL += log(pdf); // need positive f
		}
		else { cout << "WARNING -- pdf is negative!!!" << endl; }
	}
	f = -2.0 * lnL; // factor of -2 so minuit gets the errors right
}

void fillVec(TTree* t, TH1D* h)
{
	for(int v=0 ; v<(int)VVCOSTH.size() ; v++) VVCOSTH[v]->clear();
	if(t==0) return;
	
	TAxis* xaxis = h->GetXaxis();
	
	for (Long64_t l64t_jentry=0 ; l64t_jentry<t->GetEntries() ; l64t_jentry++)
	{
		t->GetEntry(l64t_jentry);
		int bin = (int)xaxis->FindBin((Double_t)IMASS);
		if(bin<=0 || bin>(int)VVCOSTH.size()) continue;
		if(CHARGE->at(0)*CHARGE->at(1)>=0)   continue;
		pa->SetPxPyPzE(PX->at(0)*MeV2TeV,PY->at(0)*MeV2TeV,PZ->at(0)*MeV2TeV,E->at(0)*MeV2TeV);
		pb->SetPxPyPzE(PX->at(1)*MeV2TeV,PY->at(1)*MeV2TeV,PZ->at(1)*MeV2TeV,E->at(1)*MeV2TeV);
		if(REFNAME=="CosThetaCS") COSTH = cosThetaCollinsSoper( pa, CHARGE->at(0), pb, CHARGE->at(1) );
		else COSTH = 0.;
		
		VVCOSTH[bin-1]->push_back(COSTH);
	
		/*	
		if(ISZMUMU)  VHIST_Z0[bin-1]->Fill(COSTH);
		if(ISZPRIME) VHIST_ZPRIME[bin-1]->Fill(COSTH);
		if(ISDATA)   VHIST_DATA[bin-1]->Fill(COSTH);
		*/
	}
}

void minimize(double guess, double& A4, double& dA4)
{
	int npar = 1;      // the number of parameters
	TMinuit minuit(npar);
	minuit.SetFCN(fcn);

	double par[npar];        // the start values
	double stepSize[npar];   // step sizes
	double minVal[npar];     // minimum bound on parameter
	double maxVal[npar];     // maximum bound on parameter
	string parName[npar];

	par[0]      = guess; 
	stepSize[0] = 1e-4;
	minVal[0]   = -1.;
	maxVal[0]   = +1.;
	parName[0]  = "A4";

	for (int i=0; i<npar; i++)
	{
		minuit.DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);
	}

	// Do the minimization!
	minuit.Migrad(); // Minuit's best minimization algorithm
	double outpar[npar], err[npar];
	for (int i=0; i<npar; i++)
	{
		minuit.GetParameter(i,outpar[i],err[i]);
	}
	
	/*
	TF1* func = new TF1("pdf", cosThetaPdf, -1., 1., npar);
	func->SetParameters(outpar);
	func->SetLineStyle(1);
	func->SetLineColor(1);
	func->SetLineWidth(1);
	func->GetXaxis()->SetTitle("cos#left(#theta^{*}#right)");
	func->GetYaxis()->SetTitle("pdf");
	*/
	
	A4  = 0.;
	dA4 = 0.;
	minuit.GetParameter(0,A4,dA4);
}


void execute(string isHistos = "")
{
	/*
	if(isHistos=="only_histos")
	{
		gROOT->ProcessLine(".x rootlogon_atlas.C");
		gROOT->SetStyle("ATLAS");
		gROOT->ForceStyle();
	}
	*/


	int minEntriesDATA = 1;//10;
	int minEntriesMC   = 1;//10;
	string refframe = REFNAME;
	bool doLogM = true;
	bool doLogx = true;
	
	
	/*
	// logarithmic boundries and bins of histograms
	const int imass_nbins = 12;
	const int ncol_pads   = 2; // = imass_nbins/nrow_pads !!!
	const int nrow_pads   = 6; // = imass_nbins/ncol_pads !!!
	double imass_min      = 70.*GeV2TeV;
	double imass_max      = 400.*GeV2TeV;
	Double_t logMmin = log10(imass_min);
	Double_t logMmax = log10(imass_max);
	Double_t imass_bins[imass_nbins+1];
	Double_t M_binwidth = (Double_t)( (logMmax-logMmin)/(Double_t)imass_nbins );
	imass_bins[0] = imass_min;
	for(Int_t i=1 ; i<=imass_nbins ; i++) imass_bins[i] = TMath::Power( 10,(logMmin + i*M_binwidth) );
	*/
	
	Double_t logMmin;
	Double_t logMmax;
	Double_t logMbinwidth;

	const double imass_min_const = 75.*GeV2TeV; //72.62*GeV2TeV
	const double imass_max_const = 400.*GeV2TeV; //381.09*GeV2TeV


	const int imass_asymmetry_nbins = 7; //14;
	double    imass_asymmetry_min   = imass_min_const;
	double    imass_asymmetry_max   = imass_max_const;
	logMmin     = log10(imass_asymmetry_min);
	logMmax     = log10(imass_asymmetry_max);
	Double_t imass_asymmetry_bins[imass_asymmetry_nbins+1];
	logMbinwidth = (Double_t)( (logMmax-logMmin)/(Double_t)imass_asymmetry_nbins );
	imass_asymmetry_bins[0] = imass_asymmetry_min;
	for(Int_t i=1 ; i<=imass_asymmetry_nbins ; i++) imass_asymmetry_bins[i] = TMath::Power( 10,(logMmin + i*logMbinwidth) );


	/*	
	const int imass_nbins = 12;
	const int ncol_pads   = 2; // = imass_nbins/nrow_pads !!!
	const int nrow_pads   = 6; // = imass_nbins/ncol_pads !!!
	double imass_min   = imass_min_const;
	double imass_max   = imass_max_const;
	Double_t imass_bins[imass_nbins+1] = {72.62*GeV2TeV, 83.37*GeV2TeV, 95.73*GeV2TeV, 109.91*GeV2TeV,
									  126.19*GeV2TeV, 144.89*GeV2TeV, 166.35*GeV2TeV, 191.00*GeV2TeV,
									  219.30*GeV2TeV, 251.79*GeV2TeV, 289.09*GeV2TeV, 331.92*GeV2TeV, 381.09*GeV2TeV};
	doLogx = false;
	//imass_bins = {2.00*GeV, 2.30*GeV, 2.64*GeV, 3.03*GeV, 3.48*GeV, 3.99*GeV, 4.58*GeV, 5.26*GeV, 6.04*GeV, 6.93*GeV, 7.96*GeV, 9.14*GeV, 10.50*GeV, 12.05*GeV, 13.84*GeV, 15.89*GeV, 18.24*GeV, 20.94*GeV, 24.05*GeV, 27.61*GeV, 31.70*GeV, 36.39*GeV, 41.79*GeV, 47.98*GeV, 55.08*GeV, 63.25*GeV, 72.62*GeV, 83.37*GeV, 95.73*GeV, 109.91*GeV, 126.19*GeV, 144.89*GeV, 166.35*GeV, 191.00*GeV, 219.30*GeV, 251.79*GeV, 289.09*GeV, 331.92*GeV, 381.09*GeV, 437.55*GeV, 502.38*GeV, 576.81*GeV, 662.26*GeV, 760.38*GeV, 873.03*GeV, 1002.37*GeV, 1150.88*GeV, 1321.39*GeV, 1517.16*GeV, 1741.93*GeV, 2000.00*GeV};
	*/

	const int imass_nbins = 50;
	double    imass_min   = imass_min_const;
	double    imass_max   = imass_max_const;
	logMmin     = log10(imass_min);
	logMmax     = log10(imass_max);
	Double_t imass_bins[imass_nbins+1];
	logMbinwidth = (Double_t)( (logMmax-logMmin)/(Double_t)imass_nbins );
	imass_bins[0] = imass_min;
	for(Int_t i=1 ; i<=imass_nbins ; i++) imass_bins[i] = TMath::Power( 10,(logMmin + i*logMbinwidth) );
	
	//////////////////////////////////////////////////////////////////////////////////
	// fill the vector with new vector<double> pointers //////////////////////////////
	
	/*
	CNV->Divide(ncol_pads,nrow_pads);
	CNV->SetFillColor(0);
	stringstream strm;
	string str;
	for(int i=0 ; i<imass_nbins ; i++)
	{
		VVCOSTH.push_back(new vector<double>); ////////
		strm.clear();
		str.clear();
		strm << imass_bins[i]+(imass_bins[i+1]-imass_bins[i])/2.;
		strm >> str;
		str = "#hat{m}_{#mu#mu} = " + str + " TeV";
		VHIST_Z0.push_back(new TH1D(("Z^{0} "+str).c_str(), ("Z^{0} "+str).c_str(), 20, -1., +1.));
		VHIST_Z0[i]->SetLineColor(kAzure-5);
		VHIST_DATA.push_back(new TH1D(("Data "+str).c_str(), ("Data "+str).c_str(), 20, -1., +1.));
		VHIST_ZPRIME.push_back(new TH1D(("Z' "+str).c_str(), ("Z' "+str).c_str(), 20, -1., +1.));
		VHIST_ZPRIME[i]->SetLineColor(kRed);
		if(i<nrow_pads) VPAD.push_back( CNV->cd(2*i+1) );
		else            VPAD.push_back( CNV->cd(2*(i-nrow_pads)+2) );
		VPAD[i]->SetFillColor(0);
		VPAD[i]->SetTicky(1);
		VPAD[i]->SetTickx(1);
	}
	*/
	
	stringstream strm;
	string str;
	for(int i=0 ; i<imass_nbins ; i++) VVCOSTH.push_back(new vector<double>);
	//////////////////////////////////////////////////////////////////////////////////
	

	string dir   = "/data/hod/D3PDfin/rel16/";
	string hDir  = "allCuts";
	string hName = "Afb";
	string xTitle = "#hat{m}_{#mu#mu} TeV";
	string yTitle= "A_{FB}";

	string m_dataAnalysisSelector = "digest";	
	//string m_muonSelector = "staco/";
	string m_muonSelector = "muid/";

	double m_miny = -1.;
	double m_maxy = +1.;

	string hNameFixed = hName;

	gStyle->SetOptStat(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	//gStyle->SetFillColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.12);
	gStyle->SetPadBottomMargin(0.16);
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
	gStyle->SetStatBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatX(0);
	gStyle->SetStatY(0);
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0);
	gStyle->SetOptStat(0);
	gStyle->SetStatW(0);
	gStyle->SetStatH(0);
	
	
	TLegend* leg = new TLegend(0.3687291,0.8108808,0.5802676,0.9326425,NULL,"brNDC");
	leg->SetFillColor(kWhite);
	
	TLegend* leg_histos = new TLegend(0.85, 0.15, 0.97, 0.45,NULL,"brNDC");
	leg_histos->SetFillColor(kWhite);
	
	string muonLabel = m_muonSelector.substr(0, m_muonSelector.length()-1);
	string lumilabel = "#intLdt~42 pb^{-1}";
	string label = "#splitline{" + lumilabel + "}{" + muonLabel + " 2010}";
	TPaveText* pvtxt = new TPaveText(0.1655518,0.2072539,0.3102007,0.3316062,"brNDC");
	pvtxt->SetFillColor(kWhite);
	TText* txt  = pvtxt->AddText( label.c_str() );
	//TText* txt  = pvtxt->AddText( lumilabel.c_str() );
	//TText* txt1 = pvtxt->AddText( muonLabel.c_str() );
	
	string cName = "cnv_" + hNameFixed;
	TCanvas* cnv = new TCanvas(cName.c_str(), cName.c_str(), 0,0,1200,800);
	
	TPad *pad_mHat = new TPad("pad_mHat","",0,0,1,1);
	pad_mHat->SetFillColor(kWhite);
	pad_mHat->SetTicky(0);
	pad_mHat->SetLogy();
	if(doLogx) pad_mHat->SetLogx();

	TPad *pad_Afb  = new TPad("pad_Afb", "",0,0,1,1);
	//pad_Afb->SetGridy();
	pad_Afb->SetTicky(0);
	pad_Afb->SetTickx(1);
	pad_Afb->SetFillStyle(0);
	pad_Afb->SetFrameFillStyle(4000); //will be transparent
	pad_Afb->SetFrameFillColor(0);
	if(doLogx) pad_Afb->SetLogx();
	
	string path = "";
	string sProc = "";
	string channel = "";
	//string analysisType = (m_dataAnalysisSelector=="digest") ? "mcDigestControl_" : "mcOfflineControl_";
	string analysisType = "mcLocalControl_";
	
	// Data
	TH1D* hDataM;
	if(doLogM) hDataM  = new TH1D("mHat_data","mHat_data", imass_nbins, imass_bins );
	else       hDataM  = new TH1D("mHat_data","mHat_data", imass_nbins, imass_min, imass_max );
	hDataM->SetTitle("");
	hDataM->SetYTitle("Events");
	hDataM->SetLineColor(kBlack);
	hDataM->SetLineWidth(1);
	TH1D* hData;
	if(doLogM) hData   = new TH1D("Afb_data","Afb_data", imass_asymmetry_nbins, imass_asymmetry_bins );
	else       hData   = new TH1D("Afb_data","Afb_data", imass_asymmetry_nbins, imass_asymmetry_min, imass_asymmetry_max );
	hData->SetMarkerStyle(20);
	hData->SetMarkerColor(kBlack);
	hData->SetMarkerSize(1.2);
	leg->AddEntry( hData, "2010 Data: A_{FB} fit", "lep");
	leg->AddEntry( hData, "2010 Data: Events", "l");
	
	//string sData = (m_dataAnalysisSelector=="digest") ? "digestControl" : "offlineControl";
	string sData = "analysisLocalControl";
	
	//path = dir + "AtoI2_ZprimeGRL/" + m_muonSelector + sData + ".root";
	path = dir + "data_A-I2_eta24/" + m_muonSelector + sData + ".root";
	cout << "path=" << path << endl;
	TFile* fData = new TFile( path.c_str(), "READ" );
	TTree* Afb_data_tree = (TTree*)fData->Get("allCuts/allCuts_tree");
	Afb_data_tree->SetBranchAddress( "Mhat",   &IMASS );
	//Afb_data_tree->SetBranchAddress( refframe.c_str(), &COSTH );
	Afb_data_tree->SetBranchAddress( "charge", &CHARGE );
	Afb_data_tree->SetBranchAddress( "px", &PX );
	Afb_data_tree->SetBranchAddress( "py", &PY );
	Afb_data_tree->SetBranchAddress( "pz", &PZ );
	Afb_data_tree->SetBranchAddress( "E", &E );
	
	// fill the imass histo
	if (Afb_data_tree==0) return;
	for (Long64_t l64t_jentry=0 ; l64t_jentry<Afb_data_tree->GetEntries() ; l64t_jentry++)
	{
		Afb_data_tree->GetEntry(l64t_jentry);
		hDataM->Fill(IMASS);
	}
	
	TAxis *xaxis = hData->GetXaxis();
	
	double a4  = 0.;
	double da4 = 0.;
	double afb  = 0.;
	double dafb = 0.;
	double guess = -0.1;
	ISZMUMU  = false;
	ISZPRIME = false;
	ISDATA   = true;
	
	cout << "### 1 ###" << endl;
	
	fillVec(Afb_data_tree, hData); // the VVCOSTH vectors are full
	
	cout << "### 2 ###" << endl;
	
	for(Int_t b=1 ; b<=hData->GetNbinsX() ; b++)
	{
		// norm to bin width
		//hDataM->SetBinContent(b, hDataM->GetBinContent(b)/hDataM->GetBinWidth(b));
		CURRENTBIN = (unsigned int)(b-1);
		if((int)VVCOSTH[CURRENTBIN]->size()<minEntriesDATA) continue;
		minimize(guess, a4, da4);
		afb = (3./8.)*a4;
		dafb = (3./8.)*da4;
		guess = a4;
		hData->SetBinContent(b,afb);
		hData->SetBinError(b,dafb);
		cout << "mHat["<< b << "/" << hData->GetNbinsX() << "]=" << xaxis->GetBinCenter(b) << ", Afb=" << afb << " +- " << dafb << "\n" << endl;
	}
	
	
	///////////////////////////////////////////////
	// Z' /////////////////////////////////////////
	channel = "0.25 TeV Z' SSM: A_{FB} fit";
	TH1D* hSignal;
	if(doLogM) hSignal = new TH1D("Afb_sig","Afb_sig", imass_asymmetry_nbins, imass_asymmetry_bins );
	else       hSignal = new TH1D("Afb_sig","Afb_sig", imass_asymmetry_nbins, imass_asymmetry_min, imass_asymmetry_max );
	hSignal->SetLineColor(kRed);
	hSignal->SetFillColor(kRed);
	hSignal->SetLineWidth(1);
	hSignal->SetMarkerSize(0);
	hSignal->SetMarkerColor(0);
	hSignal->SetTitle("");
	hSignal->SetXTitle( xTitle.c_str() );
	hSignal->SetYTitle( yTitle.c_str() );
	//leg->AddEntry( hSignal, channel.c_str(), "f");
	
	sProc = "Zprime_mumu_SSM250";
	path = dir + "Zprime_mumu/" + m_muonSelector + analysisType + sProc + ".root";
	cout << "path=" << path << endl;
	TFile* fZprime = new TFile( path.c_str(), "READ" );
	TTree* Afb_sig_tree = (TTree*)fZprime->Get("allCuts/allCuts_tree");
	Afb_sig_tree->SetBranchAddress( "Mhat",       &IMASS );
	//Afb_sig_tree->SetBranchAddress( refframe.c_str(), &COSTH );
	Afb_sig_tree->SetBranchAddress( "charge", &CHARGE );
	Afb_sig_tree->SetBranchAddress( "px", &PX );
	Afb_sig_tree->SetBranchAddress( "py", &PY );
	Afb_sig_tree->SetBranchAddress( "pz", &PZ );
	Afb_sig_tree->SetBranchAddress( "E", &E );
	
	a4   = 0.;
	da4  = 0.;
	afb  = 0.;
	dafb = 0.;
	guess = -0.1;
	ISZMUMU  = false;
	ISZPRIME = true;
	ISDATA   = false;
	fillVec(Afb_sig_tree, hSignal); // the VVCOSTH vectors are full
	for(Int_t b=1 ; b<=hSignal->GetNbinsX() ; b++)
	{
		CURRENTBIN = (unsigned int)(b-1);
		if((int)VVCOSTH[CURRENTBIN]->size()<minEntriesDATA) continue;
		minimize(guess, a4, da4);
		afb = (3./8.)*a4;
		dafb = (3./8.)*da4;
		guess = a4;
		hSignal->SetBinContent(b,afb);
		hSignal->SetBinError(b,dafb);
		cout << "mHat["<< b << "/" << hSignal->GetNbinsX() << "]=" << xaxis->GetBinCenter(b) << ", Afb=" << afb << " +- " << dafb << "\n" << endl;
	}
	
	TH1D* hSignalTmp = (TH1D*)hSignal->Clone("");
	hSignalTmp->Reset();
	for(Int_t b=1 ; b<hSignalTmp->GetNbinsX() ; b++) hSignalTmp->SetBinContent(b,hSignal->GetBinContent(b));
	hSignalTmp->SetLineColor(kRed+8);
	
	/////////////////////////////////////////////////////////////////////
	// Backgrounds //////////////////////////////////////////////////////
	channel = "MC(Z#rightarrow#mu#mu): A_{FB} fit";
	TH1D* hBGsum;
	if(doLogM) hBGsum = new TH1D("Afb_sumBG","Afb_sumBG", imass_asymmetry_nbins, imass_asymmetry_bins );
	else       hBGsum = new TH1D("Afb_sumBG","Afb_sumBG", imass_asymmetry_nbins, imass_asymmetry_min, imass_asymmetry_max );
	hBGsum->SetLineColor(kAzure-5);
	hBGsum->SetFillColor(kAzure-5);
	hBGsum->SetFillStyle(3003);
	hBGsum->SetLineWidth(1);
	hBGsum->SetMarkerSize(0);
	hBGsum->SetMarkerColor(0);
	hBGsum->SetTitle("");
	hBGsum->SetXTitle( xTitle.c_str() );
	hBGsum->SetYTitle( yTitle.c_str() );
	leg->AddEntry( hBGsum, channel.c_str(), "f");
	
	/*
	sProc = "Zmumu";
	path = dir + "Zmumu/" + m_muonSelector + analysisType + sProc + ".root";
	cout << "path=" << path << endl;
	TFile* fZmumu = new TFile( path.c_str(), "READ" );
	TTree* Afb_sumBG_tree = (TTree*)fZmumu->Get("allCuts/allCuts_tree");
	Afb_sumBG_tree->SetBranchAddress( "Mhat",       &IMASS );
	Afb_sumBG_tree->SetBranchAddress( refframe.c_str(), &COSTH );
	*/
	
	TList* Afb_sumBG_list = new TList();
	
	sProc = "Zmumu";
	path = dir + "Zmumu/" + m_muonSelector + analysisType + sProc + ".root";
	cout << "path=" << path << endl;
	TFile* fZmumu = new TFile( path.c_str(), "READ" );
	Afb_sumBG_list->Add( (TTree*)fZmumu->Get("allCuts/allCuts_tree") );
	
	sProc = "DYmumu_250M400";
	path = dir + "DYmumu/" + m_muonSelector + analysisType + sProc + ".root";
	cout << "path=" << path << endl;
	TFile* fDYmumu_250M400 = new TFile( path.c_str(), "READ" );
	Afb_sumBG_list->Add( (TTree*)fDYmumu_250M400->Get("allCuts/allCuts_tree") );

	sProc = "DYmumu_400M600";
	path = dir + "DYmumu/" + m_muonSelector + analysisType + sProc + ".root";
	cout << "path=" << path << endl;
	TFile* fDYmumu_400M600 = new TFile( path.c_str(), "READ" );
	Afb_sumBG_list->Add( (TTree*)fDYmumu_400M600->Get("allCuts/allCuts_tree") );
	
	TFile* mergedFile = new TFile("Afb_sumBG_merged.root", "RECREATE");
	cout << "Merging trees...patience..." << endl;
	TTree::MergeTrees(Afb_sumBG_list);
	mergedFile->Write();
	mergedFile->Close();
	
	TFile* fBG = new TFile("Afb_sumBG_merged.root", "READ");
	TTree* Afb_sumBG_tree = (TTree*)fBG->Get("allCuts_tree");
	Afb_sumBG_tree->SetBranchAddress( "Mhat",       &IMASS );
	//Afb_sumBG_tree->SetBranchAddress( refframe.c_str(), &COSTH );
	Afb_sumBG_tree->SetBranchAddress( "charge", &CHARGE );
	Afb_sumBG_tree->SetBranchAddress( "px", &PX );
	Afb_sumBG_tree->SetBranchAddress( "py", &PY );
	Afb_sumBG_tree->SetBranchAddress( "pz", &PZ );
	Afb_sumBG_tree->SetBranchAddress( "E", &E );
	
	
	a4   = 0.;
	da4  = 0.;
	afb  = 0.;
	dafb = 0.;
	guess = -0.1;
	ISZMUMU  = true;
	ISZPRIME = false;
	ISDATA   = false;
	fillVec(Afb_sumBG_tree, hBGsum); // the VVCOSTH vectors are full
	for(Int_t b=1 ; b<=hBGsum->GetNbinsX() ; b++)
	{
		CURRENTBIN = (unsigned int)(b-1);
		if((int)VVCOSTH[CURRENTBIN]->size()<minEntriesDATA) continue;
		minimize(guess, a4, da4);
		afb = (3./8.)*a4;
		dafb = (3./8.)*da4;
		guess = a4;
		hBGsum->SetBinContent(b,afb);
		hBGsum->SetBinError(b,dafb);
		cout << "mHat["<< b << "/" << hBGsum->GetNbinsX() << "]=" << xaxis->GetBinCenter(b) << ", Afb=" << afb << " +- " << dafb << "\n" << endl;
	}
	TH1D* hBGsumTmp = (TH1D*)hBGsum->Clone("");
	hBGsumTmp->Reset();
	for(Int_t b=1 ; b<hBGsumTmp->GetNbinsX() ; b++) hBGsumTmp->SetBinContent(b,hBGsum->GetBinContent(b));
	hBGsumTmp->SetLineColor(kAzure+8);
	
	
	/*	
	TH1D* hDataMerr = (TH1D*)hDataM->Clone("");
	hDataMerr->SetMarkerStyle(1);
	hDataMerr->SetMarkerSize(1);
	for(Int_t i=0 ; i<=hDataMerr->GetNbinsX()+1 ; i++) hDataMerr->SetBinError(i,sqrt(hDataMerr->GetBinContent(i)));
	*/
	bool isXerr = false;
	TGraphAsymmErrors* gDataMpoissonErr = GetPoissonizedGraph(hDataM,isXerr);
	gDataMpoissonErr->SetMarkerStyle(1);
	gDataMpoissonErr->SetMarkerSize(1);
	gDataMpoissonErr->SetLineColor(kBlack);
	gDataMpoissonErr->SetLineStyle(1);

	
	pad_mHat->Draw();
	pad_mHat->cd();
	if(doLogx) hDataM->GetYaxis()->SetRangeUser(0.1,1.5*hDataM->GetMaximum());
	else       hDataM->GetYaxis()->SetRangeUser(1,1.5*hDataM->GetMaximum());
	hDataM->GetXaxis()->SetMoreLogLabels(); 
	hDataM->GetXaxis()->SetNoExponent(); 
	hDataM->Draw();
	gDataMpoissonErr->GetXaxis()->SetMoreLogLabels(); 
	gDataMpoissonErr->GetXaxis()->SetNoExponent(); 
	gDataMpoissonErr->Draw("SAMES");

	cnv->cd();

	pad_Afb->Draw();
	pad_Afb->cd();
	//hSignal->GetYaxis()->SetRangeUser(m_miny,m_maxy);
	//hSignal->Draw("E5 Y+");
	//hSignalTmp->Draw("CSAMES");
	hBGsum->GetYaxis()->SetRangeUser(m_miny,m_maxy);
	hBGsum->GetXaxis()->SetMoreLogLabels(); 
	hBGsum->GetXaxis()->SetNoExponent(); 
	hBGsum->Draw("E5 Y+");
	//hBGsumTmp->Draw("CSAMES");
	hData->Draw("e1x0SAMES");
	pvtxt->Draw("SAMES");
	leg->Draw("SAMES");
	TLine* lUnit = new TLine(imass_min,0.,imass_max,0.);
	lUnit->SetLineColor(kBlack);
	lUnit->SetLineStyle(2);
	//lUnit->Draw("SAMES");
	//pad_Afb->RedrawAxis();
	
	cnv->cd();
	
	pad_mHat->cd();
	pad_mHat->RedrawAxis();
	pad_mHat->Update();

	cnv->Update();
	
	TString fName = "figures/" + (TString)hNameFixed + "_" + (TString)muonLabel;
	cnv->SaveAs(fName+".eps");
	cnv->SaveAs(fName+".C");
	cnv->SaveAs(fName+".root");
	cnv->SaveAs(fName+".png");
	
	/*
	CNV->Draw();
	for(int i=0 ; i<imass_nbins ; i++)
	{
		if(i==0)
		{
			leg_histos->AddEntry(VHIST_Z0[i], "Z#rightarrow#mu#mu");
			leg_histos->AddEntry(VHIST_Z0[i], "0.25 TeV Z' SSM#rightarrow#mu#mu");
			leg_histos->AddEntry(VHIST_Z0[i], "Data", "lep");
		}
	
		VPAD[i]->cd();
		
		Double_t NZ0     = VHIST_Z0[i]->GetEntries();
		Double_t NZprime = VHIST_ZPRIME[i]->GetEntries();
		Double_t NDATA   = VHIST_DATA[i]->GetEntries();
		
		Scale(VHIST_Z0[i], 1./NZ0);
		Scale(VHIST_ZPRIME[i], 1./NZprime);
		ScaleWerrors(VHIST_DATA[i], 1./NDATA);
		
		Double_t max = 0.;
		Double_t tmp = 0.;
		tmp = VHIST_Z0[i]->GetMaximum();
		max = (tmp>max) ? tmp : max;
		if(NZprime>=minEntriesDATA)
		{
			tmp = VHIST_ZPRIME[i]->GetMaximum();
			max = (tmp>max) ? tmp : max;
		}
		if(NDATA>=minEntriesDATA)
		{
			tmp = VHIST_DATA[i]->GetMaximum();
			max = (tmp>max) ? tmp : max;
		}
		
		
		VHIST_Z0[i]->SetMaximum(1.2*max);
		VHIST_Z0[i]->SetMinimum(0.);
		VHIST_Z0[i]->SetTitle("");
		VHIST_Z0[i]->SetXTitle("#cos#theta^{*}");
		VHIST_Z0[i]->SetYTitle("Events (normalized)");
		VHIST_Z0[i]->Draw();

		VHIST_ZPRIME[i]->SetTitle("");
		if(NZprime>=minEntriesDATA) VHIST_ZPRIME[i]->Draw("SAMES");
		
		VHIST_DATA[i]->SetTitle("");
		if(NDATA>=minEntriesDATA) VHIST_DATA[i]->Draw("e1x0SAMES");
		
		leg_histos->Draw("SAMES");
	}
	fName += ".costh";
	CNV->SaveAs(fName+".eps");
	CNV->SaveAs(fName+".C");
	CNV->SaveAs(fName+".root");
	CNV->SaveAs(fName+".png");
	*/
}
