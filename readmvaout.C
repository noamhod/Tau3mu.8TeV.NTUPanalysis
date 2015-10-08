/////////////////////////////////////////////////////////////////////////////////////////
//// root -b -l -q readmvaout.C++\(1450.,1690.,1870.,2110.,30.,0.8000,0.9990,0.001\) ////
/////////////////////////////////////////////////////////////////////////////////////////
#include "allFitSystAuto.h"
#include "postBDTcuts.h"
#include <iterator>
#include <cstdlib>

Bool_t reject;
Double_t xbmin = 1700.;
Double_t xbmax = 1860.;
Int_t NminFlat = 10000;//7;
Int_t npar0 = 6;
Int_t npar1 = 3;
Int_t npar2 = 1;
//// http://root.cern.ch/root/html534/tutorials/fit/fitExclude.C.html 
Double_t fBg0(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*exp(x*par[2]) + par[3]*exp(-pow((x-par[4])/par[5],2));
}
Double_t fBg1(Double_t* X, Double_t* par)
{
        Double_t x = X[0];
        if(reject && x>xbmin && x<xbmax)
        {
                TF1::RejectPoint();
                return 0;
        }
        return par[0] + par[1]*exp(x*par[2]);
}
Double_t fBg2(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	if(reject && x>xbmin && x<xbmax)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0];
}


// TString tstr(float x, int prcn=-1)
// {
// 	stringstream strm;
// 	string str;
// 	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
// 	else         strm << fixed << x;
// 	strm >> str;
// 	return (TString)str;
// }

TGraphAsymmErrors* poisson(TH1* h)
{
	double value = 0;
	double error_poisson_up   = 0;
	double error_poisson_down = 0;
	double y1 = 0;
	double y2 = 0;
	double d1 = 0;
	double d2 = 0;
	TGraphAsymmErrors *ga = new TGraphAsymmErrors();
	for(int i=1; i<=h->GetNbinsX(); i++)
	{
		value = h->GetBinContent(i);
		// if(value!=0)
		if(h->GetXaxis()->GetBinUpEdge(i)<xbmin || h->GetXaxis()->GetBinUpEdge(i)>xbmax)
		{
			y1 = value + 1.0;
			d1 = 1.0 - 1.0/(9.0*y1) + 1.0/(3*TMath::Sqrt(y1));
			error_poisson_up = y1*d1*d1*d1 - value;
			y2 = value;
			d2 = 1.0 - 1.0/(9.0*y2) - 1.0/(3.0*TMath::Sqrt(y2));
			error_poisson_down = value - y2*d2*d2*d2;
			ga->SetPoint(i-1, h->GetBinCenter(i), value);
			ga->SetPointError(i-1, 0, 0, error_poisson_down, error_poisson_up);
		}
	}
	ga->SetMarkerColor(kBlack);
	ga->SetMarkerStyle(20);
	ga->SetMarkerSize(0.6);
	// ga->SetLineWidth(2);
	ga->SetLineWidth(1);
	ga->SetLineColor(kBlack);
	ga->GetXaxis()->SetLimits(h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
	return ga;
}

vector<string> split_string_by_spaces(string line, bool print=true)
{
	if(print) cout << "Line to split: " << line << endl;
	istringstream buf(line);
	istream_iterator<std::string> beg(buf), end;
	vector<string> tokens(beg,end); // done!
	// for(unsigned int i=0 ; i<tokens.size() ; ++i) cout << "[" << i << "]" << '"' << tokens[i] << '"' << endl;
	return tokens;
}                              
float get_word_to_float(string word)
{
	float number = ::strtof(word.c_str(), 0);
	return number;
}

TPaveText* ptxt;
void makeAtlasLabel()
{
	// ptxt = new TPaveText(0.15,0.70,0.40,0.87,"NDC");
	ptxt = new TPaveText(0.55,0.40,0.80,0.57,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}
TLegend* legN;
void makeLegend()
{
	legN = new TLegend(0.4,0.65,0.88,0.88,NULL,"brNDC");
	legN->SetFillStyle(4000); //will be transparent
	legN->SetFillColor(0);
	legN->SetTextFont(42);
	legN->SetBorderSize(0);
}

string makeHistFitterPy(vector<string>& lines, string score, float nbkgSR, float nbkgSRerr, float systUp, float systDwn, float ndata)
{
	string fname = "python/tau3mu/tau3mu."+score+".py";
	cout << "Creating file: " << fname << endl;
	ofstream hfpyx(fname.c_str(),std::ofstream::out);
	for(unsigned int i=0 ; i<lines.size() ; ++i)
	{
		if(lines[i].find("nbkg = ")!=string::npos)
		{
			vector<string> words = split_string_by_spaces(lines[i],false);
			unsigned int nwords = words.size();
			for(unsigned int j=0 ; j<nwords ; ++j)
			{
				if     (j==2)      hfpyx << nbkgSR << " ";
				else if(j==nwords) hfpyx << endl;
				else               hfpyx << words[j] << " ";
			}
			hfpyx << endl;
		}
		else if(lines[i].find("nbkgErr = ")!=string::npos)
		{
			vector<string> words = split_string_by_spaces(lines[i],false);
			unsigned int nwords = words.size();
			for(unsigned int j=0 ; j<nwords ; ++j)
			{	
				// if  (j==2)      { float nbkgErr = (nbkgSB>0) ? (nbkgSR/nbkgSB)*sqrt(nbkgSB) : 1; hfpyx << nbkgErr << " "; }
				if     (j==2)      { float nbkgErr = nbkgSRerr; hfpyx << nbkgErr << " "; }
				else if(j==nwords) hfpyx << endl;
				else               hfpyx << words[j] << " ";
			}
			hfpyx << endl;
		}
		else if(lines[i].find("up  = ")!=string::npos)
		{
			vector<string> words = split_string_by_spaces(lines[i],false);
			unsigned int nwords = words.size();
			for(unsigned int j=0 ; j<nwords ; ++j)
			{	
				if     (j==2)      { hfpyx << systUp << " "; }
				else if(j==nwords) hfpyx << endl;
				else               hfpyx << words[j] << " ";
			}
			hfpyx << endl;
		}
		else if(lines[i].find("dwn = ")!=string::npos)
		{
			vector<string> words = split_string_by_spaces(lines[i],false);
			unsigned int nwords = words.size();
			for(unsigned int j=0 ; j<nwords ; ++j)
			{	
				if     (j==2)      { hfpyx << systDwn << " "; }
				else if(j==nwords) hfpyx << endl;
				else               hfpyx << words[j] << " ";
			}
			hfpyx << endl;
		}
		else if(lines[i].find("ndata = ")!=string::npos)
		{
			vector<string> words = split_string_by_spaces(lines[i],false);
			unsigned int nwords = words.size();
			for(unsigned int j=0 ; j<nwords ; ++j)
			{
				if     (j==2)      hfpyx << (int)ndata << " ";
				else if(j==nwords) hfpyx << endl;
				else               hfpyx << words[j] << " ";
			}
			hfpyx << endl;
		}
		else hfpyx << lines[i] << endl;
	}
	hfpyx.close();
	return fname;
}

TCanvas* makeSummaryPlot(TH1F* hEffS, TH1F* hEffSrel, TH1F* hUL90med, TH1F* hUL90min1, TH1F* hUL90pls1, TH1F* hUL90min2, TH1F* hUL90pls2, TString outfilesname)
{
	TGraphAsymmErrors* g1sigma = new TGraphAsymmErrors();
	for(Int_t i=1 ; i<=hUL90med->GetNbinsX() ; i++) 
	{
		float x = hUL90med->GetBinCenter(i);
		// float xErrLeft  = hUL90med->GetBinLowEdge(i);
		// float xErrRight = hUL90med->GetBinLowEdge(i)+hUL90med->GetBinWidth(i);
   
		float y       = hUL90med->GetBinContent(i);
		float yErrDwn = fabs(y-hUL90min1->GetBinContent(i));
		float yErrUp  = fabs(y-hUL90pls1->GetBinContent(i));

		g1sigma->SetPoint(i-1, x, y);
		g1sigma->SetPointError(i-1, 0,0, yErrDwn,yErrUp);
	}
	TGraphAsymmErrors* g2sigma = new TGraphAsymmErrors();
	for(Int_t i=1 ; i<=hUL90med->GetNbinsX() ; i++) 
	{
		float x = hUL90med->GetBinCenter(i);
		// float xErrLeft  = hUL90med->GetBinLowEdge(i);
		// float xErrRight = hUL90med->GetBinLowEdge(i)+hUL90med->GetBinWidth(i);
    
		float y       = hUL90med->GetBinContent(i);
		float yErrDwn = fabs(y-hUL90min2->GetBinContent(i));
		float yErrUp  = fabs(y-hUL90pls2->GetBinContent(i));
    
		g2sigma->SetPoint(i-1, x, y);
		g2sigma->SetPointError(i-1, 0,0, yErrDwn,yErrUp);
	}
	g1sigma->SetFillColor(kGreen);
	g1sigma->SetLineColor(kGreen);
	g1sigma->SetTitle(hUL90med->GetTitle());
	g1sigma->GetXaxis()->SetTitle(hUL90med->GetXaxis()->GetTitle());
	g1sigma->GetYaxis()->SetTitle(hUL90med->GetYaxis()->GetTitle());
	g2sigma->SetFillColor(kYellow);
	g2sigma->SetLineColor(kYellow);
	g2sigma->SetTitle(hUL90med->GetTitle());
	g2sigma->GetXaxis()->SetTitle(hUL90med->GetXaxis()->GetTitle());
	g2sigma->GetYaxis()->SetTitle(hUL90med->GetYaxis()->GetTitle());
	TMultiGraph *mgBands = new TMultiGraph();
	mgBands->Add(g2sigma);
	mgBands->Add(g1sigma);
	// mgBands->SetMaximum(1.2*hUL90med->GetMaximum());
	mgBands->SetMaximum(2.e-6);
	mgBands->SetMinimum(0.);
	
	hEffS->SetLineColor(kGray+1);
	hEffS->SetLineStyle(3);
	hEffS->SetLineWidth(2);
	hEffS->SetTitle("");
	
	hEffSrel->SetLineColor(kBlue+1);
	hEffSrel->SetLineStyle(3);
	hEffSrel->SetLineWidth(2);
	hEffSrel->SetTitle("");
	
	// fout->cd();
	// hEffS->Write();
	// // hEffB->Write();
	// hUL90med->Write();
	// hUL90med->Write();
	// g1sigma->Write();
	// g2sigma->Write();
	// fout->Write();
	// // fout->Close();
	
	float bestmedian = hUL90med->GetMinimum();
	float bestscore = -999.;
	Int_t bestbin = -1;
	for(Int_t j=1 ; j<=hUL90med->GetNbinsX() ; ++j) { if(hUL90med->GetBinContent(j)==bestmedian) { bestscore = hUL90med->GetBinCenter(j); bestbin = j; break; } }
	cout << "bestbin=" << bestbin << endl;
	TString sbestscore = tstr(bestscore,4);
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/m3body.fit."+sbestscore+".eps /tmp/hod/HistFitter/eps/"+outfilesname+"m3body.fit.optcut.eps");
	cout << "----------------------" << endl;	
	cout <<   "best median   : " << bestmedian << endl;
	cout <<   "optimal cut(f): " << bestscore  << endl;
	cout <<   "optimal cut(s): " << sbestscore << endl;
	cout << "----------------------" << endl;

	TLine* ATLASexp = new TLine(hUL90med->GetXaxis()->GetXmin(),bestmedian,bestscore,bestmedian);
	ATLASexp->SetLineColor(kRed);
	ATLASexp->SetLineWidth(2);
	
	float LHCb = 4.6e-8;
	TLine* lLHCb = new TLine(hUL90med->GetXaxis()->GetXmin(),LHCb,hUL90med->GetXaxis()->GetXmax(),LHCb);
	lLHCb->SetLineColor(kAzure+10);
	lLHCb->SetLineWidth(2);
	
	TLegend* leg = new TLegend(0.15,0.60,0.42,0.89);
	leg->SetFillStyle(4000); //will be transparent
	// leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	// leg->AddEntry(hEffSrel,"#scale[0.9]{#it{#epsilon}_{S}#times#it{A}_{S}} (rel)","l");
	leg->AddEntry(hEffS,"#scale[0.9]{#it{#epsilon}_{S}#times#it{A}_{S}}","l");
	leg->AddEntry(hUL90med,"#scale[0.9]{Median}","l");
	leg->AddEntry(g1sigma,"#scale[0.9]{Median#pm1#sigma}","f");
	leg->AddEntry(g2sigma,"#scale[0.9]{Median#pm2#sigma}","f");
	leg->AddEntry(lLHCb,"#scale[0.9]{LHCb limit}","l");
	// leg->AddEntry(lbest,"#scale[0.9]{Optimal point}","l");
	leg->AddEntry(ATLASexp,"#scale[0.9]{Optimal point}","l");
	
	TCanvas* c1 = new TCanvas("cnv","",1200,800);
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
	
	hUL90med->SetLineStyle(2);
	hUL90med->SetLineWidth(2);
	hUL90med->SetTitle("Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs");
	hUL90med->GetXaxis()->SetTitle("BDT(G) score cut");
	hUL90med->GetYaxis()->SetTitle("Br(#it{#tau#rightarrow3#mu})");
	hUL90med->Draw("L");
	mgBands->Draw("a3");
	hUL90med->Draw("L same");
	lLHCb->Draw("L same");
	ATLASexp->Draw("L same");
	
	p1->Update();
	TLine* lbest = new TLine(bestscore,p1->GetUymax(),bestscore,0.);
	lbest->SetLineColor(kRed);
	lbest->SetLineWidth(2);
	lbest->SetLineStyle(1);
	lbest->Draw("L same");
	
	mgBands->SetTitle(hUL90med->GetTitle());
	mgBands->GetXaxis()->SetTitle(hUL90med->GetXaxis()->GetTitle());
	mgBands->GetYaxis()->SetTitle(hUL90med->GetYaxis()->GetTitle());
	mgBands->GetXaxis()->SetLimits(hUL90med->GetXaxis()->GetXmin(),hUL90med->GetXaxis()->GetXmax());
	p1->Modified();
	p1->RedrawAxis();
	c1->cd();
	
	//compute the pad range with suitable margins
	Double_t ymin = 0;
	Double_t ymax = 0.04;
	Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
	Double_t xmin = hUL90med->GetXaxis()->GetXmin();
	Double_t xmax = hUL90med->GetXaxis()->GetXmax();
	Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
	p2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
	p2->Draw();
	p2->cd();
	leg->Draw("same");
	ptxt->Draw("same");
	hEffS->Draw("][L sames");
	// hEffSrel->Draw("][L sames");
	// hEffB->Draw("][hist sames");
	// p2->RedrawAxis();

	// draw axis on the right side of the pad
	TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
	axis->SetTitle("Efficiency");
	axis->SetLabelColor(kBlack);
	axis->SetLabelFont(42);
	axis->SetLabelSize(0.04);
	axis->SetTextFont(42);
	axis->SetTextSize(0.04);
	axis->Draw();
	
	return c1;
}



Int_t nBDTpars = 6;
Double_t fBDTlocal(Double_t* X, Double_t* par)
{
	Double_t x = X[0];
	return par[0]+par[1]*exp(x*par[2])+par[3]*pow(fabs(x),x*par[4])*exp(x*par[5]);
}


void readmvaout(float mMinSBleft, float mMaxSBleft, float mMinSBright, float mMaxSBright, float m3bodyBinSize, float MinBDTcut, float MaxBDTcut, float dX)
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
	
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_200k_3mu");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");
	
	TString sidebands    = tstr(mMinSBleft,0)+"-"+tstr(mMaxSBleft,0)+"-"+tstr(mMinSBright,0)+"-"+tstr(mMaxSBright,0);
	TString bdtmin       = tstr(MinBDTcut,2); bdtmin.ReplaceAll("-","neg"); if(MinBDTcut>0) bdtmin="pos"+bdtmin; if(MinBDTcut==0) bdtmin="zro"+bdtmin; bdtmin.ReplaceAll(".",""); 
	TString bdtmax       = tstr(MaxBDTcut,3); bdtmax="pos"+bdtmax; bdtmax.ReplaceAll(".","");
	TString outfilesname = "Optimization."+sidebands+"."+bdtmin+"."+bdtmax;

	makeAtlasLabel();
	makeLegend();
	
	setFiles();
	TFile* fBDTpdf = getNominal();
	TVectorF* AbsBDTmin = (TVectorF*)fBDTpdf->Get("absBDTmin");
	TVectorF* AbsBDTmax = (TVectorF*)fBDTpdf->Get("absBDTmax");
	TVectorF* NBDTbins  = (TVectorF*)fBDTpdf->Get("nBDTbins");
	float minBDTabs = ((*AbsBDTmin))[0];
	float maxBDTabs = ((*AbsBDTmax))[0];
	int   BDTbins  = (int)((*NBDTbins))[0];


	TH1F* BDTb = new TH1F("BDTb",";BDT score;Events",BDTbins,minBDTabs,maxBDTabs);
	BDTb->SetLineColor(kBlack);
	BDTb->SetLineWidth(1);
	BDTb->SetLineStyle(1);
	BDTb->SetMarkerStyle(20);
	BDTb->SetMarkerSize(0.6);
	BDTb->SetMarkerColor(kBlack);
	BDTb->Sumw2();
	TH1F* BDTs = new TH1F("BDTs",";BDT score;Events",BDTbins,minBDTabs,maxBDTabs);
	BDTs->SetLineColor(kBlue);
	BDTs->SetLineWidth(1);
	BDTs->SetLineStyle(1);
	BDTs->SetMarkerStyle(20);
	BDTs->SetMarkerSize(0.6);
	BDTs->SetMarkerColor(kBlue);
	BDTs->Sumw2();


	TCanvas* cnv;
	
	//// set the blinded area
	xbmin = mMaxSBleft;
	xbmax = mMinSBright;

	TString smMinSBleft  = tstr(mMinSBleft,0);
	TString smMaxSBleft  = tstr(mMaxSBleft,0);
	TString smMinSBright = tstr(mMinSBright,0);
	TString smMaxSBright = tstr(mMaxSBright,0);

	/*
	TString basecuts     = "m3body>"+tstr(mMinSBleft,0)+" && m3body<"+tstr(mMaxSBright,0);  // "m3body>"+smMinSBleft+".";
	TString inblinded    = "(m3body>"+tstr(mMaxSBleft,0)+" && m3body<"+tstr(mMinSBright,0)+")";
	TString basecuts_bkg = basecuts+" && !"+inblinded;  // basecuts+" && (m3body<"+tstr(xbmin,0)+". || m3body>"+tstr(xbmax,0)+".)";
	TString basecuts_sig = basecuts;
	*/

	TString basecuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","bkg");
	TString basecuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-1","-1","sig");

	cout << "basecuts_bkg=" << basecuts_bkg << endl;
	cout << "basecuts_sig=" << basecuts_sig << endl;


	Float_t ninitS = 200000;
	Float_t nallS  = tS->GetEntries(basecuts_sig);
	//Float_t nallD  = tD->GetEntries(basecuts_bkg);

	

	
	tD->SetEventList(0);
	tD->Draw(">>elistDall",basecuts_bkg);
	TEventList* elistDall = (TEventList*)gDirectory->Get("elistDall");
	tD->SetEventList(elistDall);
	tD->Draw("score>>BDTb");
	BDTb->SetMinimum(0.);
	tD->SetEventList(0);

	tS->SetEventList(0);
	tS->Draw(">>elistSall",basecuts_sig);
	TEventList* elistSall = (TEventList*)gDirectory->Get("elistSall");
	tS->SetEventList(elistSall);
	tS->Draw("score>>BDTs");
	BDTs->SetMinimum(0.);
	tS->SetEventList(0);



	TH1F* mS = new TH1F("mS", ";m_{3body} [MeV];Events", (Int_t)((mMaxSBright-mMinSBleft)/m3bodyBinSize),mMinSBleft,mMaxSBright);
	mS->Sumw2();
	tS->Draw("m3body>>mS",basecuts_sig);
	
	mS->Fit("gaus","0");
	TF1* fGauss = mS->GetFunction("gaus");
	Double_t chi2 = fGauss->GetChisquare();
	Int_t    ndf = fGauss->GetNDF();
	Double_t pvalue = TMath::Prob(chi2,ndf);
	// Double_t N  = fGauss->GetParameter(0);
	// Double_t dN = fGauss->GetParError(0);
	Double_t avg  = fGauss->GetParameter(1);
	Double_t davg = fGauss->GetParError(1);
	Double_t sigma  = fGauss->GetParameter(2);
	Double_t dsigma = fGauss->GetParError(2);
	cout << "Gauss chi2/ndf: " << chi2/ndf  << endl;
	cout << "Gauss pvalue:   " << pvalue  << endl;
	cout << "Gauss Average:  " << avg   << " +- " << davg << endl;
	cout << "Gauss Sigma:    " << sigma << " +- " << dsigma << endl;
	Double_t ximin = 1715.;//avg-2.*sigma;
	Double_t ximax = 1842.;//avg+2.*sigma;
	TLine* lmin = new TLine(ximin,mS->GetMaximum()/2.,ximin,0.);
	TLine* lmax = new TLine(ximax,mS->GetMaximum()/2.,ximax,0.);
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->SetTicks(1,1);
	mS->Draw("p");
	fGauss->Draw("same");
	lmin->Draw("same");
	lmax->Draw("same");
	cnv->Update();
	cnv->SaveAs(outfilesname+".pdf(");
	
	
	TH1F* mD = new TH1F("mD", ";m_{3body} [MeV];Events", (Int_t)((mMaxSBright-mMinSBleft)/m3bodyBinSize),mMinSBleft,mMaxSBright);
	mD->SetLineColor(kBlack);
	mD->SetLineWidth(1);
	mD->SetLineStyle(1);
	mD->SetMarkerStyle(20);
	mD->SetMarkerSize(0.6);
	mD->SetMarkerColor(kBlack);
	mD->Sumw2();



	bool hFit = true; // true for histogram fit, false for graph fit

	// remember the python script
	vector<string> hfpylines;
	string hfpyline;
	fstream hfpy;
	hfpy.open("python/MyUserAnalysis.py");
	while(getline(hfpy,hfpyline)) hfpylines.push_back(hfpyline);
	hfpy.close();
	

	gSystem->Exec("rm -f python/tau3mu/*");
	gSystem->Exec("rm -f txt/*");
	
	gSystem->Exec("rm -rf /tmp/hod/HistFitter");
	gSystem->Exec("mkdir /tmp/hod/HistFitter");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/eps");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/txt");
	
	float acceffrel_prev=-999, acceff_prev=-999., prev_upper=-999, prev_median=-999, prev_m1sigma=-999, prev_p1sigma=-999, prev_m2sigma=-999, prev_p2sigma=-999;
	
	float dx = dX;
	float xmin = MinBDTcut;
	float xmax = MaxBDTcut;
	Float_t hXmin = xmin-dx/2.;
	Float_t hXmax = xmax+dx/2.;
	Int_t hNbins  = (hXmax-hXmin)/dx;
	// if((int)((xmax-xmin)/dx)%2==0) hNbins++;
	TH1F* hEffSrel  = new TH1F("SignalEfficiencyRel", ";BDT score cut;Signal #it{A}#times#it{#epsilon} (relative)",  hNbins,hXmin,hXmax);
	TH1F* hEffS     = new TH1F("SignalEfficiency",    ";BDT score cut;Signal #it{A}#times#it{#epsilon} (absolute)",  hNbins,hXmin,hXmax);
	TH1F* hUL90min1 = new TH1F("UpperLimitMin1Sig",   "Expected upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90pls1 = new TH1F("UpperLimitPls1Sig",   "Expected upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90min2 = new TH1F("UpperLimitMin2Sig",   "Expected upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90pls2 = new TH1F("UpperLimitPls2Sig",   "Expected upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90med  = new TH1F("UpperLimitMedian",    "Expected upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);

	TH1F* hnSB1Count = new TH1F("SB1 Count",";BDT score cut;Events passing BDT cut in SB1", hNbins,hXmin,hXmax);
	TH1F* hnSB1Shape = new TH1F("SB1 Shape",";BDT score cut;Events passing BDT cut in SB1", hNbins,hXmin,hXmax);
	TH1F* hnSR1Count = new TH1F("SR1 Count",";BDT score cut;Events passing BDT cut in SR1", hNbins,hXmin,hXmax);
	TH1F* hnSR1Shape = new TH1F("SR1 Shape",";BDT score cut;Events passing BDT cut in SR1", hNbins,hXmin,hXmax);

	
	TF1* fbg0 = NULL;
	TF1* fbg1 = NULL;
	TF1* fbg2 = NULL;


	Int_t hbin = 1;
	for(float x=xmin ; x<xmax ; x+=dx, ++hbin)
	{
		cout << "\n" << endl;
	
		stringstream strm;
		string       str;
		strm << setprecision(4) << fixed << x;
		strm >> str;
		TString score = str;
		string  sscore = str;


		// //float iBkgSB1 = tf1BDT->Integral(x,maxBDTabs)/BDTbinWidth;
		// float iBkgSB1 = tf1BDTub->Integral(x,maxBDTabs)/BDTbinWidth;
		// float f01 = iBkgSB1/iBkgSB0;
		// float diBkgSB1stat = diBkgSB0stat*f01;
		// float iBkgSR1 = iBkgSR0*f01;
		// float diBkgSR1stat = diBkgSR0stat*f01;
		
		// Result* res = getExtrapolationResult(x,false);
		// float iBkgSB1      = res->nSB1;
		// float diBkgSB1stat = res->dnSB1stat;
		// float f01          = res->f01;
		// float iBkgSR1      = res->nSR1;
		// float diBkgSR1stat = res->dnSR1stat;
		// float diBkgSR1syst = res->dnSR1quad;
		// delete res;
		
		Result* res = getExtrapolationResult(x,63/*=iSR*/,false);
		// float srwidth      = res->xSRmax-res->xSRmin; 
		float iBkgSB1      = res->nSB1;
		float diBkgSB1stat = res->dnSB1stat;
		float f01          = res->f01;
		float iBkgSR1      = res->VARnSR1;
		float diBkgSR1stat = res->VARdnSR1stat;
		float diBkgSR1syst = res->VARdnSR1quad;
		TString sxSRmin    = tstr(res->xSRmin,0);
		TString sxSRmax    = tstr(res->xSRmax,0);
		ximin              = res->xSRmin;
		ximax              = res->xSRmax;
		delete res;


		TString cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9",score,"bkg");
		TString cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9",score,"sig",sxSRmin,sxSRmax);
		cout << "cuts_bkg=" << cuts_bkg << endl;
		cout << "cuts_sig=" << cuts_sig << endl;

		tS->SetEventList(0);
		tS->Draw(">>elistS",cuts_sig);
		TEventList* elistS = (TEventList*)gDirectory->Get("elistS");
		Float_t npassedS = elistS->GetN();  // number of events to pass cuts
	
		tD->SetEventList(0);
		tD->Draw(">>elistD",cuts_bkg);
		TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
		Float_t npassedD = elistD->GetN();  // number of events to pass cuts
		cout << cuts_bkg << ":  npassedS=" << npassedS << ", npassedD=" << npassedD << ", nfitD=" << iBkgSB1 << "+-" << diBkgSB1stat << "(for f01="<<f01<<")" << ", acc*eff=" << npassedS/ninitS*100 << "\%" << endl;
		
		
		//////////////////////////////////////////////
		if(npassedS/ninitS*100<0.005 && x>0.90)
		{
			cout << "####################################" << endl;
			cout << "### Skipping this point: score=" << x << ", AxE=" << npassedS/ninitS*100 << " ###" << endl;
			cout << "####################################" << endl;
			continue;
		}
		//////////////////////////////////////////////
		
		

		////////////////
		if(fbg0) delete fbg0;
		fbg0 = new TF1("Background0", fBg0, mMinSBleft,mMaxSBright, npar0);
		fbg0->SetParName(0,"Const"); fbg0->SetParLimits(0,0,100);
		fbg0->SetParName(1,"Nexp");  fbg0->SetParLimits(1,0,100);
		fbg0->SetParName(2,"Gexp");  fbg0->SetParLimits(2,-10,0);
		fbg0->SetParName(3,"Ngaus"); fbg0->SetParLimits(3,0,50);
		fbg0->SetParName(4,"Mean");  fbg0->SetParLimits(4,1550,1650);
		fbg0->SetParName(5,"Sigma"); fbg0->SetParLimits(5,1,100);
	
		if(fbg1) delete fbg1;
		fbg1 = new TF1("Background1", fBg1, mMinSBleft,mMaxSBright, npar1);
		fbg1->SetParName(0,"Const");  fbg1->SetParLimits(0,0,1000);
		fbg1->SetParName(1,"Nexp");   fbg1->SetParLimits(1,0,1000);
		fbg1->SetParName(2,"Gexp");   fbg1->SetParLimits(2,-10,0);

		if(fbg2) delete fbg2;
		fbg2 = new TF1("Background2", fBg2, mMinSBleft,mMaxSBright, npar2);
		fbg2->SetParName(0,"c0");  fbg2->SetParLimits(0,0,1000);
		////////////////



		mD->Reset();
		tD->SetEventList(elistD);
		tD->Draw("m3body>>mD", basecuts_bkg);
		for(int b=1 ; b<=mD->GetNbinsX() ; ++b) { if(mD->GetBinContent(b)<1) mD->SetBinError(b,1); }
		for(int b=1 ; b<=mD->GetNbinsX() ; ++b) { if(mD->GetBinContent(b)<1 && mD->GetXaxis()->GetBinLowEdge(b)>=xbmin && mD->GetXaxis()->GetBinUpEdge(b)<=xbmax) mD->SetBinContent(b,-10); }
		mD->SetMinimum(0.);
		if(mD->GetMaximum()<3.) mD->SetMaximum(3.);
		TGraphAsymmErrors* gD = (TGraphAsymmErrors*)poisson(mD);
		reject = kTRUE;
		if(hFit)
		{
			//if(npassedD<9) mD->Fit("Background2","L0 EMR");
			//else
			//{
				if(npassedD>NminFlat) mD->Fit("Background0","L0 EMR");
				else                  mD->Fit("Background1","L0 EMR");
			//}
		}
		else
		{
			//if(npassedD<2) gD->Fit("Background2","EX0");
			//else
			//{
				if(npassedD>NminFlat) gD->Fit("Background0","EX0");
				else                  gD->Fit("Background1","EX0");
			//}
		}
		
		/*
		*option: The second parameter is the fitting option. Here is the list of fitting options: 
		o "W" Set all weights to 1 for non empty bins; ignore error bars 
		o "WW" Set all weights to 1 including empty bins; ignore error bars 
		o "I" Use integral of function in bin instead of value at bin center 
		o "L" Use log likelihood method (default is chi-square method) 
		o "U" Use a user specified fitting algorithm 
		o "Q" Quiet mode (minimum printing) 
		o "V" Verbose mode (default is between Q and V) 
		o "E" Perform better errors estimation using the Minos technique 
		o "M" Improve fit results 
		o "R" Use the range specified in the function range 
		o "N" Do not store the graphics function, do not draw 
		o "0" Do not plot the result of the fit. By default the fitted function is drawn unless 
		 the option "N" above is specified. 
		o "+" Add this new fitted function to the list of fitted functions (by default, the 
		 previous function is deleted and only the last one is kept) 
		o "B" Use this option when you want to fix one or more parameters and the fitting 
		 function is like polN, expo, landau, gaus. 
		o LL An improved Log Likelihood fit in case of very low statistics and when bin
		 contents are not integers. Do not use this option if bin contents are large 
		 (greater than 100). 
		o C In case of linear fitting, don't calculate the chisquare (saves time). 
		o F If fitting a polN, switch to Minuit fitter (by default, polN functions are 
		 fitted by the linear fitter).
		
		Likelihood Fits
		     When using option "L" a likelihood fit is used instead of the default chi2 square fit.
		     The likelihood is built assuming a Poisson probability density function for each bin.
		     The negative log-likelihood to be minimized is
		      NLL = Sum{ log Poisson( y(i) |{ f(x(i) | p ) ) }
		     The exact likelihood used is the Poisson likelihood described in this paper:
		     S. Baker and R. D. Cousins, â€œClarification of the use of chi-square and likelihood functions in fits to histograms,â€
		     Nucl. Instrum. Meth. 221 (1984) 437.
		     This method can then be used only when the bin content represents counts (i.e. errors are sqrt(N) ).
		     The likelihood method has the advantage of treating correctly bins with low statistics. In case of high
		     statistics/bin the distribution of the bin content becomes a normal distribution and the likelihood and chi2 fit
		     give the same result.
		     The likelihood method, although a bit slower, it is therefore the recommended method in case of low
		     bin statistics, where the chi2 method may give incorrect results, in particular when there are
		     several empty bins (see also below).
		     In case of a weighted histogram, it is possible to perform a likelihood fit by using the
		     option "WL". Note a weighted histogram is an histogram which has been filled with weights and it
		     contains the sum of the weight square ( TH1::Sumw2() has been called). The bin error for a weighted
		     histogram is the square root of the sum of the weight square.

		Treatment of Empty Bins
		     Empty bins, which have the content equal to zero AND error equal to zero,
		     are excluded by default from the chisquare fit, but they are considered in the likelihood fit.
		     since they affect the likelihood if the function value in these bins is not negligible.
		     When using option "WW" these bins will be considered in the chi2 fit with an error of 1.
		     Note that if the histogram is having bins with zero content and non zero-errors they are considered as
		     any other bins in the fit. Instead bins with zero error and non-zero content are excluded in the chi2 fit.
		     A likelihood fit should also not be peformed on such an histogram, since we are assuming a wrong pdf for each bin.
		     In general, one should not fit an histogram with non-empty bins and zero errors, apart if all the bins have zero errors.
		     In this case one could use the option "w", which gives a weight=1 for each bin (unweighted least-square fit).
		*/
		
		
		reject = kFALSE;
		chi2   = (npassedD>NminFlat) ? fbg0->GetChisquare() : fbg1->GetChisquare();
		ndf    = (npassedD>NminFlat) ? fbg0->GetNDF()       : fbg1->GetNDF();
		//if(npassedD<9) chi2 = fbg2->GetChisquare();
		//if(npassedD<9) ndf = fbg2->GetNDF();
		pvalue = TMath::Prob(chi2,ndf);
		float iBkgSR = (npassedD>NminFlat) ? fbg0->Integral(ximin,ximax)/mD->GetBinWidth(1) : fbg1->Integral(ximin,ximax)/mD->GetBinWidth(1);
		//if(npassedD<9) iBkgSR = fbg2->Integral(ximin,ximax)/mD->GetBinWidth(1);
		float iBkgSB = (npassedD>NminFlat) ? (fbg0->Integral(mMinSBleft,xbmin)+fbg0->Integral(xbmax,mMaxSBright))/mD->GetBinWidth(1) : (fbg1->Integral(mMinSBleft,xbmin)+fbg1->Integral(xbmax,mMaxSBright))/mD->GetBinWidth(1);
		//if(npassedD<9) iBkgSB = (fbg2->Integral(mMinSBleft,xbmin)+fbg2->Integral(xbmax,mMaxSBright))/mD->GetBinWidth(1);
		//if(npassedD<9) cout << "Using: " << fbg2->GetName() << " as fit function" << endl;
		//else
		//{
			if(npassedD>NminFlat) cout << "Using: " << fbg0->GetName() << " as fit function" << endl;
			else                  cout << "Using: " << fbg1->GetName() << " as fit function" << endl;
		//}
		cout << "Sidebands   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
		//if(npassedD<9) cout << "Integral(" << ximin << "," << ximax << "): " << fbg2->Integral(ximin,ximax) << endl;
		//else
		//{
			if(npassedD>NminFlat) cout << "Integral(" << ximin << "," << ximax << "): " << fbg0->Integral(ximin,ximax) << endl;
			else                  cout << "Integral(" << ximin << "," << ximax << "): " << fbg1->Integral(ximin,ximax) << endl;
		//}
		cout << "Bin width: " << mD->GetBinWidth(1) << endl;
		cout << "SB expectation (" << npassedD << "): " << iBkgSB << "(cut) and " << iBkgSB1 << "(shape)" << endl;
		cout << "SR expectation [" << ximin << "," << ximax << "]: " << iBkgSR << "(cut) and " << iBkgSR1 << "(shape)" << endl;
		TF1 *fleft = (npassedD>NminFlat) ? new TF1("fleft",fBg0,mMinSBleft,xbmin,npar0) : new TF1("fleft",fBg1,mMinSBleft,xbmin,npar1);
		//if(npassedD<9)
		//{
		//	fleft = new TF1("fleft",fBg2,mMinSBleft,xbmin,npar2);
		//	fleft->SetParameters(fbg2->GetParameters());
		//}
		//else
		//{
			if(npassedD>NminFlat) fleft->SetParameters(fbg0->GetParameters());
			else                  fleft->SetParameters(fbg1->GetParameters());
		//}
		if(hFit) mD->GetListOfFunctions()->Add(fleft);
		else     gD->GetListOfFunctions()->Add(fleft);
		gROOT->GetListOfFunctions()->Remove(fleft);
		TF1 *fright = (npassedD>NminFlat) ? new TF1("fright",fBg0,xbmax,mMaxSBright,npar0) : new TF1("fright",fBg1,xbmax,mMaxSBright,npar1);
		//if(npassedD<9)
		//{
		//	fright = new TF1("fright",fBg2,mMinSBleft,xbmin,npar2);
		//	fright->SetParameters(fbg2->GetParameters());
		//}
		//else
		//{
			if(npassedD>NminFlat) fright->SetParameters(fbg0->GetParameters());
			else                  fright->SetParameters(fbg1->GetParameters());
		//}
		if(hFit) mD->GetListOfFunctions()->Add(fright);
		else     gD->GetListOfFunctions()->Add(fright);
		gROOT->GetListOfFunctions()->Remove(fright);
		
		delete cnv;
		cnv = new TCanvas("","",600,400);
		cnv->Draw();
		cnv->cd();
		cnv->SetTicks(1,1);
		if(hFit) mD->Draw("p0 e0");
		else     gD->Draw("AP");

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(iBkgSR<0. || iBkgSB<0.) { iBkgSR=0.; iBkgSB=0.; cout << "iBkgSR<0 or iBkgSB<0 ==> setting to 0" << endl; } ///////
		if(iBkgSR1<0. || iBkgSB1<0.) { iBkgSR1=0.; iBkgSB1=0.; diBkgSB1stat=0; diBkgSR1stat=0; cout << "iBkgSR1<0 or iBkgSB1<0 ==> setting to 0" << endl; } ///////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		hnSB1Shape->SetBinContent(hbin,iBkgSB1);
		hnSB1Shape->SetBinError(hbin,diBkgSB1stat);
		hnSR1Shape->SetBinContent(hbin,iBkgSR1);
		hnSR1Shape->SetBinError(hbin,diBkgSR1stat);


		//// remake the python steering script for this evaluation
		float dnSR1systUp = 1+diBkgSR1syst/iBkgSR1;
		float dnSR1systDn = 1-diBkgSR1syst/iBkgSR1;
		cout << "SR uncertainty [" << ximin << "," << ximax << "]: " << "stat=" << tstr(diBkgSR1stat/iBkgSR1*100,1) << "\%, syst=" << tstr(diBkgSR1syst/iBkgSR1*100,1) << "\%" << endl;
		
		string pyfile = makeHistFitterPy(hfpylines,sscore,iBkgSR1,diBkgSR1stat,dnSR1systUp,dnSR1systDn,floor(iBkgSR1+0.5));
		// string pyfile = makeHistFitterPy(hfpylines,sscore,iBkgSR1,diBkgSR1stat,floor(iBkgSR1+0.5));
		//string pyfile = makeHistFitterPy(hfpylines,sscore,iBkgSR,iBkgSB,floor(iBkgSR+0.5));
		gSystem->Exec("HistFitter.py -w -f -l "+pyfile+" > /tmp/hod/HistFitter/txt/histfitter"+score+".txt 2>&1");
		gSystem->Exec("/bin/cp -f results/upperlimit_cls_poi_Sig_Asym_CLs_grid_ts3.root.eps /tmp/hod/HistFitter/eps/upperlimit_cls_poi_Sig_Asym_CLs_grid_ts3."+score+".root.eps");
		gSystem->Exec("/bin/cp -f results/MyUserAnalysis_Output_upperlimit.root /tmp/hod/HistFitter/eps/MyUserAnalysis_Output_upperlimit."+score+".root");
		ifstream fileInput;
		string line;
		// open file to search
		string fname = "/tmp/hod/HistFitter/txt/histfitter"+(string)score+".txt";
		fileInput.open(fname.c_str());
		vector<string> words;
		float upper=-999, median=-999, m1sigma=-999, p1sigma=-999, m2sigma=-999, p2sigma=-999;
		bool ok = false;
		cout << "Reading file " << fname << endl;
		while(getline(fileInput,line))
		{
			/*
			<INFO> HypoTestTool: The computed upper limit is: 1.37945 +/- 0
			<INFO> HypoTestTool:  expected limit (median) 1.37945
			<INFO> HypoTestTool:  expected limit (-1 sig) 0.664826
			<INFO> HypoTestTool:  expected limit (+1 sig) 2.91938
			<INFO> HypoTestTool:  expected limit (-2 sig) 0.366556
			<INFO> HypoTestTool:  expected limit (+2 sig) 5.75149
			*/
			if(line.find("The computed upper limit")!=string::npos) { words = split_string_by_spaces(line); upper   = get_word_to_float(words[7]); }
			if(line.find("expected limit (median)")!=string::npos)  { words = split_string_by_spaces(line); median  = get_word_to_float(words[5]); }
			if(line.find("expected limit (-1 sig)")!=string::npos)  { words = split_string_by_spaces(line); m1sigma = get_word_to_float(words[6]); }
			if(line.find("expected limit (+1 sig)")!=string::npos)  { words = split_string_by_spaces(line); p1sigma = get_word_to_float(words[6]); }
			if(line.find("expected limit (-2 sig)")!=string::npos)  { words = split_string_by_spaces(line); m2sigma = get_word_to_float(words[6]); }
			if(line.find("expected limit (+2 sig)")!=string::npos)  { words = split_string_by_spaces(line); p2sigma = get_word_to_float(words[6]); }
			if(upper>=0 && median>=0 && m1sigma>=0 && p1sigma>=0 && m2sigma>=0 && p2sigma>=0) { ok = true; break; }
		}
				
		if(!ok) { _ERR(1,"Could not find one of the values: upper="<<upper<<", median="<<median<<", m1sigma="<<m1sigma<<", p1sigma="<<p1sigma<<", m2sigma="<<m2sigma<<", p2sigma"<<p2sigma); }
		else    { cout << "HistFitter summary: upper="<<upper<<", median="<<median<<", m1sigma="<<m1sigma<<", p1sigma="<<p1sigma<<", m2sigma="<<m2sigma<<", p2sigma="<<p2sigma << endl; }

		float AtimesErel = npassedS/nallS;
		float AtimesE = npassedS/ninitS;
		float Ntau = 2.41e8*(luminosity/20.28);
	
		bool isBadHistFitter = (median<1.e-10 || upper<1.e-10 || m1sigma<1.e-10 || p1sigma<1.e-10 || m2sigma<1.e-10 || p2sigma<1.e-10);
		if(isBadHistFitter) // something is wrong --> use the result from the previous iteration
		{
			cout << "#########################################" << endl;
			cout << "######### bad HistFitter result #########" << endl;
			cout << "#########################################" << endl;
			AtimesErel = acceffrel_prev;
			AtimesE = acceff_prev;
			upper   = prev_upper;
			median  = prev_median;
			m1sigma = prev_m1sigma;
			p1sigma = prev_p1sigma;
			m2sigma = prev_m2sigma;
			p2sigma = prev_p2sigma;
		}
		else // everything is OK --> update the previous result to be the current one
		{
			acceffrel_prev  = AtimesErel;
			acceff_prev  = AtimesE;
			prev_upper   = upper;
			prev_median  = median;
			prev_m1sigma = m1sigma;
			prev_p1sigma = p1sigma;
			prev_m2sigma = m2sigma;
			prev_p2sigma = p2sigma;
		}

		hEffSrel->SetBinContent(hbin,AtimesErel);
		hEffS->SetBinContent(hbin,AtimesE);
		hUL90med->SetBinContent(hbin,median/AtimesE/Ntau);
		hUL90min1->SetBinContent(hbin,m1sigma/AtimesE/Ntau);
		hUL90pls1->SetBinContent(hbin,p1sigma/AtimesE/Ntau);
		hUL90min2->SetBinContent(hbin,m2sigma/AtimesE/Ntau);
		hUL90pls2->SetBinContent(hbin,p2sigma/AtimesE/Ntau);
		hnSB1Count->SetBinContent(hbin,iBkgSB);
		hnSB1Count->SetBinError(hbin,sqrt(iBkgSB));
		hnSR1Count->SetBinContent(hbin,iBkgSR);
		hnSR1Count->SetBinError(hbin, (iBkgSB>0) ? (iBkgSR/iBkgSB)*sqrt(iBkgSB) : 1);


		cout << "Upper limit on Br(tau->3mu)=" << median/AtimesE/Ntau << " for x=" << x << endl;
		
		float iBkgSBerr = (iBkgSB>0) ? (iBkgSR/iBkgSB)*sqrt(iBkgSB) : 1;
		float effS      = npassedS/ninitS*100;
		TString nSB = "N_{SB}^{obs} = "+tstr(npassedD,0)+" (#scale[0.6]{#int}f_{SB}^{L}+#scale[0.6]{#int}f_{SB}^{R}="+tstr(iBkgSB,2)+")";
		TString nSR = "N_{SR}^{exp} = "+tstr(iBkgSR,2)+"#pm"+tstr(iBkgSBerr,2)+" (stat)";
		TString acceff = "Signal #it{A}#times#epsilon = "+tstr(effS,2)+"\%";
		
		TPaveText *pt = new TPaveText(0.1,0.6,0.9,0.9,"NDC");
		pt->SetFillStyle(4000); //will be transparent
		pt->SetFillColor(0);
		pt->SetTextFont(42);
		pt->SetBorderSize(0);
		pt->AddText(cuts_bkg);
		pt->AddText(nSB);
		pt->AddText(nSR);
		pt->AddText(acceff);
		pt->Draw("same");
		
		cnv->Update();
		cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".m3body.fit."+score+".eps");
		cnv->SaveAs(outfilesname+".pdf");
		
		tD->SetEventList(0);
	}
	
	cnv = makeSummaryPlot(hEffS,hEffSrel,hUL90med,hUL90min1,hUL90pls1,hUL90min2,hUL90pls2,outfilesname);
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization.eps");
	cnv->SaveAs(outfilesname+".pdf");


	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	BDTb->Draw("p0 e0");
	//tf1BDT->Draw("same");
	// tf1BDTub->Draw("same");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSB0.eps");
	cnv->SaveAs(outfilesname+".pdf");

	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	BDTs->Draw("p0 e0");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSignal.eps");
	cnv->SaveAs(outfilesname+".pdf");



	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	hnSB1Count->SetMinimum(0.); hnSB1Count->SetFillColor(kGray+2);  hnSB1Count->SetLineColor(kBlack); hnSB1Count->SetMarkerColor(kBlack); hnSB1Count->SetMarkerStyle(20); hnSB1Count->SetMarkerSize(0.7);  hnSB1Count->Draw("e2");
	hnSB1Shape->SetMinimum(0.); hnSB1Shape->SetFillColor(kRed);     hnSB1Shape->SetLineColor(kRed);   hnSB1Shape->SetMarkerColor(kRed); hnSB1Shape->SetMarkerStyle(20); hnSB1Shape->SetMarkerSize(0.7); hnSB1Shape->Draw("e1x0 same");
	legN->Clear();
	legN->AddEntry(hnSB1Count,"N_{SB1} Cut&count (stat. error only)", "plef");
	legN->AddEntry(hnSB1Shape,"N_{SB1} BDT&SB0 fits (stat. error only)", "ple");
	legN->Draw("smae");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSB.eps");
	cnv->SaveAs(outfilesname+".pdf");

	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	hnSR1Count->SetMinimum(0.); hnSR1Count->SetFillColor(kGray+2);  hnSR1Count->SetLineColor(kBlack); hnSR1Count->SetMarkerColor(kBlack); hnSR1Count->SetMarkerStyle(20); hnSR1Count->SetMarkerSize(0.7); hnSR1Count->Draw("e2");
	hnSR1Shape->SetMinimum(0.); hnSR1Shape->SetFillColor(kRed); hnSR1Shape->SetLineColor(kRed);  hnSR1Shape->SetMarkerColor(kRed); hnSR1Shape->SetMarkerStyle(20); hnSR1Shape->SetMarkerSize(0.7); hnSR1Shape->Draw("e1x0 same");
	legN->Clear();
	legN->AddEntry(hnSR1Count,"N_{SR1} Cut&SB1 fit (stat. error only)", "plef");
	legN->AddEntry(hnSR1Shape,"N_{SR1} BDT&SB0 fits (stat. error only)", "ple");
	legN->Draw("smae");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSR.eps");
	cnv->SaveAs(outfilesname+".pdf");

	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->Update();
	cnv->SaveAs(outfilesname+".pdf)");

	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".m3body.fit.optcut.eps /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSB0.eps       /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSignal.eps    /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSB.eps               /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSR.eps               /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f "+outfilesname+".pdf                                           /afs/cern.ch/user/h/hod/data/HistFitter/");
}
