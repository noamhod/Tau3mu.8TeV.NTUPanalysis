//////////////////////////////////////////////////////////
//// root -l -b -q fitParErr.C++\(1450,2110,-0.9,+1\) ////
//////////////////////////////////////////////////////////

#include "std.h"
#include "const.h"
#include "postBDTcuts.h"
#include "roofit.h"

bool blinded = false;

Double_t xbdtmin  = minBDTcut;
Double_t xbdtmax  = +1.0;
Double_t xbdtopt  = optBDTcut;
Int_t    nbdtbins = nBDTbins;

Double_t m3bodyMin = mSideBandLeftLowerMeVGlob;
Double_t m3bodyMax = mSideBandRightUpperMeVGlob;
Double_t xbmin  = mBlindMinGlob;
Double_t xbmax  = mBlindMaxGlob;
Double_t m3bodyBinSize = mBinSize;

TRandom* randGen;

string str(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return str;
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

TGraphErrors* inflate(TH1* h)
{
	TGraphErrors* g = new TGraphErrors(h);
	Int_t nbins       = h->GetNbinsX();
	Double_t xmin     = h->GetXaxis()->GetXmin();
	Double_t xmax     = h->GetXaxis()->GetXmax();
	Double_t binxerr  = h->GetBinWidth(1)/2;	
	g->SetPoint(0,xmin,h->GetBinContent(1));
	g->SetPointError(0,binxerr,h->GetBinError(1));
	g->SetPoint(nbins,xmax,h->GetBinContent(nbins));
	g->SetPointError(nbins,binxerr,h->GetBinError(nbins));
	return g;
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

float maxdistance(vector<float>& v, float reference, bool print=false)
{
	float maxd = -1.e20;
	for(unsigned int i=0 ; i<v.size() ; ++i)
	{
		float d = fabs(v[i]-reference);
		if(print) cout << "["<<i<<"] " << d << endl;
		maxd = (d>maxd) ? d : maxd;
	}
	if(print) cout << "max distance=" << maxd <<endl;
	return maxd;
}
float rmsd(vector<float>& v, float reference=0, bool doreference=false, bool print=false)
{
	float x = 0;
	float av = reference;
	float n = 0;
	for(unsigned int i=0 ; i<v.size() ; ++i)
	{
		if(!isnan(v[i])) n++;
	}
	if(!doreference)
	{
		av = 0;
		for(unsigned int i=0 ; i<v.size() ; ++i)
		{
			if(!isnan(v[i])) av += v[i];
		}
		av = av/n;
	}
	for(unsigned int i=0 ; i<v.size() ; ++i)
	{
		if(!isnan(v[i])) x += (av-v[i])*(av-v[i]);
	}
	x = sqrt(x/n);
	if(print) cout << "RMSD=" << x <<endl;
	return x;
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

void plotAtlasLabel()
{
	Double_t x = 0.64;
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

TLegend* leg;
TLegend* legR;
void makeLegend()
{
	leg = new TLegend(0.25,0.67,0.5,0.87,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	
	legR = new TLegend(0.65,0.60,0.90,0.67,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
}


TString expandNscaleFormula(TString formula, TString varname, RooRealVar par[], Int_t npar, float scale=1, int variation=0)
{
	TString sscale = (scale!=1) ? "*"+tstr(scale) : "";
	formula.ReplaceAll(varname,"x");
	cout << "formula=" << formula << endl;
	for(Int_t n=0 ; n<npar ; ++n)
	{
		double val0  = par[n].getVal();              TString sval0  = (val0>1.e-4)  ? tstr(val0)  : tstr(val0,12);
		double valHi = val0+par[n].getAsymErrorHi(); TString svalHi = (valHi>1.e-4) ? tstr(valHi) : tstr(valHi,12);
		double valLo = val0+par[n].getAsymErrorLo(); TString svalLo = (valLo>1.e-4) ? tstr(valLo) : tstr(valLo,12);
		
		cout << "par["<< n << "] (" << par[n].GetName() << ") -> " << par[n].getVal() << " +" << par[n].getAsymErrorHi() << " -" << par[n].getAsymErrorLo() << endl;
		if     (variation==0)  formula.ReplaceAll(par[n].GetName(),sval0);
		else if(variation==+1) formula.ReplaceAll(par[n].GetName(),svalHi);
		else if(variation==-1) formula.ReplaceAll(par[n].GetName(),svalLo);
		else _FAT("wrong variation value");
	}
	TString expanded = "("+formula+")"+sscale;
	cout << "expanded=" << expanded << endl;
	return expanded;
}


bool minuitStatus(TMinuit* m) 
{
	if(!m) return false;
	TString stat = gMinuit->fCstatu;
	_INF(1,"Minuit: "<<stat<<". ");
	if(stat.Contains("SUCCESSFUL") || stat.Contains("CONVERGED") || stat.Contains("OK") ) return true;
	return false;
}

double randomizeItialGuess(double min, double max)
{
	return min + (max-min)*randGen->Uniform(); // Uniform(x1=0) returns a uniform deviate on the interval [0,x1].
}

int readData(TTree* t, RooRealVar* score, RooRealVar* m3body, RooAbsData* data_score, RooAbsData* data_m3body,
			float xFullMin,float xBlindMin,float xBlindMax,float xFullMax, float xBDTMin,float xBDTOpt, TString type, TString CRi, bool isBlind=true, TH1* hScore=NULL, TH1* hMass=NULL)
{
	_INF(1,"========== Reading candidates =========");
	
	data_score->Clear();
	data_m3body->Clear();
		
	vector<float>* bdtscore = 0;
	vector<float>* mass     = 0;
	vector<float>* mOS1     = 0;
	vector<float>* mOS2     = 0;
	vector<float>* mSS      = 0;
	vector<float>* pt       = 0;
	vector<float>* pval     = 0;
	vector<float>* pval1    = 0;
	vector<float>* pval2    = 0;
	vector<float>* pval3    = 0;
	vector<int>*   code     = 0;
	
	vector<float>* mettrk    = 0;
	vector<float>* metcal    = 0;
	vector<float>* trkspval  = 0;
	vector<float>* maxpbalsig= 0;
	vector<float>* iso003    = 0;
	vector<float>* iso014    = 0;
	vector<float>* iso020    = 0;
	vector<float>* iso030    = 0;
	vector<float>* metsdphi  = 0;
	vector<float>* metsdptrel= 0;
	vector<float>* mtcal     = 0;
	vector<float>* dphical   = 0;
	vector<float>* mttrk     = 0;
	vector<float>* dphitrk   = 0;
	vector<float>* lxy       = 0;
	vector<float>* dlxy      = 0;
	vector<float>* lxysig    = 0;
	vector<int>*   pvntrk    = 0;
	vector<float>* a0xy      = 0;
	vector<float>* da0xy     = 0;
	vector<float>* a0xysig   = 0;
	vector<int>*   njets     = 0;
	vector<int>*   bjets     = 0;
	vector<float>* ht        = 0;
	vector<float>* mhtcal    = 0;
	vector<float>* dphihtcal = 0;
	vector<float>* mhttrk    = 0;
	vector<float>* dphihttrk = 0;
	vector<float>* drht      = 0;
	vector<float>* dptcalrel = 0;
	vector<float>* dpttrkrel = 0;
	vector<float>* dhtcalrel = 0;
	vector<float>* dhttrkrel = 0;
	
	t->SetBranchAddress("score",&bdtscore);
	t->SetBranchAddress("m3body",&mass);
	t->SetBranchAddress("mOS1",&mOS1);
	t->SetBranchAddress("mOS2",&mOS2);
	t->SetBranchAddress("mSS",&mSS);
	t->SetBranchAddress("pt3body",&pt);
	t->SetBranchAddress("pval",&pval);
	t->SetBranchAddress("pval1",&pval1);
	t->SetBranchAddress("pval2",&pval2);
	t->SetBranchAddress("pval3",&pval3);
	t->SetBranchAddress("code",&code);
	
	t->SetBranchAddress("mettrk",     &mettrk    );
	t->SetBranchAddress("metcal",     &metcal    );
	t->SetBranchAddress("trkspval",   &trkspval  );
	t->SetBranchAddress("maxpbalsig", &maxpbalsig);
	t->SetBranchAddress("iso003",     &iso003    );
	t->SetBranchAddress("iso020",     &iso020    );
	t->SetBranchAddress("iso014",     &iso014    );
	t->SetBranchAddress("iso030",     &iso030    );
	t->SetBranchAddress("metsdphi",   &metsdphi  );
	t->SetBranchAddress("mtcal",      &mtcal     );
	t->SetBranchAddress("dphical",    &dphical   );
	t->SetBranchAddress("mttrk",      &mttrk     );
	t->SetBranchAddress("dphitrk",    &dphitrk   );
	t->SetBranchAddress("lxy",        &lxy       );
	t->SetBranchAddress("dlxy",       &dlxy      );
	t->SetBranchAddress("lxysig",     &lxysig    );
	t->SetBranchAddress("pvntrk",     &pvntrk    );
	t->SetBranchAddress("a0xy",       &a0xy      );
	t->SetBranchAddress("da0xy",      &da0xy     );
	t->SetBranchAddress("a0xysig",    &a0xysig   );
	t->SetBranchAddress("njets",      &njets     );
	t->SetBranchAddress("bjets",      &bjets     );
	t->SetBranchAddress("ht",         &ht        );
	t->SetBranchAddress("mhtcal",     &mhtcal    );
	t->SetBranchAddress("dphihtcal",  &dphihtcal );
	t->SetBranchAddress("mhttrk",     &mhttrk    );
	t->SetBranchAddress("dphihttrk",  &dphihttrk );
	t->SetBranchAddress("drht",       &drht      );
	// t->SetBranchAddress("mets_dptrelcal", &dptcalrel );
	// t->SetBranchAddress("mets_dptreltrk", &dpttrkrel );
	// t->SetBranchAddress("mets_dhtrelcal", &dptcalrel );
	// t->SetBranchAddress("mets_dhtreltrk", &dpttrkrel );
	
	TString channel = t->GetName();

	int ncandidates = 0;
	for(Int_t entry=1 ; entry<=t->GetEntries() ; entry++)
	{
		t->GetEntry(entry);
		
		if(bdtscore->size()>1) continue;
		
		TMapTSP2vi vi;
		TMapTSP2vf vf;
		vector<int>* loose = new vector<int>; loose->push_back(1); // dummy - hast to be deleted later !
		vi.insert(make_pair("njets"     , njets));
		vi.insert(make_pair("pass_loose", loose));
		vf.insert(make_pair("pt3body"   , pt));
		vf.insert(make_pair("m3body"    , mass));
		vf.insert(make_pair("mSS"       , mSS));
		vf.insert(make_pair("mOS1"      , mOS1));
		vf.insert(make_pair("mOS2"      , mOS2));
		vf.insert(make_pair("score"     , bdtscore));
		vf.insert(make_pair("pval"      , pval));
		vf.insert(make_pair("iso003"    , iso003));
		vf.insert(make_pair("iso020"    , iso020));
		vf.insert(make_pair("iso030"    , iso030));
		vf.insert(make_pair("trkspval"  , trkspval));
		vf.insert(make_pair("lxysig"    , lxysig));
		vf.insert(make_pair("a0xysig"   , a0xysig));
		vf.insert(make_pair("mettrk"    , mettrk));
		vf.insert(make_pair("metcal"    , metcal));
		vf.insert(make_pair("mttrk"     , mttrk));
		vf.insert(make_pair("mtcal"     , mtcal));
		vf.insert(make_pair("dphihttrk" , dphihttrk));
		vf.insert(make_pair("dphihtcal" , dphihtcal));
		
		float xBDT1 = (CRi=="CR1") ? xBDTOpt : xBDTMin;
		bool pass = passPostBDTcut(0,vf,vi,xFullMin,xBlindMin,xBlindMax,xFullMax, xBDTMin,xBDT1, type, mSRminMeVGlob,mSRmaxMeVGlob,isBlind);
		delete loose;
		
	
		// if(pass && bdtscore->at(0)>0.85 && type=="bkg")
		// {
		// 	cout << "------- entry=" << entry << " ------" << endl;
		// 	cout << "code         " << code->at(0) << endl;
		// 	cout << "Ncandidates  " << bdtscore->size() << endl;
		// 	cout << "xBDT         " << bdtscore->at(0) << endl;
		// 	cout << "mass         " << mass->at(0) << endl;
		// 	cout << "mOS1         " << mOS1->at(0) << endl;
		// 	cout << "mOS2         " << mOS2->at(0) << endl;
		// 	cout << "mSS          " << mSS->at(0) << endl;
		// 	cout << "pval         " << pval->at(0) << endl;
		// 	cout << "trkspval     " << trkspval->at(0) << endl;
		// 	cout << "njets        " << njets->at(0) << endl;
		// 	cout << "pt           " << pt->at(0) << endl;
		// 	cout << "ht           " << ht->at(0) << endl;
		// 	cout << "dpthtrel     " << pt->at(0)/ht->at(0)-1 << endl;
		// 	cout << "mettrk       " << mettrk    ->at(0) << endl;
		// 	cout << "metcal       " << metcal    ->at(0) << endl;
		// 	cout << "dmetsrel     " << (mettrk->at(0)-pt->at(0))/(metcal->at(0)+mettrk->at(0)) << endl;
		// 	cout << "dptmetcal    " << fabs(pt->at(0)-metcal->at(0))/metcal->at(0) << endl;
		// 	cout << "dptmettrk    " << fabs(pt->at(0)-mettrk->at(0))/mettrk->at(0) << endl;
		// 	cout << "dhtmetcal    " << fabs(ht->at(0)-metcal->at(0))/metcal->at(0) << endl;
		// 	cout << "dhtmettrk    " << fabs(ht->at(0)-mettrk->at(0))/mettrk->at(0) << endl;
		// 	cout << "metsdphi     " << metsdphi  ->at(0) << endl;
		// 	cout << "dphical      " << dphical   ->at(0) << endl;
		// 	cout << "dphitrk      " << dphitrk   ->at(0) << endl;
		// 	cout << "mtcal        " << mtcal     ->at(0) << endl;
		// 	cout << "mttrk        " << mttrk     ->at(0) << endl;
		// 	
		// 	cout << "iso003       " << iso003    ->at(0) << endl;
		// 	cout << "iso014       " << iso014    ->at(0) << endl;
		// 	cout << "iso020       " << iso020    ->at(0) << endl;
		// 	cout << "iso030       " << iso030    ->at(0) << endl;
		// 	
		// 	cout << "lxy          " << lxy       ->at(0) << endl;
		// 	cout << "dlxy         " << dlxy      ->at(0) << endl;
		// 	cout << "lxysig       " << lxysig    ->at(0) << endl;
		// 	cout << "a0xy         " << a0xy      ->at(0) << endl;
		// 	cout << "da0xy        " << da0xy     ->at(0) << endl;
		// 	cout << "a0xysig      " << a0xysig   ->at(0) << endl;
		// 	
		// 	cout << "pvntrk       " << pvntrk    ->at(0) << endl;
		// }

		//// add this candidtate to the set
		if(pass)
		{
			ncandidates++;
			*score = bdtscore->at(0);
			*m3body = mass->at(0);
			data_m3body->add(RooArgSet(*m3body));
			if(hScore)
			{
				if(type=="bkg")
				{
					// the BDT for the data should be always only in the SB, even if unblinded
					if(!(mass->at(0)>xBlindMin && mass->at(0)<xBlindMax))
					{
						data_score->add(RooArgSet(*score));
						hScore->Fill(bdtscore->at(0));
					}
				}
				else
				{
					// the BDT for the signal should be always filled
					data_score->add(RooArgSet(*score));
					hScore->Fill(bdtscore->at(0));
				}
			}
			if(hMass)  hMass->Fill(mass->at(0));
		}
	}
	return ncandidates;
}

void fitParErr(Double_t xSBmin=0, Double_t xSBmax=0, Double_t xbdtcutoff=-1, Double_t xbdtmaxuser=+1)
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
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadBottomMargin(0.14);
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
	
	
	
	
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tD = (TTree*)fmvaout->Get("fltmva_Data");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_200k_3mu");
	
	makeAtlasLabel();
	makeLegend();
	
	//// for systematics (BDT range) - only for the nominal sidebands !
	if(xbdtcutoff>-1  &&  (xSBmin==mSideBandLeftLowerMeVGlob && xSBmax==mSideBandRightUpperMeVGlob))
	{
		_INF(1,"New BDT range: "<< xbdtcutoff << "-->" << xbdtmax);
		xbdtmin = xbdtcutoff;
	}
	else xbdtcutoff = minBDTcut;
	
	//// for systematics (SB range) - only with the nominal bdtcutoff
	if((xSBmin>0  && xSBmax>0) && xbdtcutoff==minBDTcut)
	{
		_INF(1,"New sidebands: "<< xSBmin << "-->" << xSBmax);
		m3bodyMin = xSBmin;
		m3bodyMax = xSBmax;
	}
	
	TString sidebands = tstr(m3bodyMin,0)+"-"+tstr(xbmin,0)+"-"+tstr(xbmax,0)+"-"+tstr(m3bodyMax,0);
	TString bdtcutoff = tstr(xbdtcutoff,2); bdtcutoff.ReplaceAll("-","neg"); if(xbdtcutoff>0) bdtcutoff="pos"+bdtcutoff; if(xbdtcutoff==0) bdtcutoff="zro"+bdtcutoff; bdtcutoff.ReplaceAll(".",""); 
	TString bdtmaxusr = tstr(xbdtmaxuser,3); bdtmaxusr="pos"+bdtmaxusr; bdtmaxusr.ReplaceAll(".",""); 
	
	
	
	
	float Sigma = SigmaRhoOmegaPhiMeV;
	Int_t nm3bodybins = (Int_t)((m3bodyMax-m3bodyMin)/m3bodyBinSize);
	// float xbdtmaxtmp = getXBDTMax(tD,m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, Sigma, "bkg","CR0",true);
	
	randGen = new TRandom();
	randGen->SetSeed(0); // Note that the machine clock is returned with a precision of 1 second.
						 // If one calls SetSeed(0) within a loop and the loop time is less than 1s,
						 // all generated numbers will be identical!
	
	// RooRealVar* scoreS = new RooRealVar("scoreS","BDT score",xbdtmin,xbdtmax);
	RooRealVar* scoreS = new RooRealVar("score","BDT score",xbdtmin,xbdtmax);
	RooRealVar* score = new RooRealVar("score","BDT score",xbdtmin,xbdtmax);
	RooRealVar* m3bodyS = new RooRealVar("m3bodyS","#it{m}_{3#mu} [MeV]",m3bodyMin,m3bodyMax);
	RooRealVar* m3body = new RooRealVar("m3body","#it{m}_{3#mu} [MeV]",m3bodyMin,m3bodyMax);
	
	
	// Roo Data holder for the BDT score
	RooAbsData*  UnbinnedDataSet_score = new RooDataSet("data_score","data_score",RooArgSet(*score));
	RooAbsData*  UnbinnedDataSet_score_noblind = new RooDataSet("data_score_noblind","data_score_noblind",RooArgSet(*score));
	RooAbsData*  UnbinnedDataSet_score_noblind_opt = new RooDataSet("data_score_noblind_opt","data_score_noblind_opt",RooArgSet(*score));
	RooAbsData*  UnbinnedDataSet_scoreS = new RooDataSet("signal_score","signal_score",RooArgSet(*scoreS));

	
	// Roo Data holder for the 3body mass
	RooAbsData*  UnbinnedDataSet_m3body = new RooDataSet("data_m3body","data_m3body",RooArgSet(*m3body));
	RooAbsData*  UnbinnedDataSet_m3body_noblind = new RooDataSet("data_m3body_noblind","data_m3body_noblind",RooArgSet(*m3body));
	RooAbsData*  UnbinnedDataSet_m3body_noblind_opt = new RooDataSet("data_m3body_noblind_opt","data_m3body_noblind_opt",RooArgSet(*m3body));
	RooAbsData*  UnbinnedDataSet_m3bodyS = new RooDataSet("signal_m3body","signal_m3body",RooArgSet(*m3bodyS));
	
	
	//// for the final estimation: xbdtmaxuser=optimal cut
	//// for the inputs for optimizations, xbdtmaxuser = +1
	//// This is relevant for the SB fit only and lot for the BDT fit.
	//// For the final estimation this has to be called agin later with xbdtmaxuser instead of xbdtmax !
	
	TString mytitle = "Events"; // "Events / "+tstr((m3bodyMax-m3bodyMin)/nm3bodybins,0)+" MeV";
	TString sytitle = "Events"; // "Events / "+tstr((xbdtmax-xbdtmin)/nbdtbins,3);
	
	TH1* hScoreSig  = new TH1F("s",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreSig->Sumw2();
	TH1* hM3bodySig = new TH1F("m",";#it{m}_{3#mu} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodySig->Sumw2();
	TH1* hScoreBkg  = new TH1F("s",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreBkg->Sumw2();
	TH1* hM3bodyBkg = new TH1F("m",";#it{m}_{3#mu} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodyBkg->Sumw2();
	TH1* hScoreBkg_noblind  = new TH1F("s_noblind",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreBkg_noblind->Sumw2();
	TH1* hM3bodyBkg_noblind = new TH1F("#it{m}_noblind",";#it{m}_{3#mu} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodyBkg_noblind->Sumw2();
	TH1* hScoreBkg_noblind_opt  = new TH1F("s_noblind_opt",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreBkg_noblind_opt->Sumw2();
	TH1* hM3bodyBkg_noblind_opt = new TH1F("m_noblind_opt",";#it{m}_{3#mu} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodyBkg_noblind_opt->Sumw2();
	
	// Int_t ncandidatesS = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
	// 							m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "full","",true,hScoreSig,hM3bodySig);
	Int_t ncandidatesS = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
								m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "full","",true,hScoreSig,hM3bodySig);
	Int_t ncandidatesDfull = readData(tD,score,m3body,UnbinnedDataSet_score,UnbinnedDataSet_m3body,
								m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "bkg","CR0",blinded,hScoreBkg,hM3bodyBkg);

	// return;

	Int_t ncandidatesDfullNoBlind = readData(tD,score,m3body,UnbinnedDataSet_score_noblind,UnbinnedDataSet_m3body_noblind,
											m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "bkg","CR0", false, hScoreBkg_noblind, hM3bodyBkg_noblind);
	Int_t ncandidatesDoptNoBlind = readData(tD,score,m3body,UnbinnedDataSet_score_noblind_opt,UnbinnedDataSet_m3body_noblind_opt,
										m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtopt,xbdtmax, "bkg","CR0", false, hScoreBkg_noblind_opt, hM3bodyBkg_noblind_opt);
										
	
	float scaleScore = Sum(hScoreBkg)/Sum(hScoreSig);
	hScoreSig->Scale(scaleScore);
	hScoreSig->SetMinimum(1e-5);
	hScoreSig->SetMaximum(25);
	
	float scaleM3body = Sum(hM3bodyBkg)/Sum(hM3bodySig);
	hM3bodySig->Scale(scaleM3body);
	hM3bodySig->SetMinimum(1e-5);
	hM3bodySig->SetMaximum(25);
										
	// hScoreSig->Scale(1./200.); hScoreSig->SetMaximum(30);
	// hM3bodySig->Scale(1./115.); hM3bodySig->SetMinimum(1e-5); hM3bodySig->SetMaximum(30); // hM3bodySig->SetMaximum(120);
	
	scoreS->setBins(nbdtbins);
	
	RooDataHist* BinnedDataSet_score; // set later
	RooDataHist* BinnedDataSet_scoreS; // set later
	RooDataHist* BinnedDataSet_m3body; // set later
	RooDataHist* BinnedDataSet_m3bodyS; // set later
	
	
	TCanvas* cnv;
	
	
	
	
	//////////////////
	//// Fit the score
	RooDataSet* rds_score = (RooDataSet*)UnbinnedDataSet_score;
	BinnedDataSet_score = rds_score->binnedClone();
	RooDataSet* rds_scoreS = (RooDataSet*)UnbinnedDataSet_scoreS;
	BinnedDataSet_scoreS = rds_scoreS->binnedClone();
	
	scoreS->setRange("range_score_signal",xbdtmin,xbdtmax);
	score->setRange("range_score",xbdtmin,xbdtmax);
	score->setRange("range_score_zoom",0.7,xbdtmax);
	score->setRange("range_score_full",xbdtmin,xbdtmax);
	score->setRange("range_score_SB0",xbdtmin,xbdtopt);
	score->setRange("range_score_SB1",xbdtopt,xbdtmax);
	score->setBins(nbdtbins);
	
	
		
	//// nominal
	float C0min=0.;    float C0max=+1000.;
	float C1min=0.;    float C1max=+1000.;
	float C2min=-1.5;  float C2max=0.;
	float C3min=0.;    float C3max=+1000.;
	float C4min=1.;    float C4max=+10.;
	// float C0min=0.;     float C0max=+2000.;
	// float C1min=0.;     float C1max=+20000.;
	// float C2min=-10.;   float C2max=0.;
	// float C3min=0.;     float C3max=+20000.;
	// float C4min=0.;     float C4max=+500.;
	
	RooRealVar* C0 = new RooRealVar("C0","C0",randomizeItialGuess(C0min,C0max),C0min,C0max); C0->setError(0.000001);
	RooRealVar* C1 = new RooRealVar("C1","C1",randomizeItialGuess(C1min,C1max),C1min,C1max); C1->setError(0.000001);
	RooRealVar* C2 = new RooRealVar("C2","C2",randomizeItialGuess(C2min,C2max),C2min,C2max); C2->setError(0.000001);
	RooRealVar* C3 = new RooRealVar("C3","C3",randomizeItialGuess(C3min,C3max),C3min,C3max); C3->setError(0.000001);
	RooRealVar* C4 = new RooRealVar("C4","C4",randomizeItialGuess(C4min,C4max),C4min,C4max); C4->setError(0.000001);
	// RooRealVar* C0 = new RooRealVar("C0","C0",randomizeItialGuess(0.,+1000.),0.,+1000.); C0->setError(0.000001);
	// RooRealVar* C1 = new RooRealVar("C1","C1",randomizeItialGuess(0.,1000),0.,1000.);    C1->setError(0.000001);
	// RooRealVar* C2 = new RooRealVar("C2","C2",randomizeItialGuess(-1.5,0.),-1.5,0.);     C2->setError(0.000001);
	// RooRealVar* C3 = new RooRealVar("C3","C3",randomizeItialGuess(0.,1000),0.,1000.);    C3->setError(0.000001);
	// RooRealVar* C4 = new RooRealVar("C4","C4",randomizeItialGuess(1.,+10.),1.,+10.);     C4->setError(0.000001);
	TString bkgBDTpdfFormula1 = "C0+C1*pow((score+1),C2)+C3*pow((score+1),C4)";
	TString bkgBDTpdfFormulaX = "A0+A1*pow((score+1),A2)+A3*pow((score+1),A4)";
	RooGenericPdf* bkgBDTpdf1 = new RooGenericPdf("bkgBDTpdf1","bkgBDTpdf1",bkgBDTpdfFormula1,RooArgSet(*score,*C0,*C1,*C2,*C3,*C4));
	RooFitResult* fitresult_score1 = bkgBDTpdf1->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	TMinuit* gFit_score1 = gMinuit;
	bool FitStatus_score1 = minuitStatus(gFit_score1);
	fitresult_score1->Print("v");
	
	RooAbsPdf* bkgBDTpdfNominal0 = bkgBDTpdf1; TString bkgBDTpdfNominalFormula0 = bkgBDTpdfFormula1; const unsigned int nparBDT0 = 5;
	
	RooAbsReal* integralSB0 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB1 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float R0 = integralSB1->getVal()/integralSB0->getVal();
	
	
	
	cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* scoreFrame = score->frame(Name("scoreFrame"),Title("Background BDT response in the sidebands: "+sidebands+" and BDT>"+tstr(xbdtcutoff,2)));
	scoreFrame->SetMaximum(20);
	UnbinnedDataSet_score->plotOn(scoreFrame,Name("BDT score"),MarkerSize(1.2),Binning(nbdtbins));
	bkgBDTpdfNominal0->plotOn(scoreFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),NormRange("range_score"));
	Double_t chi2BDTpdfNominal = scoreFrame->chiSquare("Nominal","BDT score", nparBDT0);
	TF1* pdf = (TF1*)bkgBDTpdfNominal0->asTF(*score, RooArgList(*C0,*C1,*C2,*C3,*C4), *score);
	
	// Sample dataset with parameter values according to distribution of covariance matrix of fit result
	RooRealVar* A0 = new RooRealVar("A0","A0",C0->getVal(),C0min,C0max); 
	RooRealVar* A1 = new RooRealVar("A1","A1",C1->getVal(),C1min,C1max);       
	RooRealVar* A2 = new RooRealVar("A2","A2",C2->getVal(),C2min,C2max);         
	RooRealVar* A3 = new RooRealVar("A3","A3",C3->getVal(),C3min,C3max);       
	RooRealVar* A4 = new RooRealVar("A4","A4",C4->getVal(),C4min,C4max);     
	TH1F* hChi2 = new TH1F("hChi2","hChi2",200,C0min,C0max);
	TH1F* hA0 = new TH1F("hA0",";A0;N_{toys}",200,C0min,C0max);
	TH1F* hA1 = new TH1F("hA1",";A1;N_{toys}",200,C1min,C1max);
	TH1F* hA2 = new TH1F("hA2",";A2;N_{toys}",200,C2min,C2max);
	TH1F* hA3 = new TH1F("hA3",";A3;N_{toys}",200,C3min,C3max);
	TH1F* hA4 = new TH1F("hA4",";A4;N_{toys}",200,C4min,C4max);
	vector<float> R,c0,c1,c2,c3,c4,chi2i,pdfneg09;
	for(Int_t i=0 ; i<10000 ; i++)
	{
		const RooArgList randPars = fitresult_score1->randomizePars();
		TIterator* iter = randPars.createIterator();
		RooRealVar* par;
		
		TString pdfname = "bkgBDTpdfX"+tstr(i,0);
		
		while( par = (RooRealVar*)iter->Next() )
		{
			TString parname = par->GetName();
			
			if(parname=="C0") { A0->setVal( par->getVal() ); }
			if(parname=="C1") { A1->setVal( par->getVal() ); }
			if(parname=="C2") { A2->setVal( par->getVal() ); }
			if(parname=="C3") { A3->setVal( par->getVal() ); }
			if(parname=="C4") { A4->setVal( par->getVal() ); }
		}
		delete iter;
		
		RooGenericPdf bkgBDTpdfX(pdfname.Data(),pdfname.Data(),bkgBDTpdfFormulaX,RooArgSet(*score,*A0,*A1,*A2,*A3,*A4));
		bkgBDTpdfX.plotOn(scoreFrame,Name(pdfname.Data()),LineWidth(1),LineColor(kRed),NormRange("range_score"));
		Double_t chi2BDTpdf = scoreFrame->chiSquare(pdfname.Data(),"BDT score",nparBDT0);
		
		TF1* pdfX = (TF1*)bkgBDTpdfX.asTF(*score, RooArgList(*A0,*A1,*A2,*A3,*A4), *score);
		if(pdfX->Eval(-0.9)<pdfX->Eval(+1) || pdfX->Eval(-0.9)<0.)
		{
			scoreFrame->remove(pdfname.Data());
			continue;
		}

		chi2i.push_back(chi2BDTpdf);
		c0.push_back(A0->getVal());
		c1.push_back(A1->getVal());
		c2.push_back(A2->getVal());
		c3.push_back(A3->getVal());
		c4.push_back(A4->getVal());
		
		hChi2->Fill(chi2BDTpdf);
		hA0->Fill(A0->getVal());
        hA1->Fill(A1->getVal());
		hA2->Fill(A2->getVal());
        hA3->Fill(A3->getVal());
        hA4->Fill(A4->getVal());

		

		// bool isbad = false;
		// for(float xx=-0.9 ; xx<=+1.0 ; xx+=0.01)
		// {
		// 	if(pdfX->Eval(xx)<=0.1) { isbad=true; break; }
		// }
		// if(isbad)
		// {
		// 	scoreFrame->remove(pdfname.Data());
		// 	continue;
		// }
		
		
		// if(chi2BDTpdf>0.5) //chi2BDTpdfNominal)
		// {
		// 	scoreFrame->remove(pdfname.Data());
		// 	continue;
		// }
		
		RooAbsReal* integralSB0 = bkgBDTpdfX.createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
		RooAbsReal* integralSB1 = bkgBDTpdfX.createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
		float Rx = integralSB1->getVal()/integralSB0->getVal();
		R.push_back(Rx);
	}
	
	// for(Int_t i=0 ; i<10000 ; i++)
	// {
	// // 	cout << "c0=" << tstr(c0[i],2) << ", c1=" << tstr(c1[i],3) << ", c2=" << tstr(c2[i],6) << ", c3=" << tstr(c3[i],3) << ", c4=" << tstr(c4[i],6) << ", chi2=" << tstr(chi2i[i],3) << endl;
	// }

	float dR = maxdistance(R,R0);
	float RMSD = rmsd(R,R0,true);
	float RMSDtrue = rmsd(R);
	cout << "chi^2(nom)=" << chi2BDTpdfNominal << endl;
	cout << "R0=" << R0 << endl;
	cout << "dR=" << dR << endl;
	cout << "dR/R0=" << tstr(dR/R0*100,1) << "%" << endl;
	cout << "RMSD(R,R0) =" << tstr(RMSD,5) << endl;
	cout << "RMSD(R,R)  =" << tstr(RMSDtrue,5) << endl;
	cout << "RMSD(R,R0)/R0 =" << tstr(RMSD/R0*100,1) << "%" << endl;
	cout << "RMSD(R,R)/R0  =" << tstr(RMSDtrue/R0*100,1) << "%" << endl;
	/*
	chi^2(nom)=0.439794
	R0=0.0264858
	dR=0.0216681
	dR/R0=81.8%
	RMSD(R)=0.00641
	RMSD(R)/R0=24.2%
	*/
	/*
	chi^2(nom)=0.438588
	R0=0.0266285
	dR=0.0597689
	dR/R0=224.5%
	RMSD(R)=0.00453
	RMSD(R)/R0=17.0%
	*/
	/*
	chi^2(nom)=0.438654
	R0=0.0266178
	dR=0.0270291
	dR/R0=101.5%
	RMSD(R)=0.00769
	RMSD(R)/R0=28.9%
	*/
	
	UnbinnedDataSet_score->plotOn(scoreFrame,Name("BDT score"),MarkerSize(1.2),Binning(nbdtbins));
	bkgBDTpdfNominal0->plotOn(scoreFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),NormRange("range_score"));
	
	leg->Clear();
	leg->AddEntry(scoreFrame->findObject("BDT score"),"Sidebands data","ple");
	leg->AddEntry(scoreFrame->findObject("Nominal"),"Nominal fit","l");
	leg->AddEntry(scoreFrame->findObject("bkgBDTpdfX0"),"Random pars","l");
	scoreFrame->SetTitleOffset(2,"Y");
	scoreFrame->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimateErrs.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrs.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf(");
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	hChi2->Draw();
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/BkgEstimateErrs.Chi2."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrs.Chi2."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf(");
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",1200,900);
	cnv->Divide(3,2);
	cnv->cd(1); hA0->Draw();   gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); gPad->SetLogx();
	cnv->cd(2); hA1->Draw();   gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); gPad->SetLogx();
	cnv->cd(3); hA2->Draw();   gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); //gPad->SetLogx();
	cnv->cd(4); hA3->Draw();   gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); gPad->SetLogx();
	cnv->cd(5); hA4->Draw();   gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); gPad->SetLogx();
	cnv->cd(6); hChi2->Draw(); gPad->Update(); gPad->RedrawAxis(); gPad->SetTicks(1,1); gPad->SetLogy(); gPad->SetLogx();
	cnv->Update();
	cnv->SaveAs("figures/BkgEstimateErrs.AiChi2."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrs.AiChi2."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf(");
	
	return;
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* scoreFrameZoom = score->frame(Name("scoreFrameZoom"),Title("Background BDT response in the sidebands: "+sidebands+" and BDT>"+tstr(xbdtcutoff,2)),Range("range_score_zoom"));
	UnbinnedDataSet_score->plotOn(scoreFrameZoom,Name("BDT score zoom"),MarkerSize(1.2),Binning(nbdtbins),Range("range_score_zoom"));
	bkgBDTpdfNominal0->plotOn(scoreFrameZoom,Name("Nominal"),LineWidth(2),LineColor(kBlue),Range("range_score"));
	leg->Clear();
	leg->AddEntry(scoreFrameZoom->findObject("BDT score zoom"),"Sidebands data","pl");
	leg->AddEntry(scoreFrameZoom->findObject("Nominal"),"Nominal fit","l");
	
	scoreFrameZoom->SetMaximum(10);
	scoreFrameZoom->SetTitleOffset(2,"Y");
	scoreFrameZoom->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimateErrs.BDTzoom."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrs.BDTzoom."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	
	
	
	
	TFile* fUncert = new TFile("Uncert.root","READ");
	TH1F* hdBDT = (TH1F*)((TH1F*)fUncert->Get("hdBDT"))->Clone();
	hdBDT->SetFillColor(kBlue);
	hdBDT->SetLineColor(kBlue);
	hdBDT->SetFillStyle(3344);
	TGraphErrors* grBDTfitErr = (TGraphErrors*)inflate(hdBDT)->Clone();
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* scoreFrameS = score->frame(Name("scoreFrameS"),Title("BDT response"));
	UnbinnedDataSet_score->plotOn(scoreFrameS,Name("BDT score background"),MarkerSize(1.2),Binning(nbdtbins));
	// UnbinnedDataSet_score_noblind->plotOn(scoreFrameS,Name("BDT score background noblind"),MarkerSize(1.2),MarkerColor(kRed),Binning(nbdtbins))
	// UnbinnedDataSet_scoreS->plotOn(scoreFrameS,Name("BDT score signal"),LineStyle(kDashed),LineColor(kRed),MarkerStyle(22),MarkerColor(kRed),Binning(nbdtbins));
	hScoreSig->SetLineColor(kGray); hScoreSig->SetFillColor(kGray);
	hScoreBkg_noblind->SetLineColor(kRed); hScoreBkg_noblind->SetLineStyle(2);
	bkgBDTpdfNominal0->plotOn(scoreFrameS,Name("Nominal"),LineWidth(2),LineColor(kBlue),NormRange("range_score"));
	
	// delete leg;
	leg = new TLegend(0.15,0.68,0.55,0.92,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(scoreFrameS->findObject("BDT score background"),"SB data (tight+x>x_{0} selection)","ple");
	leg->AddEntry(scoreFrameS->findObject("Nominal"),"Fit to the SB data","l");
	leg->AddEntry(hdBDT,"Fit uncertainty","F");
	leg->AddEntry(hScoreSig,"Signal (tight+x>x_{0} selection)","f");
	
	cnv->SetTicks(1,1);
	hScoreSig->Draw("hist");
	grBDTfitErr->Draw("3 same");
	scoreFrameS->Draw("same");
	plotAtlasLabel();
	leg->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	
	
	Double_t chi2BDTpdf0 = scoreFrame->chiSquare("Nominal","BDT score",      nparBDT0);
	cout << "chi2BDTpdf0=" << chi2BDTpdf0 << endl;
	
	
	_INF(1,"------------------- Done fitting the BDT shape -------------------");







	
	
	/////////////////////////////
	//// Fit the m3body sidebands
	//////////////////////////////
	
	
	//// for the final estimation: xbdtmaxuser=optimal cut
	//// for the inputs for optimizations, xbdtmaxuser = +1
	//// This is relevant for the SB fit only and lot for the BDT fit.
	//// For the final estimation this has to be called here with xbdtmaxuser (not xbdtmax !)
	Int_t ncandidatesD = readData(tD,score,m3body,UnbinnedDataSet_score,UnbinnedDataSet_m3body,
								m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmaxuser, "bkg","CR0");
	
	
	RooDataSet* rds_m3body = (RooDataSet*)UnbinnedDataSet_m3body;
	BinnedDataSet_m3body = rds_m3body->binnedClone();
	
	m3body->setRange("range_m3body",m3bodyMin,m3bodyMax);
	m3body->setBins(nm3bodybins);
	m3body->setRange("range_SBleft",  m3bodyMin,xbmin);
	m3body->setRange("range_SBright", xbmax,m3bodyMax);
	m3body->setRange("range_blinded", xbmin,xbmax);
	

	RooRealVar* k0 = new RooRealVar("k0","k0",1.,0.,100.); k0->setError(0.00001);
	RooRealVar* k1 = new RooRealVar("k1","k1",1.,0.,100.); k1->setError(0.00001);
	RooBernstein* bkgm3bodypdf11 = new RooBernstein("bkgm3bodypdf11","bkgm3bodypdf11",*m3body,RooArgSet(*k0,*k1));

	RooFitResult* fitresult_m3body11;
	fitresult_m3body11 = bkgm3bodypdf11->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body11 = gMinuit;
	bool FitStatus_m3body11 = minuitStatus(gFit_m3body11);
	fitresult_m3body11->Print("v");
	
	RooAbsPdf* bkgm3bodypdfNominal0 = bkgm3bodypdf11; Int_t nparSB0 = 2;
	
	








	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* m3bodyFrame = m3body->frame(Name("m3bodyFrame"),Title("3body mass in the sidebands of CR_{0} (BDT>"+tstr(xbdtcutoff,2)+")"));
	m3bodyFrame->SetMaximum(20);
	if(!blinded) UnbinnedDataSet_m3body_noblind->plotOn(m3bodyFrame,Name("m3body SB0"),MarkerSize(1.2),Binning(nm3bodybins));
	else         UnbinnedDataSet_m3body->plotOn(m3bodyFrame,Name("m3body SB0"),MarkerSize(1.2),Binning(nm3bodybins));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	
	leg->Clear();
	leg->AddEntry(m3bodyFrame->findObject("m3body SB0"),"Data in SB_{0}","pl");
	leg->AddEntry(m3bodyFrame->findObject("Nominal"),"Nominal fit","l");
	
	m3bodyFrame->SetTitleOffset(2,"Y");
	m3bodyFrame->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimateErrs.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrs.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");

	
	///////////
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* m3bodyFrameS = m3body->frame(Name("m3bodyFrameS"),Title("3body mass in the sidebands (BDT>"+tstr(xbdtcutoff,2)+")"));
	if(blinded) UnbinnedDataSet_m3body->plotOn(m3bodyFrameS,Name("m3body background"),MarkerSize(1.2),Binning(nm3bodybins));
	else        UnbinnedDataSet_m3body_noblind->plotOn(m3bodyFrameS,Name("m3body background"),MarkerSize(1.2),Binning(nm3bodybins));
	hM3bodySig->SetLineColor(kGray); hM3bodySig->SetFillColor(kGray);
	hM3bodyBkg_noblind_opt->SetLineColor(kRed); hM3bodyBkg_noblind_opt->SetMarkerStyle(20); hM3bodyBkg_noblind_opt->SetMarkerColor(kRed);
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_left"),LineWidth(2),LineColor(kBlue),Range("range_SBleft"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_right"),LineWidth(2),LineColor(kBlue),Range("range_SBright"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_blinded"),LineWidth(2)/*,LineStyle(2)*/,LineColor(kBlue),Range("range_blinded"),NormRange("range_SBleft,range_SBright"));
	if(!blinded) UnbinnedDataSet_m3body_noblind_opt->plotOn(m3bodyFrameS,Name("m3body background noblind opt"),MarkerSize(1.2),MarkerColor(kRed+1),LineColor(kRed+1),MarkerStyle(21),Binning(nm3bodybins));
	m3bodyFrameS->SetMinimum(1e-5);
	
	TH1F* hdSB = (TH1F*)((TH1F*)fUncert->Get("hdSB"))->Clone();
	hdSB->SetFillColor(kBlue);
	hdSB->SetLineColor(kBlue);
	hdSB->SetFillStyle(3344);
	TGraphErrors* grSBfitErr = (TGraphErrors*)inflate(hdSB)->Clone();
	
	delete leg;
	leg = new TLegend(0.15,0.64,0.55,0.92,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(m3bodyFrameS->findObject("m3body background"),"Data (tight+x>x_{0} selection)","ple");
	if(!blinded) leg->AddEntry(m3bodyFrameS->findObject("m3body background noblind opt"),"Data (tight+x>x_{1} selection)","ple");
	leg->AddEntry(m3bodyFrameS->findObject("Nominal_left"),"Fit to the SB data","l");
	leg->AddEntry(hdSB,"Fit uncertainty","F");
	leg->AddEntry(hM3bodySig,"Signal (tight+x>x_{0} selection)","f");
	
	cnv->SetTicks(1,1);
	cnv->Update();
	hM3bodySig->Draw("hist");
	grSBfitErr->Draw("3 same");
	m3bodyFrameS->Draw("same");
	plotAtlasLabel();
	leg->Draw("same");
	
	Double_t ymax = hM3bodySig->GetMaximum()*0.85;//cnv->GetUymax();
	Double_t y = -1;
	
	y = ymax/1.5;
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
	// TLine* lBlinded = new TLine(1690,y,1870,y);
	// lBlinded->SetLineColor(kBlack);
	// lBlinded->SetLineWidth(2);
	// lBlinded->SetLineStyle(3);
	// lBlinded->Draw("L same");
	TArrow* lSBleft = new TArrow(1690,y,1630,y,0.02);
	lSBleft->SetLineColor(kBlack);
	lSBleft->SetLineWidth(2);
	lSBleft->SetLineStyle(3);
	lSBleft->Draw();
	TArrow* lSBright = new TArrow(1870,y,1930,y,0.02);
	lSBright->SetLineColor(kBlack);
	lSBright->SetLineWidth(2);
	lSBright->SetLineStyle(3);
	lSBright->Draw();
	
	y = ymax/1.7;
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
	TArrow* lSRleft = new TArrow(1713,y,1750,y,0.02);
	lSRleft->SetLineColor(kRed);
	lSRleft->SetLineWidth(2);
	lSRleft->SetLineStyle(5);
	lSRleft->Draw();
	TArrow* lSRright = new TArrow(1841,y,1810,y,0.02);
	lSRright->SetLineColor(kRed);
	lSRright->SetLineWidth(2);
	lSRright->SetLineStyle(5);
	lSRright->Draw();
	
	// delete legR;
	legR = new TLegend(0.65,0.63,0.90,0.73,NULL,"brNDC");
	legR->SetFillStyle(4000); //will be transparent
	legR->SetFillColor(0);
	legR->SetTextFont(42);
	legR->SetBorderSize(0);
	legR->AddEntry(lSBright,"Sidebands (SB)","l");
	legR->AddEntry(lSRright,"Signal region","l");
	legR->Draw("same");
	
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateErrsWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	cnv->SaveAs("figures/BkgEstimateErrs."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf)");
	
	Double_t chi2SBpdf0 = m3bodyFrame->chiSquare("Nominal",      "m3body SB0",nparSB0); // last parameter is Ndof for the chi^2
	cout << "chi2SBpdf0=" << chi2SBpdf0 << endl;
	
	
	
	
	
	Float_t xSRmin = 1713;
	Float_t xSRmax = 1841;
	m3body->setRange("range_m3body_SR",xSRmin,xSRmax);
	
	

	RooAbsReal* integralSB00 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB10 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float transferFactor0 = integralSB10->getVal()/integralSB00->getVal();
	float nSB00bdt = integralSB00->getVal();

	RooAbsReal* integralSR00      = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB00left  = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB00right = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBright"));
	float nSB00 = integralSB00left->getVal()+integralSB00right->getVal();
	float nSR00 = integralSR00->getVal();
	float dnSR00 = (nSR00/nSB00)*sqrt(nSB00);
	float RooFitScale0 = ncandidatesD/nSB00;
	
	
	
	
	
	
	int ncandidatesSsr1 = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
								   m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, "sig","CR1");
	float efficiencyS = ((float)ncandidatesSsr1/200000.)*100.;
	
	int ncandidatesD_CR1 = readData(tD,score,m3body,UnbinnedDataSet_score,UnbinnedDataSet_m3body,
								    m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, "bkg","CR1");
	
	
	//// Write the TF1 object and more numbers to a root file
	TFile* fpdf = new TFile("BDTfitErrsResults."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".root","RECREATE");
	TF1* pdf0 = (TF1*)bkgBDTpdfNominal0->asTF(*score, RooArgList(*C0,*C1,*C2,*C3,*C4), *score); RooRealVar parBDT0[nparBDT0]={*C0,*C1,*C2,*C3,*C4}; TString pdfStr0 = expandNscaleFormula(bkgBDTpdfNominalFormula0, "score", parBDT0, nparBDT0, RooFitScale0, 0);
	pdf0->Write("fBDT0.0");
	
	cout << "fBDT0=" << pdfStr0 << endl;
	TNamed* bkgBDTpdfNominal0Str = new TNamed("fBDT0str.0",pdfStr0.Data());
	bkgBDTpdfNominal0Str->Write("fBDT0str.0");
	
	TF1* fSB0 = (TF1*)bkgm3bodypdfNominal0->asTF(*m3body, RooArgList(*k0,*k1), *m3body);
	fSB0->Write("fSB0.0");
	
	cout << "fSB0=" << fSB0->GetExpFormula("P") << endl;
	TString sbStr0 = fSB0->GetExpFormula("P")+"*"+tstr(RooFitScale0,10);
	TNamed* bkgm3bodypdfNominal0Str = new TNamed("fSB0str.0",sbStr0.Data());
	bkgm3bodypdfNominal0Str->Write("fSB0str.0");
	
	TVectorF chi2BDT(1);
	chi2BDT[0]=chi2BDTpdf0;
	chi2BDT.Write("chi2BDT");
	
	TVectorF NSB(1);
	NSB[0]=ncandidatesD;
	NSB.Write("nSB");
	
	TVectorF NSR(1);
	TVectorF chi2SB(1);
	NSR[0]=RooFitScale0*nSR00; chi2SB[0]=chi2SBpdf0;
	NSR.Write("nSR");
	chi2SB.Write("chi2SB");
	
	TVectorF leftSBmin(1);  leftSBmin[0]=m3bodyMin;  leftSBmin.Write("leftSBmin");
	TVectorF leftSBmax(1);  leftSBmax[0]=xbmin;      leftSBmax.Write("leftSBmax");
	TVectorF rightSBmin(1); rightSBmin[0]=xbmax;     rightSBmin.Write("rightSBmin");
	TVectorF rightSBmax(1); rightSBmax[0]=m3bodyMax; rightSBmin.Write("rightSBmax");
	
	TVectorF AbsBDTmin(1);     AbsBDTmin[0]=xbdtmin; AbsBDTmin.Write("absBDTmin");
	TVectorF AbsBDTmax(1);     AbsBDTmax[0]=xbdtmax; AbsBDTmax.Write("absBDTmax");
	TVectorF AbsBDTmaxUser(1); AbsBDTmaxUser[0]=xbdtmaxuser; AbsBDTmaxUser.Write("absBDTmaxUser");
	TVectorF NBDTbins(1);      NBDTbins[0]=nbdtbins; NBDTbins.Write("nBDTbins");
	
	// and to read the file do: TVectorF* NSB0 = (TVectorF*)fpdf->Get("nSB0"); float nSB0 = ((*NSB0))[0];
	fpdf->Write();
	fpdf->Close();
	delete fpdf;
	
	string filename = "BkgEstimateErrs."+(string)sidebands+"."+(string)bdtcutoff+"."+(string)bdtmaxusr+".txt";
	ofstream* ofstr = new ofstream(filename.c_str());
	stringstream strm; strm << "";
	string tmp = "";
	string summary = "";
	summary += "========================================================================================\n";
	summary += "m3body left sideband    ["+str(m3bodyMin,0)+","+str(xbmin,0)+"] MeV\n";
	summary += "m3body right sideband   ["+str(xbmax,0)+","+str(m3bodyMax,0)+"] MeV\n";
	summary += "Minimal BDT cut (x0)    "+str(xbdtmin,3)+"\n";
	summary += "Optimal BDT cut (x1)    "+str(xbdtopt,3)+"\n";
	summary += "Maximal BDT cut (x2)    "+str(xbdtmax,3)+"\n";
	summary += "========================================================================================\n";
	summary += "Ch2^2/DOF BDT           "+str(chi2BDTpdf0,2)+"\n";
	summary += "Ch2^2/DOF SB            "+str(chi2SBpdf0,2)+"\n";
	summary += "========================================================================================\n";
	summary += "Count:    nSB0          "+str(ncandidatesD,0)+"\n";
	summary += "Mass fit: nSB0          "+str(nSB00,3)+"\n";
	summary += "Mass fit: nSR0          "+str(nSR00,3)+"\n";
	summary += "RooFit scale            "+str(RooFitScale0,5)+"\n";
	summary += "Scaled nSR0             "+str(RooFitScale0*nSR00,3)+" +- "+str((nSR00/nSB00)*sqrt(RooFitScale0*nSB00),3)+"\n";
	summary += "Fit: f01                "+str(transferFactor0,6)+"\n";
	summary += "Cut: nSB1               "+str(ncandidatesD_CR1,3)+" +- "+str(sqrt(ncandidatesD_CR1),3)+"\n";
	summary += "Scaled: nSR1            "+str((RooFitScale0*nSR00)*transferFactor0,3)+" +- "+str((nSR00/nSB00)*sqrt(RooFitScale0*nSB00)*transferFactor0,3)+"\n";
	summary += "Acc*Eff in SR1          "+str(efficiencyS,3)+"\%\n";
	summary += "N candidates(S)         "+str(ncandidatesS,0)+"\n";
	summary += "N candidates(D)         "+str(ncandidatesDfull,0)+"\n";
	summary += "========================================================================================\n";
	cout << summary << endl;
	(*ofstr) << summary;
	ofstr->close();
}