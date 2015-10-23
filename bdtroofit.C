//////////////////////////////////////////////////////////////////////////////////////
//// root -l -b -q bdtroofit.C++\(1450,2110,-0.9,+1\) //// inputs for optimization
//// root -l -b -q bdtroofit.C++\(1450,2110,-0.9,optBDTcut\) //// for final estimation
//////////////////////////////////////////////////////////////////////////////////////

#include "std.h"
#include "const.h"
#include "postBDTcuts.h"
#include "roofit.h"

bool blinded = true;

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
		
	
		if(pass && bdtscore->at(0)>0.85 && type=="bkg")
		{
			cout << "------- entry=" << entry << " ------" << endl;
			cout << "code         " << code->at(0) << endl;
			cout << "Ncandidates  " << bdtscore->size() << endl;
			cout << "xBDT         " << bdtscore->at(0) << endl;
			cout << "mass         " << mass->at(0) << endl;
			cout << "mOS1         " << mOS1->at(0) << endl;
			cout << "mOS2         " << mOS2->at(0) << endl;
			cout << "mSS          " << mSS->at(0) << endl;
			cout << "pval         " << pval->at(0) << endl;
			cout << "trkspval     " << trkspval->at(0) << endl;
			cout << "njets        " << njets->at(0) << endl;
			cout << "pt           " << pt->at(0) << endl;
			cout << "ht           " << ht->at(0) << endl;
			cout << "dpthtrel     " << pt->at(0)/ht->at(0)-1 << endl;
			cout << "mettrk       " << mettrk    ->at(0) << endl;
			cout << "metcal       " << metcal    ->at(0) << endl;
			cout << "dmetsrel     " << (mettrk->at(0)-pt->at(0))/(metcal->at(0)+mettrk->at(0)) << endl;
			cout << "dptmetcal    " << fabs(pt->at(0)-metcal->at(0))/metcal->at(0) << endl;
			cout << "dptmettrk    " << fabs(pt->at(0)-mettrk->at(0))/mettrk->at(0) << endl;
			cout << "dhtmetcal    " << fabs(ht->at(0)-metcal->at(0))/metcal->at(0) << endl;
			cout << "dhtmettrk    " << fabs(ht->at(0)-mettrk->at(0))/mettrk->at(0) << endl;
			cout << "metsdphi     " << metsdphi  ->at(0) << endl;
			cout << "dphical      " << dphical   ->at(0) << endl;
			cout << "dphitrk      " << dphitrk   ->at(0) << endl;
			cout << "mtcal        " << mtcal     ->at(0) << endl;
			cout << "mttrk        " << mttrk     ->at(0) << endl;
			
			cout << "iso003       " << iso003    ->at(0) << endl;
			cout << "iso014       " << iso014    ->at(0) << endl;
			cout << "iso020       " << iso020    ->at(0) << endl;
			cout << "iso030       " << iso030    ->at(0) << endl;
			
			cout << "lxy          " << lxy       ->at(0) << endl;
			cout << "dlxy         " << dlxy      ->at(0) << endl;
			cout << "lxysig       " << lxysig    ->at(0) << endl;
			cout << "a0xy         " << a0xy      ->at(0) << endl;
			cout << "da0xy        " << da0xy     ->at(0) << endl;
			cout << "a0xysig      " << a0xysig   ->at(0) << endl;
			
			cout << "pvntrk       " << pvntrk    ->at(0) << endl;
		}

		//// add this candidtate to the set
		if(pass)
		{
			ncandidates++;
			*score = bdtscore->at(0);
			*m3body = mass->at(0);
			data_score->add(RooArgSet(*score));
			data_m3body->add(RooArgSet(*m3body));
			if(hScore) hScore->Fill(bdtscore->at(0));
			if(hMass)  hMass->Fill(mass->at(0));
		}
	}
	return ncandidates;
}

void bdtroofit(Double_t xSBmin=0, Double_t xSBmax=0, Double_t xbdtcutoff=-1, Double_t xbdtmaxuser=+1)
{
	// gStyle->SetFrameBorderMode(0);
	// gStyle->SetCanvasBorderMode(0);
	// gStyle->SetPadBorderMode(0);
	// gStyle->SetPadColor(0);
	// gStyle->SetCanvasColor(0);
	// gStyle->SetFrameFillColor(0);
	// gStyle->SetTitleFillColor(0);
	// gStyle->SetPaperSize(20,26);
	// // gStyle->SetPadTopMargin(0.13);
	// // gStyle->SetPadRightMargin(0.15);
	// // gStyle->SetPadBottomMargin(0.14);
	// // gStyle->SetPadLeftMargin(0.12);
	// Int_t font=42;
	// Double_t tsize=0.04;
	// gStyle->SetTextFont(font);
	// gStyle->SetTextSize(tsize);
	// gStyle->SetLabelFont(font,"x");
	// gStyle->SetTitleFont(font,"x");
	// gStyle->SetLabelFont(font,"y");
	// gStyle->SetTitleFont(font,"y");
	// gStyle->SetLabelFont(font,"z");
	// gStyle->SetTitleFont(font,"z");
	// gStyle->SetLabelSize(tsize,"x");
	// gStyle->SetTitleSize(tsize,"x");
	// gStyle->SetLabelSize(tsize,"y");
	// gStyle->SetTitleSize(tsize,"y");
	// gStyle->SetLabelSize(tsize,"z");
	// gStyle->SetTitleSize(tsize,"z");
	// gStyle->SetStatColor(0);
	// gStyle->SetStatBorderSize(0);
	// gStyle->SetStatColor(0);
	// gStyle->SetStatX(0);
	// gStyle->SetStatY(0);
	// gStyle->SetStatFont(42);
	// gStyle->SetStatFontSize(0);
	// gStyle->SetOptStat(0);
	// gStyle->SetStatW(0);
	// gStyle->SetStatH(0);
	// gStyle->SetTitleX(0.55); //title X location 
	// gStyle->SetTitleY(0.96); //title Y location 
	// gStyle->SetTitleW(0.5); //title width 
	// gStyle->SetTitleH(0.05); //title height
	// gStyle->SetTitleBorderSize(0);
	
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
	RooRealVar* m3bodyS = new RooRealVar("m3bodyS","#it{m}_{3body} [MeV]",m3bodyMin,m3bodyMax);
	RooRealVar* m3body = new RooRealVar("m3body","#it{m}_{3body} [MeV]",m3bodyMin,m3bodyMax);
	
	
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
	TH1* hM3bodySig = new TH1F("m",";#it{m}_{3body} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodySig->Sumw2();
	TH1* hScoreBkg_noblind  = new TH1F("s_noblind",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreBkg_noblind->Sumw2();
	TH1* hM3bodyBkg_noblind = new TH1F("#it{m}_noblind",";#it{m}_{3body} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodyBkg_noblind->Sumw2();
	TH1* hScoreBkg_noblind_opt  = new TH1F("s_noblind_opt",";BDT score;"+sytitle,nbdtbins,xbdtmin,xbdtmax); hScoreBkg_noblind_opt->Sumw2();
	TH1* hM3bodyBkg_noblind_opt = new TH1F("m_noblind_opt",";#it{m}_{3body} [MeV];"+mytitle,nm3bodybins,m3bodyMin,m3bodyMax); hM3bodyBkg_noblind_opt->Sumw2();
	
	// Int_t ncandidatesS = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
	// 							m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "full","",true,hScoreSig,hM3bodySig);
	Int_t ncandidatesS = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
								m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "full","",true,hScoreSig,hM3bodySig);
	Int_t ncandidatesDfull = readData(tD,score,m3body,UnbinnedDataSet_score,UnbinnedDataSet_m3body,
								m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "bkg","CR0");

	// return;

	Int_t ncandidatesDfullNoBlind = readData(tD,score,m3body,UnbinnedDataSet_score_noblind,UnbinnedDataSet_m3body_noblind,
											m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtmax, "bkg","CR0", false, hScoreBkg_noblind, hM3bodyBkg_noblind);
	Int_t ncandidatesDoptNoBlind = readData(tD,score,m3body,UnbinnedDataSet_score_noblind_opt,UnbinnedDataSet_m3body_noblind_opt,
										m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtopt,xbdtmax, "bkg","CR0", false, hScoreBkg_noblind_opt, hM3bodyBkg_noblind_opt);
	hScoreSig->Scale(1./200.); hScoreSig->SetMaximum(30);
	hM3bodySig->Scale(1./115.); hM3bodySig->SetMinimum(1e-5); hM3bodySig->SetMaximum(30); // hM3bodySig->SetMaximum(120);
	
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
	score->setRange("range_score_SB0",xbdtmin,xbdtopt);
	score->setRange("range_score_SB1",xbdtopt,xbdtmax);
	score->setBins(nbdtbins);
	
	

	//// alternative
	RooRealVar* A0 = new RooRealVar("A0","A0",randomizeItialGuess(0.,+1000.),0.,+1000.); A0->setError(0.000001);
	RooRealVar* A1 = new RooRealVar("A1","A1",randomizeItialGuess(0.,10.),0.,10.);       A1->setError(0.000001);
	RooRealVar* A2 = new RooRealVar("A2","A2",randomizeItialGuess(-10.,0.),-10.,0.);     A2->setError(0.000001);
	RooRealVar* A3 = new RooRealVar("A3","A3",randomizeItialGuess(0.,1000.),0.,1000.);   A3->setError(0.000001);
	RooRealVar* A4 = new RooRealVar("A4","A4",randomizeItialGuess(1.,10.),1.,10.);       A4->setError(0.000001);
	TString bkgBDTpdfFormula0 = "A0+A1*exp(score*A2)*pow(score,4)+A3*pow((score+1),A4*score)";
	RooGenericPdf* bkgBDTpdf0 = new RooGenericPdf("bkgBDTpdf0","bkgBDTpdf0",bkgBDTpdfFormula0,RooArgSet(*score,*A0,*A1,*A2,*A3,*A4));
	RooFitResult* fitresult_score0 = bkgBDTpdf0->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	//// alternative 
	RooRealVar* B0 = new RooRealVar("B0","B0", randomizeItialGuess(0.0, 1000.0),   0.0, 1000.0); B0->setError(0.00001);
	RooRealVar* B1 = new RooRealVar("B1","B1", randomizeItialGuess(0.0, 10.0),     0.0, 10.0); B1->setError(0.00001);
	RooRealVar* B2 = new RooRealVar("B2","B2", randomizeItialGuess(-10.0, 10.0), -10.0, 10.0); B2->setError(0.00001);
	RooRealVar* B3 = new RooRealVar("B3","B3", randomizeItialGuess(0.0, 10.0),     0.0, 10.0); B3->setError(0.00001);
	RooRealVar* B4 = new RooRealVar("B4","B4", randomizeItialGuess(0.0, 10.0),     0.0, 10.0); B4->setError(0.00001);
	TString bkgBDTpdfFormula3 = "B0+B1*TMath::Power((score+1),B2*TMath::Log(score+1)+B3*(score+1))+B4/(score+1)";
	RooGenericPdf* bkgBDTpdf3 = new RooGenericPdf("bkgBDTpdf3","bkgBDTpdf3",bkgBDTpdfFormula3,RooArgSet(*score,*B0,*B1,*B2,*B3,*B4));
	RooFitResult* fitresult_score3 = bkgBDTpdf3->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	//// nominal
	RooRealVar* C0 = new RooRealVar("C0","C0",randomizeItialGuess(0.,+1000.),0.,+1000.); C0->setError(0.000001);
	RooRealVar* C1 = new RooRealVar("C1","C1",randomizeItialGuess(0.,1000),0.,1000.);    C1->setError(0.000001);
	RooRealVar* C2 = new RooRealVar("C2","C2",randomizeItialGuess(-1.5,0.),-1.5,0.);     C2->setError(0.000001);
	RooRealVar* C3 = new RooRealVar("C3","C3",randomizeItialGuess(0.,1000),0.,1000.);    C3->setError(0.000001);
	RooRealVar* C4 = new RooRealVar("C4","C4",randomizeItialGuess(1.,+10.),1.,+10.);     C4->setError(0.000001);
	TString bkgBDTpdfFormula1 = "C0+C1*pow((score+1),C2)+C3*pow((score+1),C4)";
	RooGenericPdf* bkgBDTpdf1 = new RooGenericPdf("bkgBDTpdf1","bkgBDTpdf1",bkgBDTpdfFormula1,RooArgSet(*score,*C0,*C1,*C2,*C3,*C4));
	RooFitResult* fitresult_score1 = bkgBDTpdf1->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	//// alternative
	RooRealVar* D0 = new RooRealVar("D0","D0",randomizeItialGuess(0.,100.), 0.,100.);   D0->setError(0.000001);
	RooRealVar* D1 = new RooRealVar("D1","D1",randomizeItialGuess(1.,100.), 1.,100.);   D1->setError(0.000001);
	RooRealVar* D2 = new RooRealVar("D2","D2",randomizeItialGuess(0.,10.),  0.,10.);    D2->setError(0.000001);
	RooRealVar* D3 = new RooRealVar("D3","D3",randomizeItialGuess(0.,10.),  0.,10.);    D3->setError(0.000001);
	RooGenericPdf* bkgBDTpdf2 = new RooGenericPdf("bkgBDTpdf2","bkgBDTpdf2","D0+D1*pow(D2,score)/(score+1)+D3*pow(score,3)",RooArgSet(*score,*D0,*D1,*D2,*D3));
	RooFitResult* fitresult_score2 = bkgBDTpdf2->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	RooRealVar* E0 = new RooRealVar("E0","E0",randomizeItialGuess(0.,100.), 0.,100.);   E0->setError(0.000001);
	RooRealVar* E1 = new RooRealVar("E1","E1",randomizeItialGuess(0.,100.), 0.,100.);   E1->setError(0.000001);
	RooRealVar* E2 = new RooRealVar("E2","E2",randomizeItialGuess(-10.,0.), -10.,0.);   E2->setError(0.000001);
	RooRealVar* E3 = new RooRealVar("E3","E3",randomizeItialGuess(0.,100.), 0.,100.);   E3->setError(0.000001);
	RooRealVar* E4 = new RooRealVar("E4","E4",randomizeItialGuess(0.,+10.), 0.,+10.);   E4->setError(0.000001);
	TString bkgBDTpdfFormula4 = "E0+E1*TMath::Exp(E2*score)+E3*TMath::Exp(E4*score)";
	RooGenericPdf* bkgBDTpdf4 = new RooGenericPdf("bkgBDTpdf4","bkgBDTpdf4",bkgBDTpdfFormula4,RooArgSet(*score,*E0,*E1,*E2,*E3,*E4));
	RooFitResult* fitresult_score4 = bkgBDTpdf4->fitTo( *UnbinnedDataSet_score,Minos(kTRUE),Range("range_score"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	
	TMinuit* gFit_score0 = gMinuit;
	bool FitStatus_score0 = minuitStatus(gFit_score0);
	fitresult_score0->Print("v");
	
	TMinuit* gFit_score1 = gMinuit;
	bool FitStatus_score1 = minuitStatus(gFit_score1);
	fitresult_score1->Print("v");
	
	TMinuit* gFit_score2 = gMinuit;
	bool FitStatus_score2 = minuitStatus(gFit_score2);
	fitresult_score2->Print("v");
	
	TMinuit* gFit_score3 = gMinuit;
	bool FitStatus_score3 = minuitStatus(gFit_score3);
	fitresult_score3->Print("v");
	
	TMinuit* gFit_score4 = gMinuit;
	bool FitStatus_score4 = minuitStatus(gFit_score4);
	fitresult_score4->Print("v");
	
	
	RooAbsPdf* bkgBDTpdfNominal0 = bkgBDTpdf1; TString bkgBDTpdfNominalFormula0 = bkgBDTpdfFormula1; const unsigned int nparBDT0 = 5;
	RooAbsPdf* bkgBDTpdfNominal1 = bkgBDTpdf4; TString bkgBDTpdfNominalFormula1 = bkgBDTpdfFormula4; const unsigned int nparBDT1 = 5;
	RooAbsPdf* bkgBDTpdfNominal2 = bkgBDTpdf0; TString bkgBDTpdfNominalFormula2 = bkgBDTpdfFormula0; const unsigned int nparBDT2 = 5;
	RooAbsPdf* bkgBDTpdfNominal3 = bkgBDTpdf3; TString bkgBDTpdfNominalFormula3 = bkgBDTpdfFormula3; const unsigned int nparBDT3 = 5;
	
	
	
	
	
	
	
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* scoreFrame = score->frame(Name("scoreFrame"),Title("Background BDT response in the sidebands: "+sidebands+" and BDT>"+tstr(xbdtcutoff,2)));
	scoreFrame->SetMaximum(20);
	// scoreFrame->SetAxisRange(0,50,"Y");
	UnbinnedDataSet_score->plotOn(scoreFrame,Name("BDT score"),MarkerSize(1),Binning(nbdtbins));
	// UnbinnedDataSet_score->plotOn(scoreFrame,Name("BDT score"),XErrorSize(0),MarkerSize(1),Binning(nbdtbins));
	// BinnedDataSet_score->plotOn(scoreFrame,Name("BDT score"),XErrorSize(0),MarkerSize(0.3),Binning(nbdtbins));
	
	bkgBDTpdfNominal0->plotOn(scoreFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),NormRange("range_score"));
	bkgBDTpdfNominal1->plotOn(scoreFrame,Name("Alternative A"),LineWidth(2),LineColor(kRed+1),LineStyle(kDashed),NormRange("range_score"));
	bkgBDTpdfNominal2->plotOn(scoreFrame,Name("Alternative B"),LineWidth(2),LineColor(kGreen+1),LineStyle(kDashed),NormRange("range_score"));
	bkgBDTpdfNominal3->plotOn(scoreFrame,Name("Alternative C"),LineWidth(2),LineColor(kGray),LineStyle(kDashed),NormRange("range_score"));
	// bkgBDTpdfNominal4->plotOn(scoreFrame,Name("Alternative D"),LineWidth(2),LineColor(kMagenta-4),LineStyle(kDashed),NormRange("range_score"));
	// bkgBDTpdfNominal0->paramOn(scoreFrame,Layout(0.5,0.88,0.87), Format("NEU", AutoPrecision(3))); // X size of box is from 50% to 88% of Xaxis range, top of box is at 87% of Yaxis range)
	// scoreFrame->getAttText()->SetTextSize(0.03);
	
	leg->Clear();
	leg->AddEntry(scoreFrame->findObject("BDT score"),"Sidebands data","ple");
	leg->AddEntry(scoreFrame->findObject("Nominal"),"Nominal fit","l");
	leg->AddEntry(scoreFrame->findObject("Alternative A"),"Alternative A","l");
	leg->AddEntry(scoreFrame->findObject("Alternative B"),"Alternative B","l");
	leg->AddEntry(scoreFrame->findObject("Alternative C"),"Alternative C","l");
	// leg->AddEntry(scoreFrame->findObject("Alternative D"),"Alternative D","l");
	
	// cnv->SetLeftMargin(0.2);
	scoreFrame->SetTitleOffset(2,"Y");
	scoreFrame->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimate.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimate.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf(");
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* scoreFrameZoom = score->frame(Name("scoreFrameZoom"),Title("Background BDT response in the sidebands: "+sidebands+" and BDT>"+tstr(xbdtcutoff,2)),Range("range_score_zoom"));
	UnbinnedDataSet_score->plotOn(scoreFrameZoom,Name("BDT score zoom"),MarkerSize(1),Binning(nbdtbins),Range("range_score_zoom"));
	bkgBDTpdfNominal0->plotOn(scoreFrameZoom,Name("Nominal"),LineWidth(2),LineColor(kBlue),Range("range_score"));
	bkgBDTpdfNominal1->plotOn(scoreFrameZoom,Name("Alternative A"),LineWidth(2),LineColor(kRed+1),LineStyle(kDashed));
	bkgBDTpdfNominal2->plotOn(scoreFrameZoom,Name("Alternative B"),LineWidth(2),LineColor(kGreen+1),LineStyle(kDashed));
	bkgBDTpdfNominal3->plotOn(scoreFrameZoom,Name("Alternative C"),LineWidth(2),LineColor(kGray),LineStyle(kDashed));
	// bkgBDTpdfNominal0->paramOn(Zoom,Layout(0.5,0.88,0.87), Format("NEU", AutoPrecision(3))); // X size of box is from 50% to 88% of Xaxis range, top of box is at 87% of Yaxis range)
	// Zoom->getAttText()->SetTextSize(0.03);
	
	leg->Clear();
	leg->AddEntry(scoreFrameZoom->findObject("BDT score zoom"),"Sidebands data","pl");
	leg->AddEntry(scoreFrameZoom->findObject("Nominal"),"Nominal fit","l");
	leg->AddEntry(scoreFrameZoom->findObject("Alternative A"),"Alternative A","l");
	leg->AddEntry(scoreFrameZoom->findObject("Alternative B"),"Alternative B","l");
	leg->AddEntry(scoreFrameZoom->findObject("Alternative C"),"Alternative C","l");
	
	// // cnv->SetLeftMargin(0.2);
	scoreFrameZoom->SetMaximum(10);
	scoreFrameZoom->SetTitleOffset(2,"Y");
	scoreFrameZoom->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimate.BDTzoom."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimate.BDTzoom."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	
	
	
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* scoreFrameS = score->frame(Name("scoreFrameS"),Title("BDT response"));
	UnbinnedDataSet_score->plotOn(scoreFrameS,Name("BDT score background"),MarkerSize(1),Binning(nbdtbins));
	// UnbinnedDataSet_score_noblind->plotOn(scoreFrameS,Name("BDT score background noblind"),MarkerSize(1),MarkerColor(kRed),Binning(nbdtbins));
	// UnbinnedDataSet_scoreS->plotOn(scoreFrameS,Name("BDT score signal"),LineStyle(kDashed),LineColor(kRed),MarkerStyle(22),MarkerColor(kRed),Binning(nbdtbins));
	hScoreSig->SetLineColor(kGray); hScoreSig->SetFillColor(kGray);
	hScoreBkg_noblind->SetLineColor(kRed); hScoreBkg_noblind->SetLineStyle(2);
	bkgBDTpdfNominal0->plotOn(scoreFrameS,Name("Nominal"),LineWidth(2),LineColor(kBlue),NormRange("range_score"));
	
	// delete leg;
	leg = new TLegend(0.15,0.74,0.55,0.92,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(scoreFrameS->findObject("BDT score background"),"Sidebands data (tight+x>x_{0})","ple");
	leg->AddEntry(scoreFrameS->findObject("Nominal"),"Fit to the data","l");
	leg->AddEntry(hScoreSig,"Signal (tight+x>x_{0})","f");
	
	// // cnv->SetLeftMargin(0.2);
	cnv->SetTicks(1,1);
	hScoreSig->Draw("hist");
	// hScoreBkg_noblind->Draw("hist same");
	scoreFrameS->Draw("same");
	plotAtlasLabel(); // ptxt->Draw("same");
	leg->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/BkgEstimateWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateWithSignal.BDT."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	
	

	Double_t chi2BDTpdf0 = scoreFrame->chiSquare("Nominal","BDT score",      nparBDT0);
	Double_t chi2BDTpdf1 = scoreFrame->chiSquare("Alternative A","BDT score",nparBDT1);
	Double_t chi2BDTpdf2 = scoreFrame->chiSquare("Alternative B","BDT score",nparBDT2);
	Double_t chi2BDTpdf3 = scoreFrame->chiSquare("Alternative C","BDT score",nparBDT3);
	cout << "chi2BDTpdf0=" << chi2BDTpdf0 << endl;
	cout << "chi2BDTpdf1=" << chi2BDTpdf1 << endl;
	cout << "chi2BDTpdf2=" << chi2BDTpdf2 << endl;
	cout << "chi2BDTpdf3=" << chi2BDTpdf3 << endl;
	
	
	_INF(1,"------------------- Done fitting the BDT shape -------------------");
	
	
	
	
	
	
	
	//////////////////
	//// Fit the 3body mass for the signal
	RooDataSet* rds_m3bodyS = (RooDataSet*)UnbinnedDataSet_m3bodyS;
	BinnedDataSet_m3bodyS = rds_m3bodyS->binnedClone();
	
	m3bodyS->setRange("range_m3bodyS",m3bodyMin,m3bodyMax);
	m3bodyS->setBins(nm3bodybins);
	
	// Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
	// RooRealVar* mean   = new RooRealVar("mean","mean of gaussians",1777) ;
	// RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussians",30) ;
	// RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussians",30) ;
	// RooGaussian* sig1 = new RooGaussian("sig1","Signal component 1",*m3bodyS,*mean,*sigma1) ;  
	// RooGaussian* sig2 = new RooGaussian("sig2","Signal component 2",*m3bodyS,*mean,*sigma2) ;
	// // Sum the signal components into a composite signal p.d.f.
	// RooRealVar* sig1frac = new RooRealVar("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
	// RooAddPdf* sigm3bodypdf = new RooAddPdf("sig","Signal",RooArgList(*sig1,*sig2),*sig1frac) ;
	
	// RooRealVar* mean   = new RooRealVar("mean","mean of gaussians",1777);
	// RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussian1",30);
	// RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussian2",30);
	// RooGaussian* sig1 = new RooGaussian("sig1","Signal component 1",*m3bodyS,*mean,*sigma1);
	// RooGaussian* sig2 = new RooGaussian("sig2","Signal component 2",*m3bodyS,*mean,*sigma2);
	// RooAddPdf* sigm3bodypdf = new RooAddPdf("sigm3bodypdf","sigm3bodypdf",RooArgList(*sig1,*sig2)) ;	
	
	RooRealVar* mean   = new RooRealVar("Mean","Mean",1777,1750,1800);
	RooRealVar* sigma = new RooRealVar("#sigma","#sigma",30,0,100);
	RooGaussian* sigm3bodypdf = new RooGaussian("sigm3bodypdf","sigm3bodypdf",*m3bodyS,*mean,*sigma);
	
	// // Voigtian is convolution of Gaussian and Breit-Wigner (with numerical normalization integral)
	// RooRealVar* mass = new RooRealVar("M","M",1776.82);
	// RooRealVar* sigma = new RooRealVar("#sigma","#sigma",30,0.,100);
	// RooRealVar* width = new RooRealVar("#Gamma","#Gamma",0,0,100);
	// RooVoigtian* sigm3bodypdf = new RooVoigtian("sigm3bodypdf","sigm3bodypdf",*m3bodyS,*mass,*width,*sigma);
	
	RooFitResult* fitresult_m3bodyS;
	fitresult_m3bodyS = sigm3bodypdf->fitTo( *UnbinnedDataSet_m3bodyS,Minos(kTRUE),Range("range_m3bodyS"),NormRange("range_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	// fitresult_m3bodyS = sigm3bodypdf->fitTo( *UnbinnedDataSet_m3bodyS,Minos(kTRUE),Range("rrange_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	// fitresult_m3bodyS = sigm3bodypdf->fitTo( *BinnedDataSet_m3bodyS,Minos(kTRUE),Range("range_m3bodyS"),NormRange("range_m3bodyS"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	
	TMinuit* gFit_m3bodyS = gMinuit;
	bool FitStatus_m3bodyS = minuitStatus(gFit_m3bodyS);
	fitresult_m3bodyS->Print("v");
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* m3bodySFrame = m3bodyS->frame(Name("m3bodySFrame"),Title("Signal 3body mass"));
	UnbinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("m3bodyS"),MarkerSize(1),Binning(nm3bodybins));
	// UnbinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("BDT m3bodyS"),XErrorSize(0),MarkerSize(1),Binning(nm3bodybins));
	// BinnedDataSet_m3bodyS->plotOn(m3bodySFrame,Name("BDT m3bodyS"),XErrorSize(0),MarkerSize(0.3),Binning(nm3bodybins));
	
	// sigm3bodypdf->plotOn(m3bodySFrame,LineWidth(2),LineColor(kBlue),Range("range_SBleft,range_SBright"),NormRange("range_SBleft,range_SBright"));
	sigm3bodypdf->plotOn(m3bodySFrame,Name("Signal"),LineWidth(2),LineColor(kBlue),Range("range_m3bodyS"),NormRange("range_m3bodyS"));
	// sigm3bodypdf->plotOn(m3bodySFrame,LineWidth(2),LineColor(kBlue),Range("range_m3bodyS"));
	sigm3bodypdf->paramOn(m3bodySFrame,Layout(0.55,0.88,0.30), Format("NEU", AutoPrecision(2))); // X size of box is from 60% to 88% of Xaxis range, top of box is at 87% of Yaxis range)
	// m3bodySFrame->getAttText()->SetTextSize(0.03);
	
	Float_t xSRmin = 1713; //mean->getVal()-2.*sigma->getVal(); // 1717.97;
	Float_t xSRmax = 1841; //mean->getVal()+2.*sigma->getVal(); // 1838.06;
	m3body->setRange("range_m3body_SR",xSRmin,xSRmax);
	
	// // cnv->SetLeftMargin(0.2);
	m3bodySFrame->SetTitleOffset(2,"Y");
	m3bodySFrame->Draw();
	ptxt->Draw("same");
	cnv->SaveAs("figures/BkgEstimate.Gaus."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimate.Gaus."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");

	
	Double_t chi2Signalpdf = m3bodySFrame->chiSquare("Signal","m3bodyS",2);
	cout << "chi2Signalpdf=" << chi2Signalpdf << endl;
	
	_INF(1,"------------------- Done fitting the signal shape -------------------");
	





	//return;








	
	
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
	

	
	RooRealVar* a0 = new RooRealVar("a0","a0", randomizeItialGuess(-1.0,1.0), -1.0, 1.0); a0->setError(0.00001);
	RooChebychev* bkgm3bodypdf0 = new RooChebychev("bkgm3bodypdf0","bkgm3bodypdf0",*m3body,RooArgSet(*a0));
	
	RooRealVar* b0 = new RooRealVar("b0","b0", randomizeItialGuess(-1.0,0.0)/*-1.42803e-01*/, -1.0, 1.0); b0->setError(0.00001);
	RooRealVar* b1 = new RooRealVar("b1","b1", randomizeItialGuess( 0.0,0.1)/* 2.87161e-02*/, -1.0, 1.0); b1->setError(0.00001);
	RooChebychev* bkgm3bodypdf1 = new RooChebychev("bkgm3bodypdf1","bkgm3bodypdf1",*m3body,RooArgSet(*b0,*b1));
	
	RooRealVar* c0 = new RooRealVar("c0","c0", randomizeItialGuess(-1.00,0.0)/*-1.42803e-01*/, -1.0, 1.0); c0->setError(0.00001);
	RooRealVar* c1 = new RooRealVar("c1","c1", randomizeItialGuess( 0.00,0.1)/* 2.87161e-02*/, -1.0, 1.0); c1->setError(0.00001);
	RooRealVar* c2 = new RooRealVar("c2","c2", randomizeItialGuess(-0.01,0.0)/*-1.45196e-03*/, -1.0, 1.0); c2->setError(0.00001);
	RooChebychev* bkgm3bodypdf2 = new RooChebychev("bkgm3bodypdf2","bkgm3bodypdf2",*m3body,RooArgSet(*c0,*c1,*c2));
	
	RooRealVar* d0 = new RooRealVar("d0","d0", randomizeItialGuess(-1.00,0.0)/*-1.42803e-01*/, -1.0, 1.0); d0->setError(0.00001);
	RooRealVar* d1 = new RooRealVar("d1","d1", randomizeItialGuess( 0.00,0.1)/* 2.87161e-02*/, -1.0, 1.0); d1->setError(0.00001);
	RooRealVar* d2 = new RooRealVar("d2","d2", randomizeItialGuess(-0.01,0.0)/*-1.45196e-03*/, -1.0, 1.0); d2->setError(0.00001);
	RooRealVar* d3 = new RooRealVar("d3","d3", randomizeItialGuess(-0.01,0.0)/*-6.00993e-04*/, -1.0, 1.0); d3->setError(0.00001);
	RooChebychev* bkgm3bodypdf3 = new RooChebychev("bkgm3bodypdf3","bkgm3bodypdf3",*m3body,RooArgSet(*d0,*d1,*d2,*d3));
	
	RooRealVar* e0 = new RooRealVar("e0","e0", randomizeItialGuess(-1.000,0.0)/*-1.42803e-01*/, -1.0, 1.0); e0->setError(0.00001);
	RooRealVar* e1 = new RooRealVar("e1","e1", randomizeItialGuess( 0.000,0.1)/* 2.87161e-02*/, -1.0, 1.0); e1->setError(0.00001);
	RooRealVar* e2 = new RooRealVar("e2","e2", randomizeItialGuess(-0.010,0.0)/*-1.45196e-03*/, -1.0, 1.0); e2->setError(0.00001);
	RooRealVar* e3 = new RooRealVar("e3","e3", randomizeItialGuess(-0.010,0.0)/*-6.00993e-04*/, -1.0, 1.0); e3->setError(0.00001);
	RooRealVar* e4 = new RooRealVar("e4","e4", randomizeItialGuess(-0.001,0.0)/*-6.00993e-04*/, -1.0, 1.0); e4->setError(0.00001);
	RooChebychev* bkgm3bodypdf4 = new RooChebychev("bkgm3bodypdf4","bkgm3bodypdf4",*m3body,RooArgSet(*e0,*e1,*e2,*e3,*e4));
	
	RooRealVar* f0 = new RooRealVar("f0","f0", randomizeItialGuess(-1.0000,0.0)/*-1.42803e-01*/, -1.0, 1.0); f0->setError(0.00001);
	RooRealVar* f1 = new RooRealVar("f1","f1", randomizeItialGuess( 0.0000,0.1)/* 2.87161e-02*/, -1.0, 1.0); f1->setError(0.00001);
	RooRealVar* f2 = new RooRealVar("f2","f2", randomizeItialGuess(-0.0100,0.0)/*-1.45196e-03*/, -1.0, 1.0); f2->setError(0.00001);
	RooRealVar* f3 = new RooRealVar("f3","f3", randomizeItialGuess(-0.0100,0.0)/*-6.00993e-04*/, -1.0, 1.0); f3->setError(0.00001);
	RooRealVar* f4 = new RooRealVar("f4","f4", randomizeItialGuess(-0.0010,0.0)/*-6.00993e-04*/, -1.0, 1.0); f4->setError(0.00001);
	RooRealVar* f5 = new RooRealVar("f5","f5", randomizeItialGuess(-0.0001,0.0)/*-6.00993e-04*/, -1.0, 1.0); f5->setError(0.00001);
	RooChebychev* bkgm3bodypdf5 = new RooChebychev("bkgm3bodypdf5","bkgm3bodypdf5",*m3body,RooArgSet(*f0,*f1,*f2,*f3,*f4,*f5));
	
	RooRealVar* g0 = new RooRealVar("g0","g0",1.,0.,100.); g0->setError(0.00001);
	RooRealVar* g1 = new RooRealVar("g1","g1",1.,0.,100.); g1->setError(0.00001);
	RooRealVar* g2 = new RooRealVar("g2","g2",1.,0.,100.); g2->setError(0.00001);
	RooRealVar* g3 = new RooRealVar("g3","g3",1.,0.,100.); g3->setError(0.00001);
	RooRealVar* g4 = new RooRealVar("g4","g4",1.,0.,100.); g4->setError(0.00001);
	RooRealVar* g5 = new RooRealVar("g5","g5",1.,0.,100.); g5->setError(0.00001);
	RooRealVar* g6 = new RooRealVar("g6","g6",1.,0.,100.); g6->setError(0.00001);
	RooBernstein* bkgm3bodypdf6 = new RooBernstein("bkgm3bodypdf6","bkgm3bodypdf6",*m3body,RooArgSet(*g0,*g1,*g2,*g3));
	
	RooRealVar* h0 = new RooRealVar("h0","h0", randomizeItialGuess(0.0,100.0), 0., 100.0); h0->setError(0.00001);
	RooRealVar* h1 = new RooRealVar("h1","h1", randomizeItialGuess(0.0,100.0), 0., 100.0); h1->setError(0.00001);
	RooGenericPdf* bkgm3bodypdf8 = new RooGenericPdf("bkgm3bodypdf8","bkgm3bodypdf8","m3body/pow((h0+m3body),h1)",RooArgSet(*m3body,*h0,*h1));

	RooRealVar* i0 = new RooRealVar("i0","i0", randomizeItialGuess(0.0,100.0), 0., 100.0);   i0->setError(0.00001);
	RooRealVar* i1 = new RooRealVar("i1","i1", randomizeItialGuess(0.0,100.0), 0., 100.0);   i1->setError(0.00001);
	RooRealVar* i2 = new RooRealVar("i2","i2", randomizeItialGuess(-1000.0,1000.0), -1000., 1000.0); i2->setError(0.00001);
	RooGenericPdf* bkgm3bodypdf9 = new RooGenericPdf("bkgm3bodypdf9","bkgm3bodypdf9","i0+i1/(m3body+i2)",RooArgSet(*m3body,*i0,*i1,*i2));

	RooRealVar* j0 = new RooRealVar("j0","j0", randomizeItialGuess(0.0,100.0), 0., 100.0);     j0->setError(0.00001);
	RooRealVar* j1 = new RooRealVar("j1","j1", randomizeItialGuess(0.0,100.0), 0., 100.0);     j1->setError(0.00001);
	RooRealVar* j2 = new RooRealVar("j2","j2", randomizeItialGuess(-1000.0,0.0), -1000., 0.0); j2->setError(0.00001);
	RooGenericPdf* bkgm3bodypdf10 = new RooGenericPdf("bkgm3bodypdf10","bkgm3bodypdf10","j0+j1*TMath::Exp(m3body*j2)",RooArgSet(*m3body,*j0,*j1,*j2));
	
	RooRealVar* k0 = new RooRealVar("k0","k0",1.,0.,100.); k0->setError(0.00001);
	RooRealVar* k1 = new RooRealVar("k1","k1",1.,0.,100.); k1->setError(0.00001);
	// RooRealVar* k2 = new RooRealVar("k2","k2",1.,0.,100.); k2->setError(0.00001);
	RooBernstein* bkgm3bodypdf11 = new RooBernstein("bkgm3bodypdf11","bkgm3bodypdf11",*m3body,RooArgSet(*k0,*k1));
	
	RooRealVar* l0 = new RooRealVar("l0","l0",1.,0.,100.); l0->setError(0.00001);
	RooRealVar* l1 = new RooRealVar("l1","l1",1.,0.,100.); l1->setError(0.00001);
	RooRealVar* l2 = new RooRealVar("l2","l2",1.,0.,100.); l2->setError(0.00001);
	RooRealVar* l3 = new RooRealVar("l3","l3",1.,0.,100.); l3->setError(0.00001);
	RooBernstein* bkgm3bodypdf12 = new RooBernstein("bkgm3bodypdf12","bkgm3bodypdf12",*m3body,RooArgSet(*l0,*l1,*l2,*l3));
	
	
	RooFitResult* fitresult_m3body0;
	fitresult_m3body0 = bkgm3bodypdf0->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body0 = gMinuit;
	bool FitStatus_m3body0 = minuitStatus(gFit_m3body0);
	fitresult_m3body0->Print("v");
	
	RooFitResult* fitresult_m3body1;
	fitresult_m3body1 = bkgm3bodypdf1->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body1 = gMinuit;
	bool FitStatus_m3body1 = minuitStatus(gFit_m3body1);
	fitresult_m3body1->Print("v");
	
	RooFitResult* fitresult_m3body2;
	fitresult_m3body2 = bkgm3bodypdf2->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body2 = gMinuit;
	bool FitStatus_m3body2 = minuitStatus(gFit_m3body2);
	fitresult_m3body2->Print("v");
	
	RooFitResult* fitresult_m3body3;
	fitresult_m3body3 = bkgm3bodypdf3->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body3 = gMinuit;
	bool FitStatus_m3body3 = minuitStatus(gFit_m3body3);
	fitresult_m3body3->Print("v");
	
	RooFitResult* fitresult_m3body4;
	fitresult_m3body4 = bkgm3bodypdf4->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body4 = gMinuit;
	bool FitStatus_m3body4 = minuitStatus(gFit_m3body4);
	fitresult_m3body4->Print("v");
	
	RooFitResult* fitresult_m3body5;
	fitresult_m3body5 = bkgm3bodypdf5->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body5 = gMinuit;
	bool FitStatus_m3body5 = minuitStatus(gFit_m3body5);
	fitresult_m3body5->Print("v");
	
	RooFitResult* fitresult_m3body6;
	fitresult_m3body6 = bkgm3bodypdf6->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body6 = gMinuit;
	bool FitStatus_m3body6 = minuitStatus(gFit_m3body6);
	fitresult_m3body6->Print("v");

	RooFitResult* fitresult_m3body8;
	fitresult_m3body8 = bkgm3bodypdf8->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body8 = gMinuit;
	bool FitStatus_m3body8 = minuitStatus(gFit_m3body8);
	fitresult_m3body8->Print("v");
	
	RooFitResult* fitresult_m3body9;
	fitresult_m3body9 = bkgm3bodypdf9->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body9 = gMinuit;
	bool FitStatus_m3body9 = minuitStatus(gFit_m3body9);
	fitresult_m3body9->Print("v");
	
	RooFitResult* fitresult_m3body10;
	fitresult_m3body10 = bkgm3bodypdf10->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body10 = gMinuit;
	bool FitStatus_m3body10 = minuitStatus(gFit_m3body10);
	fitresult_m3body10->Print("v");

	RooFitResult* fitresult_m3body11;
	fitresult_m3body11 = bkgm3bodypdf11->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body11 = gMinuit;
	bool FitStatus_m3body11 = minuitStatus(gFit_m3body11);
	fitresult_m3body11->Print("v");
	
	RooFitResult* fitresult_m3body12;
	fitresult_m3body12 = bkgm3bodypdf12->fitTo( *UnbinnedDataSet_m3body,Minos(kTRUE),Range("range_SBleft,range_SBright"),Strategy(2),Save(kTRUE),Timer(kTRUE));
	TMinuit* gFit_m3body12 = gMinuit;
	bool FitStatus_m3body12 = minuitStatus(gFit_m3body12);
	fitresult_m3body12->Print("v");
	
	
	RooAbsPdf* bkgm3bodypdfNominal0 = bkgm3bodypdf11; Int_t nparSB0 = 2;
	RooAbsPdf* bkgm3bodypdfNominal1 = bkgm3bodypdf0;  Int_t nparSB1 = 1;
	RooAbsPdf* bkgm3bodypdfNominal2 = bkgm3bodypdf9;  Int_t nparSB2 = 3;
	RooAbsPdf* bkgm3bodypdfNominal3 = bkgm3bodypdf8;  Int_t nparSB3 = 2;
	RooAbsPdf* bkgm3bodypdfNominal4 = bkgm3bodypdf12; Int_t nparSB4 = 4;
	
	








	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	
	RooPlot* m3bodyFrame = m3body->frame(Name("m3bodyFrame"),Title("3body mass in the sidebands of CR_{0} (BDT>"+tstr(xbdtcutoff,2)+")"));
	m3bodyFrame->SetMaximum(20);
	if(!blinded) UnbinnedDataSet_m3body_noblind->plotOn(m3bodyFrame,Name("m3body SB0"),MarkerSize(1),Binning(nm3bodybins));
	else         UnbinnedDataSet_m3body->plotOn(m3bodyFrame,Name("m3body SB0"),MarkerSize(1),Binning(nm3bodybins));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrame,Name("Nominal"),LineWidth(2),LineColor(kBlue),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal1->plotOn(m3bodyFrame,Name("Alternative A"),LineWidth(2),LineColor(kRed),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal2->plotOn(m3bodyFrame,Name("Alternative B"),LineWidth(2),LineColor(kViolet),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal3->plotOn(m3bodyFrame,Name("Alternative C"),LineWidth(2),LineColor(kGreen+1),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	// bkgm3bodypdfNominal4->plotOn(m3bodyFrame,Name("Alternative D"),LineWidth(2),LineColor(kOrange+1),LineStyle(kDashed),Range("range_m3body"),NormRange("range_SBleft,range_SBright"));
	
	// bkgm3bodypdfNominal0->paramOn(m3bodyFrame,Layout(0.25,0.5,0.4), Format("NEU", AutoPrecision(3))); // X size of box is from 25% to 50% of Xaxis range, top of box is at 40% of Yaxis range)
	// m3bodyFrame->getAttText()->SetTextSize(0.018);
	
	leg->Clear();
	leg->AddEntry(m3bodyFrame->findObject("m3body SB0"),"Data in SB_{0}","pl");
	leg->AddEntry(m3bodyFrame->findObject("Nominal"),"Nominal fit","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative A"),"Alternative A","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative B"),"Alternative B","l");
	leg->AddEntry(m3bodyFrame->findObject("Alternative C"),"Alternative C","l");
	// leg->AddEntry(m3bodyFrame->findObject("Alternative D"),"Alternative D","l");
	
	// cnv->SetLeftMargin(0.2);
	// m3bodyFrame->SetMaximum(1.3*m3bodyFrame->GetMaximum());
	m3bodyFrame->SetTitleOffset(2,"Y");
	m3bodyFrame->Draw();
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->SaveAs("figures/BkgEstimate.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimate.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	
	
	

	
	///////////
	if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	RooPlot* m3bodyFrameS = m3body->frame(Name("m3bodyFrameS"),Title("3body mass in the sidebands (BDT>"+tstr(xbdtcutoff,2)+")"));
	if(blinded) UnbinnedDataSet_m3body->plotOn(m3bodyFrameS,Name("m3body background"),MarkerSize(1),Binning(nm3bodybins));
	else        UnbinnedDataSet_m3body_noblind->plotOn(m3bodyFrameS,Name("m3body background"),MarkerSize(1),Binning(nm3bodybins));
	hM3bodySig->SetLineColor(kGray); hM3bodySig->SetFillColor(kGray);
	hM3bodyBkg_noblind_opt->SetLineColor(kRed); hM3bodyBkg_noblind_opt->SetMarkerStyle(20); hM3bodyBkg_noblind_opt->SetMarkerColor(kRed);
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_left"),LineWidth(2),LineColor(kBlue),Range("range_SBleft"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_right"),LineWidth(2),LineColor(kBlue),Range("range_SBright"),NormRange("range_SBleft,range_SBright"));
	bkgm3bodypdfNominal0->plotOn(m3bodyFrameS,Name("Nominal_blinded"),LineWidth(2),LineStyle(2),LineColor(kBlue),Range("range_blinded"),NormRange("range_SBleft,range_SBright"));
	if(!blinded) UnbinnedDataSet_m3body_noblind_opt->plotOn(m3bodyFrameS,Name("m3body background noblind opt"),MarkerSize(1),MarkerColor(kRed),LineColor(kRed),MarkerStyle(22),Binning(nm3bodybins));
	m3bodyFrameS->SetMinimum(1e-5);
	
	
	
	delete leg;
	leg = new TLegend(0.15,0.70,0.55,0.92,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(m3bodyFrameS->findObject("m3body background"),"Sidebands data (tight+x>x_{0})","ple");
	if(!blinded) leg->AddEntry(m3bodyFrameS->findObject("m3body background noblind opt"),"Sidebands data (tight+x>x_{1})","ple");
	leg->AddEntry(m3bodyFrameS->findObject("Nominal_left"),"Sidebands fit","l");
	leg->AddEntry(m3bodyFrameS->findObject("Nominal_blinded"),"Interpolation","l");
	leg->AddEntry(hM3bodySig,"Signal (tight+x>x_{0})","f");
	
	// cnv->SetLeftMargin(0.2);
	cnv->SetTicks(1,1);
	cnv->Update();
	hM3bodySig->Draw("hist");
	m3bodyFrameS->Draw("same");
	plotAtlasLabel(); // ptxt->Draw("same");
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
	
	y = ymax/1.6;
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
	legR->AddEntry(lSBright,"Sidebands","l");
	legR->AddEntry(lSRright,"Signal region","l");
	legR->Draw("same");
	
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("figures/BkgEstimateWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".png");
	cnv->SaveAs("figures/BkgEstimateWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".eps");
	cnv->SaveAs("figures/BkgEstimateWithSignal.SB0."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf)");
	///////////
	
	
	
	
	Double_t chi2SBpdf0 = m3bodyFrame->chiSquare("Nominal",      "m3body SB0",nparSB0); // last parameter is Ndof for the chi^2
	Double_t chi2SBpdf1 = m3bodyFrame->chiSquare("Alternative A","m3body SB0",nparSB1); // last parameter is Ndof for the chi^2
	Double_t chi2SBpdf2 = m3bodyFrame->chiSquare("Alternative B","m3body SB0",nparSB2); // last parameter is Ndof for the chi^2
	Double_t chi2SBpdf3 = m3bodyFrame->chiSquare("Alternative C","m3body SB0",nparSB3); // last parameter is Ndof for the chi^2
	cout << "chi2SBpdf0=" << chi2SBpdf0 << endl;
	cout << "chi2SBpdf1=" << chi2SBpdf1 << endl;
	cout << "chi2SBpdf2=" << chi2SBpdf2 << endl;
	cout << "chi2SBpdf3=" << chi2SBpdf3 << endl;
	
	
	
	
	
	
	// if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	// tD->Draw("mOS1");
	// cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	// 
	// if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	// tD->Draw("mOS2");
	// cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf");
	// 
	// if(cnv) delete cnv; cnv = new TCanvas("cnv","",800,600);
	// tD->Draw("mSS");
	// cnv->SaveAs("figures/BkgEstimate."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".pdf)");
	
	
	
	
	
	
	
	

	RooAbsReal* integralSB00 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB10 = bkgBDTpdfNominal0->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float transferFactor0 = integralSB10->getVal()/integralSB00->getVal();
	float nSB00bdt = integralSB00->getVal();
	
	RooAbsReal* integralSB01 = bkgBDTpdfNominal1->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB11 = bkgBDTpdfNominal1->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float transferFactor1 = integralSB11->getVal()/integralSB01->getVal();
	float nSB01bdt = integralSB01->getVal();
	
	RooAbsReal* integralSB02 = bkgBDTpdfNominal2->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB12 = bkgBDTpdfNominal2->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float transferFactor2 = integralSB12->getVal()/integralSB02->getVal();
	float nSB02bdt = integralSB02->getVal();
	
	RooAbsReal* integralSB03 = bkgBDTpdfNominal3->createIntegral(*score,NormSet(*score),Range("range_score_SB0"));
	RooAbsReal* integralSB13 = bkgBDTpdfNominal3->createIntegral(*score,NormSet(*score),Range("range_score_SB1"));
	float transferFactor3 = integralSB13->getVal()/integralSB03->getVal();
	float nSB03bdt = integralSB03->getVal();
	

	
	
	RooAbsReal* integralSR00      = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB00left  = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB00right = bkgm3bodypdfNominal0->createIntegral(*m3body,Range("range_SBright"));
	float nSB00 = integralSB00left->getVal()+integralSB00right->getVal();
	float nSR00 = integralSR00->getVal();
	float dnSR00 = (nSR00/nSB00)*sqrt(nSB00);
	float RooFitScale0 = ncandidatesD/nSB00;
	
	RooAbsReal* integralSR01      = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB01left  = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB01right = bkgm3bodypdfNominal1->createIntegral(*m3body,Range("range_SBright"));
	float nSB01 = integralSB01left->getVal()+integralSB01right->getVal();
	float nSR01 = integralSR01->getVal();
	float dnSR01 = (nSR01/nSB01)*sqrt(nSB01);
	float RooFitScale1 = ncandidatesD/nSB01;
	
	RooAbsReal* integralSR02      = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB02left  = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB02right = bkgm3bodypdfNominal2->createIntegral(*m3body,Range("range_SBright"));
	float nSB02 = integralSB02left->getVal()+integralSB02right->getVal();
	float nSR02 = integralSR02->getVal();
	float dnSR02 = (nSR02/nSB02)*sqrt(nSB02);
	float RooFitScale2 = ncandidatesD/nSB02;
	
	RooAbsReal* integralSR03      = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_m3body_SR"));
	RooAbsReal* integralSB03left  = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_SBleft"));
	RooAbsReal* integralSB03right = bkgm3bodypdfNominal3->createIntegral(*m3body,Range("range_SBright"));
	float nSB03 = integralSB03left->getVal()+integralSB03right->getVal();
	float nSR03 = integralSR03->getVal();
	float dnSR03 = (nSR03/nSB03)*sqrt(nSB03);
	float RooFitScale3 = ncandidatesD/nSB03;
	
	
	vector<float> vxSRmin, vxSRmax;
	vector<float> vnSR00, vnSR01, vnSR02, vnSR03;
	float dxSR = 1; // MeV
	int d = 1;
	Float_t xSRminTmp = floor(tauMassMeV+0.5)-d*dxSR;
	Float_t xSRmaxTmp = floor(tauMassMeV+0.5)+d*dxSR;
	bool run = true;
	while(run)
	{
		TString tmp_label = tstr(d,0);
		TString SRnameTmp = "range_m3body_SR_"+tmp_label;
		m3body->setRange(SRnameTmp,xSRminTmp,xSRmaxTmp);
		vxSRmin.push_back(xSRminTmp); vxSRmax.push_back(xSRmaxTmp);
		
		RooAbsReal* integralSR00tmp = bkgm3bodypdfNominal0->createIntegral(*m3body,Range(SRnameTmp));
		vnSR00.push_back(integralSR00tmp->getVal());
		
		RooAbsReal* integralSR01tmp = bkgm3bodypdfNominal1->createIntegral(*m3body,Range(SRnameTmp));
		vnSR01.push_back(integralSR01tmp->getVal());
		
		RooAbsReal* integralSR02tmp = bkgm3bodypdfNominal2->createIntegral(*m3body,Range(SRnameTmp));
		vnSR02.push_back(integralSR02tmp->getVal());
		
		RooAbsReal* integralSR03tmp = bkgm3bodypdfNominal3->createIntegral(*m3body,Range(SRnameTmp));
		vnSR03.push_back(integralSR03tmp->getVal());
		
		/////////
		//// stop
		run = !(xSRminTmp==mBlindMinGlob && xSRmaxTmp==mBlindMaxGlob);
		
		///////////////
		//// propagate
		++d;
		xSRminTmp = (xSRminTmp!=mBlindMinGlob) ? floor(tauMassMeV+0.5)-d*dxSR : mBlindMinGlob;
		xSRmaxTmp = (xSRmaxTmp!=mBlindMaxGlob) ? floor(tauMassMeV+0.5)+d*dxSR : mBlindMaxGlob;
	}
	
	
	
	int ncandidatesSsr1 = readData(tS,scoreS,m3bodyS,UnbinnedDataSet_scoreS,UnbinnedDataSet_m3bodyS,
								   m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, "sig","CR1");
	float efficiencyS = ((float)ncandidatesSsr1/200000.)*100.;
	
	int ncandidatesD_CR1 = readData(tD,score,m3body,UnbinnedDataSet_score,UnbinnedDataSet_m3body,
								    m3bodyMin,xbmin,xbmax,m3bodyMax, xbdtmin,xbdtopt, "bkg","CR1");
	
	
	//// Write the TF1 object and more numbers to a root file
	TFile* fpdf = new TFile("BDTfitResults."+sidebands+"."+bdtcutoff+"."+bdtmaxusr+".root","RECREATE");
	TF1* pdf0 = (TF1*)bkgBDTpdfNominal0->asTF(*score, RooArgList(*C0,*C1,*C2,*C3,*C4), *score); RooRealVar parBDT0[nparBDT0]={*C0,*C1,*C2,*C3,*C4}; TString pdfStr0 = expandNscaleFormula(bkgBDTpdfNominalFormula0, "score", parBDT0, nparBDT0, RooFitScale0, 0);
	TF1* pdf1 = (TF1*)bkgBDTpdfNominal1->asTF(*score, RooArgList(*E0,*E1,*E2,*E3,*E4), *score); RooRealVar parBDT1[nparBDT1]={*E0,*E1,*E2,*E3,*E4}; TString pdfStr1 = expandNscaleFormula(bkgBDTpdfNominalFormula1, "score", parBDT1, nparBDT1, RooFitScale1, 0);
	TF1* pdf2 = (TF1*)bkgBDTpdfNominal2->asTF(*score, RooArgList(*A0,*A1,*A2,*A3,*A4), *score); RooRealVar parBDT2[nparBDT2]={*A0,*A1,*A2,*A3,*A4}; TString pdfStr2 = expandNscaleFormula(bkgBDTpdfNominalFormula2, "score", parBDT2, nparBDT2, RooFitScale2, 0);
	TF1* pdf3 = (TF1*)bkgBDTpdfNominal3->asTF(*score, RooArgList(*B0,*B2,*B1,*B3,*B4), *score); RooRealVar parBDT3[nparBDT3]={*B0,*B2,*B1,*B3,*B4}; TString pdfStr3 = expandNscaleFormula(bkgBDTpdfNominalFormula3, "score", parBDT3, nparBDT3, RooFitScale3, 0);
	pdf0->Write("fBDT0.0");
	pdf1->Write("fBDT0.1");
	pdf2->Write("fBDT0.2");
	pdf3->Write("fBDT0.3");
	
	cout << "fBDT0=" << pdfStr0 << endl;
	cout << "fBDT1=" << pdfStr1 << endl;
	cout << "fBDT2=" << pdfStr2 << endl;
	cout << "fBDT3=" << pdfStr3 << endl;
	// TString pdfStr0 = pdf0->GetExpFormula("P")+"*"+tstr(RooFitScale0,10);
	// TString pdfStr1 = pdf1->GetExpFormula("P")+"*"+tstr(RooFitScale1,10);
	// TString pdfStr2 = pdf2->GetExpFormula("P")+"*"+tstr(RooFitScale2,10);
	// TString pdfStr3 = pdf3->GetExpFormula("P")+"*"+tstr(RooFitScale3,10);
	TNamed* bkgBDTpdfNominal0Str = new TNamed("fBDT0str.0",pdfStr0.Data());
	TNamed* bkgBDTpdfNominal1Str = new TNamed("fBDT0str.1",pdfStr1.Data());
	TNamed* bkgBDTpdfNominal2Str = new TNamed("fBDT0str.2",pdfStr2.Data());
	TNamed* bkgBDTpdfNominal3Str = new TNamed("fBDT0str.3",pdfStr3.Data());
	bkgBDTpdfNominal0Str->Write("fBDT0str.0");
	bkgBDTpdfNominal1Str->Write("fBDT0str.1");
	bkgBDTpdfNominal2Str->Write("fBDT0str.2");
	bkgBDTpdfNominal3Str->Write("fBDT0str.3");
	
	TF1* fSB0 = (TF1*)bkgm3bodypdfNominal0->asTF(*m3body, RooArgList(*k0,*k1), *m3body);
	TF1* fSB1 = (TF1*)bkgm3bodypdfNominal1->asTF(*m3body, RooArgList(*a0), *m3body);
	// TF1* fSB1 = (TF1*)bkgm3bodypdfNominal1->asTF(*m3body, RooArgList(*j0,*j1,*j2), *m3body);
	// TF1* fSB1 = (TF1*)bkgm3bodypdfNominal1->asTF(*m3body, RooArgList(*b0,*b1),     *m3body);
	TF1* fSB2 = (TF1*)bkgm3bodypdfNominal2->asTF(*m3body, RooArgList(*i0,*i1,*i2), *m3body);
	// TF1* fSB2 = (TF1*)bkgm3bodypdfNominal2->asTF(*m3body, RooArgList(*c0,*c1,*c2), *m3body);
	TF1* fSB3 = (TF1*)bkgm3bodypdfNominal3->asTF(*m3body, RooArgList(*h0,*h1),     *m3body);
	fSB0->Write("fSB0.0");
	fSB1->Write("fSB0.1");
	fSB2->Write("fSB0.2");
	fSB3->Write("fSB0.3");
	
	cout << "fSB0=" << fSB0->GetExpFormula("P") << endl;
	cout << "fSB1=" << fSB1->GetExpFormula("P") << endl;
	cout << "fSB2=" << fSB2->GetExpFormula("P") << endl;
	cout << "fSB3=" << fSB3->GetExpFormula("P") << endl;
	TString sbStr0 = fSB0->GetExpFormula("P")+"*"+tstr(RooFitScale0,10);
	TString sbStr1 = fSB1->GetExpFormula("P")+"*"+tstr(RooFitScale1,10);
	TString sbStr2 = fSB2->GetExpFormula("P")+"*"+tstr(RooFitScale2,10);
	TString sbStr3 = fSB3->GetExpFormula("P")+"*"+tstr(RooFitScale3,10);
	// TString sbStr0 = bkgm3bodypdfNominal0->formula()->GetTitle()+"*"+tstr(RooFitScale0,10);
	// TString sbStr1 = bkgm3bodypdfNominal1->formula()->GetTitle()+"*"+tstr(RooFitScale1,10);
	// TString sbStr2 = bkgm3bodypdfNominal2->formula()->GetTitle()+"*"+tstr(RooFitScale2,10);
	// TString sbStr3 = bkgm3bodypdfNominal3->formula()->GetTitle()+"*"+tstr(RooFitScale3,10);
	TNamed* bkgm3bodypdfNominal0Str = new TNamed("fSB0str.0",sbStr0.Data());
	TNamed* bkgm3bodypdfNominal1Str = new TNamed("fSB0str.1",sbStr1.Data());
	TNamed* bkgm3bodypdfNominal2Str = new TNamed("fSB0str.2",sbStr2.Data());
	TNamed* bkgm3bodypdfNominal3Str = new TNamed("fSB0str.3",sbStr3.Data());
	bkgm3bodypdfNominal0Str->Write("fSB0str.0");
	bkgm3bodypdfNominal1Str->Write("fSB0str.1");
	bkgm3bodypdfNominal2Str->Write("fSB0str.2");
	bkgm3bodypdfNominal3Str->Write("fSB0str.3");
	
	
	TVectorF chi2BDT(4);
	chi2BDT[0]=chi2BDTpdf0;
	chi2BDT[1]=chi2BDTpdf1;
	chi2BDT[2]=chi2BDTpdf2;
	chi2BDT[3]=chi2BDTpdf3;
	chi2BDT.Write("chi2BDT");
	
	TVectorF NSB(1);
	NSB[0]=ncandidatesD;
	NSB.Write("nSB");
	
	TVectorF NSR(4);
	TVectorF chi2SB(4);
	NSR[0]=RooFitScale0*nSR00; chi2SB[0]=chi2SBpdf0;
	NSR[1]=RooFitScale1*nSR01; chi2SB[1]=chi2SBpdf1;
	NSR[2]=RooFitScale2*nSR02; chi2SB[2]=chi2SBpdf2;
	NSR[3]=RooFitScale3*nSR03; chi2SB[3]=chi2SBpdf3;
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

	unsigned int nNSR = vxSRmin.size();
	TVectorF vXSRmin(nNSR); TVectorF vXSRmax(nNSR);
	TVectorF vNSR00(nNSR); TVectorF vdNSR00(nNSR);
	TVectorF vNSR01(nNSR); TVectorF vdNSR01(nNSR);
	TVectorF vNSR02(nNSR); TVectorF vdNSR02(nNSR);
	TVectorF vNSR03(nNSR); TVectorF vdNSR03(nNSR);
	for(unsigned int k=0 ; k<vxSRmin.size() ; ++k)
	{
		vXSRmin[k]=vxSRmin[k]; vXSRmax[k]=vxSRmax[k];
		vNSR00[k]=RooFitScale0*vnSR00[k]; vdNSR00[k]=(vnSR00[k]/nSB00)*sqrt(RooFitScale0*nSB00);
		vNSR01[k]=RooFitScale1*vnSR01[k]; vdNSR01[k]=(vnSR01[k]/nSB01)*sqrt(RooFitScale1*nSB01);
		vNSR02[k]=RooFitScale2*vnSR02[k]; vdNSR02[k]=(vnSR02[k]/nSB02)*sqrt(RooFitScale2*nSB02);
		vNSR03[k]=RooFitScale3*vnSR03[k]; vdNSR03[k]=(vnSR03[k]/nSB03)*sqrt(RooFitScale3*nSB03);
	}
	vXSRmin.Write("xSRmin"); vXSRmax.Write("xSRmax");
	vNSR00.Write("nSR00");   vNSR01.Write("nSR01");   vNSR02.Write("nSR02");   vNSR03.Write("nSR03");
	vdNSR00.Write("dnSR00"); vdNSR01.Write("dnSR01"); vdNSR02.Write("dnSR02"); vdNSR03.Write("dnSR03");
	
	// and to read the file do: TVectorF* NSB0 = (TVectorF*)fpdf->Get("nSB0"); float nSB0 = ((*NSB0))[0];
	fpdf->Write();
	fpdf->Close();
	delete fpdf;
	
	string filename = "BkgEstimate."+(string)sidebands+"."+(string)bdtcutoff+"."+(string)bdtmaxusr+".txt";
	ofstream* ofstr = new ofstream(filename.c_str());
	stringstream strm; strm << "";
	string tmp = "";
	string summary = "";
	summary += "========================================================================================\n";
	summary += "m3body left sideband    ["+str(m3bodyMin,0)+","+str(xbmin,0)+"] MeV\n";
	summary += "m3body right sideband   ["+str(xbmax,0)+","+str(m3bodyMax,0)+"] MeV\n";
	summary += "m3body signal region    ["+str(xSRmin,0)+","+str(xSRmax,0)+"] MeV (for sigma="+str(sigma->getVal(),2)+" MeV)\n";
	summary += "Minimal BDT cut (x0)    "+str(xbdtmin,3)+"\n";
	summary += "Optimal BDT cut (x1)    "+str(xbdtopt,3)+"\n";
	summary += "Maximal BDT cut (x2)    "+str(xbdtmax,3)+"\n";
	summary += "========================================================================================\n";
	summary += "Ch2^2/DOF BDT           "+str(chi2BDTpdf0,2)+" | "+str(chi2BDTpdf1,2)+" | "+str(chi2BDTpdf2,2)+" | "+str(chi2BDTpdf3,2)+"\n";
	summary += "Ch2^2/DOF Signal        "+str(chi2Signalpdf,2)+"\n";
	summary += "Ch2^2/DOF SB            "+str(chi2SBpdf0,2)+" | "+str(chi2SBpdf1,2)+" | "+str(chi2SBpdf2,2)+" | "+str(chi2SBpdf3,2)+"\n";
	summary += "========================================================================================\n";
	summary += "Count:    nSB0          "+str(ncandidatesD,0)+"\n";
	summary += "Mass fit: nSB0          "+str(nSB00,3)+" | "+str(nSB01,3)+" | "+str(nSB02,3)+"\n";
	summary += "Mass fit: nSR0          "+str(nSR00,3)+" | "+str(nSR02,3)+" | "+str(nSR02,3)+"\n";
	summary += "RooFit scale            "+str(RooFitScale0,5)+" | "+str(RooFitScale1,5)+" | "+str(RooFitScale2,5)+"\n";
	summary += "Scaled nSR0             "+str(RooFitScale0*nSR00,3)+" +- "+str((nSR00/nSB00)*sqrt(RooFitScale0*nSB00),3)+" | "
										 +str(RooFitScale1*nSR01,3)+" +- "+str((nSR01/nSB01)*sqrt(RooFitScale1*nSB01),3)+" | "
										 +str(RooFitScale2*nSR02,3)+" +- "+str((nSR02/nSB02)*sqrt(RooFitScale2*nSB02),3)+" | "
										 +str(RooFitScale3*nSR03,3)+" +- "+str((nSR03/nSB03)*sqrt(RooFitScale3*nSB03),3)+"\n";
	summary += "Fit: f01                "+str(transferFactor0,6)+" | "+str(transferFactor1,6)+" | "+str(transferFactor2,6)+" | "+str(transferFactor3,6)+"\n";
	summary += "Cut: nSB1               "+str(ncandidatesD_CR1,3)+" +- "+str(sqrt(ncandidatesD_CR1),3)+"\n";
	summary += "Scaled: nSR1            "+str((RooFitScale0*nSR00)*transferFactor0,3)+" +- "+str((nSR00/nSB00)*sqrt(RooFitScale0*nSB00)*transferFactor0,3)+" | "
										 +str((RooFitScale1*nSR01)*transferFactor0,3)+" +- "+str((nSR01/nSB01)*sqrt(RooFitScale1*nSB01)*transferFactor0,3)+" | "
										 +str((RooFitScale2*nSR02)*transferFactor0,3)+" +- "+str((nSR02/nSB02)*sqrt(RooFitScale2*nSB02)*transferFactor0,3)+" | "
										 +str((RooFitScale3*nSR03)*transferFactor0,3)+" +- "+str((nSR03/nSB03)*sqrt(RooFitScale3*nSB03)*transferFactor0,3)+"\n";
	summary += "Acc*Eff in SR1          "+str(efficiencyS,3)+"\%\n";
	summary += "N candidates(S)         "+str(ncandidatesS,0)+"\n";
	summary += "N candidates(D)         "+str(ncandidatesDfull,0)+"\n";
	summary += "========================================================================================\n";
	cout << summary << endl;
	(*ofstr) << summary;
	ofstr->close();
}
