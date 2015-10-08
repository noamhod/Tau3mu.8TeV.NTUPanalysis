#include "all.h"
#include "roofit.h"
#include "RooAngular.h"
#include "RooCollinsSoper.h"

/*
	1. There are 2 kinds of unbinned, non-detector weights:
		a. cross section weight, wXS, because the samples are binned in mass - relevant for the three models
		b. model weight, xR, the ratio between ZP or KK and Z0 - relevant for the ZP and KK models only
	2. There is the acceptance binned weight (trigger efficiency is ignored currently)
	3. The truth PDF is A0*(1+x^2)+A4*x   ==>  Afb = 3B/8A
	4. The PDF should be corrected for the acceptance by a simple multiplication:      PDF(x) = PDFtru(x)*Acc(x)
	5. Each event should be weighted (unbinned) according to the non-detector weights: w = wXS*wR
*/

// bins for the mHat histos
const  int      nxbins = 30;
static Double_t logXmin;
static Double_t logXmax;
static Double_t logXbinwidth;
static Double_t xbins[nxbins+1];
Int_t  iMassNbins = 30; // for non log bins

double   costmin   = minCosTheta;
double   costmax   = maxCosTheta;
int      ncostbins = nCosThetaBins;

double minA0 = -10.;
double maxA0 = +10.;
double minA4 = -10.; //-2.666;
double maxA4 = +10.; //+2.666;
double _A0 = 0.;
double _A4 = 0.;
struct fitpars
{
	double A0;
	double A4;
};
TRandom* randGen;
vector<vector<fitpars> > vvInitialGuess;

TFile* file = new TFile("weights.root", "READ");
vector<TH1D*> vhMassBins;
vector<TH1D*> vhAfbBins;

vector<TTree*> vtData;
float mass_tru, mass_rec, mass_wgt;
float cost_tru, cost_rec, cost_wgt;
float xscn_wgt;

vector<TH1D*>  vhAcc;
vector<TH1D*>  vhMass;
vector<TString> vModelName;
TLegend* leg;
TLegend* leg_mHat;
TPaveText* pvtxt;
TText* txt;
vector<TGraphAsymmErrors*> gMpoissonErr;

TPaveText* pvtxt_lumi;
TPaveText* pvtxt_atlas;
TCanvas* cnvAfb;
TPad *pad_mHat;
TPad *pad_Afb;
TPad *pad_compare;

RooRealVar* cosThe; // the variable 
RooRealVar* weight; // the weight
RooRealVar* A0;    // the parameter to find
RooRealVar* A4;    // the parameter to find


vector<RooAbsData*>   vUnbinnedDataSet; // Roo Data holder
vector<RooDataHist*>  vBinnedDataSet;
vector<RooFitResult*> vFitResult;

vector<bool>          vbFitStatus;
vector<vector<bool> > vvbFitStatus;


vector<RooAbsPdf*> vModel;   // the final model pdf
vector<RooAbsPdf*> vDetAcc;  // will be taken from the acceptance histogram
vector<RooAbsPdf*> vRecEff;  // will be taken from the reconstruction efficiency histogram
vector<RooAbsPdf*> vTrigEff; // will be taken from the trigger efficiency histogram
vector<RooAbsPdf*> vEffPdf;  // will be the product of all the eff-like pdf's (=rec*trig*acc)
RooAbsPdf* sigPdf;           // will be the truth pdf

vector<RooDataHist*> vrdhAcc;
vector<RooHistPdf*>  vrhpdfAcc;

///////////////////
// for generation
RooRealVar*  A0_gen[4];
RooRealVar*  A4_gen[4];
RooAbsPdf*   sigPdf_gen[4];
TH1D*        hAcc[4];
RooDataHist* rdhAcc[4];
RooHistPdf*  rhpdfAcc[4];
RooProdPdf*  model[4];
RooAbsData*  rad[4];
///////////////////

//////////////////////////////////
// flags to config the run
bool drawAfbErrArea  = true;
bool doLumiXSweights = false; //!!! // true if want to scale binned samples to 1 smooth sample and to scale MC to data luminosity. This affects mostly the errors
bool doBinned        = false;
bool doGeneration    = false;
bool doBinomialError = false;
bool drawKKAfb       = false;
Int_t Ngen = 10000;
Int_t minEntries2Fit = 5;
void printRunConfig()
{
	cout << "//////////////////////////////////" << endl;
	cout << "// flags to config the run" << endl;
	cout << "bool drawAfbErrArea  = " << drawAfbErrArea << ";" << endl;
	cout << "bool doLumiXSweights = " << doLumiXSweights << "; // true if want to scale binned samples to 1 smooth sample and to scale MC to data luminosity. This affects mostly the errors " << endl;
	cout << "bool doBinned        = " << doBinned << ";" << endl;
	cout << "bool doGeneration    = " << doGeneration << ";" << endl;
	cout << "bool doBinomialError = " << doBinomialError << ";" << endl;
	cout << "Int_t Ngen = " << Ngen << ";" << endl;
	cout << "Int_t minEntries2Fit = " << minEntries2Fit << ";" << endl;
	cout << "//////////////////////////////////" << endl;
}
//////////////////////////////////

double randomizeItialGuess(double min, double max)
{
	_DEBUG("randomizeItialGuess");
	return min + (max-min)*randGen->Uniform(); // Uniform(x1=0) returns a uniform deviate on the interval [0,x1].
}

double getAfb(double a, double b, bool validate)
{
	if(validate)
	{
		return (a>0.) ? 3.*b/(8.*a) : -999;
	}
	else return 3.*b/(8.*a);
	
	return -999.;
}

double getAfb(double a4)
{
	return 3./8.*a4;
}

double getAfbErr(double Afb, int N)
{
	return (N>0.) ? sqrt((1-Afb*Afb)/N) : -999;
}

double getAfbErr(double dA4)
{
	return (3./8.)*dA4;
}

void setBranches(int mod)
{
	_DEBUG("setBranches");
	vtData[mod]->SetBranchAddress( "mass_tru", &mass_tru );
	vtData[mod]->SetBranchAddress( "mass_rec", &mass_rec );
	vtData[mod]->SetBranchAddress( "mass_wgt", &mass_wgt );
	vtData[mod]->SetBranchAddress( "cost_tru", &cost_tru );
	vtData[mod]->SetBranchAddress( "cost_rec", &cost_rec );
	vtData[mod]->SetBranchAddress( "cost_wgt", &cost_wgt );
	vtData[mod]->SetBranchAddress( "xscn_wgt", &xscn_wgt );
}

void setLogMassBins(Double_t xmin, Double_t xmax)
{
	_DEBUG("setLogMassBins");
	logXmin  = log10(xmin);
	logXmax  = log10(xmax);
	logXbinwidth = (Double_t)( (logXmax-logXmin)/(Double_t)nxbins );
	xbins[0] = xmin;
	for(Int_t i=1 ; i<=nxbins ; i++) xbins[i] = TMath::Power( 10,(logXmin + i*logXbinwidth) );
}

double getMassBinCenter(int massBin)
{
	return xbins[massBin-1]+(xbins[massBin]-xbins[massBin-1])/2.;
}

fitpars getTheoryAfb(int massBin, int mod)
{
	_DEBUG("getTheoryAfb");
	double imass = getMassBinCenter(massBin);
	fitpars fp;
	fp.A0 = 3./8.;
	double normflat = 1.;
	double normgausup = 0.8;
	double normgausdwn = 1.;
	double expscale = 10.;
	double expunits = 1000.;
	if(mod==Z0 || mod==KK) fp.A4 = normflat*(1.-TMath::Exp(-expscale*imass/expunits));
	// else if(mod==ZP)       fp.A4 = normflat*(1.-TMath::Exp(-expscale*imass/expunits)) - normgausdwn*TMath::Gaus(imass,520,50);
	else if(mod==ZP)       fp.A4 = normflat*(1.-TMath::Exp(-expscale*imass/expunits)) - normgausdwn*TMath::Gaus(imass,800,100);
	// else                   fp.A4 = normflat*(1.-TMath::Exp(-expscale*imass/expunits)) + normgausup*TMath::Gaus(imass,800,100);
	else                   fp.A4 = normflat*(1.-TMath::Exp(-expscale*imass/expunits));
	_INFO("A0="+tostring(fp.A0)+", A4="+tostring(fp.A4));
	return fp;
}

void generateToy(int massBin, int mod, int N)
{
	_DEBUG("generateToy");
	
	TString sMassBin = (TString)tostring(massBin);

	fitpars fp = getTheoryAfb(massBin,mod);
	A0_gen[mod] = new RooRealVar(("A0_gen"+tostring(massBin)+"_mod"+tostring(mod)).c_str(),"A0_gen",fp.A0);
	A4_gen[mod] = new RooRealVar(("A4_gen"+tostring(massBin)+"_mod"+tostring(mod)).c_str(),"A4_gen",fp.A4);

	// sigPdf_gen[mod] = new RooAngular(("SignalPdf_bin"+tostring(massBin)+"_mod"+tostring(mod)).c_str(), "SignalPdf_gen", *cosThe, *A0_gen[mod], *A4_gen[mod]);
	sigPdf_gen[mod] = new RooCollinsSoper(("SignalPdf_bin"+tostring(massBin)+"_mod"+tostring(mod)).c_str(), "SignalPdf_gen", *cosThe, *A0_gen[mod], *A4_gen[mod]);
	
	TString sName, sId;
	
	switch(mod)
	{
		case Z0:
			sName = "Z^{0}[toyMC]";
			sId = "Z0d3pd";
			break;
		case ZP:
			sName = "Z'_{SSM}[toyMC]";
			sId = "ZP";
			break;
		case KK:
			sName = "S^{1}/Z_{2} KK[toyMC]";
			sId = "KK";
			break;
		case DT:
			sName = "Data[toyMC]";
			sId = "DT";
			break;
	}

	_DEBUG("");
	
	if(mod==DT) hAcc[mod] = (TH1D*)file->Get("cosTheta_histograms/hCosThZ0d3pd_acceptance_"+sMassBin)->Clone("");      // !!!!!!!!!!!!!!!!!!!!!!!! TO CHANGE
	else        hAcc[mod] = (TH1D*)file->Get("cosTheta_histograms/hCosTh"+sId+"_acceptance_"+sMassBin)->Clone("");
	rdhAcc[mod]   = new RooDataHist("rdhAcc_gen"+sId,"rdhAcc_gen"+sId,RooArgSet(*cosThe),hAcc[mod]);
	rhpdfAcc[mod] = new RooHistPdf("rhpdfAcc_gen"+sId,"rhpdfAcc_gen"+sId,RooArgSet(*cosThe),*rdhAcc[mod],5); // last argument is the order of polinomial interpulation
	model[mod]    = new RooProdPdf("model_"+sId+"_gen","truPdf*accPdf_gen",*sigPdf_gen[mod],*rhpdfAcc[mod]);
	rad[mod]      = model[mod]->generate(*cosThe,N);
	_INFO("Data entries = "+tostring(rad[mod]->numEntries()));
}

void init(int massBin, int mod)
{
	_DEBUG("init");

	TH1D* hMassBinsDummy = (TH1D*)file->Get("all_histograms/hDummy_afb")->Clone("");
	Double_t iMassMin = hMassBinsDummy->GetBinLowEdge(massBin);
	Double_t iMassMax = iMassMin + hMassBinsDummy->GetBinWidth(massBin);
	
	TString sMassBin = (TString)tostring(massBin);
	string sMassMin = tostring((double)iMassMin);
	string sMassMax = tostring((double)iMassMax);
	TString sTitle = "Mass-bin[" + sMassBin + "] " + sMassMin + "#rightarrow" + sMassMax + " GeV";
	
	setLogMassBins(iMassMin, iMassMax);

	cosThe = new RooRealVar("cosTheta","cos#theta*",costmin,costmax);
	cosThe->setRange("range_cosThe",costmin,costmax);
	cosThe->setBins(ncostbins);
	
	weight = new RooRealVar("weight","weight",0.,1e10);
	//weight->setRange("range_weight",costmin,costmax);
	//weight->setBins(ncostbins);

	_A0 = randomizeItialGuess(minA0,maxA0);
	_A4 = randomizeItialGuess(minA4,maxA4);
	fitpars fp;
	fp.A0 = _A0;
	fp.A4 = _A4;
	vvInitialGuess[massBin-1].push_back(fp);
	A0 = new RooRealVar("A0","A0",_A0,minA0,maxA0);
	A4 = new RooRealVar("A4","A4",_A4,minA4,maxA4);
	A0->setError(0.00001);
	A4->setError(0.00001);
	
	// sigPdf = new RooAngular("SignalPdf", "SignalPdf", *cosThe, *A0, *A4);
	sigPdf = new RooCollinsSoper("SignalPdf", "SignalPdf", *cosThe, *A0, *A4);
	
	TString sName, sId, sIdShort, sChannelFit, sChannelMass;
	Int_t fillStyle = 0;
	Int_t lineStyle = 0;
	Int_t markerStyle = 0;
	
	switch(mod)
	{
		case Z0:
			sName = "Z^{0}";
			sId = "Z0d3pd";
			sChannelFit = "#gamma/Z^{0} [MC10b]: A_{FB} fit";
			sChannelMass = "#gamma/Z^{0} [MC10b]: Events";
			fillStyle = 3003;
			lineStyle = 2;
			markerStyle = 20;
			break;
		case ZP:
			sName = "Z'_{SSM}";
			sId = "ZP";
			sChannelFit = "1 TeV Z'_{SSM} [Template MC10b]: A_{FB} fit";
			sChannelMass = "1 TeV Z'_{SSM} [Template MC10b]: Events";
			fillStyle = 3017;
			lineStyle = 3;
			markerStyle = 22;
			break;
		case KK:
			sName = "S^{1}/Z_{2} KK";
			sId = "KK";
			sChannelFit = "1 TeV #gamma_{KK}/Z_{KK} [Template MC10b]: A_{FB} fit";
			sChannelMass = "1 TeV #gamma_{KK}/Z_{KK} [Template MC10b]: Events";
			fillStyle = 3018;
			lineStyle = 5;
			markerStyle = 23;
			break;
		case DT:
			sName = "Data";
			sId = "DT";
			sChannelFit = "2011 Data: A_{FB} fit";
			sChannelMass = "2011 Data: Events";
			fillStyle = 0;
			lineStyle = 1;
			markerStyle = 24;
			break;
	}
	
	if(mod==Z0) sIdShort = "Z0";
	else        sIdShort = sId;
	
	if(massBin==1) // only one copy of this vector
	{
		vhMassBins.push_back( (TH1D*)file->Get("all_histograms/hMassTmp")->Clone("") );
		vhMassBins[mod]->Reset();
		vhMassBins[mod]->SetDefaultSumw2(); // The error per bin will be computed as sqrt(sum of squares of weight) for each bin. 
		vhMassBins[mod]->SetTitle("");
		vhMassBins[mod]->SetYTitle("Events");
		vhMassBins[mod]->SetLineColor(vModelColor[mod]);
		vhMassBins[mod]->SetLineWidth(1);
		vhMassBins[mod]->SetLineStyle(lineStyle);
		if((mod==KK && drawKKAfb) || mod!=KK) leg_mHat->AddEntry( vhMassBins[mod], sChannelMass, "l");
		
		vhAfbBins.push_back( (TH1D*)file->Get("all_histograms/hDummy_afb")->Clone("") );
		vhAfbBins[mod]->Reset();
		vhAfbBins[mod]->SetTitle("");
		vhAfbBins[mod]->SetXTitle( "m_{#mu#mu} GeV" );
		vhAfbBins[mod]->SetYTitle( "A_{FB}" );
		vhAfbBins[mod]->SetLineColor(vModelColor[mod]);
		if(drawAfbErrArea  &&  mod!=DT)
		{
			vhAfbBins[mod]->SetFillColor(vModelColor[mod]);
			vhAfbBins[mod]->SetFillStyle(fillStyle);
			vhAfbBins[mod]->SetLineWidth(1);
			vhAfbBins[mod]->SetMarkerSize(0);
			vhAfbBins[mod]->SetMarkerColor(0);
			if((mod==KK && drawKKAfb) || mod!=KK) leg_mHat->AddEntry( vhAfbBins[mod], sChannelFit, "f");
		}
		else
		{	
			vhAfbBins[mod]->SetLineWidth(2);
			vhAfbBins[mod]->SetMarkerSize(2);
			vhAfbBins[mod]->SetMarkerStyle(markerStyle);
			vhAfbBins[mod]->SetMarkerColor(vModelColor[mod]);
			if((mod==KK && drawKKAfb) || mod!=KK) leg_mHat->AddEntry( vhAfbBins[mod], sChannelFit, "lep");
		}             
	}
	
	vModelName.push_back( sName );
	
	//vhMass.push_back( new TH1D("hMass_"+sId,"",iMassNbins, iMassMin, iMassMax) );
	vhMass.push_back( new TH1D("hMass_"+sId,"",nxbins,xbins) );
	vhMass[mod]->Reset();
	vhMass[mod]->SetDefaultSumw2(); // The error per bin will be computed as sqrt(sum of squares of weight) for each bin.
	vhMass[mod]->SetLineColor(vModelColor[mod]);
	vhMass[mod]->SetTitle(sTitle);
	vhMass[mod]->SetXTitle("m_{#mu#mu} GeV");
	vhMass[mod]->SetYTitle("Events");
	vhMass[mod]->GetXaxis()->SetMoreLogLabels(); 
	vhMass[mod]->GetXaxis()->SetMoreLogLabels(); 
	
	vtData.push_back( (TTree*)file->Get("ntuples/tree_"+sIdShort+"_"+sMassBin) );
	if(mod==DT) vhAcc.push_back( (TH1D*)file->Get("cosTheta_histograms/hCosThZ0d3pd_acceptance_"+sMassBin)->Clone("") ); // !!!!!!!!!!!!!!!!!!!!!!!! TO CHANGE 
	else        vhAcc.push_back( (TH1D*)file->Get("cosTheta_histograms/hCosTh"+sId+"_acceptance_"+sMassBin)->Clone("") );
	vrdhAcc.push_back( new RooDataHist("rdhAcc"+sId,"rdhAcc"+sId,RooArgSet(*cosThe),vhAcc[mod]) );
	vrhpdfAcc.push_back( new RooHistPdf("rhpdfAcc"+sId,"rhpdfAcc"+sId,RooArgSet(*cosThe),*vrdhAcc[mod],5) ); // last argument is the order of polinomial interpulation
	vDetAcc.push_back( vrhpdfAcc[mod] );
	vModel.push_back( new RooProdPdf("model_"+sId,"truPdf*accPdf",*sigPdf,*vDetAcc[mod]) );
	if(doGeneration)
	{
		Int_t nGen = Ngen;
		// if(mod==DT) nGen = (Int_t)(Ngen/TMath::Exp(10.*getMassBinCenter(massBin)/1000.));
		if(mod==DT) nGen = (Int_t)(Ngen/TMath::Power(10.,massBin-2)); // for 6 bins: Ngen/0.1, Ngen/1, Ngen/10, Ngen/100, Ngen/1000, Ngen/10000
		generateToy(massBin, mod, nGen);
		RooAbsData* r = (RooAbsData*)rad[mod]->Clone("");
		Int_t N = r->numEntries();
		vUnbinnedDataSet.push_back( r );
		_INFO("Data entries = "+tostring(N));
	}
	else
	{
		if(mod==DT || doGeneration) vUnbinnedDataSet.push_back( new RooDataSet("data_"+sId,"data_"+sId,RooArgSet(*cosThe)) ); // no weights
		else                        vUnbinnedDataSet.push_back( new RooDataSet("data_"+sId,"data_"+sId,RooArgSet(*cosThe,*weight),WeightVar(weight->GetName())) );
	}
}

void reset()
{
	_DEBUG("reset");
	
	for(int i=0 ; i<(int)vModelName.size() ; i++)
	{
		if(vhMass[i]!=NULL) delete vhMass[i];
		if(vrdhAcc[i]!=NULL) delete vrdhAcc[i];
		if(vrhpdfAcc[i]!=NULL) delete vrhpdfAcc[i];
		if(vModel[i]!=NULL) delete vModel[i];
		if(vUnbinnedDataSet[i]!=NULL) delete vUnbinnedDataSet[i];
		if(doBinned && vBinnedDataSet.size()>0 && vBinnedDataSet[i]!=NULL) delete vBinnedDataSet[i];
	}

	_DEBUG("");
	
	vModelName.clear();
	vhMass.clear();
	vtData.clear();
	vhAcc.clear();
	vrdhAcc.clear();
	vrhpdfAcc.clear();
	vDetAcc.clear();
	vModel.clear();
	vUnbinnedDataSet.clear();
	if(doBinned) vBinnedDataSet.clear();
	
	_DEBUG("");
	
	delete cosThe;
	delete weight;
	delete A0;
	delete A4;
	delete sigPdf;
	
	_DEBUG("");
	
	// if(doGeneration)
	// {
		// for(int i=0 ; i<4 ; i++)
		// {
			// _DEBUG("");
			// if(hAcc[i]!=NULL)       delete hAcc[i];
			// _DEBUG("");
			// if(rdhAcc[i]!=NULL)     delete rdhAcc[i];
			// _DEBUG("");
			// if(rhpdfAcc[i]!=NULL)   delete rhpdfAcc[i];
			// _DEBUG("");
			// if(model[i]!=NULL)      delete model[i];
			// _DEBUG("");
			// if(rad[i]!=NULL)        delete rad[i];
			// _DEBUG("");
			// if(sigPdf_gen[i]!=NULL) delete sigPdf_gen[i];
			// _DEBUG("");
			// if(A0_gen[i]!=NULL)      delete A0_gen[i];
			// _DEBUG("");
			// if(A4_gen[i]!=NULL)      delete A4_gen[i];
			// _DEBUG("");
		// }
	// }
}

Int_t loop(int mod)
{
	_DEBUG("loop");
	Int_t N = 0;
	
	if(doGeneration) // this is unweighted
	{
		N = rad[mod]->numEntries();
	}
	else
	{
		setBranches(mod);
		N = vtData[mod]->GetEntries();
		_DEBUG("N = "+tostring(N));
		for(Int_t i=0 ; i<N ; i++)
		{
			vtData[mod]->GetEntry(i);
			*cosThe = cost_rec;
			float w;
			if(mod==DT)
			{
				w = 1.;
				vhMass[mod]->Fill(mass_rec);
				vhMassBins[mod]->Fill(mass_rec);
			}
			else if(mod==Z0)
			{
				if(doLumiXSweights) w = (xscn_wgt*luminosity)*1.; // THIS IS CORRECT (cost_wgt=1) !!!!!
				else          w = 1.;
				vhMass[mod]->Fill(mass_rec,xscn_wgt*luminosity);
				vhMassBins[mod]->Fill(mass_rec,xscn_wgt*luminosity);
			}
			else
			{
				if(doLumiXSweights) w = (xscn_wgt*luminosity)*cost_wgt; // THIS IS CORRECT !!!!!
				else          w = cost_wgt;
				vhMass[mod]->Fill(mass_rec,xscn_wgt*mass_wgt*luminosity);
				vhMassBins[mod]->Fill(mass_rec,xscn_wgt*mass_wgt*luminosity);
			}
			
			if(mod==DT) vUnbinnedDataSet[mod]->add(RooArgSet(*cosThe));   // UNWEIGHTED
			else        vUnbinnedDataSet[mod]->add(RooArgSet(*cosThe),w); // WEIGHTED   
		}
	}
	
	if(doBinned)
	{
		RooDataSet* rds = (RooDataSet*)vUnbinnedDataSet[mod];
		vBinnedDataSet.push_back( rds->binnedClone() );
	}
	
	return N;
}

bool minuitStatus(TMinuit* m) 
{
	_DEBUG("minuitStatus");
	if (!m) return false;

	TString stat = gMinuit->fCstatu;
	_INFO("Minuit: "+(string)stat+". ");
	if ( stat.Contains("SUCCESSFUL")  || stat.Contains("CONVERGED")  ||  stat.Contains("OK") ) return true;

	return false;
}

RooFitResult* fit(int mod)
{
	_DEBUG("fit");
	TMinuit* gFit(0);
	
	_A0 = 0.;
	_A4 = 0.;
	
	const RooArgSet* fitParsInital;
	if(doBinned) fitParsInital = vModel[mod]->getParameters(vBinnedDataSet[mod],false);
	else         fitParsInital = vModel[mod]->getParameters(vUnbinnedDataSet[mod],false);
	RooRealVar* x = (RooRealVar*)fitParsInital->find("A0");
	RooRealVar* y = (RooRealVar*)fitParsInital->find("A4");
	if(x) *x = 0.;
	if(y) *y = 0.;
	delete fitParsInital;
	
	RooFitResult* fitresult;
	if(doBinned)
	{
		if(mod==DT || doGeneration)
		{
			if(vBinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is weighted $$$$$$$$$$$");
			fitresult = vModel[mod]->fitTo( *vBinnedDataSet[mod],Minos(kTRUE),Range("range_cosThe"),Strategy(2),Save(kTRUE),Timer(kTRUE),NumCPU(8));
		}
		else 
		{
			if(!vBinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is unweighted $$$$$$$$$$$");
			fitresult = vModel[mod]->fitTo( *vUnbinnedDataSet[mod],Minos(kTRUE),Range("range_cosThe"),Strategy(2),Save(kTRUE),Timer(kTRUE),SumW2Error(kTRUE),NumCPU(8));		
		}
	}
	else
	{
		if(mod==DT || doGeneration)
		{
			if(vUnbinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is weighted $$$$$$$$$$$");
			fitresult = vModel[mod]->fitTo( *vUnbinnedDataSet[mod],Minos(kTRUE),Range("range_cosThe"),Strategy(2),Save(kTRUE),Timer(kTRUE),NumCPU(8));
		}
		else
		{
			if(!vUnbinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is unweighted $$$$$$$$$$$");
			fitresult = vModel[mod]->fitTo( *vUnbinnedDataSet[mod],Minos(kTRUE),Range("range_cosThe"),Strategy(2),Save(kTRUE),Timer(kTRUE),SumW2Error(kTRUE),NumCPU(8));
		}
	}
	gFit = gMinuit;
	
	vbFitStatus.push_back( minuitStatus(gFit) );
	
	fitresult->Print("v");
	
	return fitresult;
}

void plot(int mod, TVirtualPad* pad)
{
	_DEBUG("plot");
	pad->cd();
	
	RooPlot* cosThetaFrame = cosThe->frame(Name("cosThetaFrame"), Title( vModelName[mod] ));
	if(doBinned)
	{
		if(mod==DT || doGeneration)
		{
			if(vBinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is weighted $$$$$$$$$$$");
			vBinnedDataSet[mod]->plotOn(cosThetaFrame,Name("cos#theta*"),XErrorSize(0),MarkerSize(0.3),Binning(ncostbins),NumCPU(8));
		}
		else
		{
			if(!vBinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is unweighted $$$$$$$$$$$");
			vBinnedDataSet[mod]->plotOn(cosThetaFrame,Name("cos#theta*"),XErrorSize(0),MarkerSize(0.3),Binning(ncostbins),DataError(RooAbsData::SumW2),NumCPU(8));
		}
	}
	else
	{
		if(mod==DT || doGeneration)
		{
			if(vUnbinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is weighted $$$$$$$$$$$");
			vUnbinnedDataSet[mod]->plotOn(cosThetaFrame,Name("cos#theta*"),XErrorSize(0),MarkerSize(0.3),Binning(ncostbins),NumCPU(8));
		}
		else
		{
			if(!vUnbinnedDataSet[mod]->isWeighted()) _WARNING("$$$$$$$$$$$ The dataset is unweighted $$$$$$$$$$$");
			vUnbinnedDataSet[mod]->plotOn(cosThetaFrame,Name("cos#theta*"),XErrorSize(0),MarkerSize(0.3),Binning(ncostbins),DataError(RooAbsData::SumW2),NumCPU(8));
		}
	}
	vDetAcc[mod]->plotOn(cosThetaFrame,LineWidth(1),LineColor(cAcceptance),NumCPU(8));
	//vDetAcc[mod]->plotOn(cosThetaFrame,LineWidth(1),LineColor(cAcceptance),NormRange("range_cosThe"));

	vModel[mod]->plotOn(cosThetaFrame,LineWidth(1),LineColor(cPdf),NormRange("range_cosThe"),NumCPU(8));
	//vModel[mod]->paramOn(cosThetaFrame,Layout(0.7,1.,0.4),Format("NEU",AutoPrecision(1)));
	vModel[mod]->paramOn(cosThetaFrame,Layout(0.4,0.88,1), Format("NEU", AutoPrecision(1)));
	cosThetaFrame->getAttText()->SetTextSize(0.05);
	//cosThetaFrame->getAttLine()->SetLineWidth(0.05);
	
	pad->SetLeftMargin(0.2);
	cosThetaFrame->SetTitleOffset(2,"Y");
	cosThetaFrame->Draw();
}

void plotempty(TVirtualPad* pad)
{
	_DEBUG("plotempty");
	pad->cd();
	
	// RooPlot* cosThetaFrame = cosThe->frame(Name("cosThetaFrame"), Title( vModelName[mod] ));
	// cosThetaFrame->getAttText()->SetTextSize(0.05);
	
	pad->SetLeftMargin(0.2);
	// cosThetaFrame->SetTitleOffset(2,"Y");
	// cosThetaFrame->Draw();
}

void getFit(int mod)
{
	_DEBUG("getFit");
	vFitResult.push_back( fit(mod) );
}


void drawAfb()
{
	_DEBUG("drawAfb");
	/*
	resetHistogramErrors(vhMassBins[Z0]);
	resetHistogramErrors(vhMassBins[ZP]);
	resetHistogramErrors(vhMassBins[KK]);
	resetHistogramErrors(vhMassBins[DT]);

	bool isXerr = false;
	//gMpoissonErr.push_back( GetPoissonizedGraph(vhMassBins[Z0], isXerr) );
	//gMpoissonErr.push_back( GetSqrtErrorsGraph(vhMassBins[Z0], isXerr) );
	gMpoissonErr.push_back( GetDefaultErrorsGraph(vhMassBins[Z0], isXerr) );
	gMpoissonErr[Z0]->SetMarkerStyle(1);
	gMpoissonErr[Z0]->SetMarkerSize(1);
	gMpoissonErr[Z0]->SetLineColor(vModelColor[Z0]);
	gMpoissonErr[Z0]->SetLineStyle(2);

	//gMpoissonErr.push_back( GetPoissonizedGraph(vhMassBins[ZP], isXerr) );
	//gMpoissonErr.push_back( GetSqrtErrorsGraph(vhMassBins[ZP], isXerr) );
	gMpoissonErr.push_back( GetDefaultErrorsGraph(vhMassBins[ZP], isXerr) );
	gMpoissonErr[ZP]->SetMarkerStyle(1);
	gMpoissonErr[ZP]->SetMarkerSize(1);
	gMpoissonErr[ZP]->SetLineColor(vModelColor[ZP]);	
	gMpoissonErr[ZP]->SetLineStyle(3);	

	//gMpoissonErr.push_back( GetPoissonizedGraph(vhMassBins[DT], isXerr) );
	//gMpoissonErr.push_back( GetSqrtErrorsGraph(vhMassBins[DT], isXerr) );
	gMpoissonErr.push_back( GetDefaultErrorsGraph(vhMassBins[DT], isXerr) );
	gMpoissonErr[DT]->SetMarkerStyle(1);
	gMpoissonErr[DT]->SetMarkerSize(1);
	gMpoissonErr[DT]->SetLineColor(vModelColor[DT]);
	gMpoissonErr[DT]->SetLineStyle(1);
	*/

	
	cnvAfb->cd();
	pad_mHat->Draw();
	pad_mHat->cd();
	if(!drawKKAfb)
	{
		vhMassBins[Z0]->SetMaximum(1.5*vhMassBins[Z0]->GetMaximum());
		vhMassBins[Z0]->SetMinimum((vhMassBins[Z0]->GetMinimum()<=0) ? 1. : 0.5*vhMassBins[Z0]->GetMinimum());
		vhMassBins[Z0]->GetXaxis()->SetMoreLogLabels();
		vhMassBins[Z0]->GetXaxis()->SetNoExponent();
		vhMassBins[Z0]->Draw();
	}
	else
	{
		vhMassBins[KK]->SetMaximum(1.5*vhMassBins[KK]->GetMaximum());
		vhMassBins[KK]->SetMinimum((vhMassBins[Z0]->GetMinimum()<=0) ? 1. : 0.5*vhMassBins[Z0]->GetMinimum());
		vhMassBins[KK]->GetXaxis()->SetMoreLogLabels();
		vhMassBins[KK]->GetXaxis()->SetNoExponent();
		vhMassBins[KK]->Draw();
		//gMpoissonErr[KK]->GetXaxis()->SetMoreLogLabels();
		//gMpoissonErr[KK]->GetXaxis()->SetNoExponent();
		//gMpoissonErr[KK]->Draw("SAMES");
	
		vhMassBins[Z0]->GetXaxis()->SetMoreLogLabels();
		vhMassBins[Z0]->GetXaxis()->SetNoExponent();
		vhMassBins[Z0]->Draw("SAMES");
		//gMpoissonErr[Z0]->GetXaxis()->SetMoreLogLabels();
		//gMpoissonErr[Z0]->GetXaxis()->SetNoExponent();
		//gMpoissonErr[Z0]->Draw("SAMES");
	}

	vhMassBins[ZP]->GetXaxis()->SetMoreLogLabels();
	vhMassBins[ZP]->GetXaxis()->SetNoExponent();
	vhMassBins[ZP]->Draw("SAMES");
	//gMpoissonErr[ZP]->GetXaxis()->SetMoreLogLabels(); 
	//gMpoissonErr[ZP]->GetXaxis()->SetNoExponent(); 
	//gMpoissonErr[ZP]->Draw("SAMES");
	
	vhMassBins[DT]->GetXaxis()->SetMoreLogLabels();
	vhMassBins[DT]->GetXaxis()->SetNoExponent();
	vhMassBins[DT]->Draw("SAMES");
	//gMpoissonErr[DT]->GetXaxis()->SetMoreLogLabels(); 
	//gMpoissonErr[DT]->GetXaxis()->SetNoExponent(); 
	//gMpoissonErr[DT]->Draw("SAMES");

	cnvAfb->cd();

	pad_Afb->Draw();
	pad_Afb->cd();
	if(!drawKKAfb)
	{
		vhAfbBins[Z0]->GetYaxis()->SetRangeUser(-1.5,+1.5);
		if(drawAfbErrArea) vhAfbBins[Z0]->Draw("E5 Y+");
		else               vhAfbBins[Z0]->Draw("Y+ e1x1");
		vhAfbBins[Z0]->GetXaxis()->SetMoreLogLabels(); 
		vhAfbBins[Z0]->GetXaxis()->SetNoExponent();
	}
	else
	{
		vhAfbBins[KK]->GetYaxis()->SetRangeUser(-1.5,+1.5);
		if(drawAfbErrArea) vhAfbBins[KK]->Draw("E5 Y+");
		else               vhAfbBins[KK]->Draw("e1x1 Y+");
		vhAfbBins[KK]->GetXaxis()->SetMoreLogLabels(); 
		vhAfbBins[KK]->GetXaxis()->SetNoExponent();

		if(drawAfbErrArea) vhAfbBins[Z0]->Draw("E5 Y+ SAMES");
		else               vhAfbBins[Z0]->Draw("Y+ e1x1 SAMES");
		vhAfbBins[Z0]->GetXaxis()->SetMoreLogLabels(); 
		vhAfbBins[Z0]->GetXaxis()->SetNoExponent();
	}

	if(drawAfbErrArea) vhAfbBins[ZP]->Draw("E5 Y+ SAMES");
	else               vhAfbBins[ZP]->Draw("Y+ e1x1 SAMES");
	vhAfbBins[ZP]->GetXaxis()->SetMoreLogLabels(); 
	vhAfbBins[ZP]->GetXaxis()->SetNoExponent();
	
	vhAfbBins[DT]->Draw("Y+ e1x1 SAMES");
	vhAfbBins[DT]->GetXaxis()->SetMoreLogLabels(); 
	vhAfbBins[DT]->GetXaxis()->SetNoExponent(); 
	
	pvtxt_lumi->Draw("SAMES");
	pvtxt_atlas->Draw("SAMES");
	leg_mHat->Draw("SAMES");
	//pad_Afb->RedrawAxis();
	
	cnvAfb->cd();
	pad_mHat->cd();
	pad_mHat->RedrawAxis();
	pad_mHat->Update();
	cnvAfb->Update();
	
	cnvAfb->SaveAs("fitplots/FitAfb.root");
	cnvAfb->SaveAs("fitplots/FitAfb.C");
	cnvAfb->SaveAs("fitplots/FitAfb.eps");
	cnvAfb->SaveAs("fitplots/FitAfb.ps");
	cnvAfb->SaveAs("fitplots/FitAfb.pdf");
	cnvAfb->SaveAs("fitplots/FitAfb.png");
}


void Afb_RooFit_weighted()
{
	setMSGlevel(VISUAL,VISUAL,VISUAL);

	_DEBUG("Afb_RooFit_weighted");
	style();
	colors();
	
	randGen = new TRandom();
	randGen->SetSeed(0); // Note that the machine clock is returned with a precision of 1 second.
						 // If one calls SetSeed(0) within a loop and the loop time is less than 1s,
						 // all generated numbers will be identical!
	
	int nCnvColumns = 6;
	
	
	TCanvas* cnvTmp;
	TVirtualPad* padTmp_cosTheta;
	TVirtualPad* padTmp_imass;
	
	
	//TCanvas* cnv = new TCanvas("fit", "fit", 1024,1280);
	TCanvas* cnv = new TCanvas("fit", "fit", 1200,800);
	cnv->Draw();
	cnv->Divide(nCnvColumns,nMassBins); // const int nMassBins is defined in constants.h
	
	vector<vector<TVirtualPad*> > vPad;
	vector<vector<TVirtualPad*> > vPadBin;
	vector<vector<fitpars> >      vvFitParsResult;
	vector<vector<fitpars> >      vvFitParsResultErr;
	vector<vector<double> >       vvAfbResult;
	vector<vector<double> >       vvAfbError;
	
	vector<TVirtualPad*> vPadTmp;
	vector<double>       vAfbResultTmp;
	vector<double>       vAfbErrorTmp;
	vector<fitpars>      vGuessTmp;
	vector<fitpars>      vResultTmp;
	vector<fitpars>      vResultErrTmp;
	
	vector<TString>      vTitles;
	//vector<TCanvas*>   vCanvases;
	vector<TVirtualPad*> vCanvases;
	vector<TH1D*>        vhMassTmp;
	
	leg = new TLegend(0.006269594,0.03457447,0.9902473,0.4069149,NULL,"brNDC");
	leg->SetFillColor(kWhite);
	
	leg_mHat = new TLegend(0.1557789,0.1818182,0.3806533,0.458042,NULL,"brNDC");
	leg_mHat->SetFillColor(kWhite);
	
	pvtxt = new TPaveText(0.006269594,0.5531915,0.9902473,0.9255319,"brNDC");
	pvtxt->SetFillColor(kWhite);
	
	pvtxt_lumi = new TPaveText(0.1620603,0.458042,0.3140704,0.5716783,"brNDC");
	pvtxt_lumi->SetFillColor(kWhite);
	TString sLumi = (TString)tostring(luminosity,2);
	txt = pvtxt_lumi->AddText( "#intLdt~"+ sLumi +" fb^{-1}" );
	
	pvtxt_atlas = new TPaveText(0.3277592,0.8056995,0.4916388,0.9002591,"brNDC");
	pvtxt_atlas->SetFillColor(0);
	pvtxt_atlas->SetTextFont(42);
	txt = pvtxt_atlas->AddText("#bf{#splitline{#it{ATLAS}}{#scale[0.42]{work in progress}}}");
	
	cnvAfb = new TCanvas("cnvAfb", "cnvAfb", 0,0,1200,800);
	pad_mHat = new TPad("padMhat","",0,0,1,1);
	pad_mHat->SetFillColor(kWhite);
	pad_mHat->SetTicky(0);
	pad_mHat->SetLogy();
	pad_mHat->SetLogx();
	pad_Afb = new TPad("padAfb", "",0,0,1,1);
	pad_Afb->SetTicky(0);
	pad_Afb->SetTickx(1);
	pad_Afb->SetFillStyle(0);
	pad_Afb->SetFrameFillStyle(4000); //will be transparent
	pad_Afb->SetFrameFillColor(0);
	pad_Afb->SetLogx();
	
	//pad_compare = new TPad("padMhat","",0,0,1,1); !!!!!!!!!!!!!!!!!!!!!!!!
	
	int padCounter = 1;
	
	for(int massBin=1 ; massBin<=nMassBins ; massBin++)
	{
		TString sMassBin = (TString)tostring(massBin);
	
		vPadTmp.clear();
		vAfbResultTmp.clear();
		vAfbErrorTmp.clear();
		vhMassTmp.clear();
		vbFitStatus.clear();
		vResultTmp.clear();
		vResultErrTmp.clear();
		
		vPad.push_back(vPadTmp);
		vPadBin.push_back(vPadTmp);
		vvAfbResult.push_back(vAfbResultTmp);
		vvAfbError.push_back(vAfbErrorTmp);
		vvFitParsResult.push_back(vResultTmp);
		vvFitParsResultErr.push_back(vResultErrTmp);
		
		_DEBUG("");
		
		//if(cnvTmp!=NULL) delete cnvTmp;
		cnvTmp = new TCanvas("tmp_"+sMassBin, "", 800,400);
		cnvTmp->Divide(2,1);
		padTmp_cosTheta = cnvTmp->cd(1);
		padTmp_cosTheta->SetPad(0, 0.1, 0.65, 0.9);
		padTmp_imass = cnvTmp->cd(2);
		padTmp_imass->SetPad(0.65, 0.1, 1, 0.9);
		
		_DEBUG("");
		
		// vCanvases.push_back( new TCanvas("tmp", "", 604,400) );
		vCanvases.push_back( cnvTmp->cd(1) );
		vCanvases[massBin-1]->SetName( sMassBin );
		vCanvases[massBin-1]->Divide(2,2);
		
		vvInitialGuess.push_back(vGuessTmp);
		
		_DEBUG("");
		
		for(int mod=Z0 ; mod<=DT ; mod++)
		{
			_INFO("\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~mass bin "+tostring(massBin)+", model "+tostring(mod)+"~~~~~~~~~~~~~~~~~~~~~~~~");
			bool skip = false;
			/////////////////////////////////////////
			init(massBin, mod); /////////////////////
			Int_t N = loop(mod); ////////////////////
			if(N<=minEntries2Fit) skip=true; ////////
			if(!skip) getFit(mod); //////////////////
			/////////////////////////////////////////
			
			vPad[massBin-1].push_back( cnv->cd( padCounter ) );
			padCounter++;
			vPadBin[massBin-1].push_back( vCanvases[massBin-1]->cd(mod+1) );
			
			if(!skip)
			{
				fitpars fp;
				fitpars dfp;
				double Afb;
				double dAfb;
				fp.A0 = A0->getVal();
				fp.A4 = A4->getVal();
				dfp.A0 = A0->getError();
				dfp.A4 = A4->getError();
				vvFitParsResult[massBin-1].push_back( fp );
				vvFitParsResultErr[massBin-1].push_back( dfp );
				Afb  = getAfb(A4->getVal());
				if(doBinomialError) dAfb = getAfbErr(Afb,N);
				else                dAfb = getAfbErr(dfp.A4);
				vvAfbResult[massBin-1].push_back( Afb );
				vvAfbError[massBin-1].push_back( dAfb );
				plot(mod,vPad[massBin-1][mod]);
				plot(mod,vPadBin[massBin-1][mod]);
				_INFO((string)vModelName[mod]+"(N="+tostring(N)+") --> Afb = "+tostring(vvAfbResult[massBin-1][mod])+" +- "+tostring(vvAfbError[massBin-1][mod]));
			}
			else
			{
				vvAfbResult[massBin-1].push_back( -999. );
				vvAfbError[massBin-1].push_back( -999. );
				plotempty(vPad[massBin-1][mod]);
				plotempty(vPadBin[massBin-1][mod]);
				_INFO((string)vModelName[mod]+" --> 0 entries, skipping.");
			}
			
			vhAfbBins[mod]->SetBinContent(massBin,vvAfbResult[massBin-1][mod]);
			vhAfbBins[mod]->SetBinError(massBin,vvAfbError[massBin-1][mod]);
		}
		
		vvbFitStatus.push_back(vbFitStatus);
		
		vTitles.push_back( vhMass[Z0]->GetTitle() );
		
		vPad[massBin-1].push_back( cnv->cd( padCounter ) );
		padCounter++;
		vPad[massBin-1][DT+1]->SetLogy();
		vPad[massBin-1][DT+1]->SetLogx();
		
		Double_t hMin = getYmin(vhMass[Z0]);
		hMin = (getYmin(vhMass[ZP]) < hMin) ? getYmin(vhMass[ZP]) : hMin;
		hMin = (getYmin(vhMass[KK]) < hMin) ? getYmin(vhMass[KK]) : hMin;
		Double_t hMax = vhMass[Z0]->GetMaximum();
		hMax = (vhMass[ZP]->GetMaximum() > hMax) ? vhMass[ZP]->GetMaximum() : hMax;
		hMax = (vhMass[KK]->GetMaximum() > hMax) ? vhMass[KK]->GetMaximum() : hMax;
		
		vhMass[Z0]->SetMinimum(0.5*hMin);
		vhMass[Z0]->SetMaximum(1.5*hMax);
		
		vhMassTmp.push_back( (TH1D*)vhMass[Z0]->Clone("") );
		vhMassTmp.push_back( (TH1D*)vhMass[ZP]->Clone("") );
		vhMassTmp.push_back( (TH1D*)vhMass[KK]->Clone("") );
		vhMassTmp.push_back( (TH1D*)vhMass[DT]->Clone("") );
		
		if(massBin==1) leg->AddEntry((TH1D*)vhMass[Z0]->Clone(""), "SM #gamma/Z^{0} (#it{ATLAS} MC10 rec')", "l");
		vhMassTmp[Z0]->Clone("")->Draw();
		if(massBin==1) leg->AddEntry((TH1D*)vhMass[ZP]->Clone(""), "1000 GeV Z' SSM Template", "l");
		vhMassTmp[ZP]->Clone("")->Draw("SAMES");
		if(massBin==1) leg->AddEntry((TH1D*)vhMass[KK]->Clone(""), "1000 GeV S^{1}/Z_{2} KK Template", "l");
		vhMassTmp[KK]->Clone("")->Draw("SAMES");
		if(massBin==1) leg->AddEntry((TH1D*)vhMass[DT]->Clone(""), "2011 Data", "l");
		vhMassTmp[DT]->Clone("")->Draw("SAMES");
		vPad[massBin-1][DT+1]->RedrawAxis();
		
		
		pvtxt->Clear();
		txt = pvtxt->AddText( vTitles[massBin-1] );
		//pvtxt->InseretTxt( vTitles[massBin-1] );
		vPad[massBin-1].push_back( cnv->cd( padCounter ) );
		padCounter++;
		vPad[massBin-1][DT+2]->cd();
		pvtxt->Clone("")->Draw();
		leg->Draw("SAMES");
		
		
		//------------------------------------------------------------
		// for the binned canvases
		// vCanvases[massBin-1]->cd();
		// vPadBin[massBin-1].push_back( vCanvases[massBin-1]->cd( DT+2 ) );
		// vPadBin[massBin-1][DT+1]->Draw();
		// vPadBin[massBin-1][DT+1]->SetLogy();
		// vPadBin[massBin-1][DT+1]->SetLogx();
		
		vPadBin[massBin-1].push_back( padTmp_imass->cd()/*cnvTmp->cd(2)*/ );
		vPadBin[massBin-1][DT+1]->Draw();
		vPadBin[massBin-1][DT+1]->SetLogy();
		vPadBin[massBin-1][DT+1]->SetLogx();
		vhMassTmp[Z0]->Clone("")->Draw();
		vhMassTmp[ZP]->Clone("")->Draw("SAMES");
		vhMassTmp[KK]->Clone("")->Draw("SAMES");
		vhMassTmp[DT]->Clone("")->Draw("SAMES");
		vPadBin[massBin-1][DT+1]->RedrawAxis();
		
		// vCanvases[massBin-1]->Update();
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".root");
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".C");
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".eps");
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".ps");
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".pdf");
		// vCanvases[massBin-1]->SaveAs("fitplots/FitMassBin_" + sMassBin + ".png");
		
		cnvTmp->Update();
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".root");
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".C");
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".eps");
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".ps");
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".pdf");
		cnvTmp->SaveAs("fitplots/FitMassBin_" + sMassBin + ".png");
		//------------------------------------------------------------
		
		
		
		
		_INFO("completted massBin="+tostring(massBin));
		/////////////
		reset(); ////
		/////////////
		_INFO("reset massBin="+tostring(massBin));
	}
	
	cnv->Update();
	
	cnv->SaveAs("fitplots/FitAllMassBins.root");
	cnv->SaveAs("fitplots/FitAllMassBins.C");
	cnv->SaveAs("fitplots/FitAllMassBins.eps");
	cnv->SaveAs("fitplots/FitAllMassBins.ps");
	cnv->SaveAs("fitplots/FitAllMassBins.pdf");
	cnv->SaveAs("fitplots/FitAllMassBins.png");
	
	drawAfb();
	
	for(unsigned int i=0 ; i<vvInitialGuess.size() ; i++)
	{
		cout << "\nmass bin " << i << ":" << endl; 
		cout << "   Z0:  STATUS=" << vvbFitStatus[i][Z0]
			 << "\t[A0,A4](guess)=[" << vvInitialGuess[i][Z0].A0 << "," << vvInitialGuess[i][Z0].A4
			 // << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][Z0].A0,vvInitialGuess[i][Z0].A4,false)
			 << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][Z0].A4)
			 << "\t-> \tA0(fit)=" << vvFitParsResult[i][Z0].A0 << "+-" << vvFitParsResultErr[i][Z0].A0
			 << "\t-> \tA4(fit)=" << vvFitParsResult[i][Z0].A4 << "+-" << vvFitParsResultErr[i][Z0].A4
			 << ",\tAfb(fit)=" << vvAfbResult[i][Z0] << "+-" << vvAfbError[i][Z0] << endl;
		
		cout << "   ZP:  STATUS=" << vvbFitStatus[i][ZP]
			 << "\t[A0,A4](guess)=[" << vvInitialGuess[i][ZP].A0 << "," << vvInitialGuess[i][ZP].A4
			 // << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][ZP].A0,vvInitialGuess[i][ZP].A4,false)
			 << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][ZP].A4)
			 << "\t-> \tA0(fit)=" << vvFitParsResult[i][ZP].A0 << "+-" << vvFitParsResultErr[i][ZP].A0
			 << "\t-> \tA4(fit)=" << vvFitParsResult[i][ZP].A4 << "+-" << vvFitParsResultErr[i][ZP].A4
			 << ",\tAfb(fit)=" << vvAfbResult[i][ZP] << "+-" << vvAfbError[i][ZP] << endl;
		
		cout << "   KK:  STATUS=" << vvbFitStatus[i][KK]
			 << "\t[A0,A4](guess)=[" << vvInitialGuess[i][KK].A0 << "," << vvInitialGuess[i][KK].A4
			 // << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][KK].A0,vvInitialGuess[i][KK].A4,false)
			 << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][KK].A4)
			 << "\t-> \tA0(fit)=" << vvFitParsResult[i][KK].A0 << "+-" << vvFitParsResultErr[i][KK].A0
			 << "\t-> \tA4(fit)=" << vvFitParsResult[i][KK].A4 << "+-" << vvFitParsResultErr[i][KK].A4
			 << ",\tAfb(fit)=" << vvAfbResult[i][KK] << "+-" << vvAfbError[i][KK] << endl;
			 
		cout << "   DT:  STATUS=" << vvbFitStatus[i][DT]
			 << "\t[A0,A4](guess)=[" << vvInitialGuess[i][DT].A0 << "," << vvInitialGuess[i][DT].A4
			 // << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][DT].A0,vvInitialGuess[i][DT].A4,false)
			 << "]\t->  Afb(guess)=" << getAfb(vvInitialGuess[i][DT].A4)
			 << "\t-> \tA0(fit)=" << vvFitParsResult[i][DT].A0 << "+-" << vvFitParsResultErr[i][DT].A0
			 << "\t-> \tA4(fit)=" << vvFitParsResult[i][DT].A4 << "+-" << vvFitParsResultErr[i][DT].A4
			 << ",\tAfb(fit)=" << vvAfbResult[i][DT] << "+-" << vvAfbError[i][DT] << endl;
	}
	
	printRunConfig();
}


