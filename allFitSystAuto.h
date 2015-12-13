#ifndef ALLFITSYSTAUTO_H
#define ALLFITSYSTAUTO_H

#include "std.h"
#include "type.h"
#include "const.h"

static unsigned int vis = 0;


float maxChi2 = 3.;
TMapTSP2TFILE files;

enum BDTs
{
	BDT0, BDT1, BDT2, BDT3
};
enum SBs
{
	SB0, SB1, SB2, SB3
};

float RooFitSF = -1;
TF1* tf1BDTnorm;
Double_t fBDT(Double_t *x, Double_t *par)
{
	return tf1BDTnorm->EvalPar(x,par)*RooFitSF;
}
TF1* tf1SBnorm;
Double_t fSB(Double_t *x, Double_t *par)
{
	return tf1SBnorm->EvalPar(x,par)*RooFitSF;
}

float BDTscale = -1;
TF1* fBDTnorm;
Double_t funcBDTscaled(Double_t *x, Double_t *par)
{
	return fBDTnorm->EvalPar(x,par)*BDTscale;
}
float SBscale = -1;
TF1* fSBnorm;
Double_t funcSBscaled(Double_t *x, Double_t *par)
{
	return fSBnorm->EvalPar(x,par)*SBscale;
}

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

TString snominal = "1450-1690-1870-2110.neg090.pos1000";
void setFiles()
{
	TString name = "";
	TString prefix = "BDTfitResults.";
	
	name="1450-1690-1870-2110.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // Nominal
	
	name="1450-1690-1870-2110.neg089.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg088.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg087.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg086.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg085.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg084.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg083.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg082.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg081.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	name="1450-1690-1870-2110.neg080.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // BDT cutoff
	                                                                                
	name="1420-1690-1870-2140.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1390-1690-1870-2170.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1360-1690-1870-2200.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1330-1690-1870-2230.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1300-1690-1870-2260.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1270-1690-1870-2290.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1540-1690-1870-2020.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1510-1690-1870-2050.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	name="1480-1690-1870-2080.neg090.pos1000"; files.insert(make_pair(name,new TFile(prefix+name+".root","READ"))); // SB range
	
	_INF(vis,"");
}
TFile* getNominal() { return files[snominal]; }

class Estimate
{
	public:
		Estimate() {};
		~Estimate()
		{
			reset();
		};
	public:
		vector<float> bdtRooFitSF, sbRooFitSF;
		vector<float> bdtIntegralLeft, bdtIntegralRight;
		float         nSB, dnSBstat;
		vector<float> nSR, dnSRstat;
		vector<float> chi2BDT, chi2SB;
		vector<float> xsrmin, xsrmax;
		vector<vector<float> > vnSR, vdnSRstat;
		float minBDTabs, maxBDTabs, BDTbinWidth;
		vector<TF1*> vtf1BDT, vtf1SB;
		vector<TNamed*> vtf1BDTstr, vtf1SBstr;
	public:
		void reset()
		{
			bdtIntegralLeft.clear();
			bdtIntegralRight.clear();
			bdtRooFitSF.clear();
			nSR.clear(); dnSRstat.clear();
			chi2BDT.clear(); chi2SB.clear();
			nSB = -1; dnSBstat = -1;
		};
		void setBDTx(float xmin, float xmax, float width) { minBDTabs=xmin; maxBDTabs=xmax; BDTbinWidth=width; };
		void addBDTchi2(float x, TF1* f, TNamed* str) { chi2BDT.push_back(x); vtf1BDT.push_back(f); vtf1BDTstr.push_back(str); };
		void addSBchi2(float x,  TF1* f, TNamed* str) { chi2SB.push_back(x);  vtf1SB.push_back(f);  vtf1SBstr.push_back(str);  };
		void setSB(float x, float dx) { nSB=x; dnSBstat=dx; };
		void addSR(float x, float dxStat)
		{
			if(nSB<0) _FAT("Cannot call addSR() before calling setSB() first");
			nSR.push_back(x);
			dnSRstat.push_back(dxStat);
		};
		void addSRs(TVectorF* xSRmin, TVectorF* xSRmax,
					TVectorF* nSR00,  TVectorF* nSR01,  TVectorF* nSR02,  TVectorF* nSR03,
					TVectorF* dnSR00, TVectorF* dnSR01, TVectorF* dnSR02, TVectorF* dnSR03)
		{
			Int_t NSRs = xSRmin->GetNoElements();
			
			_INF(vis,"NSRs="<<NSRs);
			
			for(unsigned int variation=0 ; variation<4 ; ++variation)
			{
				vector<float> tmp;
				vnSR.push_back(tmp);
				vdnSRstat.push_back(tmp);
				for(Int_t sr=0 ; sr<NSRs ; ++sr)
				{
					if     (variation==0) { xsrmin.push_back(((*xSRmin))[sr]); xsrmax.push_back(((*xSRmax))[sr]); } // NEEDED ONLY ONCE
					if     (variation==0) { vnSR[variation].push_back(((*nSR00))[sr]); vdnSRstat[variation].push_back(((*dnSR00))[sr]); }
					else if(variation==1) { vnSR[variation].push_back(((*nSR01))[sr]); vdnSRstat[variation].push_back(((*dnSR01))[sr]); }
					else if(variation==2) { vnSR[variation].push_back(((*nSR02))[sr]); vdnSRstat[variation].push_back(((*dnSR02))[sr]); }
					else if(variation==3) { vnSR[variation].push_back(((*nSR03))[sr]); vdnSRstat[variation].push_back(((*dnSR03))[sr]); }
					else _FAT("OUT OF RNGE");
				}
			}
			
			_INF(vis,"");
		};
		void addBDTintegralLeft(float x)
		{
			if(nSB<0) _FAT("Cannot call addBDTintegralLeft() before calling setSB() first");
			bdtIntegralLeft.push_back(x);
			_INF(vis,"");
		};
		float getBDTintegralLeft(unsigned int entry)
		{
			if(entry>=bdtIntegralLeft.size()) _FAT("entry="<<entry<<" is out of vector size!");
			return bdtIntegralLeft[entry];
			_INF(vis,"");
		}
		void addBDTintegralRight(float x)
		{
			if(nSB<0) _FAT("Cannot call addBDTintegralRight() before calling setSB() first");
			bdtIntegralRight.push_back(x);
			_INF(vis,"");
		};
		float getBDTintegralRight(unsigned int entry)
		{
			if(entry>=bdtIntegralRight.size()) _FAT("entry="<<entry<<" is out of vector size!");
			return bdtIntegralRight[entry];
			_INF(vis,"");
		}
		
		void addBDTsf(float x)
		{
			if(nSB<0) _FAT("Cannot call addBDTsf() before calling setSB() first");
			bdtRooFitSF.push_back(x);
			_INF(vis,"");
		};
		float getBDTRooFitSF(unsigned int entry)
		{
			if(entry>=bdtRooFitSF.size()) _FAT("entry="<<entry<<" is out of vector size!");
			return bdtRooFitSF[entry];
			_INF(vis,"");
		}
		void addSBsf(float x)
		{
			if(nSB<0) _FAT("Cannot call addSBsf() before calling setSB() first");
			sbRooFitSF.push_back(x);
			_INF(vis,"");
		};
		float getSBRooFitSF(unsigned int entry)
		{
			if(entry>=sbRooFitSF.size()) _FAT("entry="<<entry<<" is out of vector size!");
			return sbRooFitSF[entry];
			_INF(vis,"");
		}
		string getSummary()
		{
			string summary = "";
			summary += "========================================================================================\n";
			summary += "BDT RooFitSF    ";
			for(unsigned int i=0 ; i<bdtRooFitSF.size() ; ++i) { summary+=str(bdtRooFitSF[i],3)+" | "; }
			summary += "\n";
			summary += "Count: nSB      "+str(nSB,0)+"\n";
			summary += "Scaled nSR      ";
			for(unsigned int i=0 ; i<nSR.size() ; ++i) { summary+=str(nSR[i],3)+"+-"+str(dnSRstat[i],3)+" | "; }
			summary += "\n";
			summary += "========================================================================================\n";
			_INF(vis,"");
			return summary;
		};
};
typedef map<TString, Estimate*> TMapTSP2Estimate;


Estimate* setEstimate(TFile* f, float xBDTcut=optBDTcut)
{
	Estimate* est = new Estimate();
	
	//// some header details
	TVectorF* AbsBDTmin     = (TVectorF*)f->Get("absBDTmin");
	TVectorF* AbsBDTmax     = (TVectorF*)f->Get("absBDTmax");
	TVectorF* AbsBDTmaxUser = (TVectorF*)f->Get("absBDTmaxUser");
	TVectorF* NBDTbins      = (TVectorF*)f->Get("nBDTbins");
	TVectorF* NSB           = (TVectorF*)f->Get("nSB");
	TVectorF* NSR           = (TVectorF*)f->Get("nSR");
	TVectorF* chi2BDT       = (TVectorF*)f->Get("chi2BDT");
	TVectorF* chi2SB        = (TVectorF*)f->Get("chi2SB");
	float minBDTabs   = ((*AbsBDTmin))[0];
	float maxBDTabs   = ((*AbsBDTmax))[0];
	float maxBDTusr   = ((*AbsBDTmaxUser))[0];
	int   BDTbins     = (int)((*NBDTbins))[0];
	float BDTbinWidth = ((maxBDTabs-minBDTabs)/BDTbins);
	TVectorF* xSRmin = (TVectorF*)f->Get("xSRmin"); TVectorF* xSRmax = (TVectorF*)f->Get("xSRmax");
	TVectorF* nSR00  = (TVectorF*)f->Get("nSR00");  TVectorF* dnSR00 = (TVectorF*)f->Get("dnSR00");
	TVectorF* nSR01  = (TVectorF*)f->Get("nSR01");  TVectorF* dnSR01 = (TVectorF*)f->Get("dnSR01");
	TVectorF* nSR02  = (TVectorF*)f->Get("nSR02");  TVectorF* dnSR02 = (TVectorF*)f->Get("dnSR02");
	TVectorF* nSR03  = (TVectorF*)f->Get("nSR03");  TVectorF* dnSR03 = (TVectorF*)f->Get("dnSR03");
	
	_INF(vis,"");
	
	//// For the optimization we have to do the SB fit in the full range, CR0+CR1, not just in CR0 !
	if(maxBDTusr!=maxBDTabs) _FAT("Running with final-estimation mode: maxBDTusr="<<maxBDTusr);
	
	//// set BDT Chi2
	est->addBDTchi2( ((*chi2BDT))[0], (TF1*)f->Get("fBDT0.0")->Clone((TString)"fBDT0.0."+(TString)f->GetName()), (TNamed*)f->Get("fBDT0str.0")->Clone((TString)"fBDT0str.0."+(TString)f->GetName()) );
	est->addBDTchi2( ((*chi2BDT))[1], (TF1*)f->Get("fBDT0.1")->Clone((TString)"fBDT0.1."+(TString)f->GetName()), (TNamed*)f->Get("fBDT0str.1")->Clone((TString)"fBDT0str.1."+(TString)f->GetName()) );
	est->addBDTchi2( ((*chi2BDT))[2], (TF1*)f->Get("fBDT0.2")->Clone((TString)"fBDT0.2."+(TString)f->GetName()), (TNamed*)f->Get("fBDT0str.2")->Clone((TString)"fBDT0str.2."+(TString)f->GetName()) );
	est->addBDTchi2( ((*chi2BDT))[3], (TF1*)f->Get("fBDT0.3")->Clone((TString)"fBDT0.3."+(TString)f->GetName()), (TNamed*)f->Get("fBDT0str.3")->Clone((TString)"fBDT0str.3."+(TString)f->GetName()) );
	
	//// set SB Chi2
	est->addSBchi2( ((*chi2SB))[0], (TF1*)f->Get("fSB0.0")->Clone((TString)"fSB0.0."+(TString)f->GetName()), (TNamed*)f->Get("fSB0str.0")->Clone((TString)"fSB0str.0."+(TString)f->GetName()) );
	est->addSBchi2( ((*chi2SB))[1], (TF1*)f->Get("fSB0.1")->Clone((TString)"fSB0.1."+(TString)f->GetName()), (TNamed*)f->Get("fSB0str.1")->Clone((TString)"fSB0str.1."+(TString)f->GetName()) );
	est->addSBchi2( ((*chi2SB))[2], (TF1*)f->Get("fSB0.2")->Clone((TString)"fSB0.2."+(TString)f->GetName()), (TNamed*)f->Get("fSB0str.2")->Clone((TString)"fSB0str.2."+(TString)f->GetName()) );
	est->addSBchi2( ((*chi2SB))[3], (TF1*)f->Get("fSB0.3")->Clone((TString)"fSB0.3."+(TString)f->GetName()), (TNamed*)f->Get("fSB0str.3")->Clone((TString)"fSB0str.3."+(TString)f->GetName()) );
	
	//// set the BDT boundaries for this set of functions
	est->setBDTx(minBDTabs,maxBDTabs,BDTbinWidth);
	
	//// SB total events and stat error
	est->setSB( ((*NSB))[0],sqrt(((*NSB))[0]) );
	
	//// SR estimations [0]=nominal, no need to remember the functions
	est->addSR( ((*NSR))[0], (((*NSR))[0]/est->nSB)*est->dnSBstat );
	est->addSR( ((*NSR))[1], (((*NSR))[1]/est->nSB)*est->dnSBstat );
	est->addSR( ((*NSR))[2], (((*NSR))[2]/est->nSB)*est->dnSBstat );
	est->addSR( ((*NSR))[3], (((*NSR))[3]/est->nSB)*est->dnSBstat );
	
	_INF(vis,"");
	
	est->addSRs(xSRmin,xSRmax,nSR00,nSR01,nSR02,nSR03,dnSR00,dnSR01,dnSR02,dnSR03);
	
	_INF(vis,"");
	
	tf1BDTnorm = (TF1*)f->Get("fBDT0.0"); est->addBDTintegralLeft(tf1BDTnorm->Integral(minBDTabs,xBDTcut));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.1"); est->addBDTintegralLeft(tf1BDTnorm->Integral(minBDTabs,xBDTcut));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.2"); est->addBDTintegralLeft(tf1BDTnorm->Integral(minBDTabs,xBDTcut));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.3"); est->addBDTintegralLeft(tf1BDTnorm->Integral(minBDTabs,xBDTcut));
	
	_INF(vis,"");
	
	tf1BDTnorm = (TF1*)f->Get("fBDT0.0"); est->addBDTintegralRight(tf1BDTnorm->Integral(xBDTcut,maxBDTabs));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.1"); est->addBDTintegralRight(tf1BDTnorm->Integral(xBDTcut,maxBDTabs));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.2"); est->addBDTintegralRight(tf1BDTnorm->Integral(xBDTcut,maxBDTabs));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.3"); est->addBDTintegralRight(tf1BDTnorm->Integral(xBDTcut,maxBDTabs));
	
	_INF(vis,"");
	
	//// BDT functions [0]=nominal, need to remember the RooFitSF per function !
	tf1BDTnorm = (TF1*)f->Get("fBDT0.0"); est->addBDTsf(est->nSB/(tf1BDTnorm->Integral(minBDTabs,maxBDTabs)/BDTbinWidth));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.1"); est->addBDTsf(est->nSB/(tf1BDTnorm->Integral(minBDTabs,maxBDTabs)/BDTbinWidth));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.2"); est->addBDTsf(est->nSB/(tf1BDTnorm->Integral(minBDTabs,maxBDTabs)/BDTbinWidth));
	tf1BDTnorm = (TF1*)f->Get("fBDT0.3"); est->addBDTsf(est->nSB/(tf1BDTnorm->Integral(minBDTabs,maxBDTabs)/BDTbinWidth));
	
	_INF(vis,"");
	
	//// SB functions [0]=nominal, need to remember the RooFitSF per function !
	tf1SBnorm = (TF1*)f->Get("fSB0.0"); est->addSBsf(est->nSB/(tf1SBnorm->Integral(mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob)/mBinSize));
	tf1SBnorm = (TF1*)f->Get("fSB0.1"); est->addSBsf(est->nSB/(tf1SBnorm->Integral(mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob)/mBinSize));
	tf1SBnorm = (TF1*)f->Get("fSB0.2"); est->addSBsf(est->nSB/(tf1SBnorm->Integral(mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob)/mBinSize));
	tf1SBnorm = (TF1*)f->Get("fSB0.3"); est->addSBsf(est->nSB/(tf1SBnorm->Integral(mSideBandLeftLowerMeVGlob,mSideBandRightUpperMeVGlob)/mBinSize));
	
	_INF(vis,"");
	
	return est;
}


struct Extrapolation
{
	float nSB0, dnSB0stat;
	float chi2BDT;
	float chi2SB;
	float BDTcut;
	float correction;
	float f01;
	float LeftInt;
	float RightInt;
	float nSB, dnSBstat;
	float nSR, dnSRstat;
	float nSR0, dnSR0stat;
	float nSB1, dnSB1stat;
	
	float xSRmin, xSRmax;
	float VARnSR0, VARdnSR0stat;
	float VARnSR, VARdnSRstat;
};
// Extrapolation* extrapulate(unsigned int iBDT, unsigned int iSB, float xBDTcut, TFile* f, Estimate* est, bool print=false)
Extrapolation* extrapulate(unsigned int iBDT, unsigned int iSB, unsigned int iSR, float xBDTcut, TFile* f, Estimate* est, bool print=false)
{
	Extrapolation* CR1 = new Extrapolation;
	
	stringstream strm;
	string str;
	strm << setprecision(0) << fixed << iBDT;
	strm >> str;
 	TString sIBDT = str;
	
	float integralLeft  = est->getBDTintegralLeft(iBDT);
	float integralRight = est->getBDTintegralRight(iBDT);
	
	RooFitSF = est->getBDTRooFitSF(iBDT);
	tf1BDTnorm = (TF1*)f->Get("fBDT0."+tstr(iBDT,0)); // this global TF1 is part of fBDT function pointer
	TF1* tf1BDT = new TF1("fBDT",fBDT,est->minBDTabs,est->minBDTabs);
	float nSB0 = tf1BDT->Integral(est->minBDTabs,xBDTcut)/est->BDTbinWidth;
	float nSB1 = tf1BDT->Integral(xBDTcut,est->maxBDTabs)/est->BDTbinWidth;
	if(fabs((nSB0+nSB1)-est->nSB)/est->nSB*100.>1e-2) _FAT("Non closure! nSB0="<<nSB0<<", nSB1="<<nSB1<<", nSB="<<est->nSB);
	float f01  = nSB1/nSB0;
	float correction = nSB0/(nSB0+nSB1);
	float dnSB0stat = sqrt(nSB0);
	
	_INF(vis,"");
	
	CR1->nSB0         = nSB0;
	CR1->dnSB0stat    = dnSB0stat;
	CR1->chi2BDT      = est->chi2BDT[iBDT];
	CR1->chi2SB       = est->chi2SB[iSB];
	CR1->BDTcut       = xBDTcut;
	CR1->correction   = correction;
	CR1->f01          = f01;
	CR1->LeftInt      = integralLeft;
	CR1->RightInt     = integralRight;
	CR1->nSB          = est->nSB*correction*f01;
	CR1->dnSBstat     = est->dnSBstat*correction*f01;
	CR1->nSR0         = est->nSR[iSB];
	CR1->dnSR0stat    = est->dnSRstat[iSB];
	CR1->nSR          = est->nSR[iSB]*correction*f01;
	CR1->dnSRstat     = est->dnSRstat[iSB]*correction*f01;
	CR1->xSRmin       = est->xsrmin[iSR];
	CR1->xSRmax       = est->xsrmax[iSR];
	CR1->VARnSR0      = est->vnSR[iSB][iSR];
	CR1->VARdnSR0stat = est->vdnSRstat[iSB][iSR];
	CR1->VARnSR       = est->vnSR[iSB][iSR]*correction*f01;
	CR1->VARdnSRstat  = est->vdnSRstat[iSB][iSR]*correction*f01;
	
	_INF(vis,"");
	
	if(print)
	{
		cout << "<<<<<<<<< Extrapolation iBDT=" << iBDT << ", iSB=" << iSB << ", iSR=" << iSR << " >>>>>>>>>" << endl;
		cout << "   CR1->nSB0       = " << CR1->nSB0       << endl;
		cout << "   CR1->dnSB0stat  = " << CR1->dnSB0stat  << endl;
		cout << "   CR1->chi2BDT    = " << CR1->chi2BDT    << endl;
		cout << "   CR1->chi2SB     = " << CR1->chi2SB     << endl;
		cout << "   CR1->BDTcut     = " << CR1->BDTcut     << endl;
		cout << "   CR1->correction = " << CR1->correction << endl;
		cout << "   CR1->f01        = " << CR1->f01        << endl;
		cout << "   CR1->nSB        = " << CR1->nSB        << endl;
		cout << "   CR1->dnSBstat   = " << CR1->dnSBstat   << endl;
		cout << "   CR1->nSR        = " << CR1->nSR        << endl;
		cout << "   CR1->dnSRstat   = " << CR1->dnSRstat   << endl;
		cout << "   -----------------"                         << endl;
		cout << "   CR1->xSRmin       = " << CR1->xSRmin       << endl;
		cout << "   CR1->xSRmax       = " << CR1->xSRmax       << endl;
		cout << "   CR1->VARnSR0      = " << CR1->VARnSR0      << endl;
		cout << "   CR1->VARdnSR0stat = " << CR1->VARdnSR0stat << endl;
		cout << "   CR1->VARnSR       = " << CR1->VARnSR       << endl;
		cout << "   CR1->VARdnSRstat  = " << CR1->VARdnSRstat  << endl;
		cout << "<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>" << endl;
	}
	
	delete tf1BDT;
	
	_INF(vis,"");
	
	return CR1;
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
	_INF(vis,"");
	return maxd;
}

struct Result
{
	float correction;
	float nSB0, dnSB0stat;
	float nSB1, dnSB1stat;
	float f01, LeftInt, RightInt, df01shape,  df01range,  df01cutoff;
	float nSR0, dnSR0stat, dnSR0shape, dnSR0range, dnSR0cutoff;
	float nSR1, dnSR1stat, dnSR1shape, dnSR1range, dnSR1cutoff, dnSR1quad;
	//
	float VARnSR0, VARdnSR0stat, VARdnSR0shape, VARdnSR0range, VARdnSR0cutoff;
	float VARnSR1, VARdnSR1stat, VARdnSR1shape, VARdnSR1range, VARdnSR1cutoff, VARdnSR1quad;
	float xSRmin, xSRmax;
};
Result* getExtrapolationResult(float xBDTcut, unsigned int iSR, bool print=false, bool write=false)
{
	ofstream ftxt;
	if(write) ftxt.open("extrapolation_result.txt");
	
	vector<float> vf01_BDTcutoff, vf01_SBrange;
	vector<float> vnSR0_BDTcutoff, vnSR0_SBrange;
	vector<float> VARvnSR0_BDTcutoff, VARvnSR0_SBrange;
	float correction       = -1;
	float chi2BDT          = -1;
	float chi2SB           = -1;
	float nSB0             = -1;
	float nSR0             = -1;
	float nSR1             = -1;
	float nSB1             = -1;
	float f01              = -1;
	float LeftInt          = -1;
	float RightInt         = -1;
	float dnSB0stat        = -1;
	float dnSR0stat        = -1;
	float dnSR1stat        = -1;
	float dnSB1stat        = -1;
	float df01_BDTshape    = -1;
	float dnSR0_SBshape    = -1;
	float VARnSR0          = -1;
	float VARnSR1          = -1;
	float VARdnSR0stat     = -1;
	float VARdnSR1stat     = -1;
	float VARdnSR0_SBshape = -1;
	float xSRmin           = -1;
	float xSRmax           = -1;
	float wSR              = -1;
	
	float nominal_correction       = -1;
	float nominal_chi2BDT          = -1;
	float nominal_chi2SB           = -1;
	float nominal_nSB0             = -1;
	float nominal_nSR0             = -1;
	float nominal_nSR1             = -1;
	float nominal_nSB1             = -1;
	float nominal_f01              = -1;
	float nominal_LeftInt          = -1;
	float nominal_RightInt         = -1;
	float nominal_dnSB0stat        = -1;
	float nominal_dnSR0stat        = -1;
	float nominal_dnSR1stat        = -1;
	float nominal_dnSB1stat        = -1;
	float nominal_df01_BDTshape    = -1;
	float nominal_dnSR0_SBshape    = -1;
	float nominal_VARnSR0          = -1;
	float nominal_VARnSR1          = -1;
	float nominal_VARdnSR0stat     = -1;
	float nominal_VARdnSR1stat     = -1;
	float nominal_VARdnSR0_SBshape = -1;
	float nominal_xSRmin           = -1;
	float nominal_xSRmax           = -1;
	float nominal_wSR              = -1;
	
	TMapTSP2Estimate estimates;
	for(TMapTSP2TFILE::iterator it=files.begin() ; it!=files.end() ; ++it)
	{
		estimates.insert(make_pair(it->first,setEstimate(it->second,xBDTcut)));
		string summary = estimates[it->first]->getSummary();
		if(print) cout << summary << endl;
		
		//// get the nominal extrapolation 
		Extrapolation* nominal_fits = extrapulate(BDT0,SB0,iSR,xBDTcut,it->second,estimates[it->first],print);
		
		//// Test if the nominal chi2's are bad
		if(nominal_fits->chi2BDT>maxChi2) { _WRN(1,"nominal_fits->chi2BDT="<<nominal_fits->chi2BDT<<" --> skipping!"); continue; }
		if(nominal_fits->chi2SB>maxChi2)  { _WRN(1,"nominal_fits->chi2SB="<<nominal_fits->chi2SB<<" --> skipping!");   continue; }
		
		if(write && it->first==snominal) ftxt << "Nominals: R=" << nominal_fits->f01 << ", Nb(x0)=" << nominal_fits->nSR0 << ", correction=" << nominal_fits->correction << endl;
		if(write && it->first==snominal) ftxt << "------------------------------------------------------" << endl;
		
		//// Keep the nominal values to calculate the distances
		if(it->first.Contains("neg090"))
		{
			if(write) ftxt << "[" << it->first << "] BDT range is fixed and the SB range is shifted (using only nominal BDT function): R=" << tstr(nominal_fits->f01,7) << ", Left=" << tstr(nominal_fits->LeftInt,5) << ", Right=" << tstr(nominal_fits->RightInt,5) << endl;
			
			//// when the BDT range is fixed and the SB range is shifted
			vf01_SBrange.push_back(nominal_fits->f01);
			vnSR0_SBrange.push_back(nominal_fits->nSR0);
			VARvnSR0_SBrange.push_back(nominal_fits->VARnSR0);
		}
		else
		{
			if(write) ftxt << "[" << it->first << "] SB range is fixed and the BDT range is shifted (using only nominal SB function): R=" << tstr(nominal_fits->f01,7) << ", Left=" << tstr(nominal_fits->LeftInt,5) << ", Right=" << tstr(nominal_fits->RightInt,5) << endl;
			
			//// when the SB range is fixed and the BDT range is shifted
			vf01_BDTcutoff.push_back(nominal_fits->f01);
			vnSR0_BDTcutoff.push_back(nominal_fits->nSR0);
			VARvnSR0_BDTcutoff.push_back(nominal_fits->VARnSR0);
		}

		
		///////////////////////////////////
		//// below code is relevant
		//// only for nominal_fits
		if(it->first!=snominal) continue;
		///////////////////////////////////
		
		nSB0         = nominal_fits->nSB0;
		dnSB0stat    = nominal_fits->dnSB0stat;
		chi2BDT      = nominal_fits->chi2BDT; if(0) cout << "chi2BDT=" << chi2BDT << endl;
		chi2SB       = nominal_fits->chi2SB;  if(0) cout << "chi2SB=" << chi2SB << endl;
		f01          = nominal_fits->f01;
		nSR0         = nominal_fits->nSR0;
		dnSR0stat    = nominal_fits->dnSR0stat;
		nSR1         = nominal_fits->nSR;
		dnSR1stat    = nominal_fits->dnSRstat;
		nSB1         = nominal_fits->nSB;
		dnSB1stat    = nominal_fits->dnSBstat;
		correction   = nominal_fits->correction;
		VARnSR0      = nominal_fits->VARnSR0;
		VARnSR1      = nominal_fits->VARnSR;
		VARdnSR0stat = nominal_fits->VARdnSR0stat;
		VARdnSR1stat = nominal_fits->VARdnSRstat;
		xSRmin       = nominal_fits->xSRmin;
		xSRmax       = nominal_fits->xSRmax;
		wSR          = nominal_fits->xSRmax-nominal_fits->xSRmin;

	
		//// get the variation extrapolations
		Extrapolation* variant_BDTfit_A = extrapulate(BDT1,SB0,iSR,xBDTcut,it->second,estimates[it->first],print);
		Extrapolation* variant_BDTfit_B = extrapulate(BDT2,SB0,iSR,xBDTcut,it->second,estimates[it->first],print);
		Extrapolation* variant_BDTfit_C = extrapulate(BDT3,SB0,iSR,xBDTcut,it->second,estimates[it->first],print);
		Extrapolation* variant_SBfit_A  = extrapulate(BDT0,SB1,iSR,xBDTcut,it->second,estimates[it->first],print);
		Extrapolation* variant_SBfit_B  = extrapulate(BDT0,SB2,iSR,xBDTcut,it->second,estimates[it->first],print);
		Extrapolation* variant_SBfit_C  = extrapulate(BDT0,SB3,iSR,xBDTcut,it->second,estimates[it->first],print);

		
		//// Test if the variation chi2's are bad
		bool BDTfitAflag = 1;
		bool BDTfitBflag = 1;
		bool BDTfitCflag = 1;
		bool SBfitAflag  = 1;
		bool SBfitBflag  = 1;
		bool SBfitCflag  = 1;
		if(variant_BDTfit_A->chi2BDT>maxChi2) { _WRN(1,"variant_BDTfit_A->chi2BDT="<<variant_BDTfit_A->chi2BDT); BDTfitAflag = 0; }
		if(variant_BDTfit_B->chi2BDT>maxChi2) { _WRN(1,"variant_BDTfit_B->chi2BDT="<<variant_BDTfit_B->chi2BDT); BDTfitBflag = 0; }
		if(variant_BDTfit_C->chi2BDT>maxChi2) { _WRN(1,"variant_BDTfit_C->chi2BDT="<<variant_BDTfit_C->chi2BDT); BDTfitCflag = 0; }
		if(variant_SBfit_A->chi2SB>maxChi2)   { _WRN(1,"variant_SBfit_A->chi2SB="<<variant_SBfit_A->chi2SB);     SBfitAflag  = 0; }
		if(variant_SBfit_B->chi2SB>maxChi2)   { _WRN(1,"variant_SBfit_B->chi2SB="<<variant_SBfit_B->chi2SB);     SBfitBflag  = 0; }
		if(variant_SBfit_C->chi2SB>maxChi2)   { _WRN(1,"variant_SBfit_C->chi2SB="<<variant_SBfit_C->chi2SB);     SBfitCflag  = 0; }
		

		
		//// BDT shape fit choice systematic
		float df01A = (BDTfitAflag) ? fabs(nominal_fits->f01-variant_BDTfit_A->f01) : 0.;
		float df01B = (BDTfitBflag) ? fabs(nominal_fits->f01-variant_BDTfit_B->f01) : 0.;
		float df01C = (BDTfitCflag) ? fabs(nominal_fits->f01-variant_BDTfit_C->f01) : 0.;
		float df01  = (df01A>df01B) ? df01A : df01B;
		df01 = (df01C>df01) ? df01C : df01;
		df01_BDTshape = df01;
		// cout << "df01=" << df01 << " (" << df01/nominal_fits->f01*100 << "\%)" << endl;

		if(write) ftxt << "[" << it->first << "] SB and BDT ranges are fixed and the BDT function is varied (a): R=" << tstr(variant_BDTfit_A->f01,7) << ", Left=" << tstr(variant_BDTfit_A->LeftInt,5) << ", Right=" << tstr(variant_BDTfit_A->RightInt,5) << endl;
		if(write) ftxt << "[" << it->first << "] SB and BDT ranges are fixed and the BDT function is varied (b): R=" << tstr(variant_BDTfit_B->f01,7) << ", Left=" << tstr(variant_BDTfit_B->LeftInt,5) << ", Right=" << tstr(variant_BDTfit_B->RightInt,5) << endl;
		if(write) ftxt << "[" << it->first << "] SB and BDT ranges are fixed and the BDT function is varied (c): R=" << tstr(variant_BDTfit_C->f01,7) << ", Left=" << tstr(variant_BDTfit_C->LeftInt,5) << ", Right=" << tstr(variant_BDTfit_C->RightInt,5) << endl;

		
		//// SR shape fit choice systematic
		float dnSR0A = (SBfitAflag) ? fabs(nominal_fits->nSR0-variant_SBfit_A->nSR0) : 0.;
		float dnSR0B = (SBfitBflag) ? fabs(nominal_fits->nSR0-variant_SBfit_B->nSR0) : 0.;
		float dnSR0C = (SBfitCflag) ? fabs(nominal_fits->nSR0-variant_SBfit_C->nSR0) : 0.;
		float dnSR0 = 0.;
		dnSR0 = (dnSR0A>dnSR0 && SBfitAflag) ? dnSR0A : dnSR0;
		dnSR0 = (dnSR0B>dnSR0 && SBfitBflag) ? dnSR0B : dnSR0;
		dnSR0 = (dnSR0C>dnSR0 && SBfitCflag) ? dnSR0C : dnSR0;
		dnSR0_SBshape = dnSR0;
		// cout << "dnSR0=" << dnSR0 << " (" << dnSR0/nominal_fits->nSR0*100 << "\%)" << endl;

		
		//// SRi shape fit choice systematic
		float VARdnSR0A = (SBfitAflag) ? fabs(nominal_fits->VARnSR0-variant_SBfit_A->VARnSR0) : 0.;
		float VARdnSR0B = (SBfitBflag) ? fabs(nominal_fits->VARnSR0-variant_SBfit_B->VARnSR0) : 0.;
		float VARdnSR0C = (SBfitCflag) ? fabs(nominal_fits->VARnSR0-variant_SBfit_C->VARnSR0) : 0.;
		float VARdnSR0 = 0.;
		VARdnSR0 = (VARdnSR0A>VARdnSR0 && SBfitAflag) ? VARdnSR0A : VARdnSR0;
		VARdnSR0 = (VARdnSR0B>VARdnSR0 && SBfitBflag) ? VARdnSR0B : VARdnSR0;
		VARdnSR0 = (VARdnSR0C>VARdnSR0 && SBfitCflag) ? VARdnSR0C : VARdnSR0;
		VARdnSR0_SBshape = VARdnSR0;
		// cout << "VARdnSR0=" << VARdnSR0 << " (" << VARdnSR0/nominal_fits->VARnSR0*100 << "\%)" << endl;
		
		delete variant_BDTfit_A;
		delete variant_BDTfit_B;
		delete variant_BDTfit_C;
		delete variant_SBfit_A;
		delete variant_SBfit_B;
		delete variant_SBfit_C;
		
		if(it->first==snominal) 
		{
			nominal_correction       = correction      ;
			nominal_chi2BDT          = chi2BDT         ;
			nominal_chi2SB           = chi2SB          ;
			nominal_nSB0             = nSB0            ;
			nominal_nSR0             = nSR0            ;
			nominal_nSR1             = nSR1            ;
			nominal_nSB1             = nSB1            ;
			nominal_f01              = f01             ;
			nominal_dnSB0stat        = dnSB0stat       ;
			nominal_dnSR0stat        = dnSR0stat       ;
			nominal_dnSR1stat        = dnSR1stat       ;
			nominal_dnSB1stat        = dnSB1stat       ;
			nominal_df01_BDTshape    = df01_BDTshape   ;
			nominal_dnSR0_SBshape    = dnSR0_SBshape   ;
			nominal_VARnSR0          = VARnSR0         ;
			nominal_VARnSR1          = VARnSR1         ;
			nominal_VARdnSR0stat     = VARdnSR0stat    ;
			nominal_VARdnSR1stat     = VARdnSR1stat    ;
			nominal_VARdnSR0_SBshape = VARdnSR0_SBshape;
			nominal_xSRmin           = xSRmin          ;
			nominal_xSRmax           = xSRmax          ;
			nominal_wSR              = wSR             ;
		}
	}
	
	float df01_BDTcutoff     = maxdistance(vf01_BDTcutoff, f01);         vf01_BDTcutoff.clear();
	float df01_SBrange       = maxdistance(vf01_SBrange, f01);           vf01_SBrange.clear();
	float dnSR0_BDTcutoff    = maxdistance(vnSR0_BDTcutoff, nSR0);       vnSR0_BDTcutoff.clear();
	float dnSR0_SBrange      = maxdistance(vnSR0_SBrange, nSR0);         vnSR0_SBrange.clear();
	float VARdnSR0_BDTcutoff = maxdistance(VARvnSR0_BDTcutoff, VARnSR0); VARvnSR0_BDTcutoff.clear();
	float VARdnSR0_SBrange   = maxdistance(VARvnSR0_SBrange, VARnSR0);   VARvnSR0_SBrange.clear();
	
	float dnSR1shape  = nSR1*sqrt(pow(df01_BDTshape/f01,2)  + pow(dnSR0_SBshape/nSR0,2));
	float dnSR1range  = nSR1*sqrt(pow(df01_SBrange/f01,2)   + pow(dnSR0_SBrange/nSR0,2));
	float dnSR1cutoff = nSR1*sqrt(pow(df01_BDTcutoff/f01,2) + pow(dnSR0_BDTcutoff/nSR0,2));
	
	float VARdnSR1shape  = VARnSR1*sqrt(pow(df01_BDTshape/f01,2)  + pow(VARdnSR0_SBshape/VARnSR0,2));
	float VARdnSR1range  = VARnSR1*sqrt(pow(df01_SBrange/f01,2)   + pow(VARdnSR0_SBrange/VARnSR0,2));
	float VARdnSR1cutoff = VARnSR1*sqrt(pow(df01_BDTcutoff/f01,2) + pow(VARdnSR0_BDTcutoff/VARnSR0,2));
	
	float quad = nSR1*sqrt(pow(dnSR0_SBshape/nSR0,2)   +
						   pow(dnSR0_SBrange/nSR0,2)   +
						   pow(dnSR0_BDTcutoff/nSR0,2) +
						   pow(df01_BDTshape/f01,2)    +
						   pow(df01_SBrange/f01,2)     +
						   pow(df01_BDTcutoff/f01,2));
						
	float VARquad = VARnSR1*sqrt(pow(VARdnSR0_SBshape/VARnSR0,2)   +
								 pow(VARdnSR0_SBrange/VARnSR0,2)   +
								 pow(VARdnSR0_BDTcutoff/VARnSR0,2) +
								 pow(df01_BDTshape/f01,2)    +
								 pow(df01_SBrange/f01,2)     +
								 pow(df01_BDTcutoff/f01,2));
	
	if(print)
	{
		cout << "<<<<<<<<<<<< BDT cut at " << tstr(xBDTcut,4) << " >>>>>>>>>>" << endl;
		cout << "iSR = " << iSR << ": [" << xSRmin << "," << xSRmax << "] --> " << wSR << "MeV" << endl;
		cout << "f01 = " << f01 << endl;
		cout << "  df01[BDT shape]  = " << tstr(df01_BDTshape,6) << " (" << tstr(df01_BDTshape/f01*100,1) << "\%)" << endl;
		cout << "  df01[BDT cutoff] = " << tstr(df01_BDTcutoff,6) << " (" << tstr(df01_BDTcutoff/f01*100,1) << "\%)" << endl;
		cout << "  df01[SB range]   = " << tstr(df01_SBrange ,6) << " (" << tstr(df01_SBrange/f01*100 ,1) << "\%)" << endl;
		
		cout << "nSR0 = " << nSR0 << endl;
		cout << "  dnSR0[SB shape]   = " << tstr(dnSR0_SBshape ,3)  << " (" << tstr(dnSR0_SBshape/nSR0*100 ,1)  << "\%)" << endl;
		cout << "  dnSR0[BDT cutoff] = " << tstr(dnSR0_BDTcutoff,3) << " (" << tstr(dnSR0_BDTcutoff/nSR0*100,1) << "\%)" << endl;
		cout << "  dnSR0[SB range]   = " << tstr(dnSR0_SBrange ,3)  << " (" << tstr(dnSR0_SBrange/nSR0*100 ,1)  << "\%)" << endl;
		
		cout << "Systematic quadrature (nSR1) = " << tstr(quad ,3) << " (" << tstr(quad/(nSR0*f01)*100 ,1) << "\%)" << endl;
		cout << "Statistic uncertainty (nSR1) = " << tstr((dnSR0stat*f01) ,3) << " (" << tstr((dnSR0stat*f01)/(nSR0*f01)*100 ,1) << "\%)" << endl;
		
		cout << "VARnSR0 = " << VARnSR0 << endl;
		cout << "  VARdnSR0[SB shape]   = " << tstr(VARdnSR0_SBshape ,3)  << " (" << tstr(VARdnSR0_SBshape/VARnSR0*100 ,1)  << "\%)" << endl;
		cout << "  VARdnSR0[BDT cutoff] = " << tstr(VARdnSR0_BDTcutoff,3) << " (" << tstr(VARdnSR0_BDTcutoff/VARnSR0*100,1) << "\%)" << endl;
		cout << "  VARdnSR0[SB range]   = " << tstr(VARdnSR0_SBrange ,3)  << " (" << tstr(VARdnSR0_SBrange/VARnSR0*100 ,1)  << "\%)" << endl;
		
		cout << "Systematic quadrature (VARnSR1) = " << tstr(VARquad ,3) << " (" << tstr(VARquad/(VARnSR0*f01)*100 ,1) << "\%)" << endl;
		cout << "Statistic uncertainty (VARnSR1) = " << tstr((VARdnSR0stat*f01) ,3) << " (" << tstr((VARdnSR0stat*f01)/(VARnSR0*f01)*100 ,1) << "\%)" << endl;
	}
	
	_INF(vis,"");
	
	Result* res = new Result;
	res->correction = correction;
	res->nSB0 = nSB0; res->dnSB0stat = dnSB0stat;
	res->nSB1 = nSB1; res->dnSB1stat = dnSB1stat;
	res->f01  = f01;                                          res->df01shape     = df01_BDTshape;    res->df01range     = df01_SBrange;     res->df01cutoff     = df01_BDTcutoff;
	res->nSR0 = nSR0; res->dnSR0stat = dnSR0stat;             res->dnSR0shape    = dnSR0_SBshape;    res->dnSR0range    = dnSR0_SBrange;    res->dnSR0cutoff    = dnSR0_BDTcutoff;
	res->nSR1 = nSR1; res->dnSR1stat = dnSR1stat;             res->dnSR1shape    = dnSR1shape;       res->dnSR1range    = dnSR1range;       res->dnSR1cutoff    = dnSR1cutoff;        res->dnSR1quad = quad;
	res->VARnSR0 = VARnSR0; res->VARdnSR0stat = VARdnSR0stat; res->VARdnSR0shape = VARdnSR0_SBshape; res->VARdnSR0range = VARdnSR0_SBrange; res->VARdnSR0cutoff = VARdnSR0_BDTcutoff;
	res->VARnSR1 = VARnSR1; res->VARdnSR1stat = VARdnSR1stat; res->VARdnSR1shape = VARdnSR1shape;    res->VARdnSR1range = VARdnSR1range;    res->VARdnSR1cutoff = VARdnSR1cutoff;     res->VARdnSR1quad = VARquad;
	res->xSRmin = xSRmin; res->xSRmax = xSRmax;
	
	_INF(vis,"");
	
	
	
	
	
	// float df01_BDTcutoff     = maxdistance(vf01_BDTcutoff, nominal_f01);         vf01_BDTcutoff.clear();
	// float df01_SBrange       = maxdistance(vf01_SBrange, nominal_f01);           vf01_SBrange.clear();
	// float dnSR0_BDTcutoff    = maxdistance(vnSR0_BDTcutoff, nominal_nSR0);       vnSR0_BDTcutoff.clear();
	// float dnSR0_SBrange      = maxdistance(vnSR0_SBrange, nominal_nSR0);         vnSR0_SBrange.clear();
	// float VARdnSR0_BDTcutoff = maxdistance(VARvnSR0_BDTcutoff, nominal_VARnSR0); VARvnSR0_BDTcutoff.clear();
	// float VARdnSR0_SBrange   = maxdistance(VARvnSR0_SBrange, nominal_VARnSR0);   VARvnSR0_SBrange.clear();
	// 
	// float dnSR1shape  = nominal_nSR1*sqrt(pow(df01_BDTshape/nominal_f01,2)  + pow(dnSR0_SBshape/nominal_nSR0,2));
	// float dnSR1range  = nominal_nSR1*sqrt(pow(df01_SBrange/nominal_f01,2)   + pow(dnSR0_SBrange/nominal_nSR0,2));
	// float dnSR1cutoff = nominal_nSR1*sqrt(pow(df01_BDTcutoff/nominal_f01,2) + pow(dnSR0_BDTcutoff/nominal_nSR0,2));
	// 
	// float VARdnSR1shape  = nominal_VARnSR1*sqrt(pow(df01_BDTshape/nominal_f01,2)  + pow(VARdnSR0_SBshape/nominal_VARnSR0,2));
	// float VARdnSR1range  = nominal_VARnSR1*sqrt(pow(df01_SBrange/nominal_f01,2)   + pow(VARdnSR0_SBrange/nominal_VARnSR0,2));
	// float VARdnSR1cutoff = nominal_VARnSR1*sqrt(pow(df01_BDTcutoff/nominal_f01,2) + pow(VARdnSR0_BDTcutoff/nominal_VARnSR0,2));
	// 
	// float quad = nominal_nSR1*sqrt(pow(dnSR0_SBshape/nominal_nSR0,2)   +
	// 					   pow(dnSR0_SBrange/nominal_nSR0,2)   +
	// 					   pow(dnSR0_BDTcutoff/nominal_nSR0,2) +
	// 					   pow(df01_BDTshape/nominal_f01,2)    +
	// 					   pow(df01_SBrange/nominal_f01,2)     +
	// 					   pow(df01_BDTcutoff/nominal_f01,2));
	// 					
	// float VARquad = nominal_VARnSR1*sqrt(pow(VARdnSR0_SBshape/nominal_VARnSR0,2)   +
	// 							 pow(VARdnSR0_SBrange/nominal_VARnSR0,2)   +
	// 							 pow(VARdnSR0_BDTcutoff/nominal_VARnSR0,2) +
	// 							 pow(df01_BDTshape/nominal_f01,2)    +
	// 							 pow(df01_SBrange/nominal_f01,2)     +
	// 							 pow(df01_BDTcutoff/nominal_f01,2));
	// 
	// if(print)
	// {
	// 	cout << "<<<<<<<<<<<< BDT cut at " << tstr(xBDTcut,4) << " >>>>>>>>>>" << endl;
	// 	cout << "iSR = " << iSR << ": [" << xSRmin << "," << xSRmax << "] --> " << wSR << "MeV" << endl;
	// 	cout << "f01 = " << nominal_f01 << endl;
	// 	cout << "  df01[BDT shape]  = " << tstr(df01_BDTshape,6) << " (" << tstr(df01_BDTshape/nominal_f01*100,1) << "\%)" << endl;
	// 	cout << "  df01[BDT cutoff] = " << tstr(df01_BDTcutoff,6) << " (" << tstr(df01_BDTcutoff/nominal_f01*100,1) << "\%)" << endl;
	// 	cout << "  df01[SB range]   = " << tstr(df01_SBrange ,6) << " (" << tstr(df01_SBrange/nominal_f01*100 ,1) << "\%)" << endl;
	// 	
	// 	cout << "nSR0 = " << nominal_nSR0 << endl;
	// 	cout << "  dnSR0[SB shape]   = " << tstr(dnSR0_SBshape ,3)  << " (" << tstr(dnSR0_SBshape/nominal_nSR0*100 ,1)  << "\%)" << endl;
	// 	cout << "  dnSR0[BDT cutoff] = " << tstr(dnSR0_BDTcutoff,3) << " (" << tstr(dnSR0_BDTcutoff/nominal_nSR0*100,1) << "\%)" << endl;
	// 	cout << "  dnSR0[SB range]   = " << tstr(dnSR0_SBrange ,3)  << " (" << tstr(dnSR0_SBrange/nominal_nSR0*100 ,1)  << "\%)" << endl;
	// 	
	// 	cout << "Systematic quadrature (nSR1) = " << tstr(quad ,3) << " (" << tstr(quad/(nominal_nSR0*nominal_f01)*100 ,1) << "\%)" << endl;
	// 	cout << "Statistic uncertainty (nSR1) = " << tstr((dnSR0stat*nominal_f01) ,3) << " (" << tstr((dnSR0stat*nominal_f01)/(nominal_nSR0*nominal_f01)*100 ,1) << "\%)" << endl;
	// 	
	// 	cout << "VARnSR0 = " << nominal_VARnSR0 << endl;
	// 	cout << "  VARdnSR0[SB shape]   = " << tstr(VARdnSR0_SBshape ,3)  << " (" << tstr(VARdnSR0_SBshape/nominal_VARnSR0*100 ,1)  << "\%)" << endl;
	// 	cout << "  VARdnSR0[BDT cutoff] = " << tstr(VARdnSR0_BDTcutoff,3) << " (" << tstr(VARdnSR0_BDTcutoff/nominal_VARnSR0*100,1) << "\%)" << endl;
	// 	cout << "  VARdnSR0[SB range]   = " << tstr(VARdnSR0_SBrange ,3)  << " (" << tstr(VARdnSR0_SBrange/nominal_VARnSR0*100 ,1)  << "\%)" << endl;
	// 	
	// 	cout << "Systematic quadrature (VARnSR1) = " << tstr(VARquad ,3) << " (" << tstr(VARquad/(nominal_VARnSR0*nominal_f01)*100 ,1) << "\%)" << endl;
	// 	cout << "Statistic uncertainty (VARnSR1) = " << tstr((nominal_VARdnSR0stat*nominal_f01) ,3) << " (" << tstr((VARdnSR0stat*nominal_f01)/(nominal_VARnSR0*nominal_f01)*100 ,1) << "\%)" << endl;
	// }
	// 
	// _INF(vis,"");
	// 
	// Result* res = new Result;
	// res->correction = nominal_correction;
	// res->nSB0 = nominal_nSB0;       res->dnSB0stat    = nominal_dnSB0stat;
	// res->nSB1 = nominal_nSB1;       res->dnSB1stat    = nominal_dnSB1stat;
	// res->f01  = nominal_f01;                                                  res->df01shape     = df01_BDTshape;         res->df01range     = df01_SBrange;     res->df01cutoff     = df01_BDTcutoff;
	// res->nSR0 = nominal_nSR0;       res->dnSR0stat    = nominal_dnSR0stat;    res->dnSR0shape    = dnSR0_SBshape;         res->dnSR0range    = dnSR0_SBrange;    res->dnSR0cutoff    = dnSR0_BDTcutoff;
	// res->nSR1 = nominal_nSR1;       res->dnSR1stat    = nominal_dnSR1stat;    res->dnSR1shape    = dnSR1shape;            res->dnSR1range    = dnSR1range;       res->dnSR1cutoff    = dnSR1cutoff;        res->dnSR1quad = quad;
	// res->VARnSR0 = nominal_VARnSR0; res->VARdnSR0stat = nominal_VARdnSR0stat; res->VARdnSR0shape = VARdnSR0_SBshape;      res->VARdnSR0range = VARdnSR0_SBrange; res->VARdnSR0cutoff = VARdnSR0_BDTcutoff;
	// res->VARnSR1 = nominal_VARnSR1; res->VARdnSR1stat = nominal_VARdnSR1stat; res->VARdnSR1shape = VARdnSR1shape;         res->VARdnSR1range = VARdnSR1range;    res->VARdnSR1cutoff = VARdnSR1cutoff;     res->VARdnSR1quad = VARquad;
	// res->xSRmin = xSRmin; res->xSRmax = xSRmax;
	// 
	// _INF(vis,"");
	
	
	if(write) ftxt.close();
	
	return res;
}


void getEnvelopes(TH1* hdBDT, TH1* hdSB, TH1* hBDT, TH1* hSB, bool print=false)
{	
	TMapTSf BDTsf, SBsf;
	TMapTSP2Estimate estimates;
	
	TPaveText* ptxt = new TPaveText(0.15,0.70,0.40,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
	
	TCanvas* cnvBDT = new TCanvas("cnvBDT","",600,400);
	TCanvas* cnvSB  = new TCanvas("cnvSB","",600,400);
	
	cnvBDT->cd();
	cnvBDT->Draw();
	cnvBDT->SetTicks(1,1);
	
	cnvSB->cd();
	cnvSB->Draw();
	cnvSB->SetTicks(1,1);
	
	for(TMapTSP2TFILE::iterator it=files.begin() ; it!=files.end() ; ++it)
	{
		TString name = it->first;
		estimates.insert(make_pair(name,setEstimate(it->second)));
		
		TString tf1BDTstr0 = (TString)estimates[name]->vtf1BDTstr[0]->GetTitle(); // cout << "fBDT0=" << tf1BDTstr0 << endl;
		TString tf1BDTstr1 = (TString)estimates[name]->vtf1BDTstr[1]->GetTitle(); // cout << "fBDT1=" << tf1BDTstr1 << endl;
		TString tf1BDTstr2 = (TString)estimates[name]->vtf1BDTstr[2]->GetTitle(); // cout << "fBDT2=" << tf1BDTstr2 << endl;
		TString tf1BDTstr3 = (TString)estimates[name]->vtf1BDTstr[3]->GetTitle(); // cout << "fBDT3=" << tf1BDTstr3 << endl;
		
		
		float sfbdt0 = estimates[name]->nSB/(estimates[name]->vtf1BDT[0]->Integral(estimates[name]->minBDTabs,estimates[name]->maxBDTabs)/estimates[name]->BDTbinWidth);
		float sfbdt1 = estimates[name]->nSB/(estimates[name]->vtf1BDT[1]->Integral(estimates[name]->minBDTabs,estimates[name]->maxBDTabs)/estimates[name]->BDTbinWidth);
		float sfbdt2 = estimates[name]->nSB/(estimates[name]->vtf1BDT[2]->Integral(estimates[name]->minBDTabs,estimates[name]->maxBDTabs)/estimates[name]->BDTbinWidth);
		float sfbdt3 = estimates[name]->nSB/(estimates[name]->vtf1BDT[3]->Integral(estimates[name]->minBDTabs,estimates[name]->maxBDTabs)/estimates[name]->BDTbinWidth);
		BDTsf.insert(make_pair(it->first+".0", sfbdt0 ));
		BDTsf.insert(make_pair(it->first+".1", sfbdt1 ));
		BDTsf.insert(make_pair(it->first+".2", sfbdt2 ));
		BDTsf.insert(make_pair(it->first+".3", sfbdt3 ));
		
		// e.g. 1450-1690-1870-2110.neg090.pos1000
		if(print) cout << "name=" << name << endl;
		TString ranges = name;
		ranges.ReplaceAll(".",":");
		ranges.ReplaceAll("-",":");
		TPMERegexp regexp(":");
		regexp.Split(ranges);
		if(print) regexp.Print("all");
		Double_t mFullMin  = (Double_t)regexp[0].Atof();
		Double_t mBlindMin = (Double_t)regexp[1].Atof();
		Double_t mBlindMax = (Double_t)regexp[2].Atof();
		Double_t mFullMax  = (Double_t)regexp[3].Atof();
		if(print) cout << mFullMin << " -> " << mBlindMin << " -> " << mBlindMax << " -> " << mFullMax << endl;
		
		TString tf1SBstr0 = (TString)estimates[name]->vtf1SBstr[0]->GetTitle(); //cout << "fSB0=" << tf1SBstr0 << endl;
		TString tf1SBstr1 = (TString)estimates[name]->vtf1SBstr[1]->GetTitle(); //cout << "fSB1=" << tf1SBstr1 << endl;
		TString tf1SBstr2 = (TString)estimates[name]->vtf1SBstr[2]->GetTitle(); //cout << "fSB2=" << tf1SBstr2 << endl;
		TString tf1SBstr3 = (TString)estimates[name]->vtf1SBstr[3]->GetTitle(); //cout << "fSB3=" << tf1SBstr3 << endl;
		
		float sfsb0 = estimates[name]->nSB/(estimates[name]->vtf1SB[0]->Integral(mFullMin,mBlindMin)/mBinSize + estimates[name]->vtf1SB[0]->Integral(mBlindMax,mFullMax)/mBinSize);
		float sfsb1 = estimates[name]->nSB/(estimates[name]->vtf1SB[1]->Integral(mFullMin,mBlindMin)/mBinSize + estimates[name]->vtf1SB[1]->Integral(mBlindMax,mFullMax)/mBinSize);
		float sfsb2 = estimates[name]->nSB/(estimates[name]->vtf1SB[2]->Integral(mFullMin,mBlindMin)/mBinSize + estimates[name]->vtf1SB[2]->Integral(mBlindMax,mFullMax)/mBinSize);
		float sfsb3 = estimates[name]->nSB/(estimates[name]->vtf1SB[3]->Integral(mFullMin,mBlindMin)/mBinSize + estimates[name]->vtf1SB[3]->Integral(mBlindMax,mFullMax)/mBinSize);
		SBsf.insert(make_pair(it->first+".0", sfsb0 ));
		SBsf.insert(make_pair(it->first+".1", sfsb1 ));
		SBsf.insert(make_pair(it->first+".2", sfsb2 ));
		SBsf.insert(make_pair(it->first+".3", sfsb3 ));
	}
	
	float xBDTmidle = (estimates[snominal]->vtf1BDT[0]->GetXaxis()->GetXmax()-estimates[snominal]->vtf1BDT[0]->GetXaxis()->GetXmin())/2;
	float dxBDT = (hdBDT->GetXaxis()->GetXmax()-hdBDT->GetXaxis()->GetXmin())/hdBDT->GetNbinsX();
	for(Int_t b=1 ; b<=hdBDT->GetNbinsX() ; ++b)
	{
		Double_t x = hdBDT->GetBinCenter(b);
		float nomBDT = estimates[snominal]->vtf1BDT[0]->Eval(x)*BDTsf[snominal+".0"]; // if(b==2) cout << snominal << ", BDT0: x=" << x << ": " << estimates[snominal]->vtf1BDT[0]->Eval(x)*BDTsf[snominal+".0"] << endl;
		vector<float> vBDTs;		
		for(TMapTSP2Estimate::iterator it=estimates.begin() ; it!=estimates.end() ; ++it)
		{
			TString name = it->first;
						
			// //// for truncating the functions to the edge value where they are not defined, e.g. if minBDTcut=-0.8
			// float val0 = estimates[name]->vtf1BDT[0]->Eval(x)*BDTsf[name+".0"]; val0 = (val0>1e-5) ? val0 : estimates[name]->vtf1BDT[0]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[0]->GetXaxis()->GetXmin() : estimates[name]->vtf1BDT[0]->GetXaxis()->GetXmax())*BDTsf[name+".0"];
			// float val1 = estimates[name]->vtf1BDT[1]->Eval(x)*BDTsf[name+".1"]; val1 = (val1>1e-5) ? val1 : estimates[name]->vtf1BDT[1]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[1]->GetXaxis()->GetXmin() : estimates[name]->vtf1BDT[1]->GetXaxis()->GetXmax())*BDTsf[name+".1"];
			// float val2 = estimates[name]->vtf1BDT[2]->Eval(x)*BDTsf[name+".2"]; val2 = (val2>1e-5) ? val2 : estimates[name]->vtf1BDT[2]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[2]->GetXaxis()->GetXmin() : estimates[name]->vtf1BDT[2]->GetXaxis()->GetXmax())*BDTsf[name+".2"];
			// float val3 = estimates[name]->vtf1BDT[3]->Eval(x)*BDTsf[name+".3"]; val3 = (val3>1e-5) ? val3 : estimates[name]->vtf1BDT[3]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[3]->GetXaxis()->GetXmin() : estimates[name]->vtf1BDT[3]->GetXaxis()->GetXmax())*BDTsf[name+".3"];
			// val0 = ((val0/nomBDT)<5) ? val0 : estimates[name]->vtf1BDT[0]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[0]->GetXaxis()->GetXmin()+dxBDT : estimates[name]->vtf1BDT[0]->GetXaxis()->GetXmax()-dxBDT )*BDTsf[name+".0"];
			// val1 = ((val1/nomBDT)<5) ? val1 : estimates[name]->vtf1BDT[1]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[1]->GetXaxis()->GetXmin()+dxBDT : estimates[name]->vtf1BDT[1]->GetXaxis()->GetXmax()-dxBDT )*BDTsf[name+".1"];
			// val2 = ((val2/nomBDT)<5) ? val2 : estimates[name]->vtf1BDT[2]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[2]->GetXaxis()->GetXmin()+dxBDT : estimates[name]->vtf1BDT[2]->GetXaxis()->GetXmax()-dxBDT )*BDTsf[name+".2"];
			// val3 = ((val3/nomBDT)<5) ? val3 : estimates[name]->vtf1BDT[3]->Eval((x<xBDTmidle) ? estimates[name]->vtf1BDT[3]->GetXaxis()->GetXmin()+dxBDT : estimates[name]->vtf1BDT[3]->GetXaxis()->GetXmax()-dxBDT )*BDTsf[name+".3"];
			// vBDTs.push_back(val0);
			// vBDTs.push_back(val1);
			// vBDTs.push_back(val2);
			// vBDTs.push_back(val3);
						
			vBDTs.push_back(estimates[name]->vtf1BDT[0]->Eval(x)*BDTsf[name+".0"]); 
			if(name!=snominal) continue;                                            
			vBDTs.push_back(estimates[name]->vtf1BDT[1]->Eval(x)*BDTsf[name+".1"]); 
			vBDTs.push_back(estimates[name]->vtf1BDT[2]->Eval(x)*BDTsf[name+".2"]); 
			vBDTs.push_back(estimates[name]->vtf1BDT[3]->Eval(x)*BDTsf[name+".3"]); 
		}
		float distance = maxdistance(vBDTs,nomBDT,false);
		hdBDT->SetBinContent(b,nomBDT);
		hdBDT->SetBinError(b,distance);
	}
	
	
	float xSBmidle = (estimates[snominal]->vtf1SB[0]->GetXaxis()->GetXmax()-estimates[snominal]->vtf1SB[0]->GetXaxis()->GetXmin())/2;
	float dxSB = (hdSB->GetXaxis()->GetXmax()-hdSB->GetXaxis()->GetXmin())/hdSB->GetNbinsX();
	for(Int_t b=1 ; b<=hdSB->GetNbinsX() ; ++b)
	{
		Double_t x = hdSB->GetBinCenter(b);
		float nomSB = estimates[snominal]->vtf1SB[0]->Eval(x)*SBsf[snominal+".0"];

		vector<float> vSBs;
		for(TMapTSP2Estimate::iterator it=estimates.begin() ; it!=estimates.end() ; ++it)
		{
			TString name = it->first;
			
			// //// for truncating the functions to the edge value where they are not defined, e.g. if the SB are: 1540-1690-1870-2020
			// float val0 = estimates[name]->vtf1SB[0]->Eval(x)*SBsf[name+".0"]; val0 = (val0>1e-5) ? val0 : estimates[name]->vtf1SB[0]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[0]->GetXaxis()->GetXmin() : estimates[name]->vtf1SB[0]->GetXaxis()->GetXmax())*SBsf[name+".0"];
			// float val1 = estimates[name]->vtf1SB[1]->Eval(x)*SBsf[name+".1"]; val1 = (val1>1e-5) ? val1 : estimates[name]->vtf1SB[1]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[1]->GetXaxis()->GetXmin() : estimates[name]->vtf1SB[1]->GetXaxis()->GetXmax())*SBsf[name+".1"];
			// float val2 = estimates[name]->vtf1SB[2]->Eval(x)*SBsf[name+".2"]; val2 = (val2>1e-5) ? val2 : estimates[name]->vtf1SB[2]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[2]->GetXaxis()->GetXmin() : estimates[name]->vtf1SB[2]->GetXaxis()->GetXmax())*SBsf[name+".2"];
			// float val3 = estimates[name]->vtf1SB[3]->Eval(x)*SBsf[name+".3"]; val3 = (val3>1e-5) ? val3 : estimates[name]->vtf1SB[3]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[3]->GetXaxis()->GetXmin() : estimates[name]->vtf1SB[3]->GetXaxis()->GetXmax())*SBsf[name+".3"];
			// val0 = ((val0/nomSB)<5) ? val0 : estimates[name]->vtf1SB[0]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[0]->GetXaxis()->GetXmin()+dxSB : estimates[name]->vtf1SB[0]->GetXaxis()->GetXmax()-dxSB )*SBsf[name+".0"];
			// val1 = ((val1/nomSB)<5) ? val1 : estimates[name]->vtf1SB[1]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[1]->GetXaxis()->GetXmin()+dxSB : estimates[name]->vtf1SB[1]->GetXaxis()->GetXmax()-dxSB )*SBsf[name+".1"];
			// val2 = ((val2/nomSB)<5) ? val2 : estimates[name]->vtf1SB[2]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[2]->GetXaxis()->GetXmin()+dxSB : estimates[name]->vtf1SB[2]->GetXaxis()->GetXmax()-dxSB )*SBsf[name+".2"];
			// val3 = ((val3/nomSB)<5) ? val3 : estimates[name]->vtf1SB[3]->Eval((x<xSBmidle) ? estimates[name]->vtf1SB[3]->GetXaxis()->GetXmin()+dxSB : estimates[name]->vtf1SB[3]->GetXaxis()->GetXmax()-dxSB )*SBsf[name+".3"];
			// vSBs.push_back(val0);
			// vSBs.push_back(val1);
			// vSBs.push_back(val2);
			// vSBs.push_back(val3);
			
			vSBs.push_back(estimates[name]->vtf1SB[0]->Eval(x)*SBsf[name+".0"]);
			if(name!=snominal) continue;
			vSBs.push_back(estimates[name]->vtf1SB[1]->Eval(x)*SBsf[name+".1"]);
			vSBs.push_back(estimates[name]->vtf1SB[2]->Eval(x)*SBsf[name+".2"]);
			vSBs.push_back(estimates[name]->vtf1SB[3]->Eval(x)*SBsf[name+".3"]);
		}
		float distance = maxdistance(vSBs,nomSB,false);
		if(print) cout << "SB bin[" << b << "]: yield=" << nomSB << ", error=" << distance << endl;
		hdSB->SetBinContent(b,nomSB);
		hdSB->SetBinError(b,distance);
	}
	
	
	
	
	cnvBDT->cd();
	hdBDT->SetMinimum(0);
	hdBDT->SetMaximum(20);
	hdBDT->Draw("e2");
	
	Double_t xhBDTmin = hBDT->GetXaxis()->GetXmin();
	Double_t xhBDTmax = hBDT->GetXaxis()->GetXmax();
	
	// for(TMapTSP2Estimate::iterator it=estimates.begin() ; it!=estimates.end() ; ++it)
	// {	
	// 	TString name = it->first;
	// 	if(name==snominal) continue;
	// 	
	// 	if(fBDTnorm) delete fBDTnorm; fBDTnorm = (TF1*)estimates[name]->vtf1BDT[0]->Clone("fBDTnorm_"+name); BDTscale = BDTsf[name+".0"]; fBDTnorm->SetRange(-0.9,+1);
	// 	TF1* fBDTscaled = new TF1("fBDTscaled_"+name,funcBDTscaled,xhBDTmin,xhBDTmax);
	// 	fBDTscaled->SetLineColor(kYellow); fBDTscaled->SetLineWidth(1); fBDTscaled->DrawCopy("same"); cnvBDT->Update(); cnvBDT->Draw();
	// }
	// 
	// if(fBDTnorm) delete fBDTnorm; fBDTnorm = (TF1*)estimates[snominal]->vtf1BDT[3]->Clone("fBDTnorm3"); BDTscale = BDTsf[snominal+".3"]; fBDTnorm->SetRange(-0.9,+1);
	// TF1* fBDTscaled3 = new TF1("fBDTscaled3",funcBDTscaled,xhBDTmin,xhBDTmax);
	// fBDTscaled3->SetLineColor(kBlue); fBDTscaled3->SetLineWidth(1); fBDTscaled3->DrawCopy("same"); cnvBDT->Update(); cnvBDT->Draw();
	// 
	// if(fBDTnorm) delete fBDTnorm; fBDTnorm = (TF1*)estimates[snominal]->vtf1BDT[2]->Clone("fBDTnorm2"); BDTscale = BDTsf[snominal+".2"]; fBDTnorm->SetRange(-0.9,+1);
	// TF1* fBDTscaled2 = new TF1("fBDTscaled2",funcBDTscaled,xhBDTmin,xhBDTmax);
	// fBDTscaled2->SetLineColor(kRed+1); fBDTscaled2->SetLineWidth(1); fBDTscaled2->DrawCopy("same"); cnvBDT->Update(); cnvBDT->Draw();
	// 
	// if(fBDTnorm) delete fBDTnorm; fBDTnorm = (TF1*)estimates[snominal]->vtf1BDT[1]->Clone("fBDTnorm1"); BDTscale = BDTsf[snominal+".1"]; fBDTnorm->SetRange(-0.9,+1);
	// TF1* fBDTscaled1 = new TF1("fBDTscaled1",funcBDTscaled,xhBDTmin,xhBDTmax);
	// fBDTscaled1->SetLineColor(kGray+1); fBDTscaled1->SetLineWidth(1); fBDTscaled1->DrawCopy("same"); cnvBDT->Update(); cnvBDT->Draw();
	
	if(fBDTnorm) delete fBDTnorm; fBDTnorm = (TF1*)estimates[snominal]->vtf1BDT[0]->Clone("fBDTnorm0"); BDTscale = BDTsf[snominal+".0"]; fBDTnorm->SetRange(-0.9,+1);
	TF1* fBDTscaled0 = new TF1("fBDTscaled0",funcBDTscaled,xhBDTmin,xhBDTmax);
	fBDTscaled0->SetLineColor(kPink); fBDTscaled0->DrawCopy("same"); cnvBDT->Update(); cnvBDT->Draw();
	
	hBDT->Draw("p0 e same");
	ptxt->Draw("same");
	cnvBDT->Update();
	cnvBDT->RedrawAxis();
	cnvBDT->SaveAs("figures/syst.BDTenvelope.pdf");
	cnvBDT->SaveAs("figures/syst.BDTenvelope.png");
	cnvBDT->SaveAs("figures/syst.BDTenvelope.eps");
	cnvBDT->SaveAs("figures/syst.pdf(");
	
	
	
	
	cnvSB->cd();
	hdSB->SetMinimum(0);
	hdSB->SetMaximum(20);
	hdSB->Draw("e2");
	
	Double_t xhSBmin = hSB->GetXaxis()->GetXmin();
	Double_t xhSBmax = hSB->GetXaxis()->GetXmax();

	// for(TMapTSP2Estimate::iterator it=estimates.begin() ; it!=estimates.end() ; ++it)
	// {	
	// 	TString name = it->first;
	// 	if(name==snominal) continue;
	// 	
	// 	if(fSBnorm) delete fSBnorm; fSBnorm = (TF1*)estimates[name]->vtf1SB[0]->Clone("fSBnorm_"+name); SBscale = SBsf[name+".0"]; fSBnorm->SetRange(1450,2110);
	// 	TF1* fSBscaled = new TF1("fSBscaled_"+name,funcSBscaled,xhSBmin,xhSBmax);
	// 	fSBscaled->SetLineColor(kYellow); fSBscaled->SetLineWidth(1); fSBscaled->DrawCopy("same"); cnvSB->Update(); cnvSB->Draw();
	// }
	// if(fSBnorm) delete fSBnorm; fSBnorm = (TF1*)estimates[snominal]->vtf1SB[3]->Clone("fSBnorm3"); SBscale = SBsf[snominal+".3"]; fSBnorm->SetRange(1450,2110);
	// TF1* fSBscaled3 = new TF1("fSBscaled3",funcSBscaled,xhSBmin,xhSBmax);
	// fSBscaled3->SetLineColor(kBlue); fSBscaled3->SetLineWidth(1); fSBscaled3->DrawCopy("same"); cnvSB->Update(); cnvSB->Draw();
	// 
	// if(fSBnorm) delete fSBnorm; fSBnorm = (TF1*)estimates[snominal]->vtf1SB[2]->Clone("fSBnorm2"); SBscale = SBsf[snominal+".2"]; fSBnorm->SetRange(1450,2110);
	// TF1* fSBscaled2 = new TF1("fSBscaled2",funcSBscaled,xhSBmin,xhSBmax);
	// fSBscaled2->SetLineColor(kRed+1); fSBscaled2->SetLineWidth(1); fSBscaled2->DrawCopy("same"); cnvSB->Update(); cnvSB->Draw();
	// 
	// if(fSBnorm) delete fSBnorm; fSBnorm = (TF1*)estimates[snominal]->vtf1SB[1]->Clone("fSBnorm1"); SBscale = SBsf[snominal+".1"]; fSBnorm->SetRange(1450,2110);
	// TF1* fSBscaled1 = new TF1("fSBscaled1",funcSBscaled,xhSBmin,xhSBmax);
	// fSBscaled1->SetLineColor(kGray+1); fSBscaled1->SetLineWidth(1); fSBscaled1->DrawCopy("same"); cnvSB->Update(); cnvSB->Draw();
	
	if(fSBnorm) delete fSBnorm; fSBnorm = (TF1*)estimates[snominal]->vtf1SB[0]->Clone("fSBnorm0"); SBscale = SBsf[snominal+".0"]; fSBnorm->SetRange(1450,2110);
	TF1* fSBscaled0 = new TF1("fSBscaled0",funcSBscaled,xhSBmin,xhSBmax);
	fSBscaled0->SetLineColor(kPink); fSBscaled0->DrawCopy("same"); cnvSB->Update(); cnvSB->Draw();

	hSB->Draw("p0 e same");
	ptxt->Draw("same");
	cnvSB->Update();
	cnvSB->RedrawAxis();
	cnvSB->SaveAs("figures/syst.SBenvelope.pdf");
	cnvSB->SaveAs("figures/syst.SBenvelope.png");
	cnvSB->SaveAs("figures/syst.SBenvelope.eps");
	cnvSB->SaveAs("figures/syst.pdf)");	
}


#endif
