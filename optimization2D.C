//////////////////////////////////////////////////////////////////////////////////////////////////////////
//// root -b -l -q optimization2D.C++\(1450.,1690.,1870.,2110.,50.,0.8000,0.9900,0.01,20,180,10,-1\)
//// root -b -l -q optimization2D.C++\(1450.,1690.,1870.,2110.,50.,0.8400,0.9900,0.01,60,180,4,-1\)
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "allFitSystAuto.h"
#include "postBDTcuts.h"
#include <iterator>
#include <cstdlib>

TDirectory* olddir = gDirectory;

Bool_t reject;
Double_t xbmin = 1700.;
Double_t xbmax = 1860.;
Int_t NminFlat = 10000;
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

TBox* boxHighlight(TH2F* h)
{
	Int_t bin = h->GetMinimumBin();
	Int_t binx, biny, binz;
	h->GetBinXYZ(bin, binx, biny, binz);
	Double_t x1 = h->GetXaxis()->GetBinLowEdge(binx);
	Double_t x2 = h->GetXaxis()->GetBinUpEdge(binx);
	Double_t y1 = h->GetYaxis()->GetBinLowEdge(biny);
	Double_t y2 = h->GetYaxis()->GetBinUpEdge(biny);
	TBox* box = new TBox(x1,y1,x2,y2);
	box->SetFillStyle(0);
	box->SetLineWidth(4);
	box->SetLineColor(kBlack);
	return box;
}
TH2F* subtract(TH2F* hNominal, TH2F* hShited)
{
	TH2F* h = (TH2F*)hNominal->Clone();
	h->Reset();
	
	Int_t binx = h->GetNbinsX();
	Int_t biny = h->GetNbinsY();
	for(Int_t x=0 ; x<=binx ; ++x)
	{
		for(Int_t y=0 ; y<=biny ; ++y)
		{
			float d = fabs(hNominal->GetBinContent(x,y)-hShited->GetBinContent(x,y));
			h->SetBinContent(x,y,d);
		}
	}
	return h;
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
	ptxt = new TPaveText(0.15,0.73,0.40,0.90,"NDC");
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

string makeHistFitterPy(vector<string>& lines, string score, string srsize, string nobs, float nbkgSR, float nbkgSRerr, float systUp, float systDwn, float ndata)
{
	string fname = "python/tau3mu/tau3mu.BDTcut"+score+".SRsize"+srsize+".Nobs"+nobs+".py";
	string dir = "BDTcut"+score+".SRsize"+srsize+".Nobs"+nobs;
	
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
				if     (j==2)      hfpyx << " " << systUp << " ";
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
				if     (j==2)      hfpyx << systDwn << " ";
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
		// else if(lines[i].find("configMgr.analysisName = ")!=string::npos)
		// {
		// 	vector<string> words = split_string_by_spaces(lines[i],false);
		// 	unsigned int nwords = words.size();
		// 	for(unsigned int j=0 ; j<nwords ; ++j)
		// 	{
		// 		if     (j==2)      hfpyx << "\"Optimizaion2D\"" << " ";
		// 		else if(j==nwords) hfpyx << endl;
		// 		else               hfpyx << words[j] << " ";
		// 	}
		// 	hfpyx << endl;
		// }
		// else if(lines[i].find("configMgr.outputFileName = ")!=string::npos)
		// {
		// 	vector<string> words = split_string_by_spaces(lines[i],false);
		// 	unsigned int nwords = words.size();
		// 	for(unsigned int j=0 ; j<nwords ; ++j)
		// 	{
		// 		if     (j==2)      hfpyx << "\"/tmp/hod/HistFitter/results/"+dir+"/\%s_Output.root\"" << "\%configMgr.analysisName";
		// 		else if(j==nwords) hfpyx << endl;
		// 		else               hfpyx << words[j] << " ";
		// 	}
		// 	hfpyx << endl;
		// }
		// else if(lines[i].find("    if os.path.isfile(")!=string::npos)
		// {
		// 	hfpyx << "    if os.path.isfile(\"/tmp/hod/HistFitter/"+dir+"/data/\%s.root\"\%configMgr.analysisName):" << endl;
		// 	hfpyx << endl;
		// }
		// else if(lines[i].find("        os.remove(")!=string::npos)
		// {
		// 	hfpyx << "        os.remove(\"/tmp/hod/HistFitter/"+dir+"/data/\%s.root\"\%configMgr.analysisName)" << endl;
		// 	hfpyx << endl;
		// }
		else hfpyx << lines[i] << endl;
	}
	hfpyx.close();
	return fname;
}

// TString postBDTcut(TString xFullMin, TString xBlindMin, TString xBlindMax, TString xFullMax, TString xBDT, TString sigma, TString type, TString xSRmin="", TString xSRmax="")
// {
// 	TString sidebandscuts = "(m3body>"+xFullMin+". && m3body<"+xFullMax+".)";
// 	TString blindedcuts   = "(m3body<"+xBlindMin+". || m3body>"+xBlindMax+".)";
// 	
// 	TString srcuts        = (xSRmin!="" && xSRmax!="") ? "(m3body>"+xSRmin+". && m3body<"+xSRmax+".)" : "(1)";
// 	TString bdtcuts       = "score>"+xBDT;
// 	
// 	TString ncandcuts     = "(@m3body->size()==1)";
// 	
// 	TString iso           = "(iso003<0.1 && iso020<0.3)";
// 	TString trkspval      = "(trkspval>1.e-9)";
// 	TString distances     = "(lxysig>-10 && lxysig<50 && a0xysig<25)";
// 	TString kinematics    = "(pt3body>10000 && (mettrk>10000 && mettrk<250000 && mttrk>20000) && (metcal>10000 && metcal<250000 && mtcal>20000))";
// 	TString lowmass       = "(mOS1>300. && mOS2>300. && mSS>300.)";
// 	TString preTraining   = "("+iso+" && "+trkspval+" && "+distances+" && "+kinematics+" && "+lowmass+")";
// 	
// 	TString HT            = "(mttrk>60000 && mtcal>60000 && dphihttrk>2 && dphihtcal>2)";
// 	// TString HT            = "(ht>20000 && dphihttrk>2. && mhttrk>60000 && mhtcal>60000)";
// 	// TString resonances    = "!(((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50 || TMath::Abs(mSS-1020)<50) || (TMath::Abs(mOS1-782)<50 || TMath::Abs(mOS2-782)<50 || TMath::Abs(mSS-782)<50)) && (mhttrk<60000 || mhtcal<60000))";
// 	TString resonances    = "!( ((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50) || (TMath::Abs(mOS1-782)<50 || TMath::Abs(mOS2-782)<50)) && (mettrk<35000 || metcal<35000 || pt3body<35000) )";
// 	// TString low2body      = "!((mOS1<700 || mOS2<700 || mSS<700) && (mhttrk<60000 || mhtcal<60000))";
// 	// TString postTraining  = "("+HT+" && "+resonances+" && "+low2body+")";
// 	TString postTraining  = HT+" && "+resonances;
// 	
// 	TString basecuts      = ncandcuts+" && "+sidebandscuts+" && "+preTraining+" && "+postTraining;
// 	
// 	TString cut = "";
// 	if     (type=="bkg") cut = basecuts+" && "+bdtcuts+" && "+blindedcuts;
// 	else if(type=="sig") cut = basecuts+" && "+bdtcuts+" && "+srcuts;
// 	else _FAT("unknown type: "<<type);
// 	return cut;
// }

struct Limit
{
	bool  isBad;
	bool  ok;
	float upper;  
	float median; 
	float m1sigma;
	float p1sigma;
	float m2sigma;
	float p2sigma;	
};

Limit* getLimit(string pyfile, TString score, TString srsize, TString nobs)
{
	Limit* lim = new Limit;
	
	string fname = "/tmp/hod/HistFitter/txt/histfitter.BDTcut"+(string)score+".SRsize"+(string)srsize+".Nobs"+(string)nobs+".txt";

	string dir = "/tmp/hod/HistFitter/results/BDTcut"+(string)score+".SRsize"+(string)srsize+".Nobs"+(string)nobs;
	string mkdir = "mkdir "+dir;
	string rmdir = "rm -rf "+dir;
	cout << "Creating results directory " << dir << endl;
	gSystem->Exec(mkdir.c_str());
	string command = "HistFitter.py -w -f -l "+pyfile+" > "+fname+" 2>&1";
	cout << "Running python file " << pyfile << endl;
	gSystem->Exec(command.c_str());
	cout << "Cleaning results directory " << dir << endl;
	gSystem->Exec(rmdir.c_str());
	cout << "Output was written to " << fname << endl;
	
	// gSystem->Exec("/bin/cp -f results/upperlimit_cls_poi_Sig_Asym_CLs_grid_ts3.root.eps /tmp/hod/HistFitter/eps/upperlimit_cls_poi_Sig_Asym_CLs_grid_ts3.BDTcut"+score+".SRsize"+srsize+".Nobs"+nobs+".root.eps");
	// gSystem->Exec("/bin/cp -f results/MyUserAnalysis_Output_upperlimit.root /tmp/hod/HistFitter/eps/MyUserAnalysis_Output_upperlimit.BDRcut"+score+".SRsize"+srsize+".Nobs"+nobs+".root");
	ifstream fileInput;
	string line;
	// open file to search
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
		if(upper>=0 && upper<50 && median>=0 && m1sigma>=0 && p1sigma>=0 && m2sigma>=0 && p2sigma>=0) { ok = true; break; }
	}
	
	bool isBad = (median<1.e-10 || upper<1.e-10 || upper>=50. || m1sigma<1.e-10 || p1sigma<1.e-10 || m2sigma<1.e-10 || p2sigma<1.e-10);
	
	lim->upper   = upper  ;  
	lim->median  = median ; 
	lim->m1sigma = m1sigma;
	lim->p1sigma = p1sigma;
	lim->m2sigma = m2sigma;
	lim->p2sigma = p2sigma;
	lim->ok      = ok;
	lim->isBad   = isBad;
	
	if(!ok) { _ERR(1,"Could not find one of the values: upper="<<upper<<", median="<<median<<", m1sigma="<<m1sigma<<", p1sigma="<<p1sigma<<", m2sigma="<<m2sigma<<", p2sigma"<<p2sigma); }
	else    { cout << "HistFitter summary: upper="<<upper<<", median="<<median<<", m1sigma="<<m1sigma<<", p1sigma="<<p1sigma<<", m2sigma="<<m2sigma<<", p2sigma="<<p2sigma << endl; }
	
	if(isBad) // something is wrong --> use the result from the previous iteration
	{
		cout << "#########################################" << endl;
		cout << "######### bad HistFitter result #########" << endl;
		cout << "#########################################" << endl;
	}
	
	return lim;
}


float getAverage(TH2F* h, int xcenter, int ycenter, int radius)
{
	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();
	int xmin = ((xcenter-radius)>=1)  ? xcenter-radius : 1;
	int xmax = ((xcenter+radius)<=nx) ? xcenter+radius : nx;
	int ymin = ((ycenter-radius)>=1)  ? ycenter-radius : 1;
	int ymax = ((ycenter+radius)<=ny) ? ycenter+radius : ny;
	float nbinsused = 0;
	float av = 0;
	for(int x=xmin ; x<=xmax ; ++x)
	{
		for(int y=ymin ; y<=ymax ; ++y)
		{
			float content = h->GetBinContent(x,y);
			if(x==xcenter && y==ycenter) continue;
			if(content<=0.)              continue;
			av += content;
			nbinsused++;
		}
	}
	float result = (nbinsused>0) ? av/nbinsused : 0;
	cout << "Using histogram: " << h->GetName() << " [" << xcenter << "," << ycenter << "] --> setting to average: " << result << endl;
	return result;
}


void optimization2D(float mMinSBleft, float mMaxSBleft, float mMinSBright, float mMaxSBright, float m3bodyBinSize, float MinBDTcut, float MaxBDTcut, float dX, float Ymin=20, float Ymax=180, int dY=10, int Nobs=-1)
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
	TString srmin        = tstr(Ymin,0);
	TString srmax        = tstr(Ymax,0);
	TString srdelta      = tstr(dY,0);
	TString nobs         = (Nobs>=0) ? tstr(Nobs,0) : "X";
	TString outfilesname = "Optimization2D."+sidebands+"."+bdtmin+"-"+bdtmax+"."+srmin+"-"+srmax+".Nobs"+nobs;

	
	TFile* fOut = new TFile(outfilesname+".root","RECREATE");


	makeAtlasLabel();
	makeLegend();
	
	setFiles();
	TFile* fBDTpdf = getNominal();
	TVectorF* AbsBDTmin = (TVectorF*)fBDTpdf->Get("absBDTmin");
	TVectorF* AbsBDTmax = (TVectorF*)fBDTpdf->Get("absBDTmax");
	TVectorF* NBDTbins  = (TVectorF*)fBDTpdf->Get("nBDTbins");
	float minBDTabs = ((*AbsBDTmin))[0];
	float maxBDTabs = ((*AbsBDTmax))[0];
	int   BDTbins   = (int)((*NBDTbins))[0];

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
	Double_t ximin = avg-3.*sigma;
	Double_t ximax = avg+3.*sigma;
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
	
	
	// TH1F* mD = new TH1F("mD", ";m_{3body} [MeV];Events", (Int_t)((mMaxSBright-mMinSBleft)/m3bodyBinSize),mMinSBleft,mMaxSBright);
	// mD->SetLineColor(kBlack);
	// mD->SetLineWidth(1);
	// mD->SetLineStyle(1);
	// mD->SetMarkerStyle(20);
	// mD->SetMarkerSize(0.6);
	// mD->SetMarkerColor(kBlack);
	// mD->Sumw2();



	bool hFit = true; // true for histogram fit, false for graph fit

	// remember the python script
	vector<string> hfpylines;
	string hfpyline;
	fstream hfpy;
	hfpy.open("python/MyUserAnalysis.py");
	while(getline(hfpy,hfpyline)) hfpylines.push_back(hfpyline);
	hfpy.close();
	

	if(Nobs<=0)
	{
		// gSystem->Exec("rm -f python/tau3mu/*");
		// gSystem->Exec("rm -f txt/*");
		// gSystem->Exec("rm -rf /tmp/hod/HistFitter");
	}
	gSystem->Exec("mkdir /tmp/hod/HistFitter");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/eps");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/txt");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/results");
	gSystem->Exec("mkdir /tmp/hod/HistFitter/data");
	
	float dx = dX;
	float xmin = MinBDTcut;
	float xmax = MaxBDTcut;
	Float_t hXmin = xmin-dx/2.;
	Float_t hXmax = xmax+dx/2.;
	Int_t hNXbins  = (hXmax-hXmin)/dx;
	
	TFile* fTmp = new TFile("BDTfitResults.1450-1690-1870-2290.neg090.pos1000.root","READ");
	TVectorF* xsrmin = (TVectorF*)fTmp->Get("xSRmin");
	TVectorF* xsrmax = (TVectorF*)fTmp->Get("xSRmax");
	vector<int> srwidths, srindices;
	for(int i=0 ; i<xsrmin->GetNoElements() ; ++i)
	{
		int d = (int)(((*xsrmax))[i]-((*xsrmin))[i]);
		if(d%dY!=0 || d<Ymin || d>Ymax) continue;
		srwidths.push_back(d);
		srindices.push_back(i);
	}
	Int_t Ny = srwidths.size();
	cout << "xsrmin->GetNoElements()=" << xsrmin->GetNoElements() << endl;
	cout << "Ny=" << Ny << endl;
	
	Float_t hYmin = 0.;
	Float_t hYmax = Ny;
	Int_t hNYbins = Ny;

	TH2F* hEffSrel  = new TH2F("SignalEfficiencyRel", "Signal #it{A}#times#it{#epsilon};BDT score cut;SR size [MeV];Signal #it{A}#times#it{#epsilon} [\%] (relative)", hNXbins,hXmin,hXmax,  hNYbins,hYmin,hYmax);
	TH2F* hEffS     = new TH2F("SignalEfficiency",    "Signal #it{A}#times#it{#epsilon};BDT score cut;SR size [MeV];Signal #it{A}#times#it{#epsilon} [\%] (absolute)", hNXbins,hXmin,hXmax,  hNYbins,hYmin,hYmax);
	TH2F* hUL90min1 = new TH2F("UpperLimitMin1Sig",   "Expected -1#sigma Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hUL90pls1 = new TH2F("UpperLimitPls1Sig",   "Expected +1#sigma Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hUL90min2 = new TH2F("UpperLimitMin2Sig",   "Expected -2#sigma Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hUL90pls2 = new TH2F("UpperLimitPls2Sig",   "Expected +2#sigma Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hUL90med  = new TH2F("UpperLimitMedian",    "Expected (median) Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hUL90obs  = new TH2F("UpperLimitObserved",  "Observed (N_{obs}="+nobs+") Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs;BDT score cut;SR size [MeV];Br(#it{#tau#rightarrow3#mu})/10^{-7} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hNobs     = new TH2F("Nobs",                "N_{obs} (obs events in SR1 used for the fit);BDT score cut;SR size [MeV];N_{obs} (obs events in SR1 used for the fit)", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hNdat     = new TH2F("Ndat",                "N_{dat}=floor(N_{SR1}+0.5) (obs events in SR1 used for the fit);BDT score cut;SR size [MeV];N_{dat}=floor(N_{SR1}+0.5) (obs events in SR1 used for the fit)", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSB1Count = new TH2F("SB1 Count",          "SB1 Count;BDT score cut;SR size [MeV];Events passing BDT cut in SB1", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSB1Shape = new TH2F("SB1 Shape",          "SB1 Shape;BDT score cut;SR size [MeV];Events passing BDT cut in SB1", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1Count = new TH2F("SR1 Count",          "SR1 Count;BDT score cut;SR size [MeV];Events passing BDT cut in SR1", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1Shape = new TH2F("SR1 Shape",          "SR1 Shape;BDT score cut;SR size [MeV];Events passing BDT cut in SR1", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsObs      = new TH2F("SR1 NsObs",     "Observed (N_{obs}="+nobs+") N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsObsExp   = new TH2F("SR1 NsObsExp",  "Observed (N_{obs}=exp) N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} exp at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsMed      = new TH2F("SR1 NsMed",     "Expected (median) N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsMedMin1  = new TH2F("SR1 NsMedMin1", "Expected -1#sigma N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsMedPls1  = new TH2F("SR1 NsMedPls1", "Expected +1#sigma N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsMedMin2  = new TH2F("SR1 NsMedMin2", "Expected -2#sigma N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	TH2F* hnSR1NsMedPls2  = new TH2F("SR1 NsMedPls2", "Expected +2#sigma N_{s} at 90\% CLs;BDT score cut;SR size [MeV];N_{s}^{UL} at 90\% CLs", hNXbins,hXmin,hXmax, hNYbins,hYmin,hYmax);
	for(int i=0 ; i<Ny ; ++i)
	{
		TString ssrwidth = tstr(srwidths[i],0);
		hEffSrel  ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hEffS     ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90min1 ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90pls1 ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90min2 ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90pls2 ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90med  ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hUL90obs  ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hNobs     ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hNdat     ->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSB1Count->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSB1Shape->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1Count->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1Shape->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsMed->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsObs->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsObsExp->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsMedMin1->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsMedMin2->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsMedPls1->GetYaxis()->SetBinLabel(i+1,ssrwidth);
		hnSR1NsMedPls2->GetYaxis()->SetBinLabel(i+1,ssrwidth);
	}
	
	TF1* fbg0 = NULL;
	TF1* fbg1 = NULL;
	TF1* fbg2 = NULL;

	stringstream strm;
	string       str;
	
	Int_t hbinx = 1;
	for(float x=xmin ; x<xmax ; x+=dx, ++hbinx)
	{
		cout << "\n" << endl;
	
		strm.clear();
		str.clear();
		
		strm << setprecision(4) << fixed << x;
		strm >> str;
		TString score = str;

		for(Int_t hbiny=1 ; hbiny<=Ny ; ++hbiny)
		{
			cout << "\n" << endl;
			
			unsigned int y = srindices[hbiny-1];
			cout << "srindices["<<hbiny-1<<"]=" << y << endl;
			
			Result* res = getExtrapolationResult(x,y,false);
			float srwidth      = res->xSRmax-res->xSRmin; 
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
			
			cout << "srwidth=" << srwidth << endl;
			
			strm.clear();
			str.clear();
			if     (srwidth<10)  strm << "00" << setprecision(0) << fixed << srwidth << "MeV";
			else if(srwidth<100) strm << "0"  << setprecision(0) << fixed << srwidth << "MeV";
			else                 strm <<         setprecision(0) << fixed << srwidth << "MeV";
			strm >> str;
			TString srsize = str;
	    	
			TString cuts_bkg      = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9",score,"bkg");
			TString cuts_sig      = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9",score,"sig",sxSRmin,sxSRmax);
			TString cuts_sig_noSR = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,"-0.9",score,"sig");
			// TString cuts_bkg      = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,score,"30","bkg");
			// TString cuts_sig      = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,score,"30","sig",sxSRmin,sxSRmax); // with SR cut
			// TString cuts_sig_noSR = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,score,"30","sig"); // without SR cut
			cout << "cuts_bkg      = " << cuts_bkg << endl;
			cout << "cuts_sig      = " << cuts_sig << endl;
			cout << "cuts_sig_noSR = " << cuts_sig_noSR << endl;
        	
			tS->SetEventList(0);
			tS->Draw(">>elistSnoSR",cuts_sig_noSR);
			TEventList* elistSnoSR = (TEventList*)gDirectory->Get("elistSnoSR");
			Float_t npassedSnoSR = elistSnoSR->GetN();  // number of events to pass cuts
	    	
			tS->SetEventList(0);
			tS->Draw(">>elistS",cuts_sig);
			TEventList* elistS = (TEventList*)gDirectory->Get("elistS");
			Float_t npassedS = elistS->GetN();  // number of events to pass cuts
	
			tD->SetEventList(0);
			tD->Draw(">>elistD",cuts_bkg);
			TEventList* elistD = (TEventList*)gDirectory->Get("elistD");
			Float_t npassedD = elistD->GetN();  // number of events to pass cuts
			cout << cuts_bkg << ":  npassedS=" << npassedS << ", npassedD=" << npassedD << ", nfitD=" << iBkgSB1 << "+-" << diBkgSB1stat << "(for f01="<<f01<<")" << ", acc*eff(SR)=" << npassedS/ninitS*100 << "\% (inclusive: " << npassedSnoSR/ninitS*100 << "%)" << endl;
			
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
			


			TH1F* mD = new TH1F("mD", ";m_{3body} [MeV];Events", (Int_t)((mMaxSBright-mMinSBleft)/m3bodyBinSize),mMinSBleft,mMaxSBright);
			mD->SetLineColor(kBlack);
			mD->SetLineWidth(1);
			mD->SetLineStyle(1);
			mD->SetMarkerStyle(20);
			mD->SetMarkerSize(0.6);
			mD->SetMarkerColor(kBlack);
			mD->Sumw2();

			tD->SetEventList(elistD);
			tD->Draw("m3body>>mD", cuts_bkg);	
			for(int b=1 ; b<=mD->GetNbinsX() ; ++b) { if(mD->GetBinContent(b)<1) mD->SetBinError(b,1); }
			for(int b=1 ; b<=mD->GetNbinsX() ; ++b) { if(mD->GetBinContent(b)<1 && mD->GetXaxis()->GetBinLowEdge(b)>=xbmin && mD->GetXaxis()->GetBinUpEdge(b)<=xbmax) mD->SetBinContent(b,-10); }
			mD->SetMinimum(0.);
			if(mD->GetMaximum()<3.) mD->SetMaximum(3.);
			TGraphAsymmErrors* gD = (TGraphAsymmErrors*)poisson(mD);
			reject = kTRUE;
			if(hFit)
			{
				if(npassedD>NminFlat) { mD->Fit("Background0","L0 EMR"); _INF(1,"fitted: Background0 (TH1)"); }
				else                  { mD->Fit("Background1","L0 EMR"); _INF(1,"fitted: Background1 (TH1)"); }
			}
			else
			{
				if(npassedD>NminFlat) { gD->Fit("Background0","EX0"); _INF(1,"fitted: Background0 (TGraph)"); }
				else                  { gD->Fit("Background1","EX0"); _INF(1,"fitted: Background1 (TGraph)"); }
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
			pvalue = TMath::Prob(chi2,ndf);
			float iBkgSR = (npassedD>NminFlat) ? fbg0->Integral(ximin,ximax)/mD->GetBinWidth(1) : fbg1->Integral(ximin,ximax)/mD->GetBinWidth(1);
			float iBkgSB = (npassedD>NminFlat) ? (fbg0->Integral(mMinSBleft,xbmin)+fbg0->Integral(xbmax,mMaxSBright))/mD->GetBinWidth(1) : (fbg1->Integral(mMinSBleft,xbmin)+fbg1->Integral(xbmax,mMaxSBright))/mD->GetBinWidth(1);
			if(npassedD>NminFlat) cout << "Using: " << fbg0->GetName() << " as fit function" << endl;
			else                  cout << "Using: " << fbg1->GetName() << " as fit function" << endl;
			cout << "Sidebands   fit-chi2(" << chi2  << ")/ndf("<<ndf<<")=" << chi2/ndf << " -> p-value=" << pvalue << endl;
			if(npassedD>NminFlat) cout << "SR integral [" << ximin << "," << ximax << "] ("<< ximax-ximin <<" MeV): " << fbg0->Integral(ximin,ximax) << endl;
			else                  cout << "SR integral [" << ximin << "," << ximax << "] ("<< ximax-ximin <<" MeV): " << fbg1->Integral(ximin,ximax) << endl;
			cout << "m3body bin width: " << mD->GetBinWidth(1) << endl;
			cout << "SB expectation: " << npassedD << "(count): " << iBkgSB << "(cut&fitSB) and " << iBkgSB1 << "(BDTshape)" << endl;
			cout << "SR expectation: " << iBkgSR << "(cut&fitSR) and " << iBkgSR1 << "(BDTshape)" << endl;
			cout << "HistFitter inputs:" << endl;
			cout << "  * NSR1 = " << iBkgSR1 << " +- " << diBkgSR1stat << "(stat.abs) +- " << diBkgSR1syst/iBkgSR1*100 << "\%(syst.rel)" << endl;
			cout << "  * Nobs = " << floor(iBkgSR1+0.5) << " Floor(NSR1+0.5)" << endl;
			TF1 *fleft = (npassedD>NminFlat) ? new TF1("fleft",fBg0,mMinSBleft,xbmin,npar0) : new TF1("fleft",fBg1,mMinSBleft,xbmin,npar1);
			if(npassedD>NminFlat) fleft->SetParameters(fbg0->GetParameters());
			else                  fleft->SetParameters(fbg1->GetParameters());
			if(hFit) mD->GetListOfFunctions()->Add(fleft);
			else     gD->GetListOfFunctions()->Add(fleft);
			gROOT->GetListOfFunctions()->Remove(fleft);
			TF1 *fright = (npassedD>NminFlat) ? new TF1("fright",fBg0,xbmax,mMaxSBright,npar0) : new TF1("fright",fBg1,xbmax,mMaxSBright,npar1);
			if(npassedD>NminFlat) fright->SetParameters(fbg0->GetParameters());
			else                  fright->SetParameters(fbg1->GetParameters());
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
        	
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if(iBkgSR<0.  || iBkgSB<0.)  { iBkgSR=0.; iBkgSB=0.; cout << "iBkgSR<0 or iBkgSB<0 ==> setting to 0" << endl; } ///////
			if(iBkgSR1<0. || iBkgSB1<0.) { iBkgSR1=0.; iBkgSB1=0.; diBkgSB1stat=0; diBkgSR1stat=0; cout << "iBkgSR1<0 or iBkgSB1<0 ==> setting to 0" << endl; } ///////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        	
			hnSB1Shape->SetBinContent(hbinx,hbiny,iBkgSB1);
			hnSB1Shape->SetBinError(hbinx,hbiny,diBkgSB1stat);
			hnSR1Shape->SetBinContent(hbinx,hbiny,iBkgSR1);
			hnSR1Shape->SetBinError(hbinx,hbiny,diBkgSR1stat);
	    	
	
			float AtimesErel = npassedS/nallS;
			float AtimesE = npassedS/ninitS;
			float Ntau = 2.41e8*(luminosity/20.28);
        	
			//// remake the python steering script for this evaluation
			Limit* limit = NULL;
			float dnSR1systUp = 1.+diBkgSR1syst/iBkgSR1;
			float dnSR1systDn = 1.-diBkgSR1syst/iBkgSR1;
			float ndata = floor(iBkgSR1+0.5);
			string pyfile = makeHistFitterPy(hfpylines,(string)score,(string)srsize,(string)nobs,iBkgSR1,diBkgSR1stat,dnSR1systUp,dnSR1systDn,ndata);
			limit  = getLimit(pyfile,score,srsize,nobs);
			float upper   = limit->upper;
			float median  = limit->median;
			float m1sigma = limit->m1sigma;
			float p1sigma = limit->p1sigma;
			float m2sigma = limit->m2sigma;
			float p2sigma = limit->p2sigma;
			bool ok       = limit->ok;
			bool isBad    = limit->isBad;
			delete limit;
			if(!ok || isBad) // something is wrong --> use the result from the previous iteration
			{
				int niterations = 0;
				float fraction  = 1;
				while((!ok || isBad) && niterations<11)
				{
					fraction = 1.-niterations*0.01;
					pyfile = makeHistFitterPy(hfpylines,(string)score,(string)srsize,(string)nobs,fraction*iBkgSR1,diBkgSR1stat,dnSR1systUp,dnSR1systDn,ndata);
					limit  = getLimit(pyfile,score,srsize,nobs);
					upper   = limit->upper;
					median  = limit->median;
					m1sigma = limit->m1sigma;
					p1sigma = limit->p1sigma;
					m2sigma = limit->m2sigma;
					p2sigma = limit->p2sigma;
					ok       = limit->ok;
					isBad    = limit->isBad;
					delete limit;
					niterations++;
					cout << "niterations=" << niterations << ": " << iBkgSR1 << "-->" << fraction*iBkgSR1 << endl;
				}
				if(ok && !isBad) cout << "converged after " << niterations << ": " << iBkgSR1 << "-->" << fraction*iBkgSR1 << endl;
			}
			
			
			float ndataObs = -1;
			float upperObs = -1;
			float AtimesEObs = -1;
			bool okObs = true;
			bool isBadObs = false;
			if(nobs!="X")
			{
				ndataObs = atoi(((string)nobs).c_str());
				string pyfileObs = makeHistFitterPy(hfpylines,(string)score,(string)srsize,(string)nobs,iBkgSR1,diBkgSR1stat,dnSR1systUp,dnSR1systDn,ndataObs);
				limit = getLimit(pyfileObs,score,srsize,nobs);
				upperObs   = limit->upper;
				AtimesEObs = AtimesE;
				okObs      = limit->ok;
				isBadObs   = limit->isBad;
				delete limit;
				if(!okObs || isBadObs)
				{
					int niterations = 0;
					float fraction  = 1; 
					while((!ok || isBad) && niterations<11)
					{
						fraction = 1.-niterations*0.01;
						ndataObs = atoi(((string)nobs).c_str());
						pyfileObs = makeHistFitterPy(hfpylines,(string)score,(string)srsize,(string)nobs,fraction*iBkgSR1,diBkgSR1stat,dnSR1systUp,dnSR1systDn,ndataObs);
						limit = getLimit(pyfileObs,score,srsize,nobs);
						upperObs   = limit->upper;
						AtimesEObs = AtimesE;
						okObs         = limit->ok;
						isBadObs      = limit->isBad;
						delete limit;
						niterations++;
						cout << "niterations=" << niterations << ": " << iBkgSR1 << "-->" << fraction*iBkgSR1 << endl;
					}
					if(okObs && !isBadObs) cout << "converged after " << niterations << ": " << iBkgSR1 << "-->" << fraction*iBkgSR1 << endl;
				}
			}
			cout << "upper    = " << upper    << " (for nObs=" << ndata << " <--floor(" << iBkgSR1 << "+0.5))" << endl;
			cout << "upperObs = " << upperObs << " (for nObs=" << ndataObs << " <--fixed)" << endl;

			
			if(ok && !isBad)
			{
				cout << "Expected (median) Upper limit on Br(tau->3mu)=" << median/AtimesE/Ntau << endl;
				
				hEffSrel->SetBinContent(hbinx,hbiny,AtimesErel);
				hEffS->SetBinContent(hbinx,hbiny,AtimesE);
				hNdat->SetBinContent(hbinx,hbiny,ndata);
			
				hnSR1NsObsExp->SetBinContent(hbinx,hbiny,upper);
				hnSR1NsMed->SetBinContent(hbinx,hbiny,median);
				hnSR1NsMedMin1->SetBinContent(hbinx,hbiny,m1sigma);
				hnSR1NsMedPls1->SetBinContent(hbinx,hbiny,p1sigma);
				hnSR1NsMedMin2->SetBinContent(hbinx,hbiny,m2sigma);
				hnSR1NsMedPls2->SetBinContent(hbinx,hbiny,p2sigma);
				
				hUL90med->SetBinContent(hbinx,hbiny,median/AtimesE/Ntau);
				hUL90min1->SetBinContent(hbinx,hbiny,m1sigma/AtimesE/Ntau);
				hUL90pls1->SetBinContent(hbinx,hbiny,p1sigma/AtimesE/Ntau);
				hUL90min2->SetBinContent(hbinx,hbiny,m2sigma/AtimesE/Ntau);
				hUL90pls2->SetBinContent(hbinx,hbiny,p2sigma/AtimesE/Ntau);
			}
			else cout << "Expected point: score=" << score << ", srsize=" << srsize << " was skipped" << endl;
			
			if(okObs && !isBadObs)
			{
				cout << "Observed (Nobs="<< Nobs << ") Upper limit on Br(tau->3mu)=" << upperObs/AtimesEObs/Ntau << endl;
				
				hnSR1NsObs->SetBinContent(hbinx,hbiny,upperObs);
				hUL90obs->SetBinContent(hbinx,hbiny,upperObs/AtimesEObs/Ntau);
				hNobs->SetBinContent(hbinx,hbiny,ndataObs);
			}
			else cout << "Observed point: score=" << score << ", srsize=" << srsize << " was skipped" << endl;
			
			
			//////////////////////////////////////////////////////
			//// regardless of a failed HistFitter procedure  ////
			//////////////////////////////////////////////////////
			
			hnSB1Count->SetBinContent(hbinx,hbiny,iBkgSB);
			hnSB1Count->SetBinError(hbinx,hbiny,sqrt(iBkgSB));
			hnSR1Count->SetBinContent(hbinx,hbiny,iBkgSR);
			hnSR1Count->SetBinError(hbinx,hbiny, (iBkgSB>0) ? (iBkgSR/iBkgSB)*sqrt(iBkgSB) : 1);
 
			
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
			// cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".m3body.fit."+score+"."+srsize+".eps");
			cnv->SaveAs(outfilesname+".pdf");
			
			tD->SetEventList(0);
			fOut->cd();
			mD->SetName("m3body_"+score+"."+srsize);
			mD->Write();
			olddir->cd();
			delete mD;
		}
	}
	
	

	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	BDTb->Draw("p0 e0");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSB0.eps");
	fOut->cd(); BDTb->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	BDTs->Draw("p0 e0");
	ptxt->Draw("same");
	cnv->Update();
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSig.eps");
	fOut->cd(); BDTs->Write(); olddir->cd();


	TBox* box;
	gStyle->SetPaintTextFormat("4.2f");


	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hNdat->SetMinimum(1.); hNdat->SetMarkerSize(0.7); hNdat->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NdatUsed.eps");
	fOut->cd(); hNdat->Write(); olddir->cd();

	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hNobs->SetMinimum(1.); hNobs->SetMarkerSize(0.7); hNobs->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NobsUsed.eps");
	fOut->cd(); hNobs->Write(); olddir->cd();


	////
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsObsExp->SetMinimum(1.); hnSR1NsObsExp->SetMarkerSize(0.7); hnSR1NsObsExp->Draw("colz text45");
	box = boxHighlight(hnSR1NsObsExp); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsULexp.observed.eps");
	fOut->cd(); hnSR1NsObsExp->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsObs->SetMinimum(1.); hnSR1NsObs->SetMarkerSize(0.7); hnSR1NsObs->Draw("colz text45");
	box = boxHighlight(hnSR1NsObs); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.observed.eps");
	fOut->cd(); hnSR1NsObs->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsMed->SetMinimum(1.); hnSR1NsMed->SetMarkerSize(0.7); hnSR1NsMed->Draw("colz text45");
	box = boxHighlight(hnSR1NsMed); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.median.eps");
	fOut->cd(); hnSR1NsMed->Write(); olddir->cd();

	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsMedMin1->SetMinimum(1.); hnSR1NsMedMin1->SetMarkerSize(0.7); hnSR1NsMedMin1->Draw("colz text45");
	box = boxHighlight(hnSR1NsMedMin1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.min1sigma.eps");
	fOut->cd(); hnSR1NsMedMin1->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsMedPls1->SetMinimum(1.); hnSR1NsMedPls1->SetMarkerSize(0.7); hnSR1NsMedPls1->Draw("colz text45");
	box = boxHighlight(hnSR1NsMedPls1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.pls1sigma.eps");
	fOut->cd(); hnSR1NsMedPls1->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsMedMin2->SetMinimum(1.); hnSR1NsMedMin2->SetMarkerSize(0.7); hnSR1NsMedMin2->Draw("colz text45");
	box = boxHighlight(hnSR1NsMedMin2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.min2sigma.eps");
	fOut->cd(); hnSR1NsMedMin2->Write(); olddir->cd();

	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1NsMedPls2->SetMinimum(1.); hnSR1NsMedPls2->SetMarkerSize(0.7); hnSR1NsMedPls2->Draw("colz text45");
	box = boxHighlight(hnSR1NsMedPls2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.pls2sigma.eps");
	fOut->cd(); hnSR1NsMedPls2->Write(); olddir->cd();
	

	////
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90obs->Scale(1./1.e-7); hUL90obs->SetMinimum(1.); hUL90obs->SetMarkerSize(0.7); hUL90obs->Draw("colz text45");
	box = boxHighlight(hUL90obs); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.observed.eps");
	fOut->cd(); hUL90obs->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90med->Scale(1./1.e-7); hUL90med->SetMinimum(1.); hUL90med->SetMarkerSize(0.7); hUL90med->Draw("colz text45");
	box = boxHighlight(hUL90med); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median.eps");
	fOut->cd(); hUL90med->Write(); olddir->cd();

	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90min1->Scale(1./1.e-7); hUL90min1->SetMinimum(1.); hUL90min1->SetMarkerSize(0.7);  hUL90min1->Draw("colz text45");
	box = boxHighlight(hUL90min1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.min1sigma.eps");
	fOut->cd(); hUL90min1->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90pls1->Scale(1./1.e-7); hUL90pls1->SetMinimum(1.); hUL90pls1->SetMarkerSize(0.7);  hUL90pls1->Draw("colz text45");
	box = boxHighlight(hUL90pls1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.pls1sigma.eps");
	fOut->cd(); hUL90pls1->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90min2->Scale(1./1.e-7); hUL90min2->SetMinimum(1.); hUL90min2->SetMarkerSize(0.7);  hUL90min2->Draw("colz text45");
	box = boxHighlight(hUL90min2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.min2sigma.eps");
	fOut->cd(); hUL90min2->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hUL90pls2->Scale(1./1.e-7); hUL90pls2->SetMinimum(1.); hUL90pls2->SetMarkerSize(0.7);  hUL90pls2->Draw("colz text45");
	box = boxHighlight(hUL90pls2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.pls2sigma.eps");
	fOut->cd(); hUL90pls2->Write(); olddir->cd();
	
	
	
	
	////////////
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hUL90medmin1 = subtract(hUL90med,hUL90min1); hUL90medmin1->SetName("-1sigmaBand");
	hUL90medmin1->SetMarkerSize(0.7); hUL90medmin1->SetTitle("-1#sigma band magnitude/10^{-7}"); hUL90medmin1->Draw("colz text45");
	box = boxHighlight(hUL90medmin1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median-1sigma.eps");
	fOut->cd(); hUL90medmin1->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hUL90medpls1 = subtract(hUL90med,hUL90pls1); hUL90medpls1->SetName("+1sigmaBand");
	hUL90medpls1->SetMarkerSize(0.7); hUL90medpls1->SetTitle("+1#sigma band magnitude/10^{-7}"); hUL90medpls1->Draw("colz text45");
	box = boxHighlight(hUL90medpls1); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median+1sigma.eps");
	fOut->cd(); hUL90medpls1->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hUL90medmin2 = subtract(hUL90med,hUL90min2); hUL90medmin2->SetName("-2sigmaBand");
	hUL90medmin2->SetMarkerSize(0.7); hUL90medmin2->SetTitle("-2#sigma band magnitude/10^{-7}"); hUL90medmin2->Draw("colz text45");
	box = boxHighlight(hUL90medmin2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median-2sigma.eps");
	fOut->cd(); hUL90medmin2->Write(); olddir->cd();
	
	delete box;
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hUL90medpls2 = subtract(hUL90med,hUL90pls2); hUL90medpls2->SetName("+2sigmaBand");
	hUL90medpls2->SetMarkerSize(0.7); hUL90medpls2->SetTitle("+2#sigma band magnitude/10^{-7}"); hUL90medpls2->Draw("colz text45");
	box = boxHighlight(hUL90medpls2); box->Draw("same");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median+2sigma.eps");
	fOut->cd(); hUL90medpls2->Write(); olddir->cd();
	
	
	
	//////////
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hEffS->Scale(100.); hEffS->SetMarkerSize(0.7); hEffS->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.efficiency.eps");
	fOut->cd(); hEffS->Write(); olddir->cd();
	
	

	/////////////
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSB1Shape->SetMarkerSize(0.7); hnSB1Shape->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSB1shape.eps");
	fOut->cd(); hnSB1Shape->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSB1Count->SetMarkerSize(0.7); hnSB1Count->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSB1count.eps");
	fOut->cd(); hnSB1Count->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hnSB1Ratio = (TH2F*)hnSB1Shape->Clone("SB1 shape/cut&count");
	hnSB1Ratio->Divide(hnSB1Count); hnSB1Ratio->SetTitle("SB1 shape/cut&count");
	hnSB1Ratio->SetMarkerSize(0.7); hnSB1Ratio->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSB1ratio.eps");
	fOut->cd(); hnSB1Ratio->Write(); olddir->cd();
	
	
	////////////
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1Shape->SetMarkerSize(0.7); hnSR1Shape->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSR1shape.eps");
	fOut->cd(); hnSR1Shape->Write(); olddir->cd();

	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	hnSR1Count->SetMarkerSize(0.7); hnSR1Count->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSR1count.eps");
	fOut->cd(); hnSR1Count->Write(); olddir->cd();
	
	delete cnv;
	cnv = new TCanvas("","",600,600);
	cnv->Draw();
	cnv->cd();
	TH2F* hnSR1Ratio = (TH2F*)hnSR1Shape->Clone("SR1 shape/cut&count");
	hnSR1Ratio->Divide(hnSR1Count); hnSR1Ratio->SetTitle("SR1 shape/cut&count");
	hnSR1Ratio->SetMarkerSize(0.7); hnSR1Ratio->Draw("colz text45");
	ptxt->Draw("same");
	cnv->SaveAs(outfilesname+".pdf");
	cnv->SaveAs("/tmp/hod/HistFitter/eps/"+outfilesname+".nSR1ratio.eps");
	fOut->cd(); hnSR1Ratio->Write(); olddir->cd();
	
	fOut->Write();
	fOut->Close();
	
	
	///////////
	delete cnv;
	cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->cd();
	cnv->Update();
	cnv->SaveAs(outfilesname+".pdf)");

	gSystem->Exec("/bin/cp -f "+outfilesname+".pdf                                                           /afs/cern.ch/user/h/hod/data/HistFitter/");
	gSystem->Exec("/bin/cp -f "+outfilesname+".root                                                          /afs/cern.ch/user/h/hod/data/HistFitter/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSB0.eps                       /afs/cern.ch/user/h/hod/data/HistFitter/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".BDTshapeSig.eps                       /afs/cern.ch/user/h/hod/data/HistFitter/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NobsUsed.eps           /afs/cern.ch/user/h/hod/data/HistFitter/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsULexp.observed.eps   /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.observed.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.median.eps        /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.min1sigma.eps     /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.pls1sigma.eps     /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.min2sigma.eps     /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.NsUL.pls2sigma.eps     /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.observed.eps           /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median.eps             /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.min1sigma.eps          /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.pls1sigma.eps          /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.min2sigma.eps          /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.pls2sigma.eps          /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median-1sigma.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median+1sigma.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median-2sigma.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.median+2sigma.eps      /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".optimization2D.efficiency.eps         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSB1shape.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSB1count.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSB1ratio.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSR1shape.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSR1count.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
	gSystem->Exec("/bin/cp -f /tmp/hod/HistFitter/eps/"+outfilesname+".nSR1ratio.eps                         /afs/cern.ch/user/h/hod/data/HistFitter/eps/");
}
