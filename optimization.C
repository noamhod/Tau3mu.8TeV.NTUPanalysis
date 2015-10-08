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
	
	
	
	// Example macro of using the TFeldmanCousins class in root.
	// get a FeldmanCousins calculation object with the default limits
	// of calculating a 90% CL with the minimum signal value scanned
	// = 0.0 and the maximum signal value scanned of 50.0
	//Author : Adrian John Bevan <bevan@SLAC.Stanford.EDU>

	if (!gROOT->GetClass("TFeldmanCousins")) gSystem->Load("libPhysics");
	TFeldmanCousins fc;

	// // calculate either the upper or lower limit for 10 observerd
	// // events with an estimated background of 3.  The calculation of
	// // either upper or lower limit will return that limit and fill
	// // data members with both the upper and lower limit for you.
	// Double_t Nobserved1   = 10.0;
	// Double_t Nbackground1 = 3.0;
	// Double_t ul = fc.CalculateUpperLimit(Nobserved1, Nbackground1);
	// Double_t ll = fc.GetLowerLimit();
	// cout << "For " <<  Nobserved1 << " data observed with and estimated background"<<endl;
	// cout << "of " << Nbackground1 << " candidates, the Feldman-Cousins method of "<<endl;
	// cout << "calculating confidence limits gives:"<<endl;
	// cout << "\tUpper Limit = " <<  ul << endl;
	// cout << "\tLower Limit = " <<  ll << endl;
	// cout << "at the 90% CL"<< endl;
 
	
	TFile* fout = new TFile("optimization.root","RECREATE");
	
	Int_t hNbins  = 50;
	Float_t hXmin = 0.899;
	Float_t hXmax = 0.999;
	
	TH1F* hEffS     = new TH1F("SignalEfficiency",  ";BDT score cut;Signal #it{A}#times#it{#epsilon}",  hNbins,hXmin,hXmax);
	TH1F* hUL90min1 = new TH1F("UpperLimitMin1Sig", "Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90pls1 = new TH1F("UpperLimitPls1Sig", "Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90min2 = new TH1F("UpperLimitMin2Sig", "Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90pls2 = new TH1F("UpperLimitPls2Sig", "Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);
	TH1F* hUL90med  = new TH1F("UpperLimitMedian",  "Upper limit on Br(#it{#tau#rightarrow3#mu}) at 90\% CLs;BDT score cut;Br(#it{#tau#rightarrow3#mu}) at 90\% CLs", hNbins,hXmin,hXmax);	
	
	hEffS->SetBinContent(1,  0.05782);  hUL90med->SetBinContent(1,  6.10E-07); hUL90min1->SetBinContent(1,  4.08E-07);	hUL90pls1->SetBinContent(1,  9.62E-07);	hUL90min2->SetBinContent(1,  2.91E-07);  hUL90pls2->SetBinContent(1,  1.50E-06);
	hEffS->SetBinContent(2,  0.05753);  hUL90med->SetBinContent(2,  6.04E-07); hUL90min1->SetBinContent(2,  4.03E-07);	hUL90pls1->SetBinContent(2,  9.50E-07);	hUL90min2->SetBinContent(2,  2.84E-07);  hUL90pls2->SetBinContent(2,  1.48E-06);
	hEffS->SetBinContent(3,  0.05735);  hUL90med->SetBinContent(3,  6.06E-07); hUL90min1->SetBinContent(3,  4.04E-07);	hUL90pls1->SetBinContent(3,  9.54E-07);	hUL90min2->SetBinContent(3,  2.85E-07);  hUL90pls2->SetBinContent(3,  1.48E-06);
	hEffS->SetBinContent(4,  0.05711);  hUL90med->SetBinContent(4,  6.05E-07); hUL90min1->SetBinContent(4,  4.03E-07);	hUL90pls1->SetBinContent(4,  9.52E-07);	hUL90min2->SetBinContent(4,  2.84E-07);  hUL90pls2->SetBinContent(4,  1.48E-06);
	hEffS->SetBinContent(5,  0.05691);  hUL90med->SetBinContent(5,  6.07E-07); hUL90min1->SetBinContent(5,  4.04E-07);	hUL90pls1->SetBinContent(5,  9.55E-07);	hUL90min2->SetBinContent(5,  2.85E-07);  hUL90pls2->SetBinContent(5,  1.48E-06);
	hEffS->SetBinContent(6,  0.05667);  hUL90med->SetBinContent(6,  6.03E-07); hUL90min1->SetBinContent(6,  3.98E-07);	hUL90pls1->SetBinContent(6,  9.53E-07);	hUL90min2->SetBinContent(6,  2.86E-07);  hUL90pls2->SetBinContent(6,  1.48E-06);
	hEffS->SetBinContent(7,  0.05639);  hUL90med->SetBinContent(7,  5.97E-07); hUL90min1->SetBinContent(7,  3.97E-07);	hUL90pls1->SetBinContent(7,  9.40E-07);	hUL90min2->SetBinContent(7,  2.79E-07);  hUL90pls2->SetBinContent(7,  1.46E-06);
	hEffS->SetBinContent(8,  0.05613);  hUL90med->SetBinContent(8,  5.89E-07); hUL90min1->SetBinContent(8,  3.93E-07);	hUL90pls1->SetBinContent(8,  9.30E-07);	hUL90min2->SetBinContent(8,  2.79E-07);  hUL90pls2->SetBinContent(8,  1.45E-06);
	hEffS->SetBinContent(9,  0.05591);  hUL90med->SetBinContent(9,  5.92E-07); hUL90min1->SetBinContent(9,  3.95E-07);	hUL90pls1->SetBinContent(9,  9.34E-07);	hUL90min2->SetBinContent(9,  2.79E-07);  hUL90pls2->SetBinContent(9,  1.46E-06);
	hEffS->SetBinContent(10, 0.05564);  hUL90med->SetBinContent(10, 5.89E-07); hUL90min1->SetBinContent(10, 3.90E-07);	hUL90pls1->SetBinContent(10, 9.26E-07);	hUL90min2->SetBinContent(10, 2.81E-07);  hUL90pls2->SetBinContent(10, 1.45E-06);
	hEffS->SetBinContent(11, 0.05532);  hUL90med->SetBinContent(11, 5.83E-07); hUL90min1->SetBinContent(11, 3.88E-07);	hUL90pls1->SetBinContent(11, 9.21E-07);	hUL90min2->SetBinContent(11, 2.73E-07);  hUL90pls2->SetBinContent(11, 1.44E-06);
	hEffS->SetBinContent(12, 0.05506);  hUL90med->SetBinContent(12, 5.83E-07); hUL90min1->SetBinContent(12, 3.87E-07);	hUL90pls1->SetBinContent(12, 9.19E-07);	hUL90min2->SetBinContent(12, 2.70E-07);  hUL90pls2->SetBinContent(12, 1.43E-06);
	hEffS->SetBinContent(13, 0.05470);  hUL90med->SetBinContent(13, 5.87E-07); hUL90min1->SetBinContent(13, 3.89E-07);	hUL90pls1->SetBinContent(13, 9.25E-07);	hUL90min2->SetBinContent(13, 2.71E-07);  hUL90pls2->SetBinContent(13, 1.44E-06);
	hEffS->SetBinContent(14, 0.05439);  hUL90med->SetBinContent(14, 5.72E-07); hUL90min1->SetBinContent(14, 3.80E-07);	hUL90pls1->SetBinContent(14, 9.05E-07);	hUL90min2->SetBinContent(14, 2.67E-07);  hUL90pls2->SetBinContent(14, 1.42E-06);
	hEffS->SetBinContent(15, 0.05415);  hUL90med->SetBinContent(15, 5.67E-07); hUL90min1->SetBinContent(15, 3.75E-07);	hUL90pls1->SetBinContent(15, 8.96E-07);	hUL90min2->SetBinContent(15, 2.62E-07);  hUL90pls2->SetBinContent(15, 1.40E-06);
	hEffS->SetBinContent(16, 0.05378);  hUL90med->SetBinContent(16, 5.67E-07); hUL90min1->SetBinContent(16, 3.74E-07);	hUL90pls1->SetBinContent(16, 8.95E-07);	hUL90min2->SetBinContent(16, 2.60E-07);  hUL90pls2->SetBinContent(16, 1.40E-06);
	hEffS->SetBinContent(17, 0.05335);  hUL90med->SetBinContent(17, 5.71E-07); hUL90min1->SetBinContent(17, 3.77E-07);	hUL90pls1->SetBinContent(17, 9.02E-07);	hUL90min2->SetBinContent(17, 2.63E-07);  hUL90pls2->SetBinContent(17, 1.41E-06);
	hEffS->SetBinContent(18, 0.05289);  hUL90med->SetBinContent(18, 5.63E-07); hUL90min1->SetBinContent(18, 3.73E-07);	hUL90pls1->SetBinContent(18, 8.87E-07);	hUL90min2->SetBinContent(18, 2.68E-07);  hUL90pls2->SetBinContent(18, 1.39E-06);
	hEffS->SetBinContent(19, 0.05244);  hUL90med->SetBinContent(19, 5.57E-07); hUL90min1->SetBinContent(19, 3.69E-07);	hUL90pls1->SetBinContent(19, 8.84E-07);	hUL90min2->SetBinContent(19, 2.60E-07);  hUL90pls2->SetBinContent(19, 1.39E-06);
	hEffS->SetBinContent(20, 0.05197);  hUL90med->SetBinContent(20, 5.59E-07); hUL90min1->SetBinContent(20, 3.69E-07);	hUL90pls1->SetBinContent(20, 8.85E-07);	hUL90min2->SetBinContent(20, 2.56E-07);  hUL90pls2->SetBinContent(20, 1.39E-06);
	hEffS->SetBinContent(21, 0.05157);  hUL90med->SetBinContent(21, 5.53E-07); hUL90min1->SetBinContent(21, 3.63E-07);	hUL90pls1->SetBinContent(21, 8.74E-07);	hUL90min2->SetBinContent(21, 2.56E-07);  hUL90pls2->SetBinContent(21, 1.38E-06);
	hEffS->SetBinContent(22, 0.05116);  hUL90med->SetBinContent(22, 5.42E-07); hUL90min1->SetBinContent(22, 3.58E-07);	hUL90pls1->SetBinContent(22, 8.62E-07);	hUL90min2->SetBinContent(22, 2.51E-07);  hUL90pls2->SetBinContent(22, 1.36E-06);
	hEffS->SetBinContent(23, 0.05062);  hUL90med->SetBinContent(23, 5.44E-07); hUL90min1->SetBinContent(23, 3.58E-07);	hUL90pls1->SetBinContent(23, 8.64E-07);	hUL90min2->SetBinContent(23, 2.48E-07);  hUL90pls2->SetBinContent(23, 1.36E-06);
	hEffS->SetBinContent(24, 0.05013);  hUL90med->SetBinContent(24, 5.28E-07); hUL90min1->SetBinContent(24, 3.46E-07);	hUL90pls1->SetBinContent(24, 8.40E-07);	hUL90min2->SetBinContent(24, 2.46E-07);  hUL90pls2->SetBinContent(24, 1.33E-06);
	hEffS->SetBinContent(25, 0.04956);  hUL90med->SetBinContent(25, 5.23E-07); hUL90min1->SetBinContent(25, 3.45E-07);	hUL90pls1->SetBinContent(25, 8.35E-07);	hUL90min2->SetBinContent(25, 2.40E-07);  hUL90pls2->SetBinContent(25, 1.32E-06);
	hEffS->SetBinContent(26, 0.04905);  hUL90med->SetBinContent(26, 5.24E-07); hUL90min1->SetBinContent(26, 3.43E-07);	hUL90pls1->SetBinContent(26, 8.36E-07);	hUL90min2->SetBinContent(26, 2.40E-07);  hUL90pls2->SetBinContent(26, 1.32E-06);
	hEffS->SetBinContent(27, 0.04850);  hUL90med->SetBinContent(27, 5.08E-07); hUL90min1->SetBinContent(27, 3.34E-07);	hUL90pls1->SetBinContent(27, 8.14E-07);	hUL90min2->SetBinContent(27, 2.31E-07);  hUL90pls2->SetBinContent(27, 1.29E-06);
	hEffS->SetBinContent(28, 0.04785);  hUL90med->SetBinContent(28, 5.15E-07); hUL90min1->SetBinContent(28, 3.38E-07);	hUL90pls1->SetBinContent(28, 8.24E-07);	hUL90min2->SetBinContent(28, 2.33E-07);  hUL90pls2->SetBinContent(28, 1.31E-06);
	hEffS->SetBinContent(29, 0.04735);  hUL90med->SetBinContent(29, 5.15E-07); hUL90min1->SetBinContent(29, 3.37E-07);	hUL90pls1->SetBinContent(29, 8.25E-07);	hUL90min2->SetBinContent(29, 2.33E-07);  hUL90pls2->SetBinContent(29, 1.31E-06);
	hEffS->SetBinContent(30, 0.04671);  hUL90med->SetBinContent(30, 5.23E-07); hUL90min1->SetBinContent(30, 3.42E-07);	hUL90pls1->SetBinContent(30, 8.37E-07);	hUL90min2->SetBinContent(30, 2.37E-07);  hUL90pls2->SetBinContent(30, 1.33E-06);
	hEffS->SetBinContent(31, 0.04587);  hUL90med->SetBinContent(31, 5.14E-07); hUL90min1->SetBinContent(31, 3.35E-07);	hUL90pls1->SetBinContent(31, 8.24E-07);	hUL90min2->SetBinContent(31, 2.33E-07);  hUL90pls2->SetBinContent(31, 1.31E-06);
	hEffS->SetBinContent(32, 0.04518);  hUL90med->SetBinContent(32, 5.16E-07); hUL90min1->SetBinContent(32, 3.39E-07);	hUL90pls1->SetBinContent(32, 8.25E-07);	hUL90min2->SetBinContent(32, 2.42E-07);  hUL90pls2->SetBinContent(32, 1.32E-06);
	hEffS->SetBinContent(33, 0.04435);  hUL90med->SetBinContent(33, 5.19E-07); hUL90min1->SetBinContent(33, 3.39E-07);	hUL90pls1->SetBinContent(33, 8.35E-07);	hUL90min2->SetBinContent(33, 2.34E-07);  hUL90pls2->SetBinContent(33, 1.33E-06);
	hEffS->SetBinContent(34, 0.04345);  hUL90med->SetBinContent(34, 5.10E-07); hUL90min1->SetBinContent(34, 3.30E-07);	hUL90pls1->SetBinContent(34, 8.22E-07);	hUL90min2->SetBinContent(34, 2.31E-07);  hUL90pls2->SetBinContent(34, 1.31E-06);
	hEffS->SetBinContent(35, 0.04256);  hUL90med->SetBinContent(35, 5.03E-07); hUL90min1->SetBinContent(35, 3.28E-07);	hUL90pls1->SetBinContent(35, 8.18E-07);	hUL90min2->SetBinContent(35, 2.33E-07);  hUL90pls2->SetBinContent(35, 1.31E-06);
	hEffS->SetBinContent(36, 0.04148);  hUL90med->SetBinContent(36, 4.90E-07); hUL90min1->SetBinContent(36, 3.16E-07);	hUL90pls1->SetBinContent(36, 7.96E-07);	hUL90min2->SetBinContent(36, 2.20E-07);  hUL90pls2->SetBinContent(36, 1.28E-06);
	hEffS->SetBinContent(37, 0.04033);  hUL90med->SetBinContent(37, 4.80E-07); hUL90min1->SetBinContent(37, 3.10E-07);	hUL90pls1->SetBinContent(37, 7.82E-07);	hUL90min2->SetBinContent(37, 2.14E-07);  hUL90pls2->SetBinContent(37, 1.26E-06);
	hEffS->SetBinContent(38, 0.03929);  hUL90med->SetBinContent(38, 4.69E-07); hUL90min1->SetBinContent(38, 2.99E-07);	hUL90pls1->SetBinContent(38, 7.69E-07);	hUL90min2->SetBinContent(38, 2.10E-07);  hUL90pls2->SetBinContent(38, 1.25E-06);
	hEffS->SetBinContent(39, 0.03804);  hUL90med->SetBinContent(39, 4.48E-07); hUL90min1->SetBinContent(39, 2.89E-07);	hUL90pls1->SetBinContent(39, 7.36E-07);	hUL90min2->SetBinContent(39, 2.03E-07);  hUL90pls2->SetBinContent(39, 1.21E-06);
	hEffS->SetBinContent(40, 0.03655);  hUL90med->SetBinContent(40, 4.35E-07); hUL90min1->SetBinContent(40, 2.74E-07);	hUL90pls1->SetBinContent(40, 7.24E-07);	hUL90min2->SetBinContent(40, 1.92E-07);  hUL90pls2->SetBinContent(40, 1.19E-06);
	hEffS->SetBinContent(41, 0.03491);  hUL90med->SetBinContent(41, 4.19E-07); hUL90min1->SetBinContent(41, 2.66E-07);	hUL90pls1->SetBinContent(41, 7.07E-07);	hUL90min2->SetBinContent(41, 1.85E-07);  hUL90pls2->SetBinContent(41, 1.17E-06);
	hEffS->SetBinContent(42, 0.03311);  hUL90med->SetBinContent(42, 4.30E-07); hUL90min1->SetBinContent(42, 2.71E-07);	hUL90pls1->SetBinContent(42, 7.29E-07);	hUL90min2->SetBinContent(42, 1.89E-07);  hUL90pls2->SetBinContent(42, 1.22E-06);
	hEffS->SetBinContent(43, 0.03110);  hUL90med->SetBinContent(43, 4.28E-07); hUL90min1->SetBinContent(43, 2.69E-07);	hUL90pls1->SetBinContent(43, 7.33E-07);	hUL90min2->SetBinContent(43, 1.87E-07);  hUL90pls2->SetBinContent(43, 1.23E-06);
	hEffS->SetBinContent(44, 0.02878);  hUL90med->SetBinContent(44, 4.13E-07); hUL90min1->SetBinContent(44, 2.55E-07);	hUL90pls1->SetBinContent(44, 7.14E-07);	hUL90min2->SetBinContent(44, 1.74E-07);  hUL90pls2->SetBinContent(44, 1.22E-06);
	hEffS->SetBinContent(45, 0.02585);  hUL90med->SetBinContent(45, 4.40E-07); hUL90min1->SetBinContent(45, 2.70E-07);	hUL90pls1->SetBinContent(45, 7.67E-07);	hUL90min2->SetBinContent(45, 1.83E-07);  hUL90pls2->SetBinContent(45, 1.33E-06);
	hEffS->SetBinContent(46, 0.02265);  hUL90med->SetBinContent(46, 3.75E-07); hUL90min1->SetBinContent(46, 2.16E-07);	hUL90pls1->SetBinContent(46, 7.00E-07);	hUL90min2->SetBinContent(46, 1.41E-07);  hUL90pls2->SetBinContent(46, 1.27E-06);
	hEffS->SetBinContent(47, 0.01886);  hUL90med->SetBinContent(47, 4.50E-07); hUL90min1->SetBinContent(47, 2.59E-07);	hUL90pls1->SetBinContent(47, 8.41E-07);	hUL90min2->SetBinContent(47, 1.69E-07);  hUL90pls2->SetBinContent(47, 1.53E-06);
	hEffS->SetBinContent(48, 0.01409);  hUL90med->SetBinContent(48, 5.99E-07); hUL90min1->SetBinContent(48, 3.45E-07);	hUL90pls1->SetBinContent(48, 1.12E-06);	hUL90min2->SetBinContent(48, 2.24E-07);  hUL90pls2->SetBinContent(48, 2.03E-06);
	hEffS->SetBinContent(49, 0.00889);  hUL90med->SetBinContent(49, 9.42E-07); hUL90min1->SetBinContent(49, 5.42E-07);	hUL90pls1->SetBinContent(49, 1.77E-06);	hUL90min2->SetBinContent(49, 3.53E-07);  hUL90pls2->SetBinContent(49, 3.21E-06);
	hEffS->SetBinContent(50, 0.00266);  hUL90med->SetBinContent(50, 2.18E-06); hUL90min1->SetBinContent(50, 1.05E-06);	hUL90pls1->SetBinContent(50, 4.61E-06);	hUL90min2->SetBinContent(50, 5.79E-07);  hUL90pls2->SetBinContent(50, 9.08E-06);
	
	
	TGraphAsymmErrors* g1sigma = new TGraphAsymmErrors();
	for(Int_t i=1 ; i<=hUL90med->GetNbinsX() ; i++) 
	{
		float x         = hUL90med->GetBinCenter(i);
		float xErrLeft  = hUL90med->GetBinLowEdge(i);
		float xErrRight = hUL90med->GetBinLowEdge(i)+hUL90med->GetBinWidth(i);
    
		float y       = hUL90med->GetBinContent(i);
		float yErrDwn = fabs(y-hUL90min1->GetBinContent(i));
		float yErrUp  = fabs(y-hUL90pls1->GetBinContent(i));

		g1sigma->SetPoint(i-1, x, y);
		// g1sigma->SetPointError(i-1, xErrLeft,xErrRight, yErrDwn,yErrUp);
		g1sigma->SetPointError(i-1, 0,0, yErrDwn,yErrUp);
	}
	TGraphAsymmErrors* g2sigma = new TGraphAsymmErrors();
	for(Int_t i=1 ; i<=hUL90med->GetNbinsX() ; i++) 
	{
		float x         = hUL90med->GetBinCenter(i);
		float xErrLeft  = hUL90med->GetBinLowEdge(i);
		float xErrRight = hUL90med->GetBinLowEdge(i)+hUL90med->GetBinWidth(i);
    
		float y       = hUL90med->GetBinContent(i);
		float yErrDwn = fabs(y-hUL90min2->GetBinContent(i));
		float yErrUp  = fabs(y-hUL90pls2->GetBinContent(i));
    
		g2sigma->SetPoint(i-1, x, y);
		// g2sigma->SetPointError(i-1, xErrLeft,xErrRight, yErrDwn,yErrUp);
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
	mgBands->SetMaximum(0.7*hUL90med->GetMaximum());
	
	hEffS->SetLineColor(kRed);
	// hEffB->SetLineColor(kRed);
	hEffS->SetLineStyle(1);
	hEffS->SetLineWidth(2);
	// hEffB->SetLineStyle(2);
	hEffS->SetTitle("");
	// hEffB->SetTitle("");

	Double_t Ymin = +1.e20;
	Double_t Xmin = -999;
	
	fout->cd();
	hEffS->Write();
	// hEffB->Write();
	hUL90med->Write();
	g1sigma->Write();
	g2sigma->Write();
	fout->Write();
	// fout->Close();
	
	
	TLegend* leg = new TLegend(0.11,0.48,0.36,0.72);
	leg->SetFillStyle(4000); //will be transparent
	// leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(hUL90med,"#scale[0.7]{#frac{N_{UL}^{med}/#it{#epsilon}#times#it{A}}{#it{k}_{F}#timesN_{#it{#tau}}}}","l");
	leg->AddEntry(g1sigma,"1#sigma","f");
	leg->AddEntry(g2sigma,"2#sigma","f");
	leg->AddEntry(hEffS,"#it{#epsilon}_{S}#times#it{A}_{S}","l");
	
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
	mgBands->SetTitle(hUL90med->GetTitle());
	mgBands->GetXaxis()->SetTitle(hUL90med->GetXaxis()->GetTitle());
	mgBands->GetYaxis()->SetTitle(hUL90med->GetYaxis()->GetTitle());
	mgBands->GetXaxis()->SetLimits(hXmin,hXmax);
	p1->Modified();
	p1->RedrawAxis();
	c1->cd();
	
	//compute the pad range with suitable margins
	Double_t ymin = 0;
	Double_t ymax = 0.07;
	Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
	Double_t xmin = hUL90med->GetXaxis()->GetXmin();
	Double_t xmax = hUL90med->GetXaxis()->GetXmax();
	Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
	p2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
	p2->Draw();
	p2->cd();
	leg->Draw("same");
	hEffS->Draw("][L sames");
	// hEffB->Draw("][hist sames");
	// p2->RedrawAxis();

	// draw axis on the right side of the pad
	TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
	axis->SetTitle("Efficiency");
	axis->SetLabelColor(kBlack);
	axis->SetLabelFont(font);
	axis->SetLabelSize(tsize);
	axis->SetTextFont(font);
	axis->SetTextSize(tsize);
	axis->Draw();
	
	c1->cd();
	// c1->Update();
	c1->SaveAs("figures.optimization.pdf");
	c1->SaveAs("figures.optimization.eps");
	c1->SaveAs("figures.optimization.png");
	c1->SaveAs("figures.optimization.root");
	// cout << "Optimal BDT score is " << Xmin << " for Sqrt(effB)/effS=" << Ymin << endl;
	// cout << "and efficiency of " << hEffS->GetBinContent(hEffS->FindBin(Xmin)) << " and rejection of " << 1-hEffB->GetBinContent(hEffB->FindBin(Xmin)) << endl;
}