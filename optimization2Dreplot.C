#include "std.h"
#include "const.h"

TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}

TH2F* Multiply(TH2F* h1, TH2F* h2, float c)
{
	TH2F* h = (TH2F*)h1->Clone("tmp");
	h->Reset();
	for(Int_t bx=1 ; bx<=h->GetNbinsX() ; bx++)
	{
		for(Int_t by=1 ; by<=h->GetNbinsY() ; by++)
		{
			h->SetBinContent(bx,by,h1->GetBinContent(bx,by)*h2->GetBinContent(bx,by)*c);
		}
	}
	return h;
}

TPaveText* ptxt;
void makeAtlasLabel()
{
	ptxt = new TPaveText(0.12,0.73,0.37,0.90,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

TLegend* legN;
TLegend* legP;
TLegend* legF;
void makeLegend()
{
	legN = new TLegend(0.675,0.75,0.845,0.89,NULL,"brNDC");
	// legN->SetFillStyle(4000); //will be transparent
	legN->SetFillColor(0);
	legN->SetTextFont(42);
	legN->SetBorderSize(0);
	
	legP = new TLegend(0.4,0.75,0.58,0.89,NULL,"brNDC");
	// legP->SetFillStyle(4000); //will be transparent
	legP->SetFillColor(0);
	legP->SetTextFont(42);
	legP->SetBorderSize(0);
	
	legF = new TLegend(0.4,0.82,0.57,0.89,NULL,"brNDC");
	// legF->SetFillStyle(4000); //will be transparent
	legF->SetFillColor(0);
	legF->SetTextFont(42);
	legF->SetBorderSize(0);
}

TBox* boxHighlight(Color_t col, Int_t linestyle, TH2F* h, Double_t xcenter=-1e10, Double_t ycenter=-1e10)
{
	Int_t binx, biny, binz, bin;
	if(xcenter<0 || ycenter<0) bin = h->GetMinimumBin();
	else                       bin = h->FindBin(xcenter,ycenter);
	h->GetBinXYZ(bin, binx, biny, binz);
	Double_t x1 = h->GetXaxis()->GetBinLowEdge(binx);
	Double_t x2 = h->GetXaxis()->GetBinUpEdge(binx);
	Double_t y1 = h->GetYaxis()->GetBinLowEdge(biny);
	Double_t y2 = h->GetYaxis()->GetBinUpEdge(biny);
	TBox* box = new TBox(x1,y1,x2,y2);
	box->SetFillStyle(0);
	box->SetLineWidth(3);
	box->SetLineStyle(linestyle);
	box->SetLineColor(col);
	return box;
}

TPaletteAxis* palette;
void fixPalette(TH2F* h)
{
	palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.855); palette->SetX2NDC(0.905); palette->SetY1NDC(0.1); palette->SetY2NDC(0.9);
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

void fixFailedFits(TH2F* h, int radius)
{
	for(Int_t bx=1 ; bx<=h->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=h->GetNbinsY() ; ++by)
		{
			float z = h->GetBinContent(bx,by);
			if(z>0.) continue;
			z = getAverage(h,bx,by,radius);
			h->SetBinContent(bx,by,z);
		}
	}
}


int factorial(int n) { return (n==1 || n==0) ? 1 : factorial(n-1)*n; }
TH2F* poissonize(TH2F* h, TH2F* eff)
{
	TH2F* hPoisson = (TH2F*)h->Clone("poisson");
	hPoisson->Reset();
	if(eff!=NULL) hPoisson->SetTitle("Poisson(N_{obs}=nearest_integer(N_{SR1}) | N_{SR1})/#it{A}#times#epsilon");
	else          hPoisson->SetTitle("Poisson(N_{obs}=nearest_integer(N_{SR1}) | N_{SR1})");
	if(eff!=NULL) hPoisson->GetZaxis()->SetTitle("Poisson/#it{A}#times#epsilon");
	else          hPoisson->GetZaxis()->SetTitle("Poisson");
	for(Int_t bx=1 ; bx<=h->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=h->GetNbinsY() ; ++by)
		{
			float z = h->GetBinContent(bx,by);
			float p = 0;
			int n = nearbyint(z);
			// for(int i=0 ; i<=n ; ++i) p += exp(-z)*pow(z,i)/factorial(i);
			p = exp(-z)*pow(z,n)/factorial(n);
			hPoisson->SetBinContent(bx,by,(eff!=NULL) ? p/eff->GetBinContent(bx,by) : p);
		}
	}
	return hPoisson;
}

TH2F* getGrid(TH2F* h, float maxreldiff, float zsize, Color_t linecolor=kBlack, Int_t linestyle=1)
{
	Int_t bin = h->GetMinimumBin();
	Int_t binx, biny, binz;
	h->GetBinXYZ(bin, binx, biny, binz);
	float minimum = h->GetBinContent(binx,biny);
	TH2F* hGrid = (TH2F*)h->Clone("tmp");
	hGrid->Reset();
	hGrid->SetLineStyle(linestyle);
	hGrid->SetLineColor(linecolor);
	hGrid->SetLineWidth(1);
	for(Int_t bx=1 ; bx<=h->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=h->GetNbinsY() ; ++by)
		{
			float z = h->GetBinContent(bx,by);
			float reldif = (z-minimum)/minimum;
			if(reldif<=maxreldiff) { hGrid->SetBinContent(bx,by,zsize); }
		}
	}
	return hGrid;
}

TH2F* getSigma(TH2F* hPls, TH2F* hMin)
{
	TH2F* hGrid = (TH2F*)hPls->Clone("tmp");
	hGrid->Reset();
	hGrid->SetLineColor(kBlack);
	hGrid->SetLineWidth(1);
	for(Int_t bx=1 ; bx<=hPls->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=hPls->GetNbinsY() ; ++by)
		{	
			float pls = hPls->GetBinContent(bx,by);
			float min = hMin->GetBinContent(bx,by);
			hGrid->SetBinContent(bx,by,pls-min);
		}
	}
	return hGrid;
}

void optimization2Dreplot()
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
	
	gStyle->SetPaintTextFormat("4.2f");
	
	makeAtlasLabel();
	makeLegend();


	TFile* f0 = new TFile("Optimization2D.1450-1690-1870-2290.pos094-pos0990.80-130.Nobs0.root", "READ");
	TFile* f1 = new TFile("Optimization2D.1450-1690-1870-2290.pos094-pos0990.80-130.Nobs1.root", "READ");
	TFile* f2 = new TFile("Optimization2D.1450-1690-1870-2290.pos094-pos0990.80-130.Nobs2.root", "READ");
	TFile* f3 = new TFile("Optimization2D.1450-1690-1870-2290.pos094-pos0990.80-130.Nobs3.root", "READ");
	TFile* fX = new TFile("Optimization2D.1450-1690-1870-2290.pos094-pos0995.44-154.NobsX.root", "READ");
	
	// TString inname = "Optimization2D.1450-1870-1870-2290.pos090-pos0995.20-180.root"; // OLD
	// TString inname = "Optimization2D.1450-1690-1870-2290.pos094-pos0995.50-150.root"; // OLD
	TString inname = "Optimization2D.1450-1690-1870-2290.pos094-pos0990.80-130.Nobs0.root";

	TFile* f = new TFile(inname, "READ");
	TString outname = "figures/"+inname;
	outname.ReplaceAll(".root","");

	TH2F* hNSR1 = (TH2F*)f0->Get("SR1 Shape");         TH2F* hNSR1_full = (TH2F*)fX->Get("SR1 Shape");
	TH2F* hNSB1 = (TH2F*)f0->Get("SB1 Shape");         TH2F* hNSB1_full = (TH2F*)fX->Get("SB1 Shape");
	TH2F* hEff  = (TH2F*)f0->Get("SignalEfficiency");  TH2F* hEff_full  = (TH2F*)fX->Get("SignalEfficiency");
	TH2F* hDat  = (TH2F*)f0->Get("Ndat");              TH2F* hDat_full  = (TH2F*)fX->Get("Ndat");
	
	TH2F* hNsUL_obsEQ0 = (TH2F*)f0->Get("SR1 NsObs");  TH2F* hObs_obsEQ0 = (TH2F*)f0->Get("UpperLimitObserved");
	TH2F* hNsUL_obsEQ1 = (TH2F*)f1->Get("SR1 NsObs");  TH2F* hObs_obsEQ1 = (TH2F*)f1->Get("UpperLimitObserved");
	TH2F* hNsUL_obsEQ2 = (TH2F*)f2->Get("SR1 NsObs");  TH2F* hObs_obsEQ2 = (TH2F*)f2->Get("UpperLimitObserved");
	TH2F* hNsUL_obsEQ3 = (TH2F*)f3->Get("SR1 NsObs");  TH2F* hObs_obsEQ3 = (TH2F*)f3->Get("UpperLimitObserved");

	TH2F* hNsUL_exp  = (TH2F*)f0->Get("SR1 NsObsExp"); TH2F* hNsUL_exp_full = (TH2F*)fX->Get("SR1 NsObsExp");
	TH2F* hNsUL_med  = (TH2F*)f0->Get("SR1 NsMed");    TH2F* hNsUL_med_full = (TH2F*)fX->Get("SR1 NsMed");
	TH2F* hNsUL_min1 = (TH2F*)f0->Get("SR1 NsMedMin1");
	TH2F* hNsUL_pls1 = (TH2F*)f0->Get("SR1 NsMedPls1");
	TH2F* hNsUL_min2 = (TH2F*)f0->Get("SR1 NsMedMin2");
	TH2F* hNsUL_pls2 = (TH2F*)f0->Get("SR1 NsMedPls2");
	
	TH2F* hMed  = (TH2F*)f0->Get("UpperLimitMedian");  TH2F* hMed_full = (TH2F*)fX->Get("UpperLimitMedian");
	TH2F* hMin1 = (TH2F*)f0->Get("UpperLimitMin1Sig");
	TH2F* hPls1 = (TH2F*)f0->Get("UpperLimitPls1Sig");
	TH2F* hMin2 = (TH2F*)f0->Get("UpperLimitMin2Sig");
	TH2F* hPls2 = (TH2F*)f0->Get("UpperLimitPls2Sig");
	
	TString title  = hMed->GetTitle();
	TString ztitle = hMed->GetZaxis()->GetTitle();
	title.ReplaceAll("Expected (median) ","#Delta["); title.ReplaceAll(" at 90\% CLs","]");
	ztitle.ReplaceAll("Br","#Delta[Br");              ztitle.ReplaceAll(" at 90\% CLs","]");
	
	hNsUL_exp     ->SetTitle("Observed (N_{obs}=floor(N_{SR1}+0.5)) N_{s} at 90\% CLs");
	hNsUL_exp_full->SetTitle("Observed (N_{obs}=floor(N_{SR1}+0.5)) N_{s} at 90\% CLs");
	
	hNSB1     ->SetTitle("N_{SB} (exp. sidebands from BDT shape)");
	hNSB1     ->GetZaxis()->SetTitle("N_{SB} (exp. sidebands from BDT shape)");
	hNSB1_full->SetTitle("N_{SB} (exp. sidebands from BDT shape)");
	hNSB1_full->GetZaxis()->SetTitle("N_{SB} (exp. sidebands from BDT shape)");
	
	hNSR1     ->SetTitle("N_{SR} (exp. bkg. from BDT shape)");
	hNSR1     ->GetZaxis()->SetTitle("N_{SR} (exp. bkg. from BDT shape)");
	hNSR1_full->SetTitle("N_{SR} (exp. bkg. from BDT shape)");
	hNSR1_full->GetZaxis()->SetTitle("N_{SR} (exp. bkg. from BDT shape)");
	
	TH2F* hMin1SigmaBand = (TH2F*)f0->Get("-1sigmaBand");  hMin1SigmaBand->SetTitle("+1#sigma band width "+title); hMin1SigmaBand->GetZaxis()->SetTitle(ztitle);
	TH2F* hPls1SigmaBand = (TH2F*)f0->Get("+1sigmaBand");  hPls1SigmaBand->SetTitle("-1#sigma band width "+title); hPls1SigmaBand->GetZaxis()->SetTitle(ztitle);
	TH2F* hMin2SigmaBand = (TH2F*)f0->Get("-2sigmaBand");  hMin2SigmaBand->SetTitle("+2#sigma band width "+title); hMin2SigmaBand->GetZaxis()->SetTitle(ztitle);
	TH2F* hPls2SigmaBand = (TH2F*)f0->Get("+2sigmaBand");  hPls2SigmaBand->SetTitle("-2#sigma band width "+title); hPls2SigmaBand->GetZaxis()->SetTitle(ztitle);
	
	
	
	TH2F* h1Ssigma = (TH2F*)getSigma(hPls1,hMin1);   h1Ssigma->SetName("1SigmaBand");    h1Ssigma->SetTitle("1#sigma band width "+title);    h1Ssigma->GetZaxis()->SetTitle(ztitle);
	TH2F* h2Ssigma = (TH2F*)getSigma(hPls2,hMin2);   h2Ssigma->SetName("2SigmaBand");    h2Ssigma->SetTitle("2#sigma band width "+title);    h2Ssigma->GetZaxis()->SetTitle(ztitle);
	
	
	TH2F* hPoisson    = (TH2F*)poissonize(hNSR1,NULL);
	TH2F* hEffPoisson = (TH2F*)poissonize(hNSR1,hEff);
	
	
	float markersize = 0.35;
	hNSR1_full->SetMarkerSize(markersize);
	hNSB1_full->SetMarkerSize(markersize);
	hEff_full->SetMarkerSize(markersize);
	hDat_full->SetMarkerSize(markersize);
	hNsUL_exp_full->SetMarkerSize(markersize);
	hNsUL_med_full->SetMarkerSize(markersize);
	hMed_full->SetMarkerSize(markersize);

	
	// float medMax = 4;
	// hMed ->SetMaximum(medMax);
	// hMin1->SetMaximum(medMax-1);
	//     hPls1->SetMaximum(medMax+1);
	//     hMin2->SetMaximum(medMax-2);
	//     // hPls2->SetMaximum(medMax+2);
	// 
	// float bandMax = 5;
	// h1Ssigma  ->SetMaximum(bandMax);
	//     hPls1SigmaBand->SetMaximum(bandMax-2);
	//     hMin1SigmaBand->SetMaximum(bandMax-2);
	//     // h2Ssigma  ->SetMaximum(5);
	//     // hPls2SigmaBand->SetMaximum(5);
	//     // hMin2SigmaBand->SetMaximum(5);
	// 
	// float medMax_full = 8;
	// hMed_full->SetMaximum(medMax_full);


	int radius = 2;
	fixFailedFits(hNsUL_obsEQ0, radius);
	fixFailedFits(hObs_obsEQ0,  radius);


	Float_t zsize      = hMed->GetMaximum()*1.25;
	Float_t zsize_full = hMed_full->GetMaximum()*1.25;
	
	float maxreldiff = 0.02;
	TString smaxreldiff = tstr(maxreldiff*100,0);
	TH2F* hGrid      = (TH2F*)getGrid(hMed,maxreldiff,zsize,kBlack,1);
	TH2F* hGrid_full = (TH2F*)getGrid(hMed_full,maxreldiff,zsize_full);
	
	TCanvas* cnv;
	TVirtualPad* p1;
	TVirtualPad* p2;
	
	
	TH2F* hGrid1 = (TH2F*)hNSR1->Clone("tmp1");
	hGrid1->Reset();
	hGrid1->SetLineColor(17);
	hGrid1->SetLineStyle(3);
	for(Int_t bx=1 ; bx<=hNSR1->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=hNSR1->GetNbinsY() ; ++by)
		{
			float z = hNSR1->GetBinContent(bx,by);
			if(z<=0.5) hGrid1->SetBinContent(bx,by,10);
		}
	}
	TH2F* hGrid1_full = (TH2F*)hNSR1_full->Clone("tmp1");
	hGrid1_full->Reset();
	hGrid1_full->SetLineColor(17);
	hGrid1_full->SetLineWidth(1);
	for(Int_t bx=1 ; bx<=hNSR1_full->GetNbinsX() ; ++bx)
	{
		for(Int_t by=1 ; by<=hNSR1_full->GetNbinsY() ; ++by)
		{
			float z = hNSR1_full->GetBinContent(bx,by);
			if(z<=0.5) hGrid1_full->SetBinContent(bx,by,10);
		}
	}
	
	
	
	TBox* boxMedian    = boxHighlight(kAzure,   1, hMed);
	// TBox* boxMin1sig   = boxHighlight(kBlack,   1, hMin1);
	// TBox* boxMin2sig   = boxHighlight(kGray+1,  1, hMin2);
	// TBox* box1sigbnd   = boxHighlight(kGreen+2, 1, h1Ssigma);
	// TBox* box2sigbnd   = boxHighlight(kPink+1,  1, h2Ssigma);
	// TBox* boxRatLoose  = boxHighlight(kCyan+4,  2, hMed, 0.946,8);
	// TBox* boxRatMedium = boxHighlight(kCyan+3,  2, hMed, 0.961,8);
	// TBox* boxRatTight  = boxHighlight(kCyan+2,  2, hMed, 0.975,8);
	TBox* boxRatLoose  = boxHighlight(kGray+1,  1, hMed, 0.975,8);
	TBox* boxRatMedium = boxHighlight(kGreen+2, 1, hMed, 0.980,8);
	TBox* boxRatTight  = boxHighlight(kBlack,   1, hMed, 0.984,8);
	
	TBox* boxMedian_full    = boxHighlight(kAzure,   1, hMed_full);
	TBox* boxRatLoose_full  = boxHighlight(kCyan+4,  2, hMed_full, 0.975,26);
	TBox* boxRatMedium_full = boxHighlight(kCyan+3,  2, hMed_full, 0.980,26);
	TBox* boxRatTight_full  = boxHighlight(kCyan+2,  2, hMed_full, 0.984,26);
	
	legN->Clear();
	legN->AddEntry(hGrid,"<"+smaxreldiff+"\%#timesMedian(min)", "l");
	legN->AddEntry(boxMedian, "Median(min)", "l");
	legN->AddEntry(boxRatLoose,"Loose(rational)", "l");
	legN->AddEntry(boxRatMedium,"Medium(rational)", "l");
	legN->AddEntry(boxRatTight,"Tight(rational)", "l");
	
	legF->Clear();
	legF->AddEntry(hGrid_full,"<"+smaxreldiff+"\%#timesMedian(min)", "l");
	legF->AddEntry(boxMedian_full, "Median(min)", "l");
	
	legP->Clear();
	legP->AddEntry(hGrid1,"N_{SR1}<0.5", "l");
	legP->AddEntry(hGrid,"<"+smaxreldiff+"\%#timesMedian(min)", "l");
	legP->AddEntry(boxMedian, "Median(min)", "l");
	legP->AddEntry(boxRatLoose,"Loose(rational)", "l");
	legP->AddEntry(boxRatMedium,"Medium(rational)", "l");
	legP->AddEntry(boxRatTight,"Tight(rational)", "l");
	

	

	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hEff_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hEff_full);
	cnv->Modified();
	ptxt->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Efficiency.Full.eps");
	cnv->SaveAs(outname+".Efficiency.Full.pdf");
	cnv->SaveAs(outname+".pdf(");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hEff->Draw("colz text45");
	cnv->Update();
	fixPalette(hEff);
	cnv->Modified();
	ptxt->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Efficiency.eps");
	cnv->SaveAs(outname+".Efficiency.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNSB1_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hNSB1_full);
	cnv->Modified();
	ptxt->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NSB1.Full.eps");
	cnv->SaveAs(outname+".NSB1.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNSB1->Draw("colz text45");
	cnv->Update();
	fixPalette(hNSB1);
	cnv->Modified();
	ptxt->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NSB1.eps");
	cnv->SaveAs(outname+".NSB1.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNSR1_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hNSR1_full);
	hGrid1_full->Draw("box same");
	cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NSR1.Full.eps");
	cnv->SaveAs(outname+".NSR1.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNSR1->Draw("colz text45");
	cnv->Update();
	fixPalette(hNSR1);
	hGrid1->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NSR1.eps");
	cnv->SaveAs(outname+".NSR1.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;


	
	cnv = new TCanvas("c","",1400,600); cnv->Divide(2,1); cnv->Draw(); cnv->cd();
	p1 = (TVirtualPad*)cnv->cd(1); p1->cd(); p1->SetRightMargin(0.15); p1->SetTicks(1,1);
	hPoisson->Draw("colz text45");
	p1->Update();
	fixPalette(hPoisson);
	hGrid->Draw("box same");
	hGrid1->Draw("box same");
	p1->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	p1->RedrawAxis();
	p1->Update();
	delete palette;
	/////////////////
	p2 = (TVirtualPad*)cnv->cd(2); p2->cd(); p2->SetRightMargin(0.15); p2->SetTicks(1,1);
	hEffPoisson->Draw("colz text45");
	p2->Update();
	fixPalette(hEffPoisson);
	hGrid->Draw("box same");
	hGrid1->Draw("box same");
	p2->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	p2->RedrawAxis();
	p2->Update();
	cnv->SaveAs(outname+".PoissonTest.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;


	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hPoisson->Draw("colz text45");
	cnv->Update();
	fixPalette(hPoisson);
	hGrid->Draw("box same");
	hGrid1->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".PoissonProb.eps");
	cnv->SaveAs(outname+".PoissonProb.pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hEffPoisson->Draw("colz text45");
	cnv->Update();
	fixPalette(hEffPoisson);
	hGrid->Draw("box same");
	hGrid1->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legP->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".EffOverPoisson.eps");
	cnv->SaveAs(outname+".EffOverPoisson.pdf");
	delete palette;
	delete cnv;
	
	
	
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hDat_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hDat_full);
	cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Nobs.fromNSR1.Full.eps");
	cnv->SaveAs(outname+".Nobs.fromNSR1.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hDat->Draw("colz text45");
	cnv->Update();
	fixPalette(hDat);
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Nobs.fromNSR1.eps");
	cnv->SaveAs(outname+".Nobs.fromNSR1.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	
	cnv = new TCanvas("c","",1400,600); cnv->Divide(2,1); cnv->Draw(); cnv->cd();
	p1 = (TVirtualPad*)cnv->cd(1); p1->cd(); p1->SetRightMargin(0.15); p1->SetTicks(1,1);
	// hNsUL_obsEQ0->SetMaximum(3);
	hNsUL_obsEQ0->Draw("colz text45");
	p1->Update();
	fixPalette(hNsUL_obsEQ0);
	// hGrid->Draw("box same");
	p1->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p1->RedrawAxis();
	p1->Update();
	delete palette;
	/////////////////
	p2 = (TVirtualPad*)cnv->cd(2); p2->cd(); p2->SetRightMargin(0.15); p2->SetTicks(1,1);
	// hObs_obsEQ0->SetMaximum(3);
	hObs_obsEQ0->Draw("colz text45");
	p2->Update();
	fixPalette(hObs_obsEQ0);
	// hGrid->Draw("box same");
	p2->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p2->RedrawAxis();
	p2->Update();
	cnv->SaveAs(outname+".UL.obs0.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_obsEQ0->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_obsEQ0);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.obs0.eps");
	cnv->SaveAs(outname+".NsUL.obs0.pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hObs_obsEQ0->Draw("colz text45");
	cnv->Update();
	fixPalette(hObs_obsEQ0);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".BRUL.obs0.eps");
	cnv->SaveAs(outname+".BRUL.obs0.pdf");
	delete palette;
	delete cnv;
	
	
	
	cnv = new TCanvas("c","",1400,600); cnv->Divide(2,1); cnv->Draw(); cnv->cd();
	p1 = (TVirtualPad*)cnv->cd(1); p1->cd(); p1->SetRightMargin(0.15); p1->SetTicks(1,1);
	// hNsUL_obsEQ1->SetMaximum(3);
	hNsUL_obsEQ1->Draw("colz text45");
	p1->Update();
	fixPalette(hNsUL_obsEQ1);
	// hGrid->Draw("box same");
	p1->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p1->RedrawAxis();
	p1->Update();
	delete palette;
	/////////////////
	p2 = (TVirtualPad*)cnv->cd(2); p2->cd(); p2->SetRightMargin(0.15); p2->SetTicks(1,1);
	// hObs_obsEQ1->SetMaximum(3);
	hObs_obsEQ1->Draw("colz text45");
	p2->Update();
	fixPalette(hObs_obsEQ1);
	// hGrid->Draw("box same");
	p2->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p2->RedrawAxis();
	p2->Update();
	cnv->SaveAs(outname+".UL.obs1.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_obsEQ1->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_obsEQ1);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.obs1.eps");
	cnv->SaveAs(outname+".NsUL.obs1.pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hObs_obsEQ1->Draw("colz text45");
	cnv->Update();
	fixPalette(hObs_obsEQ1);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".BRUL.obs1.eps");
	cnv->SaveAs(outname+".BRUL.obs1.pdf");
	delete palette;
	delete cnv;
	
	
	cnv = new TCanvas("c","",1400,600); cnv->Divide(2,1); cnv->Draw(); cnv->cd();
	p1 = (TVirtualPad*)cnv->cd(1); p1->cd(); p1->SetRightMargin(0.15); p1->SetTicks(1,1);
	// hNsUL_obsEQ2->SetMaximum(3);
	hNsUL_obsEQ2->Draw("colz text45");
	p1->Update();
	fixPalette(hNsUL_obsEQ2);
	// hGrid->Draw("box same");
	p1->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p1->RedrawAxis();
	p1->Update();
	delete palette;
	/////////////////
	p2 = (TVirtualPad*)cnv->cd(2); p2->cd(); p2->SetRightMargin(0.15); p2->SetTicks(1,1);
	// hObs_obsEQ2->SetMaximum(3);
	hObs_obsEQ2->Draw("colz text45");
	p2->Update();
	fixPalette(hObs_obsEQ2);
	// hGrid->Draw("box same");
	p2->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p2->RedrawAxis();
	p2->Update();
	cnv->SaveAs(outname+".UL.obs2.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_obsEQ2->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_obsEQ2);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.obs2.eps");
	cnv->SaveAs(outname+".NsUL.obs2.pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hObs_obsEQ2->Draw("colz text45");
	cnv->Update();
	fixPalette(hObs_obsEQ2);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".BRUL.obs2.eps");
	cnv->SaveAs(outname+".BRUL.obs2.pdf");
	delete palette;
	delete cnv;
	
	

	cnv = new TCanvas("c","",1400,600); cnv->Divide(2,1); cnv->Draw(); cnv->cd();
	p1 = (TVirtualPad*)cnv->cd(1); p1->cd(); p1->SetRightMargin(0.15); p1->SetTicks(1,1);
	// hNsUL_obsEQ3->SetMaximum(3);
	hNsUL_obsEQ3->Draw("colz text45");
	p1->Update();
	fixPalette(hNsUL_obsEQ3);
	// hGrid->Draw("box same");
	p1->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p1->RedrawAxis();
	p1->Update();
	delete palette;
	/////////////////
	p2 = (TVirtualPad*)cnv->cd(2); p2->cd(); p2->SetRightMargin(0.15); p2->SetTicks(1,1);
	// hObs_obsEQ3->SetMaximum(3);
	hObs_obsEQ3->Draw("colz text45");
	p2->Update();
	fixPalette(hObs_obsEQ3);
	// hGrid->Draw("box same");
	p2->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	p2->RedrawAxis();
	p2->Update();
	cnv->SaveAs(outname+".UL.obs3.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_obsEQ3->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_obsEQ3);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.obs3.eps");
	cnv->SaveAs(outname+".NsUL.obs3.pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hObs_obsEQ3->Draw("colz text45");
	cnv->Update();
	fixPalette(hObs_obsEQ3);
	cnv->Modified();
	boxMedian->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".BRUL.obs3.eps");
	cnv->SaveAs(outname+".BRUL.obs3.pdf");
	delete palette;
	delete cnv;


	
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_exp_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_exp_full);
	hGrid_full->Draw("box same");
	cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.Observed.fromExp.Full.eps");
	cnv->SaveAs(outname+".NsUL.Observed.fromExp.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_exp->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_exp);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.Observed.fromExp.eps");
	cnv->SaveAs(outname+".NsUL.Observed.fromExp.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_med_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_med_full);
	hGrid_full->Draw("box same");
	cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.Median.Full.eps");
	cnv->SaveAs(outname+".NsUL.Median.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hNsUL_med->Draw("colz text45");
	cnv->Update();
	fixPalette(hNsUL_med);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".NsUL.Median.eps");
	cnv->SaveAs(outname+".NsUL.Median.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hMed_full->Draw("colz text45");
	cnv->Update();
	fixPalette(hMed_full);
	hGrid_full->Draw("box same");
	cnv->Modified();
	boxMedian_full->Draw("same");
	ptxt->Draw("same");
	legF->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Median.Full.eps");
	cnv->SaveAs(outname+".Median.Full.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;

	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hMed->Draw("colz text45");
	cnv->Update();
	fixPalette(hMed);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// box2sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".Median.eps");
	cnv->SaveAs(outname+".Median.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hMin1->Draw("colz text45");
	cnv->Update();
	fixPalette(hMin1);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".MedianMin1sig.eps");
	cnv->SaveAs(outname+".MedianMin1sig.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hMin2->Draw("colz text45");
	cnv->Update();
	fixPalette(hMin2);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".MedianMin2sig.eps");
	cnv->SaveAs(outname+".MedianMin2sig.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;

	// cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	// hPls1->Draw("colz text45");
	// cnv->Update();
	// fixPalette(hPls1);
	// hGrid->Draw("box same");
	// cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	// ptxt->Draw("same");
	// legN->Draw("smae");
	// cnv->RedrawAxis();
	// cnv->Update();
	// cnv->SaveAs(outname+".MedianPls1sig.eps");
	// cnv->SaveAs(outname+".MedianPls1sig.pdf");
	// cnv->SaveAs(outname+".pdf");
	// delete palette;
	// delete cnv;
	// 
	// cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	// hPls2->Draw("colz text45");
	// cnv->Update();
	// fixPalette(hPls2);
	// hGrid->Draw("box same");
	// cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	// ptxt->Draw("same");
	// legN->Draw("smae");
	// cnv->RedrawAxis();
	// cnv->Update();
	// cnv->SaveAs(outname+".MedianPls2sig.eps");
	// cnv->SaveAs(outname+".MedianPls2sig.pdf");
	// cnv->SaveAs(outname+".pdf)");
	// delete palette;
	// delete cnv;
	

	
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	h1Ssigma->Draw("colz text45");
	cnv->Update();
	fixPalette(h1Ssigma);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".1SigBand.eps");
	cnv->SaveAs(outname+".1SigBand.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	// cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	// h2Ssigma->Draw("colz text45");
	// cnv->Update();
	// fixPalette(h2Ssigma);
	// hGrid->Draw("box same");
	// cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// box2sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	// ptxt->Draw("same");
	// legN->Draw("smae");
	// cnv->RedrawAxis();
	// cnv->Update();
	// cnv->SaveAs(outname+".2SigBand.eps");
	// cnv->SaveAs(outname+".2SigBand.pdf");
	// cnv->SaveAs(outname+".pdf");
	// delete palette;
	// delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hPls1SigmaBand->Draw("colz text45");
	cnv->Update();
	fixPalette(hPls1SigmaBand);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".1SigBandUp.eps");
	cnv->SaveAs(outname+".1SigBandUp.pdf");
	cnv->SaveAs(outname+".pdf");
	delete palette;
	delete cnv;
	
	cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	hMin1SigmaBand->Draw("colz text45");
	cnv->Update();
	fixPalette(hMin1SigmaBand);
	hGrid->Draw("box same");
	cnv->Modified();
	boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	boxRatLoose->Draw("same");
	boxRatMedium->Draw("same");
	boxRatTight->Draw("same");
	ptxt->Draw("same");
	legN->Draw("smae");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(outname+".1SigBandDn.eps");
	cnv->SaveAs(outname+".1SigBandDn.pdf");
	cnv->SaveAs(outname+".pdf)");
	delete palette;
	delete cnv;
	
	// cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	// hPls2SigmaBand->Draw("colz text45");
	// cnv->Update();
	// fixPalette(hPls2SigmaBand);
	// hGrid->Draw("box same");
	// cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// box2sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	// ptxt->Draw("same");
	// legN->Draw("smae");
	// cnv->RedrawAxis();
	// cnv->Update();
	// cnv->SaveAs(outname+".2SigBandUp.eps");
	// cnv->SaveAs(outname+".2SigBandUp.pdf");
	// cnv->SaveAs(outname+".pdf");
	// delete palette;
	// delete cnv;
	// 
	// cnv = new TCanvas("c","",700,600); cnv->Draw(); cnv->cd(); cnv->SetRightMargin(0.15); cnv->SetTicks(1,1);
	// hMin2SigmaBand->Draw("colz text45");
	// cnv->Update();
	// fixPalette(hMin2SigmaBand);
	// hGrid->Draw("box same");
	// cnv->Modified();
	// boxMedian->Draw("same");
	// boxMin1sig->Draw("same");
	// boxMin2sig->Draw("same");
	// box1sigbnd->Draw("same");
	// box2sigbnd->Draw("same");
	// boxRatLoose->Draw("same");
	// boxRatMedium->Draw("same");
	// boxRatTight->Draw("same");
	// ptxt->Draw("same");
	// legN->Draw("smae");
	// cnv->RedrawAxis();
	// cnv->Update();
	// cnv->SaveAs(outname+".2SigBandDn.eps");
	// cnv->SaveAs(outname+".2SigBandDn.pdf");
	// cnv->SaveAs(outname+".pdf)");
	// delete palette;
	// delete cnv;
}
