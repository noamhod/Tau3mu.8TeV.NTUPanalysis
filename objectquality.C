/////////////////////////////////////////
//// root -b -l -q objectquality.C++ ////
/////////////////////////////////////////
#include "std.h"
#include "type.h"

TString oname = "";

TString tstr(float x, int prcn=-1)
{
        stringstream strm;
        string str;
        if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
        else         strm << fixed << x;
        strm >> str;
        return (TString)str;
}
TString strip(TString sIn)
{
	TString sOut = sIn;
	sOut.ReplaceAll("|logx","");
	sOut.ReplaceAll("|logy","");
	sOut.ReplaceAll("|logz","");
	return sOut;
}
bool isLog(TString sIn, TString axis)
{
	if(axis=="x" && sIn.Contains("|logx")) return true;
	if(axis=="y" && sIn.Contains("|logy")) return true;
	if(axis=="z" && sIn.Contains("|logz")) return true;
	return false;
}
void setMax(TMapTSP2TH1& histos1, TString name)
{
	TString name1=name;
	TString name2=name;
	float max = histos1[name1]->GetMaximum();
	TString toreplace   = (name.Contains("onJpsi")) ? "onJpsi"  : "offJpsi";
	TString replacewith = (toreplace=="onJpsi")     ? "offJpsi" : "onJpsi";
	name2.ReplaceAll(toreplace,replacewith);
	max = (histos1[name2]->GetMaximum()>max) ? histos1[name2]->GetMaximum() : max;
	histos1[name1]->SetMaximum(1.2*max);
	histos1[name2]->SetMaximum(1.2*max);
}

void draw(TMapTSP2TH1& histos1, vector<TString>& prefix, TString varname)
{
	TString barename = strip(varname);

	TString name1 = prefix[0]+barename+"_offJpsi";
	TString name2 = prefix[0]+barename+"_onJpsi";
	// TString name3 = prefix[1]+barename+"_offJpsi";
	TString name4 = prefix[1]+barename+"_onJpsi";
	// TString name5 = prefix[2]+barename+"_offJpsi";
	// TString name6 = prefix[2]+barename+"_onJpsi";
	

	float max = (histos1[name1]->GetMaximum()>histos1[name2]->GetMaximum()) ? histos1[name1]->GetMaximum() : histos1[name2]->GetMaximum();
	histos1[name1]->SetMaximum(1.7*max);
	histos1[name2]->SetMaximum(1.7*max);
	
	TLegend* leg = new TLegend(0.5,0.6,0.9,0.85,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry(histos1[name4], "#it{bb#rightarrow J}/#psi#it{#rightarrow#mu4#mu4}","f");
	leg->AddEntry(histos1[name1], "Data (3body) OFF #it{J}/#psi","f");
	leg->AddEntry(histos1[name2], "Data (3body) ON #it{J}/#psi","p");
	
	TCanvas* cnv = new TCanvas("","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	if(isLog(varname,"y")) cnv->SetLogy();
	else histos1[name1]->SetMinimum(0.);
	cnv->cd();
	histos1[name1]->Draw("hist");
	histos1[name4]->Draw("hist same");
	histos1[name2]->Draw("psame");
	// histos1[name3]->Draw("hist same");
	// histos1[name5]->Draw("hist same");
	// histos1[name6]->Draw("hist same");
	leg->Draw("same");
	
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs(oname+"objqa.pdf");
	cnv->SaveAs(oname+"objqa."+barename+".eps");
	
	delete leg;
	delete cnv;
}

void draw(TMapTSP2TH2& histos2, vector<TString>& prefix, TString varname)
{
	TString barename = strip(varname);

	TString name1 = prefix[0]+barename+"_offJpsi";
	TString name2 = prefix[0]+barename+"_onJpsi";
	TString name3 = prefix[1]+barename+"_offJpsi";
	TString name4 = prefix[1]+barename+"_onJpsi";
	// TString name5 = prefix[2]+barename+"_offJpsi";
	// TString name6 = prefix[2]+barename+"_onJpsi";
	
	// float max = (histos2[name1]->GetMaximum()>histos2[name2]->GetMaximum()) ? histos2[name1]->GetMaximum() : histos2[name2]->GetMaximum();
	// histos2[name1]->SetMaximum(1.05*max);
	// histos2[name2]->SetMaximum(1.05*max);
	
	TCanvas* cnv1 = new TCanvas("1","",600,400);
	cnv1->Draw();
	cnv1->SetTicks(1,1);
	if(isLog(varname,"z")) cnv1->SetLogz();
	cnv1->cd();
	histos2[name1]->Draw("col");
	cnv1->RedrawAxis();
	cnv1->Update();
	cnv1->SaveAs(oname+"objqa.pdf");
	cnv1->SaveAs(oname+"objqa."+barename+"_offJpsi.eps");
	
	TCanvas* cnv2 = new TCanvas("2","",600,400);
	cnv2->Draw();
	cnv2->SetTicks(1,1);
	if(isLog(varname,"z")) cnv2->SetLogz();
	cnv2->cd();
	histos2[name2]->Draw("col");
	cnv2->RedrawAxis();
	cnv2->Update();
	cnv2->SaveAs(oname+"objqa.pdf");
	cnv2->SaveAs(oname+"objqa."+barename+"_onJpsi.eps");
	
	delete cnv1;
	delete cnv2;
}

void objectquality(TString master = "muons", TString method = "cuts")
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
	
	
	TFile* fIn = new TFile("histos."+master+"."+method+".root", "READ");
	TMapTSP2TH1 histos1;
	TMapTSP2TH2 histos2;
	TMapiTS     order1;
	TMapiTS     order2;
	
	oname = "figures/"+method+"."+master+".";
	vector<TString> pref;
	pref.push_back("Data_");
	// pref.push_back("bb_mu4mu4_");
	pref.push_back("bb_Jpsimu4mu4_");
	// pref.push_back("JZxW_");
	
	vector<TString> suff;
	suff.push_back("_onJpsi");
	suff.push_back("_offJpsi");
	
	vector<TString> var1names;
	vector<TString> var2names;
	
	var1names.push_back("m2body");
	
	var1names.push_back("trks_hitstudy_fitprob|logy");
	var1names.push_back("trks_hitstudy_pixdEdx");
	
	var1names.push_back("trk_hitstudy_eta1");
	var1names.push_back("trk_hitstudy_eta2");
	
	var1names.push_back("mu_hitstudy_eta");
	var1names.push_back("mu_sctangsig|logy");
	var1names.push_back("mu_sctngbsig|logy");
	var1names.push_back("mu_pbalsig|logy");
	var1names.push_back("mu_qoverp|logy");
	var1names.push_back("mu_chi2ndf|logy");
	var1names.push_back("mu_pvalue|logy");

	var1names.push_back("TPa_hitstudy_eta");
	var1names.push_back("TPa_qoverp|logy");
	var1names.push_back("TPa_chi2ndf|logy");
	var1names.push_back("TPa_pvalue|logy");
	
	var1names.push_back("mu_nHighThresholdTRTHits|logy");
	var1names.push_back("mu_nHighThresholdTRTHitsFraction|logy");
	var1names.push_back("mu_nHighThresholdTRTHitsRatio|logy");
	var1names.push_back("TPa_nHighThresholdTRTHits|logy");
	var1names.push_back("TPa_nHighThresholdTRTHitsFraction|logy");
	var1names.push_back("TPa_nHighThresholdTRTHitsRatio|logy");

	var1names.push_back("mu_outliersOnTrack1");
	var1names.push_back("mu_outliersOnTrack2");
	var1names.push_back("TPa_outliersOnTrack1");
	var1names.push_back("TPa_outliersOnTrack2");
	
	var1names.push_back("mu_MDTHits1");
	var1names.push_back("mu_MDTHits2");
	var1names.push_back("TPa_MDTHits1");
	var1names.push_back("TPa_MDTHits2");
	
	var1names.push_back("mu_precisionHits1");
	var1names.push_back("mu_precisionHits2");
	var1names.push_back("TPa_precisionHits1");
	var1names.push_back("TPa_precisionHits2");
	
	var1names.push_back("mu_precisionHoles1");
	var1names.push_back("mu_precisionHoles2");
	var1names.push_back("TPa_precisionHoles1");
	var1names.push_back("TPa_precisionHoles2");
	
	var1names.push_back("mu_precisionOutliers1");
	var1names.push_back("mu_precisionOutliers2");
	var1names.push_back("TPa_precisionOutliers1");
	var1names.push_back("TPa_precisionOutliers2");
	
	var1names.push_back("mu_phiLayers1");
	var1names.push_back("mu_phiLayers2");
	var1names.push_back("TPa_phiLayers1");
	var1names.push_back("TPa_phiLayers2");
	
	var1names.push_back("mu_etaphiLayers1");
	var1names.push_back("mu_etaphiLayers2");
	var1names.push_back("TPa_etaphiLayers1");
	var1names.push_back("TPa_etaphiLayers2");
	
	var2names.push_back("mu_RPC_vs_TGC_phiHits1|logz");
	var2names.push_back("mu_RPC_vs_TGC_phiHits2|logz");
	var2names.push_back("TPa_RPC_vs_TGC_phiHits1|logz");
	var2names.push_back("TPa_RPC_vs_TGC_phiHits2|logz");
	
	var2names.push_back("mu_RPC_vs_TGC_etaHits1|logz");
	var2names.push_back("mu_RPC_vs_TGC_etaHits2|logz");
	var2names.push_back("TPa_RPC_vs_TGC_etaHits1|logz");
	var2names.push_back("TPa_RPC_vs_TGC_etaHits2|logz");
	
	// Get the histos
	
	for(unsigned int p=0 ; p<pref.size() ; ++p)
	{
		for(unsigned int s=0 ; s<suff.size() ; ++s)
		{
			for(unsigned int v=0 ; v<var1names.size() ; ++v)
			{
				TString barename = strip(var1names[v]);
				TString name = pref[p]+barename+suff[s];
				cout << "Getting: " << name << endl;
				histos1.insert(make_pair(name, (TH1F*)fIn->Get(name)));
				if(s==0)
				{
					histos1[name]->SetMarkerStyle(20);
					histos1[name]->SetMarkerSize(0.8);
					histos1[name]->SetMarkerColor(kBlack);
				}
				if(s==1)
				{
					// histos1[name]->SetLineColor(kRed-3);
					// histos1[name]->SetFillColor(kRed-3);
					histos1[name]->SetFillStyle(3005);
				}
				TString xtitle = histos1[name]->GetXaxis()->GetTitle();
				TString ytitle = histos1[name]->GetYaxis()->GetTitle();
				TString title  = histos1[name]->GetTitle();
				xtitle.ReplaceAll(" (on the J/#psi)", ""); xtitle.ReplaceAll(" (off the J/#psi)", ""); histos1[name]->GetXaxis()->SetTitle(xtitle);
				ytitle.ReplaceAll(" (on the J/#psi)", ""); ytitle.ReplaceAll(" (off the J/#psi)", ""); histos1[name]->GetYaxis()->SetTitle(ytitle);
				title.ReplaceAll(" (on the J/#psi)", "");  title.ReplaceAll(" (off the J/#psi)", "");  histos1[name]->SetTitle(title);
			}
			for(unsigned int v=0 ; v<var2names.size() ; ++v)
			{
				TString barename = strip(var2names[v]);
				TString name = pref[p]+barename+suff[s];
				cout << "Getting: " << name << endl;
				histos2.insert(make_pair(name, (TH2F*)fIn->Get(name)));
			}
		}
	}
	
	TCanvas* cnv = new TCanvas("c","c",600,400);
	cnv->SaveAs(oname+"objqa.pdf(");
	for(unsigned int v=0 ; v<var1names.size() ; ++v) draw(histos1,pref,var1names[v]);
	for(unsigned int v=0 ; v<var2names.size() ; ++v) draw(histos2,pref,var2names[v]);
	cnv->SaveAs(oname+"objqa.pdf)");
}