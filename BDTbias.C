////////////////////////////////////////////////
//// root -b -l -q BDTbias.C++\(\"loose\"\) ////
////////////////////////////////////////////////

#include "postBDTcuts.h"
#include "type.h"
#include "const.h"
// #include "jets.h"

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
	// ptxt = new TPaveText(0.62,0.70,0.87,0.87,"NDC");
	ptxt = new TPaveText(0.15,0.70,0.40,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

TLegend* leg;
vector<TLegend*> legS;
void makeLegend()
{
	// leg = new TLegend(0.15,0.55,0.45,0.88,NULL,"brNDC");
	leg = new TLegend(0.65,0.55,0.95,0.87,NULL,"brNDC");
	leg->SetFillStyle(4000); //will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	
	for(unsigned int i=0 ; i<9 ; ++i)
	{
		legS.push_back(new TLegend(0.45,0.80,0.91,0.85,NULL,"brNDC"));
		legS[i]->SetFillStyle(4000); //will be transparent
		legS[i]->SetFillColor(0);
		legS[i]->SetTextFont(42);
		legS[i]->SetBorderSize(0);
	}
}


void BDTbias(TString selection)
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
	
	makeAtlasLabel();
	makeLegend();
	
	//// calculate the signal efficiency
	// TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	// TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	
	TString ntuple = "flatout"; // or "mvaout"
	TFile* f = (ntuple=="flatout") ? new TFile("flatout.periodall+MC.muons.mva.analysis.n0.j0.loose.root","READ") : new TFile("mvaout.muons.root","READ");;
	TTree* tS = (TTree*)f->Get((ntuple=="flatout") ? "flatout_Wtaunu_3mu" : "fltmva_Wtaunu_3mu");
	TTree* tD = (TTree*)f->Get((ntuple=="flatout") ? "flatout_Data"       : "fltmva_Data");
	TString score  = (ntuple=="mvaout") ? "score"  : "mva_score";
	TString m3body = (ntuple=="mvaout") ? "m3body" : "vtx_mass";

	TString smMinTraining  = tstr(mTrainingMinGlob,0);
	TString smMaxTraining  = tstr(mTrainingMaxGlob,0);
	TString smMinSBleft    = tstr(mSideBandLeftLowerMeVGlob,0);
	TString smMaxSBleft    = tstr(mBlindMinGlob,0);
	TString smMinSBright   = tstr(mBlindMaxGlob,0);
	TString smMaxSBright   = tstr(mSideBandRightUpperMeVGlob,0);
	TString sxSRmin        = tstr(mSRminMeVGlob,0);
	TString sxSRmax        = tstr(mSRmaxMeVGlob,0);
	TString soptBDTcut     = tstr(optBDTcut,4);
	TString sminBDTcut     = tstr(minBDTcut,4);
	
	cout << endl;
	

	TString cuts_sig = "";
	TString cuts_bkg = "";
	Float_t npassedS = -1;
	Float_t npassedD = -1;
	bool blinded = true; bool unblinded = false;
	bool doBDT   = true; bool noBDT     = false;
	
	const int     nmbins = 11;
	const Double_t mbins[nmbins+1] = {mSideBandLeftLowerMeVGlob,1510,1570,1630,1690,1750,1810,1870,1930,1990,2050,mSideBandRightUpperMeVGlob};

	const int     nmbins1 = 6;
	const Double_t mbins1[nmbins1+1] = {mSideBandLeftLowerMeVGlob,1570,1690,1810,1930,2050,mSideBandRightUpperMeVGlob};
	
	TString mode = (selection=="loose") ? "preTraining" : "postTraining";
	
	//// tight+x>x1 selection
	cuts_sig = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"sig", sxSRmin,sxSRmax, unblinded,doBDT,mode,ntuple);
	cuts_bkg = postBDTcut(smMinSBleft,smMaxSBleft,smMinSBright,smMaxSBright,sminBDTcut,sminBDTcut,"bkg", sxSRmin,sxSRmax, unblinded,doBDT,mode,ntuple);
	cuts_bkg += " && ("+score+"<0.8)";
	
	cout << "cuts_bkg=" << cuts_bkg << endl;
	
	TH2* BDTvsMsig = new TH2F("BDTvsMsig",";#it{m}_{3body} [MeV];BDT score;Events",nmbins,mbins,10,-0.9,+0.8);
	TH2* BDTvsMbkg = new TH2F("BDTvsMbkg",";#it{m}_{3body} [MeV];BDT score;Events",nmbins,mbins,10,-0.9,+0.8);

	TProfile* BDTvsMsig_profile = new TProfile("BDTvsMsig_profile",";#it{m}_{3body} [MeV];BDT score;Events",nmbins,mbins,-0.9,+0.8);
	BDTvsMsig_profile->SetLineColor(kBlack);
	BDTvsMsig_profile->SetMarkerColor(kBlack);
	BDTvsMsig_profile->SetMarkerStyle(20);
	
	TH1* m1 = new TH1F("m1",";#it{m}_{3body} [MeV]",nmbins,mbins); m1->SetLineColor(19); m1->SetFillColor(19); m1->Sumw2();
	TH1* m2 = new TH1F("m2",";#it{m}_{3body} [MeV]",nmbins,mbins); m2->SetLineColor(18); m2->SetFillColor(18); m2->Sumw2();
	TH1* m3 = new TH1F("m3",";#it{m}_{3body} [MeV]",nmbins,mbins); m3->SetLineColor(17); m3->SetFillColor(17); m3->Sumw2();
	TH1* m4 = new TH1F("m4",";#it{m}_{3body} [MeV]",nmbins,mbins); m4->SetLineColor(16); m4->SetFillColor(16); m4->Sumw2();
	TH1* m5 = new TH1F("m5",";#it{m}_{3body} [MeV]",nmbins,mbins); m5->SetLineColor(15); m5->SetFillColor(15); m5->Sumw2();
	TH1* m6 = new TH1F("m6",";#it{m}_{3body} [MeV]",nmbins,mbins); m6->SetLineColor(14); m6->SetFillColor(14); m6->Sumw2();
	TH1* m7 = new TH1F("m7",";#it{m}_{3body} [MeV]",nmbins,mbins); m7->SetLineColor(13); m7->SetFillColor(13); m7->Sumw2();
	TH1* m8 = new TH1F("m8",";#it{m}_{3body} [MeV]",nmbins,mbins); m8->SetLineColor(12); m8->SetFillColor(12); m8->Sumw2();
	TH1* m9 = new TH1F("m9",";#it{m}_{3body} [MeV]",nmbins,mbins); m9->SetLineColor(1);  m9->SetFillColor(1);  m9->Sumw2();
	
	TH1* m1R = new TH1F("m1R",";#it{m}_{3body} [MeV]",nmbins,mbins); m1R->SetLineColor(19); m1R->SetFillColor(19); m1R->Sumw2();
	TH1* m2R = new TH1F("m2R",";#it{m}_{3body} [MeV]",nmbins,mbins); m2R->SetLineColor(18); m2R->SetFillColor(18); m2R->Sumw2();
	TH1* m3R = new TH1F("m3R",";#it{m}_{3body} [MeV]",nmbins,mbins); m3R->SetLineColor(17); m3R->SetFillColor(17); m3R->Sumw2();
	TH1* m4R = new TH1F("m4R",";#it{m}_{3body} [MeV]",nmbins,mbins); m4R->SetLineColor(16); m4R->SetFillColor(16); m4R->Sumw2();
	TH1* m5R = new TH1F("m5R",";#it{m}_{3body} [MeV]",nmbins,mbins); m5R->SetLineColor(15); m5R->SetFillColor(15); m5R->Sumw2();
	TH1* m6R = new TH1F("m6R",";#it{m}_{3body} [MeV]",nmbins,mbins); m6R->SetLineColor(14); m6R->SetFillColor(14); m6R->Sumw2();
	TH1* m7R = new TH1F("m7R",";#it{m}_{3body} [MeV]",nmbins,mbins); m7R->SetLineColor(13); m7R->SetFillColor(13); m7R->Sumw2();
	TH1* m8R = new TH1F("m8R",";#it{m}_{3body} [MeV]",nmbins,mbins); m8R->SetLineColor(12); m8R->SetFillColor(12); m8R->Sumw2();
	TH1* m9R = new TH1F("m9R",";#it{m}_{3body} [MeV]",nmbins,mbins); m9R->SetLineColor(1);  m9R->SetFillColor(1);  m9R->Sumw2();
	
	
	tD->Draw(m3body+">>m1",cuts_bkg+" && ("+score+">-0.9)"); cout << "done m1" << endl;
	tD->Draw(m3body+">>m2",cuts_bkg+" && ("+score+">-0.7)"); cout << "done m2" << endl;
	tD->Draw(m3body+">>m3",cuts_bkg+" && ("+score+">-0.5)"); cout << "done m3" << endl;
	tD->Draw(m3body+">>m4",cuts_bkg+" && ("+score+">-0.3)"); cout << "done m4" << endl;
	tD->Draw(m3body+">>m5",cuts_bkg+" && ("+score+">-0.1)"); cout << "done m5" << endl;
	tD->Draw(m3body+">>m6",cuts_bkg+" && ("+score+">+0.1)"); cout << "done m6" << endl;
	tD->Draw(m3body+">>m7",cuts_bkg+" && ("+score+">+0.3)"); cout << "done m7" << endl;
	tD->Draw(m3body+">>m8",cuts_bkg+" && ("+score+">+0.5)"); cout << "done m8" << endl;
	tD->Draw(m3body+">>m9",cuts_bkg+" && ("+score+">+0.7)"); cout << "done m9" << endl;
	
	tD->Draw(m3body+">>m1R",cuts_bkg+" && ("+score+">-0.9 && "+score+"<=-0.7)");
	tD->Draw(m3body+">>m2R",cuts_bkg+" && ("+score+">-0.7 && "+score+"<=-0.5)");
	tD->Draw(m3body+">>m3R",cuts_bkg+" && ("+score+">-0.5 && "+score+"<=-0.3)");
	tD->Draw(m3body+">>m4R",cuts_bkg+" && ("+score+">-0.3 && "+score+"<=-0.1)");
	tD->Draw(m3body+">>m5R",cuts_bkg+" && ("+score+">-0.1 && "+score+"<=+0.1)");
	tD->Draw(m3body+">>m6R",cuts_bkg+" && ("+score+">+0.1 && "+score+"<=+0.3)");
	tD->Draw(m3body+">>m7R",cuts_bkg+" && ("+score+">+0.3 && "+score+"<=+0.5)");
	tD->Draw(m3body+">>m8R",cuts_bkg+" && ("+score+">+0.5 && "+score+"<=+0.7)");
	tD->Draw(m3body+">>m9R",cuts_bkg+" && ("+score+">+0.7 && "+score+"<+0.8)");
	
	
	TString title = (selection=="loose") ? "Loose" : "Tight";

	// npassedS = tS->Draw(m3body+":"+score+">>BDTvsMSIG",cuts_sig); cout << title+"+x<0.8 only (signal): " << npassedS << endl;
	npassedD = tD->Draw(score+":"+m3body+">>BDTvsMbkg",cuts_bkg); cout << title+"+x<0.8 only (SB data): " << npassedD << endl;
	tD->Draw(score+":"+m3body+">>BDTvsMsig_profile",cuts_bkg);
	
	TCanvas* cnv = new TCanvas("cnv","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	BDTvsMbkg->Draw("col");
	BDTvsMsig_profile->Draw("same");
	ptxt->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf(");
	cnv->SaveAs("figures/BDTbias2D."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias2D."+selection+".png");
	cnv->SaveAs("figures/BDTbias2D."+selection+".eps");
	
	
	if(cnv) delete cnv;
	cnv = new TCanvas("cnv","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	m1->SetMinimum(0);
	m1->SetMaximum(m1->GetMaximum()*1.3);
	m1->Draw("hist");      leg->AddEntry(m1,title+" && -0.9<#it{x}<0.8","f");
	m2->Draw("hist same"); leg->AddEntry(m2,title+" && -0.7<#it{x}<0.8","f");
	m3->Draw("hist same"); leg->AddEntry(m3,title+" && -0.5<#it{x}<0.8","f");
	m4->Draw("hist same"); leg->AddEntry(m4,title+" && -0.3<#it{x}<0.8","f");
	m5->Draw("hist same"); leg->AddEntry(m5,title+" && -0.1<#it{x}<0.8","f");
	m6->Draw("hist same"); leg->AddEntry(m6,title+" && +0.1<#it{x}<0.8","f");
	m7->Draw("hist same"); leg->AddEntry(m7,title+" && +0.3<#it{x}<0.8","f");
	m8->Draw("hist same"); leg->AddEntry(m8,title+" && +0.5<#it{x}<0.8","f");
	m9->Draw("hist same"); leg->AddEntry(m9,title+" && +0.7<#it{x}<0.8","f");
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dinclusive."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dinclusive."+selection+".png");
	cnv->SaveAs("figures/BDTbias1Dinclusive."+selection+".eps");
	
	if(cnv) delete cnv; leg->Clear();
	cnv = new TCanvas("cnv","",600,400);
	cnv->Draw();
	cnv->SetTicks(1,1);
	m1R->SetMinimum(0);
	m1R->SetMaximum(m1R->GetMaximum()*1.3);
	m1R->Draw("hist");      leg->AddEntry(m1R,title+" && -0.9<#it{x}<-0.7","f");
	m2R->Draw("hist same"); leg->AddEntry(m2R,title+" && -0.7<#it{x}<-0.5","f");
	m3R->Draw("hist same"); leg->AddEntry(m3R,title+" && -0.5<#it{x}<-0.3","f");
	m4R->Draw("hist same"); leg->AddEntry(m4R,title+" && -0.3<#it{x}<-0.1","f");
	m5R->Draw("hist same"); leg->AddEntry(m5R,title+" && -0.1<#it{x}<+0.1","f");
	m6R->Draw("hist same"); leg->AddEntry(m6R,title+" && +0.1<#it{x}<+0.3","f");
	m7R->Draw("hist same"); leg->AddEntry(m7R,title+" && +0.3<#it{x}<+0.5","f");
	m8R->Draw("hist same"); leg->AddEntry(m8R,title+" && +0.5<#it{x}<+0.7","f");
	m9R->Draw("hist same"); leg->AddEntry(m9R,title+" && +0.7<#it{x}<+0.8","f");
	ptxt->Draw("same");
	leg->Draw("same");
	cnv->Update();
	cnv->RedrawAxis();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dbinned."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dbinned."+selection+".png");
	cnv->SaveAs("figures/BDTbias1Dbinned."+selection+".eps");
	
	
	
	
	
	m1R->SetLineColor(kBlack); m1R->SetMarkerColor(kBlack); m1R->SetMarkerStyle(20);
	m2R->SetLineColor(kBlack); m2R->SetMarkerColor(kBlack); m2R->SetMarkerStyle(20);
	m3R->SetLineColor(kBlack); m3R->SetMarkerColor(kBlack); m3R->SetMarkerStyle(20);
	m4R->SetLineColor(kBlack); m4R->SetMarkerColor(kBlack); m4R->SetMarkerStyle(20);
	m5R->SetLineColor(kBlack); m5R->SetMarkerColor(kBlack); m5R->SetMarkerStyle(20);
	m6R->SetLineColor(kBlack); m6R->SetMarkerColor(kBlack); m6R->SetMarkerStyle(20);
	m7R->SetLineColor(kBlack); m7R->SetMarkerColor(kBlack); m7R->SetMarkerStyle(20);
	m8R->SetLineColor(kBlack); m8R->SetMarkerColor(kBlack); m8R->SetMarkerStyle(20);
	m9R->SetLineColor(kBlack); m9R->SetMarkerColor(kBlack); m9R->SetMarkerStyle(20);
	
	
	
	
	if(cnv) delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(3,3);
	cnv->Draw();
	cnv->cd(1); m1->SetLineColor(kAzure-3); m1->SetFillColor(kAzure-3); m1->SetMinimum(0); m1->SetMaximum(m1->GetMaximum()); m1->Draw("hist"); ptxt->Draw("same"); legS[0]->AddEntry(m1,title+" && -0.9<#it{x}<0.8","f"); legS[0]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(2); m2->SetLineColor(kAzure-3); m2->SetFillColor(kAzure-3); m2->SetMinimum(0); m2->SetMaximum(m1->GetMaximum()); m2->Draw("hist"); ptxt->Draw("same"); legS[1]->AddEntry(m2,title+" && -0.7<#it{x}<0.8","f"); legS[1]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(3); m3->SetLineColor(kAzure-3); m3->SetFillColor(kAzure-3); m3->SetMinimum(0); m3->SetMaximum(m1->GetMaximum()); m3->Draw("hist"); ptxt->Draw("same"); legS[2]->AddEntry(m3,title+" && -0.5<#it{x}<0.8","f"); legS[2]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(4); m4->SetLineColor(kAzure-3); m4->SetFillColor(kAzure-3); m4->SetMinimum(0); m4->SetMaximum(m1->GetMaximum()); m4->Draw("hist"); ptxt->Draw("same"); legS[3]->AddEntry(m4,title+" && -0.3<#it{x}<0.8","f"); legS[3]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(5); m5->SetLineColor(kAzure-3); m5->SetFillColor(kAzure-3); m5->SetMinimum(0); m5->SetMaximum(m1->GetMaximum()); m5->Draw("hist"); ptxt->Draw("same"); legS[4]->AddEntry(m5,title+" && -0.1<#it{x}<0.8","f"); legS[4]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(6); m6->SetLineColor(kAzure-3); m6->SetFillColor(kAzure-3); m6->SetMinimum(0); m6->SetMaximum(m1->GetMaximum()); m6->Draw("hist"); ptxt->Draw("same"); legS[5]->AddEntry(m6,title+" && +0.1<#it{x}<0.8","f"); legS[5]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(7); m7->SetLineColor(kAzure-3); m7->SetFillColor(kAzure-3); m7->SetMinimum(0); m7->SetMaximum(m1->GetMaximum()); m7->Draw("hist"); ptxt->Draw("same"); legS[6]->AddEntry(m7,title+" && +0.3<#it{x}<0.8","f"); legS[6]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(8); m8->SetLineColor(kAzure-3); m8->SetFillColor(kAzure-3); m8->SetMinimum(0); m8->SetMaximum(m1->GetMaximum()); m8->Draw("hist"); ptxt->Draw("same"); legS[7]->AddEntry(m8,title+" && +0.5<#it{x}<0.8","f"); legS[7]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(9); m9->SetLineColor(kAzure-3); m9->SetFillColor(kAzure-3); m9->SetMinimum(0); m9->SetMaximum(m1->GetMaximum()); m9->Draw("hist"); ptxt->Draw("same"); legS[8]->AddEntry(m9,title+" && +0.7<#it{x}<0.8","f"); legS[8]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->Update();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dmultiple."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1Dmultiple."+selection+".png");
	cnv->SaveAs("figures/BDTbias1Dmultiple."+selection+".eps");
	
	
	
	
	TH1* m1Rcopy = (TH1*)m1R->Clone((TString)m1R->GetName()+"copy");
	TH1* m2Rcopy = (TH1*)m2R->Clone((TString)m2R->GetName()+"copy");
	TH1* m3Rcopy = (TH1*)m3R->Clone((TString)m3R->GetName()+"copy");
	TH1* m4Rcopy = (TH1*)m4R->Clone((TString)m4R->GetName()+"copy");
	TH1* m5Rcopy = (TH1*)m5R->Clone((TString)m5R->GetName()+"copy");
	TH1* m6Rcopy = (TH1*)m6R->Clone((TString)m6R->GetName()+"copy");
	TH1* m7Rcopy = (TH1*)m7R->Clone((TString)m7R->GetName()+"copy");
	TH1* m8Rcopy = (TH1*)m8R->Clone((TString)m8R->GetName()+"copy");
	TH1* m9Rcopy = (TH1*)m9R->Clone((TString)m9R->GetName()+"copy");
	
	if(cnv) delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(3,3);
	cnv->Draw();
	cnv->cd(1); m1Rcopy->SetLineColor(kAzure-3); m1Rcopy->SetFillColor(kAzure-3); m1Rcopy->SetMinimum(0); m1Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m1Rcopy->Draw("hist"); ptxt->Draw("same"); legS[0]->Clear(); legS[0]->AddEntry(m1Rcopy,title+" && -0.9<#it{x}<-0.7","f"); legS[0]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(2); m2Rcopy->SetLineColor(kAzure-3); m2Rcopy->SetFillColor(kAzure-3); m2Rcopy->SetMinimum(0); m2Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m2Rcopy->Draw("hist"); ptxt->Draw("same"); legS[1]->Clear(); legS[1]->AddEntry(m2Rcopy,title+" && -0.7<#it{x}<-0.5","f"); legS[1]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(3); m3Rcopy->SetLineColor(kAzure-3); m3Rcopy->SetFillColor(kAzure-3); m3Rcopy->SetMinimum(0); m3Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m3Rcopy->Draw("hist"); ptxt->Draw("same"); legS[2]->Clear(); legS[2]->AddEntry(m3Rcopy,title+" && -0.5<#it{x}<-0.3","f"); legS[2]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(4); m4Rcopy->SetLineColor(kAzure-3); m4Rcopy->SetFillColor(kAzure-3); m4Rcopy->SetMinimum(0); m4Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m4Rcopy->Draw("hist"); ptxt->Draw("same"); legS[3]->Clear(); legS[3]->AddEntry(m4Rcopy,title+" && -0.3<#it{x}<-0.1","f"); legS[3]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(5); m5Rcopy->SetLineColor(kAzure-3); m5Rcopy->SetFillColor(kAzure-3); m5Rcopy->SetMinimum(0); m5Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m5Rcopy->Draw("hist"); ptxt->Draw("same"); legS[4]->Clear(); legS[4]->AddEntry(m5Rcopy,title+" && -0.1<#it{x}<+0.1","f"); legS[4]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(6); m6Rcopy->SetLineColor(kAzure-3); m6Rcopy->SetFillColor(kAzure-3); m6Rcopy->SetMinimum(0); m6Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m6Rcopy->Draw("hist"); ptxt->Draw("same"); legS[5]->Clear(); legS[5]->AddEntry(m6Rcopy,title+" && +0.1<#it{x}<+0.3","f"); legS[5]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(7); m7Rcopy->SetLineColor(kAzure-3); m7Rcopy->SetFillColor(kAzure-3); m7Rcopy->SetMinimum(0); m7Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m7Rcopy->Draw("hist"); ptxt->Draw("same"); legS[6]->Clear(); legS[6]->AddEntry(m7Rcopy,title+" && +0.3<#it{x}<+0.5","f"); legS[6]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(8); m8Rcopy->SetLineColor(kAzure-3); m8Rcopy->SetFillColor(kAzure-3); m8Rcopy->SetMinimum(0); m8Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m8Rcopy->Draw("hist"); ptxt->Draw("same"); legS[7]->Clear(); legS[7]->AddEntry(m8Rcopy,title+" && +0.5<#it{x}<+0.7","f"); legS[7]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(9); m9Rcopy->SetLineColor(kAzure-3); m9Rcopy->SetFillColor(kAzure-3); m9Rcopy->SetMinimum(0); m9Rcopy->SetMaximum(m1Rcopy->GetMaximum()); m9Rcopy->Draw("hist"); ptxt->Draw("same"); legS[8]->Clear(); legS[8]->AddEntry(m9Rcopy,title+" && +0.7<#it{x}<+0.8","f"); legS[8]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->Update();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1DmultipleBDTbins."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1DmultipleBDTbins."+selection+".png");
	cnv->SaveAs("figures/BDTbias1DmultipleBDTbins."+selection+".eps");
	
	
	
	
	
	TH1* m1Rcopy0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"copy0");
	TH1* m2Rcopy0 = (TH1*)m2R->Clone((TString)m2R->GetName()+"copy0");
	TH1* m3Rcopy0 = (TH1*)m3R->Clone((TString)m3R->GetName()+"copy0");
	TH1* m4Rcopy0 = (TH1*)m4R->Clone((TString)m4R->GetName()+"copy0");
	TH1* m5Rcopy0 = (TH1*)m5R->Clone((TString)m5R->GetName()+"copy0");
	TH1* m6Rcopy0 = (TH1*)m6R->Clone((TString)m6R->GetName()+"copy0");
	TH1* m7Rcopy0 = (TH1*)m7R->Clone((TString)m7R->GetName()+"copy0");
	TH1* m8Rcopy0 = (TH1*)m8R->Clone((TString)m8R->GetName()+"copy0");
	TH1* m9Rcopy0 = (TH1*)m9R->Clone((TString)m9R->GetName()+"copy0");
	
	TH1* m1Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m1Rdenominator0->Reset();
	TH1* m2Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m2Rdenominator0->Reset();
	TH1* m3Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m3Rdenominator0->Reset();
	TH1* m4Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m4Rdenominator0->Reset();
	TH1* m5Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m5Rdenominator0->Reset();
	TH1* m6Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m6Rdenominator0->Reset();
	TH1* m7Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m7Rdenominator0->Reset();
	TH1* m8Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m8Rdenominator0->Reset();
	TH1* m9Rdenominator0 = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator0"); m9Rdenominator0->Reset();
	
	m1Rdenominator0->Add(m1R);
	m2Rdenominator0->Add(m1Rdenominator0);
	m3Rdenominator0->Add(m2Rdenominator0); m3Rdenominator0->Add(m2R);
	m4Rdenominator0->Add(m3Rdenominator0); m4Rdenominator0->Add(m3R);
	m5Rdenominator0->Add(m4Rdenominator0); m5Rdenominator0->Add(m4R);
	m6Rdenominator0->Add(m5Rdenominator0); m6Rdenominator0->Add(m5R);
	m7Rdenominator0->Add(m6Rdenominator0); m7Rdenominator0->Add(m6R);
	m8Rdenominator0->Add(m7Rdenominator0); m8Rdenominator0->Add(m7R);
	m9Rdenominator0->Add(m8Rdenominator0); m9Rdenominator0->Add(m8R);
	
	if(cnv) delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(3,3);
	cnv->Draw();
	cnv->cd(1); m1Rcopy0->Divide(m1Rdenominator0); m1Rcopy0->SetMinimum(0); m1Rcopy0->SetMaximum(1.1); m1Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[0]->Clear(); legS[0]->AddEntry(m1Rcopy0,"#frac{-0.9<#it{x}<-0.7}{-0.9<#it{x}<-0.7}","ple"); legS[0]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(2); m2Rcopy0->Divide(m2Rdenominator0); m2Rcopy0->SetMinimum(0); m2Rcopy0->SetMaximum(1.1); m2Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[1]->Clear(); legS[1]->AddEntry(m2Rcopy0,"#frac{-0.7<#it{x}<-0.5}{-0.9<#it{x}<-0.7}","ple"); legS[1]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(3); m3Rcopy0->Divide(m3Rdenominator0); m3Rcopy0->SetMinimum(0); m3Rcopy0->SetMaximum(1.1); m3Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[2]->Clear(); legS[2]->AddEntry(m3Rcopy0,"#frac{-0.5<#it{x}<-0.3}{-0.9<#it{x}<-0.5}","ple"); legS[2]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(4); m4Rcopy0->Divide(m4Rdenominator0); m4Rcopy0->SetMinimum(0); m4Rcopy0->SetMaximum(1.1); m4Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[3]->Clear(); legS[3]->AddEntry(m4Rcopy0,"#frac{-0.3<#it{x}<-0.1}{-0.9<#it{x}<-0.3}","ple"); legS[3]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(5); m5Rcopy0->Divide(m5Rdenominator0); m5Rcopy0->SetMinimum(0); m5Rcopy0->SetMaximum(1.1); m5Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[4]->Clear(); legS[4]->AddEntry(m5Rcopy0,"#frac{-0.1<#it{x}<+0.1}{-0.9<#it{x}<-0.1}","ple"); legS[4]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(6); m6Rcopy0->Divide(m6Rdenominator0); m6Rcopy0->SetMinimum(0); m6Rcopy0->SetMaximum(1.1); m6Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[5]->Clear(); legS[5]->AddEntry(m6Rcopy0,"#frac{+0.1<#it{x}<+0.3}{-0.9<#it{x}<+0.1}","ple"); legS[5]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(7); m7Rcopy0->Divide(m7Rdenominator0); m7Rcopy0->SetMinimum(0); m7Rcopy0->SetMaximum(1.1); m7Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[6]->Clear(); legS[6]->AddEntry(m7Rcopy0,"#frac{+0.3<#it{x}<+0.5}{+0.9<#it{x}<+0.3}","ple"); legS[6]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(8); m8Rcopy0->Divide(m8Rdenominator0); m8Rcopy0->SetMinimum(0); m8Rcopy0->SetMaximum(1.1); m8Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[7]->Clear(); legS[7]->AddEntry(m8Rcopy0,"#frac{+0.5<#it{x}<+0.7}{+0.9<#it{x}<+0.5}","ple"); legS[7]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(9); m9Rcopy0->Divide(m9Rdenominator0); m9Rcopy0->SetMinimum(0); m9Rcopy0->SetMaximum(1.1); m9Rcopy0->Draw("p0 e0"); ptxt->Draw("same"); legS[8]->Clear(); legS[8]->AddEntry(m9Rcopy0,"#frac{+0.7<#it{x}<+0.8}{+0.9<#it{x}<+0.7}","ple"); legS[8]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->Update();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosSteps."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosSteps."+selection+".png");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosSteps."+selection+".eps");
	
	
	
	
	TH1* m1Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m1Rnumerator->Reset();
	TH1* m2Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m2Rnumerator->Reset();
	TH1* m3Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m3Rnumerator->Reset();
	TH1* m4Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m4Rnumerator->Reset();
	TH1* m5Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m5Rnumerator->Reset();
	TH1* m6Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m6Rnumerator->Reset();
	TH1* m7Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m7Rnumerator->Reset();
	TH1* m8Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m8Rnumerator->Reset();
	// TH1* m9Rnumerator = (TH1*)m1R->Clone((TString)m1R->GetName()+"numerator"); m9Rnumerator->Reset();
	
	TH1* m1Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m1Rdenominator->Reset();
	TH1* m2Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m2Rdenominator->Reset();
	TH1* m3Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m3Rdenominator->Reset();
	TH1* m4Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m4Rdenominator->Reset();
	TH1* m5Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m5Rdenominator->Reset();
	TH1* m6Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m6Rdenominator->Reset();
	TH1* m7Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m7Rdenominator->Reset();
	TH1* m8Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m8Rdenominator->Reset();
	// TH1* m9Rdenominator = (TH1*)m1R->Clone((TString)m1R->GetName()+"denominator"); m9Rdenominator->Reset();
	
	m1Rnumerator->Add(m2R); m1Rnumerator->Add(m3R); m1Rnumerator->Add(m4R); m1Rnumerator->Add(m5R); m1Rnumerator->Add(m6R); m1Rnumerator->Add(m7R); m1Rnumerator->Add(m8R); m1Rnumerator->Add(m9R);
	m2Rnumerator->Add(m3R); m2Rnumerator->Add(m4R); m2Rnumerator->Add(m5R); m2Rnumerator->Add(m6R); m2Rnumerator->Add(m7R); m2Rnumerator->Add(m8R); m2Rnumerator->Add(m9R);
	m3Rnumerator->Add(m4R); m3Rnumerator->Add(m5R); m3Rnumerator->Add(m6R); m3Rnumerator->Add(m7R); m3Rnumerator->Add(m8R); m3Rnumerator->Add(m9R);
	m4Rnumerator->Add(m5R); m4Rnumerator->Add(m6R); m4Rnumerator->Add(m7R); m4Rnumerator->Add(m8R); m4Rnumerator->Add(m9R);
	m5Rnumerator->Add(m6R); m5Rnumerator->Add(m7R); m5Rnumerator->Add(m8R); m5Rnumerator->Add(m9R);
	m6Rnumerator->Add(m7R); m6Rnumerator->Add(m8R); m6Rnumerator->Add(m9R);
	m7Rnumerator->Add(m8R); m7Rnumerator->Add(m9R);
	m8Rnumerator->Add(m9R);
	
	m1Rdenominator->Add(m1R); // numerator = 2+3+4+5+6+7+8+9
	m2Rdenominator->Add(m1R); m2Rdenominator->Add(m2R); // numerator = 3+4+5+6+7+8+9
	m3Rdenominator->Add(m1R); m3Rdenominator->Add(m2R); m3Rdenominator->Add(m3R); // numerator = 4+5+6+7+8+9
	m4Rdenominator->Add(m1R); m4Rdenominator->Add(m2R); m4Rdenominator->Add(m3R); m4Rdenominator->Add(m4R); // numerator = 5+6+7+8+9
	m5Rdenominator->Add(m1R); m5Rdenominator->Add(m2R); m5Rdenominator->Add(m3R); m5Rdenominator->Add(m4R); m5Rdenominator->Add(m5R); // numerator = 6+7+8+9
	m6Rdenominator->Add(m1R); m6Rdenominator->Add(m2R); m6Rdenominator->Add(m3R); m6Rdenominator->Add(m4R); m6Rdenominator->Add(m5R); m6Rdenominator->Add(m6R); // numerator = 7+8+9
	m7Rdenominator->Add(m1R); m7Rdenominator->Add(m2R); m7Rdenominator->Add(m3R); m7Rdenominator->Add(m4R); m7Rdenominator->Add(m5R); m7Rdenominator->Add(m6R); m7Rdenominator->Add(m7R); // numerator = 8+9
	m8Rdenominator->Add(m1R); m8Rdenominator->Add(m2R); m8Rdenominator->Add(m3R); m8Rdenominator->Add(m4R); m8Rdenominator->Add(m5R); m8Rdenominator->Add(m6R); m8Rdenominator->Add(m7R); m8Rdenominator->Add(m8R);	// numerator = 9
	
	if(cnv) delete cnv;
	cnv = new TCanvas("cnv","",1200,1200);
	cnv->Divide(3,3);
	cnv->Draw();
	cnv->cd(1); m1Rnumerator->Divide(m1Rdenominator); m1Rnumerator->SetMinimum(0); m1Rnumerator->SetMaximum(m1Rnumerator->GetMaximum()*2); m1Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[0]->Clear(); legS[0]->AddEntry(m1Rnumerator,"#frac{-0.7<#it{x}<+0.8}{-0.9<#it{x}<-0.7}","ple"); legS[0]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(2); m2Rnumerator->Divide(m2Rdenominator); m2Rnumerator->SetMinimum(0); m2Rnumerator->SetMaximum(m2Rnumerator->GetMaximum()*2); m2Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[1]->Clear(); legS[1]->AddEntry(m2Rnumerator,"#frac{-0.5<#it{x}<+0.8}{-0.9<#it{x}<-0.5}","ple"); legS[1]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(3); m3Rnumerator->Divide(m2Rdenominator); m3Rnumerator->SetMinimum(0); m3Rnumerator->SetMaximum(m3Rnumerator->GetMaximum()*2); m3Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[2]->Clear(); legS[2]->AddEntry(m3Rnumerator,"#frac{-0.3<#it{x}<+0.8}{-0.9<#it{x}<-0.3}","ple"); legS[2]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(4); m4Rnumerator->Divide(m3Rdenominator); m4Rnumerator->SetMinimum(0); m4Rnumerator->SetMaximum(m4Rnumerator->GetMaximum()*2); m4Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[3]->Clear(); legS[3]->AddEntry(m4Rnumerator,"#frac{-0.1<#it{x}<+0.8}{-0.9<#it{x}<-0.1}","ple"); legS[3]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(5); m5Rnumerator->Divide(m4Rdenominator); m5Rnumerator->SetMinimum(0); m5Rnumerator->SetMaximum(m5Rnumerator->GetMaximum()*2); m5Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[4]->Clear(); legS[4]->AddEntry(m5Rnumerator,"#frac{+0.1<#it{x}<+0.8}{-0.9<#it{x}<+0.1}","ple"); legS[4]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(6); m6Rnumerator->Divide(m5Rdenominator); m6Rnumerator->SetMinimum(0); m6Rnumerator->SetMaximum(m6Rnumerator->GetMaximum()*2); m6Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[5]->Clear(); legS[5]->AddEntry(m6Rnumerator,"#frac{+0.3<#it{x}<+0.8}{-0.9<#it{x}<+0.3}","ple"); legS[5]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(7); m7Rnumerator->Divide(m6Rdenominator); m7Rnumerator->SetMinimum(0); m7Rnumerator->SetMaximum(m7Rnumerator->GetMaximum()*2); m7Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[6]->Clear(); legS[6]->AddEntry(m7Rnumerator,"#frac{+0.5<#it{x}<+0.8}{-0.9<#it{x}<+0.5}","ple"); legS[6]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(8); m8Rnumerator->Divide(m7Rdenominator); m8Rnumerator->SetMinimum(0); m8Rnumerator->SetMaximum(m8Rnumerator->GetMaximum()*2); m8Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[7]->Clear(); legS[7]->AddEntry(m8Rnumerator,"#frac{+0.7<#it{x}<+0.8}{-0.9<#it{x}<+0.7}","ple"); legS[7]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	// cnv->cd(9); m9Rnumerator->Divide(m8Rdenominator); m9Rnumerator->SetMinimum(0); m9Rnumerator->SetMaximum(1.1); m9Rnumerator->Draw("p0 e0"); ptxt->Draw("same"); legS[8]->Clear(); legS[8]->AddEntry(m9Rnumerator,"#frac{+0.7<#it{x}<+0.8}{-0.9<#it{x}<+0.7}","ple"); legS[8]->Draw("same"); gPad->SetTicks(1,1); gPad->RedrawAxis(); gPad->Update();
	cnv->Update();
	cnv->SaveAs("figures/BDTbiasAll."+selection+".pdf)");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosFull."+selection+".pdf");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosFull."+selection+".png");
	cnv->SaveAs("figures/BDTbias1DmultipleRatiosFull."+selection+".eps");
}