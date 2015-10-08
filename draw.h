#ifndef DRAW_H
#define DRAW_H

#include "hist.h"

TLegend* leg;
TPaveText* ptxt;
TCanvas* c1;
TCanvas* cnv;
vector<TVirtualPad*> pads;
int divx=-1;
int divy=-1;
TString pdffilename;
static int visibility = 1;

void makeAtlasLabel()
{	
	ptxt = new TPaveText(0.12,0.65,0.37,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt="+slumi);
	ptxt->AddText("#sqrt{s}=8 TeV");
}

void makeLegend()
{
	leg = new TLegend(0.55,0.6,0.83,0.83,NULL,"brNDC");
	setLegendDefaults(leg);
}

void makeSingleCnv(bool is2d=false)
{
	c1 = (is2d) ? new TCanvas("","",400,400) : new TCanvas("","",600,400);
}
void closeSingleCnv(TString varname, bool is2d=false, TString channel="")
{
	c1->SetTicks(1,1);
	c1->RedrawAxis();
	c1->Update();
	TString epsfilename = pdffilename;
	TString newsuffix = (channel=="") ? "."+varname+".eps"  :  "."+varname+"."+channel+".eps";
	epsfilename.ReplaceAll(".pdf",newsuffix);
	c1->SaveAs(epsfilename);
	if(!is2d)
	{
		// TH1
		epsfilename.ReplaceAll(".eps",".logy.eps");
		c1->SetLogy();
		c1->SaveAs(epsfilename);
	}
	else
	{
		// TH2
		epsfilename.ReplaceAll(".eps",".logz.eps");
		c1->SetLogz();
		c1->SaveAs(epsfilename);
	}
	delete c1;
}


void drawData(TString varname, TMapTSP2TH1& histos, bool isfirst, TString drawopt="pex0", Bool_t doymin=0, Bool_t doymax=0, Double_t ymin=-1, Double_t ymax=-1)
{	
	TMapTSP2TH1::iterator it=histos.begin();
	TString hname = "";
	for(; it!=histos.end() ; it++)
	{
		hname = it->first;
		if(!isData(hname))           continue;
		if(!hname.EndsWith(varname)) continue;
		break; // if got to here, break
	}
	if(it==histos.end()) return;
	
	it->second->SetLineColor(kBlack);
	it->second->SetLineWidth(2);
	it->second->SetLineStyle(1);
	it->second->SetMarkerStyle(20);
	it->second->SetMarkerSize(1);
	it->second->SetMarkerColor(kBlack);
	
	double min = getYmin(it->second);
	double max = getYmax(it->second);
	
	if(min==max && min==0.)
	{
		min = 1.e-2;
		max = 1.;
	}
	
	if(doymin) it->second->SetMinimum(min*0.9);
	if(doymax) it->second->SetMaximum(max*2.5);
	if(ymin!=-1) it->second->SetMinimum(ymin);
	if(ymax!=-1) it->second->SetMaximum(ymax);
	
	if(isfirst) it->second->Draw(drawopt);
	else        it->second->Draw(drawopt+" same");
}
void drawData(TString varname, TMapTSP2TH2& histos, bool isfirst, TString drawopt)
{	
	TMapTSP2TH2::iterator it=histos.begin();
	for(; it!=histos.end() ; it++)
	{
		TString hname = it->first;
		if(!isData(hname))           continue;
		if(!hname.EndsWith(varname)) continue;
		break; // if got to here, break
	}
	if(it==histos.end()) return;
	
	it->second->SetLineColor(kBlack);
	it->second->SetFillColor(kBlack);
	it->second->SetMarkerColor(kBlack);
	it->second->SetLineStyle(1);
	
	if(drawopt=="")
	{
		it->second->SetContour(50);
		if(isfirst) it->second->Draw("col");
		else        it->second->Draw("col same");
		exeGray->Draw();
		it->second->Draw("col same");
	}
	else
	{
		if(isfirst) it->second->Draw(drawopt);
		else        it->second->Draw(drawopt+" same");
	}
}

void drawPadSingle(TString channel, TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH1& histos, TString drawopt="hist")
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		if(drawchannels[pit->second]==0) continue;
		if(isData(pit->second))                                    continue;
		if(pit->second!=channel)                                   continue; // draw only one channel
		
		TString name = pit->second;
		TString hname = name+"_"+varname;
		
		double min = getYmin(histos[hname]);
		double max = getYmax(histos[hname]);
		
		if(min==max && min==0.)
		{
			min = 1e-2;
			max = 1.;
		}
		
		histos[hname]->SetMinimum(min*0.9);
		histos[hname]->SetMaximum(max*2.5);
		
		c1->cd();
		histos[hname]->Draw(drawopt);
		pad->cd();
		histos[hname]->Draw(drawopt);
	}
	// the data
	c1->cd();
	if(isData(channel)) drawData(varname,histos,false,drawopt,true,true);
	pad->cd();
	if(isData(channel)) drawData(varname,histos,false,drawopt,true,true);
	
	c1->cd();
	ptxt->Draw("same");
	pad->cd();
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,false,channel);
}
void drawPadSingle(TString channel, TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH2& histos, TString drawopt="box")
{
	makeSingleCnv(true);
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	float max = -1.e20;
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		if(drawchannels[pit->second]==0) continue;
		
		TString name = pit->second;
		TString hname = name+"_"+varname;
		max = (histos[hname]->GetMaximum()>max) ? histos[hname]->GetMaximum() : max;
	}
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		if(drawchannels[pit->second]==0) continue;
		if(isData(pit->second))                                    continue;
		if(pit->second!=channel)                                   continue; // draw only one channel
		
		TString name = pit->second;
		TString hname = name+"_"+varname;
		
		c1->cd();
		histos[hname]->Draw(drawopt);
		pad->cd();
		histos[hname]->Draw(drawopt);
	}
	// the data
	if(isData(channel))
	{
		c1->cd();
		drawData(varname,histos,true,drawopt);
		pad->cd();
		drawData(varname,histos,true,drawopt);
	}
	
	c1->cd();
	ptxt->Draw("same");
	pad->cd();
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,true,channel);
}

void drawPadAll(TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH1& histos/*, TLegend* leg*/, Double_t ymin=-1, Double_t ymax=-1)
{
	makeSingleCnv();
	pad->cd();
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	unsigned int counter = 0;
	double max = -1.e20;
	double min = +1.e20;
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		TString name = pit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(isWsignal(name) || isData(name))
		{
			Double_t hmin = getYmin(histos[hname]);
			Double_t hmax = getYmax(histos[hname]);
			max = (hmax>max) ? hmax : max;
			min = (hmin<min) ? hmin : min;
		}
	}
	if(ymin==-1) ymin = min*0.9;
	if(ymax==-1) ymax = max*2.5;
	if(min==max && min==0.)
	{
		min = 1e-2;
		max = 1.;
	}
	
	// legend...
	leg->Clear();
	
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		if(drawchannels[pit->second]==0) continue;
		if(isData(pit->second))          continue;
		
		TString name = pit->second;
		TString hname = name+"_"+varname;
		
		if(isData(name))   continue;
				
		histos[hname]->SetMinimum(ymin);
		histos[hname]->SetMaximum(ymax);
		
		TString ytitle = histos[hname]->GetYaxis()->GetTitle();
		TString label = (isSignal(name) && ytitle.Contains("Events")) ? labels[name]+"#times"+sBRf : labels[name];
		leg->AddEntry(histos1[hname],label,legoptions[name]);
		
		c1->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("hist same");
		
		pad->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("hist same");
		
		counter++;
	}
	leg->AddEntry(histos1["Data_"+varname],labels["Data"],legoptions["Data"]);
	
	
	// the data
	c1->cd();
	drawData(varname,histos,false,"pex0",false,false,ymin,ymax);
	leg->Draw("same");
	ptxt->Draw("same");
	
	pad->cd();
	drawData(varname,histos,false,"pex0",false,false,ymin,ymax);
	leg->Draw("same");
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}
void drawPadAll(TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH2& histos, TString drawopt)
{
	makeSingleCnv(true);
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	
	// the MSs
	TMapuiTS::reverse_iterator rpit=procs.rbegin();
	TString name = rpit->second;
	TString hname = name+"_"+varname;
	TH2* hbkg     = NULL;
	TH2* hsignal  = NULL;
	TH2* hWsignal = NULL;
	TH2* hdata    = NULL;
	for( ; rpit!=procs.rend() ; rpit++)
	{		
		name = rpit->second;
		if(drawchannels[name]==0) continue;

		hname = name+"_"+varname;
		
		if(!isData(name)) continue;
		
		if(isData(name))    hdata    = histos[hname];
		if(isWsignal(name)) hWsignal = histos[hname];
	}
	
	double max = -1.e20;
	max = (hdata->GetMaximum()>max)    ? hdata->GetMaximum()    : max;
	max = (hWsignal->GetMaximum()>max) ? hWsignal->GetMaximum() : max;
	// max = (hsignal->GetMaximum()>max) ? hsignal->GetMaximum() : max;
	// max = (hbkg->GetMaximum()>max)    ? hbkg->GetMaximum()    : max;
	
	hdata->SetMaximum(max*2.5);
	hsignal->SetMaximum(max*2.5);
	hbkg->SetMaximum(max*2.5);
	
	// the data
	drawData(varname,histos,true,drawopt);
	
	if(drawopt=="")
	{
		c1->cd();
		hbkg->SetContour(50);
		hbkg->Draw("col same");
		exeRed->Draw();
		hbkg->Draw("col same");
		pad->cd();
		hbkg->SetContour(50);
		hbkg->Draw("col same");
		exeRed->Draw();
		hbkg->Draw("col same");
		
		c1->cd();
		hsignal->SetContour(50);
		hsignal->Draw("col same");
		exeGreen->Draw();
		hsignal->Draw("col same");
		pad->cd();
		hsignal->SetContour(50);
		hsignal->Draw("col same");
		exeGreen->Draw();
		hsignal->Draw("col same");
	}
	else
	{
		c1->cd();
		hbkg->Draw(drawopt+" same");
		hsignal->Draw(drawopt+" same");
		pad->cd();
		hbkg->Draw(drawopt+" same");
		hsignal->Draw(drawopt+" same");
	}
	
	c1->cd();
	ptxt->Draw("same");
	pad->cd();
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,true);
}
void drawPadStack(TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH1& histos, /*TLegend* leg,*/ Double_t ymin=-1, Double_t ymax=-1)
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	
	TMapuiTS::reverse_iterator rit=procs.rbegin();
	TH1* hsumBkg = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hWsig   = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hsumSig = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hData   = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TString t  = histos[rit->second+"_"+varname]->GetTitle();
	TString tx = histos[rit->second+"_"+varname]->GetXaxis()->GetTitle();
	TString ty = histos[rit->second+"_"+varname]->GetYaxis()->GetTitle();
	TString titles = t+";"+tx+";"+ty;
	THStack* hsBkg = new THStack(varname+"_bkg",titles);
	THStack* hsSig = new THStack(varname+"_sig",titles);
	hsumBkg->Reset();
	hsumSig->Reset();
	hData->Reset();
	hWsig->Reset();
	
	// legend...
	leg->Clear();
	
	for(rit=procs.rbegin() ; rit!=procs.rend() ; rit++)
	{
		TString name = rit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(!isSignal(name)) continue;
		if(isWsignal(name))
		{
			hWsig->Add(histos[hname]);
			
			TString ytitle = histos[hname]->GetYaxis()->GetTitle();
			TString label = (isSignal(name) && ytitle.Contains("Events")) ? labels[name]+"#times"+sBRf : labels[name];
			leg->AddEntry(histos1[hname],label,legoptions[name]);
		}
		
		hsSig->Add(histos[hname],"hist");
		hsumSig->Add(histos[hname]);
	}
	for(rit=procs.rbegin() ; rit!=procs.rend() ; rit++)
	{
		TString name = rit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(isData(rit->second))
		{
			hData->Add(histos[hname]);
			
			TString ytitle = histos[hname]->GetYaxis()->GetTitle();
			TString label = (isSignal(name) && ytitle.Contains("Events")) ? labels[name]+"#times"+sBRf : labels[name];
			leg->AddEntry(histos1[hname],label,legoptions[name]);
			
			continue;
		}
		if(isSignal(rit->second)) continue;
		
		hsBkg->Add(histos[hname],"hist");
		hsumBkg->Add(histos[hname]);
	}
	
	// Double_t maxBkg = getYmax(hsumBkg);
	// Double_t maxSig = getYmax(hsumSig);
	Double_t maxSig = getYmax(hWsig);
	Double_t maxDat = getYmax(hData);
	// Double_t max = (maxBkg>maxSig) ? maxBkg : maxSig;
	// max = (maxDat>max) ? maxDat : max;
	Double_t max = (maxDat>maxSig) ? maxDat : maxSig;
	
	TH1* htmp;
	TIter nextBhist((TList*)hsBkg->GetHists());
	while( (htmp=(TH1*)nextBhist())!=NULL )
	{
		if(ymin<0.) htmp->SetMinimum(0.);
		else        htmp->SetMinimum(ymin);
		if(ymax<0.) htmp->SetMaximum(max*2.5);
		else        htmp->SetMaximum(ymax);
	}
	TIter nextShist((TList*)hsSig->GetHists());
	while( (htmp=(TH1*)nextShist())!=NULL )
	{
		if(ymin<0.) htmp->SetMinimum(0.);
		else        htmp->SetMinimum(ymin);
		if(ymax<0.) htmp->SetMaximum(max*2.5);
		else        htmp->SetMaximum(ymax);
	}
	
	if(ymin<0.) hsBkg->SetMinimum(0.);
	else        hsBkg->SetMinimum(ymin);
	if(ymax<0.) hsBkg->SetMaximum(max*2.5);
	else        hsBkg->SetMaximum(ymax);
	if(ymin<0.) hsSig->SetMinimum(0.);
	else        hsSig->SetMinimum(ymin);
	if(ymax<0.) hsSig->SetMaximum(max*2.5);
	else        hsSig->SetMaximum(ymax);
	
	// the data
	c1->cd();
	hsBkg->Draw();
	hsSig->Draw("hist same");
	drawData(varname,histos,false);
	leg->Draw("same");
	ptxt->Draw("same");
	pad->cd();
	hsBkg->Draw();
	hsSig->Draw("hist same");
	drawData(varname,histos,false);
	leg->Draw("same");
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}
void drawPadSignal(TString varname, TVirtualPad* pad, TMapuiTS& procs, TMapTSP2TH1& histos/*, TLegend* leg*/)
{
	makeSingleCnv();
	pad->cd();
	
	leg->Clear();
	
	int counter = 0;
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		TString name = pit->second;
		if(drawchannels[name]==0) continue;
		
		if(!isSignal(name)) continue;
		TString hname = name+"_"+varname;
		histos[hname]->SetStats(1);
		
		TString ytitle = histos[hname]->GetYaxis()->GetTitle();
		TString label = (isSignal(name) && ytitle.Contains("Events")) ? labels[name]+"#times"+sBRf : labels[name];
		leg->AddEntry(histos1[hname],label,legoptions[name]);
		
		c1->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("same hist");
		
		pad->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("same hist");
		
		counter++;
	}
	c1->cd();
	leg->Draw("same");
	ptxt->Draw("same");
	pad->cd();
	leg->Draw("same");
	ptxt->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}

void divideCuts(TString varname_rat, TString varname_num, TString varname_den, TMapuiTS& procs, TMapTSP2TH1& histos)
{
	for(TMapuiTS::iterator pit=procs.begin() ; pit!=procs.end() ; pit++)
	{
		TString name         = pit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname_den = name+"_"+varname_den;
		TString hname_num = name+"_"+varname_num;
		TString hname_rat = name+"_"+varname_rat;
		histos[hname_rat]->Divide(histos[hname_num],histos[hname_den]);
	}
}
int getsizex(int ndivx)
{
	if     (ndivx==1) return 600;
	else if(ndivx==2) return 1200;
	else if(ndivx==3) return 1800;
	else _FAT("ndivx=4 is not supported");
	return 1000;
}
int getsizey(int ndivy)
{
	if     (ndivy==1) return 400;
	else if(ndivy==2) return 800;
	else if(ndivy==3) return 1200;
	else if(ndivy==4) return 1600;
	else _FAT("ndivy=5 is not supported");
	return 800;
}
void makeCnv(int ndivx, int ndivy, bool logy=false, TString parity="")
{
	cnv = new TCanvas("c","c",getsizex(ndivx),getsizey(ndivy));
	cnv->Divide(divx,divy);	
	pads.clear();
	for(int i=0 ; i<ndivx*ndivy ; i++)
	{
		pads.push_back(cnv->cd(i+1));
		pads[i]->SetTicks(1,1);
		if(logy)
		{
			if(parity=="" || (parity!="even" && parity!="odd")) pads[i]->SetLogy();
			else
			{
				if     (parity=="even" && i%2==0) pads[i]->SetLogy();
				else if(parity=="odd"  && i%2!=0) pads[i]->SetLogy();
			}
		}
	}
	cnv->Draw();
}
void closeCnv(TString fname)
{
	cnv->Update();
	cnv->SaveAs(fname);
	delete cnv;
}

void finalizeFigures(TString master, bool doMVA)
{
	// pdf file with all the plots
	pdffilename = (doMVA) ? "figures/mva."+master+".pdf" : "figures/cuts."+master+".pdf";
	const Int_t nlines = 20;
	TString lines[nlines];
	resetLines(nlines,lines);
	int pcounter = -1;
	
	
	// // legend...
	// TLegend* leg     = new TLegend(0.55,0.6,0.83,0.83,NULL,"brNDC");
	// // TLegend* legleft = new TLegend(0.15,0.6,0.43,0.83,NULL,"brNDC");
	// // TLegend* legsig  = new TLegend(0.15,0.6,0.43,0.83,NULL,"brNDC");
	// setLegendDefaults(leg);
	// // setLegendDefaults(legleft);
	// // setLegendDefaults(legsig);
	// for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	// {	
	// 	TString name = cit->second;
	// 	if(drawchannels[name]==0) continue;
	// 	
	// 	TString hname = name+"_m3body";
	// 	leg->AddEntry(histos1[hname],labels[name],legoptions[name]);
	// 	// legleft->AddEntry(histos1[hname],labels[name],legoptions[name]);
	// 	// if(isSignal(name)) legsig->AddEntry(histos1[hname],labels[name],legoptions[name]);
	// }
	
	
	makeAtlasLabel();
	makeLegend();

	// blank page
	_INF(visibility,"");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename+"(");
	/////////////////////////////////////////////////////////

	divx=2;
	divy=2;
	makeCnv(divx,divy);
	pcounter = -1;
	drawPadAll("triggers_absolute_after_vertexing", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("triggers_unique_after_vertexing",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("triggers_absolute_endOfSelection",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("triggers_unique_endOfSelection",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_norm_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_norm_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_norm_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_norm_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_normAB_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_normAB_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Wtaunu_3mu", "triggers_correlation_normAB_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("Data",       "triggers_correlation_normAB_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("bb_mu4mu4", "triggers_correlation_normAB_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("ccTomu15",  "triggers_correlation_normAB_after_vertexing", pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("bb_mu4mu4", "triggers_correlation_normAB_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	drawPadSingle("ccTomu15",  "triggers_correlation_normAB_endOfSelection",  pads[increment(pcounter)], channels, histos2, "col text0");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	divx=3;
	divy=2;
	makeCnv(divx,divy);
	pcounter = -1;
	drawPadAll("tripletCategories_afterVertexing",      pads[increment(pcounter)], channels, histos1, /*leg,*/0.);
	drawPadAll("tripletCategories_after_muons",         pads[increment(pcounter)], channels, histos1, /*leg,*/0.);
	drawPadAll("tripletCategories_after_triplet",       pads[increment(pcounter)], channels, histos1, /*leg,*/0.);
	drawPadAll("tripletCategories_after_ht",            pads[increment(pcounter)], channels, histos1, /*leg,*/0.);
	drawPadAll("tripletCategories_endOfSelection",      pads[increment(pcounter)], channels, histos1, /*leg,*/0.);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INF(visibility,"");
	divx=3;
	divy=2;
	makeCnv(divx,divy);
	pcounter = -1;
	drawPadAll("tripletCategories_norm_afterVertexing",      pads[increment(pcounter)], channels, histos1, /*leg,*/0);
	drawPadAll("tripletCategories_norm_after_muons",         pads[increment(pcounter)], channels, histos1, /*leg,*/0);
	drawPadAll("tripletCategories_norm_after_triplet",       pads[increment(pcounter)], channels, histos1, /*leg,*/0);
	drawPadAll("tripletCategories_norm_after_ht",            pads[increment(pcounter)], channels, histos1, /*leg,*/0);
	drawPadAll("tripletCategories_norm_endOfSelection",      pads[increment(pcounter)], channels, histos1, /*leg,*/0);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////


	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("pvalue",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("",        pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("Lxy",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("dLxy",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("a0xy",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("da0xy",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("PVNtrk_after_vertexing",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("PVNtrk_before_isolation", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("PVNtrk_after_isolation",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("PVNtrk_after_ht",         pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("trk_pt1_before_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta1_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi1_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_pt2_before_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta2_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi2_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_pt3_before_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta3_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi3_before_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("trk_pt1_after_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta1_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi1_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_pt2_after_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta2_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi2_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_pt3_after_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_eta3_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trk_phi3_after_muons", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "trk_hitstudy_eta1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trk_hitstudy_eta1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trk_hitstudy_eta2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trk_hitstudy_eta2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_pt_before_muons",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_eta_before_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_pt_before_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_eta_before_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_pt_before_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_eta_before_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_pt_after_muons",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_eta_after_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_pt_after_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_eta_after_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_pt_after_muons",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_eta_after_muons",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadAll("mu_pt",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("mu_eta",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_pt",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_eta",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPb_pt",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPb_eta",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "mu_hitstudy_pt_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_hitstudy_pt_offJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_hitstudy_eta_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_hitstudy_eta_offJpsi",  pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "TPa_hitstudy_pt_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_hitstudy_pt_offJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_hitstudy_eta_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_hitstudy_eta_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "TPb_hitstudy_pt_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_hitstudy_pt_offJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_hitstudy_eta_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_hitstudy_eta_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_isMedium1", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_isMedium2", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_isMedium3", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_isMedium",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "mu_isMedium_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_isMedium_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////


	
	
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_sctangsig", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_sctngbsig", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_pbalsig",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_pvalue",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_chi2ndf",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_qoverp",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "mu_sctangsig_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_sctangsig_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_pbalsig_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_pbalsig_offJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_qoverp_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_qoverp_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "mu_pvalue_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_pvalue_offJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_chi2ndf_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_chi2ndf_offJpsi",  pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	// divx=1;
	// divy=1;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadAll("mu_matchchi2ndf", pads[increment(pcounter)], channels, histos1/*, leg*/);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("TPa_chi2ndf",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_chi2ndf",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_pvalue",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_pvalue",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_qoverp",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_qoverp",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "TPa_chi2ndf_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_chi2ndf_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_pvalue_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_pvalue_offJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_qoverp_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_qoverp_offJpsi",  pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "TPb_chi2ndf_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_chi2ndf_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_pvalue_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_pvalue_offJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_qoverp_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPb_qoverp_offJpsi",  pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	divx=2;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	// drawPadAll("m2body",         pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("",               pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadSingle("Data", "m2body_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "m2body_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_MDTHits1",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_MDTHits1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_MDTHits2",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_MDTHits2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_MDTHits3",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_MDTHits3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_RPCHits1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_RPCHits2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_RPCHits3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_TGCHits1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_TGCHits2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("TPa_TGCHits3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_MDTHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_MDTHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_MDTHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_MDTHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_MDTHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_MDTHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_MDTHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_MDTHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,false);
	// drawPadSingle("Data", "mu_RPCHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_RPCHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_RPCHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_RPCHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,false);
	// drawPadSingle("Data", "TPa_RPCHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_RPCHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_RPCHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_RPCHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data", "mu_TGCHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_TGCHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_TGCHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "mu_TGCHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data", "TPa_TGCHits1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_TGCHits1_offJpsi", pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_TGCHits2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPa_TGCHits2_offJpsi", pads[increment(pcounter)], channels, histos1);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////




	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_outliersOnTrack1",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_outliersOnTrack1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_outliersOnTrack2",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_outliersOnTrack2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_outliersOnTrack3",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_outliersOnTrack3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_outliersOnTrack1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_outliersOnTrack1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_outliersOnTrack2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_outliersOnTrack2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_outliersOnTrack1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_outliersOnTrack1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_outliersOnTrack2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_outliersOnTrack2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////




	divx=2;
	divy=3;
	makeCnv(divx,divy,false,"even");
	pcounter = -1;
	drawPadAll("trks_fitprob_after_muons",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trks_pixdEdx_after_muons",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trks_fitprob_after_triplet",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trks_pixdEdx_after_triplet",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trks_fitprob_after_ht",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("trks_pixdEdx_after_ht",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "trks_hitstudy_fitprob_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trks_hitstudy_fitprob_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trks_hitstudy_pixdEdx_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "trks_hitstudy_pixdEdx_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////







	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_phiHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_phiHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_phiHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_phiHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_phiHits3",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_phiHits3",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_phiHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_phiHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_phiHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_phiHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_phiHits3",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_phiHits3",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",  "mu_RPC_vs_TGC_phiHits1_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_phiHits1_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_phiHits2_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_phiHits2_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_phiHits1_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_phiHits1_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_phiHits2_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_phiHits2_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_etaHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_etaHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_etaHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_etaHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_RPC_vs_TGC_etaHits3",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_RPC_vs_TGC_etaHits3",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_etaHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_etaHits1",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_etaHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_etaHits2",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "TPa_RPC_vs_TGC_etaHits3",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "TPa_RPC_vs_TGC_etaHits3",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",  "mu_RPC_vs_TGC_etaHits1_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_etaHits1_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_etaHits2_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "mu_RPC_vs_TGC_etaHits2_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_etaHits1_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_etaHits1_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_etaHits2_onJpsi",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",  "TPa_RPC_vs_TGC_etaHits2_offJpsi",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_precisionHits1",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionHits2",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionHits3",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionHoles1",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionHoles2",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionHoles3",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionOutliers1", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionOutliers2", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_precisionOutliers3", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("TPa_precisionHits1",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionHits2",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionHits3",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionHoles1",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionHoles2",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionHoles3",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionOutliers1", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionOutliers2", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_precisionOutliers3", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_precisionHits1_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHits1_offJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHits2_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHits2_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_precisionHits1_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHits1_offJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHits2_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHits2_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_precisionHoles1_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHoles1_offJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHoles2_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionHoles2_offJpsi",   pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_precisionHoles1_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHoles1_offJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHoles2_onJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionHoles2_offJpsi",   pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_precisionOutliers1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionOutliers1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionOutliers2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_precisionOutliers2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_precisionOutliers1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionOutliers1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionOutliers2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_precisionOutliers2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_phiLayers1",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_etaphiLayers1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_phiLayers2",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_etaphiLayers2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_phiLayers3",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_etaphiLayers3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("TPa_phiLayers1",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_etaphiLayers1",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_phiLayers2",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_etaphiLayers2",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_phiLayers3",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_etaphiLayers3",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_phiLayers1_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_phiLayers1_offJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_phiLayers2_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_phiLayers2_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_phiLayers1_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_phiLayers1_offJpsi",    pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_phiLayers2_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_phiLayers2_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "mu_etaphiLayers1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_etaphiLayers1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_etaphiLayers2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_etaphiLayers2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data", "TPa_etaphiLayers1_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_etaphiLayers1_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_etaphiLayers2_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_etaphiLayers2_offJpsi", pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_nHighThresholdTRTHits",           pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_nHighThresholdTRTHitsFraction",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mu_nHighThresholdTRTHitsRatio",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_nHighThresholdTRTHits",          pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_nHighThresholdTRTHitsFraction",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPa_nHighThresholdTRTHitsRatio",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_nHighThresholdTRTHits",          pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_nHighThresholdTRTHitsFraction",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("TPb_nHighThresholdTRTHitsRatio",     pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "mu_nHighThresholdTRTHits_onJpsi",           pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_nHighThresholdTRTHits_offJpsi",          pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_nHighThresholdTRTHitsFraction_onJpsi",   pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_nHighThresholdTRTHitsFraction_offJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_nHighThresholdTRTHitsRatio_onJpsi",      pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "mu_nHighThresholdTRTHitsRatio_offJpsi",     pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data", "TPa_nHighThresholdTRTHits_onJpsi",          pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_nHighThresholdTRTHits_offJpsi",         pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_nHighThresholdTRTHitsFraction_onJpsi",  pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_nHighThresholdTRTHitsFraction_offJpsi", pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_nHighThresholdTRTHitsRatio_onJpsi",     pads[increment(pcounter)], channels, histos1);
	drawPadSingle("Data", "TPa_nHighThresholdTRTHitsRatio_offJpsi",    pads[increment(pcounter)], channels, histos1);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////	
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHits_onJpsi",          pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHits_offJpsi",         pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHitsFraction_onJpsi",  pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHitsFraction_offJpsi", pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHitsRatio_onJpsi",     pads[increment(pcounter)], channels, histos1);
	// drawPadSingle("Data", "TPb_nHighThresholdTRTHitsRatio_offJpsi",    pads[increment(pcounter)], channels, histos1);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	
	
	
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadAll("met_et_after_hadclean",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_mt_et3mu_after_hadclean", pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_dphi3mu_after_hadclean",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_et_after_met",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_mt_et3mu_after_met", pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_dphi3mu_after_met",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_et_after_ht",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_mt_et3mu_after_ht", pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_dphi3mu_after_ht",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// divx=3;
	// divy=1;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadAll("met_et",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_mt_et3mu", pads[increment(pcounter)], channels, histos1/*, leg*/);
	// drawPadAll("met_dphi3mu",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////

	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("m3body_1600_beforeObjectcuts",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("",                              pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_1600_beforeWcuts",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("",                              pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_1600_afterLoosecuts",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("",                              pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("m3body_sigregion_onPhi_kaon_hypothesis_beforeObjectcuts",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_pion_hypothesis_beforeObjectcuts",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_kaon_hypothesis_beforeWcuts",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_pion_hypothesis_beforeWcuts",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_kaon_hypothesis_afterLoosecuts",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_pion_hypothesis_afterLoosecuts",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("m3body_lin_zoom_onPhi_beforeObjectcuts",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_lin_zoom_onPhi_beforeWcuts",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_lin_zoom_onPhi_afterLoosecuts",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onPhi_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("m3body_lin_zoom_onRhoOmega_beforeObjectcuts",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onRhoOmega_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_lin_zoom_onRhoOmega_beforeWcuts",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onRhoOmega_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_lin_zoom_onRhoOmega_afterLoosecuts",    pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("m3body_sigregion_onRhoOmega_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mT_mu1_onPhi_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onPhi_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onPhi_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu1_onPhi_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onPhi_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onPhi_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu1_onPhi_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onPhi_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onPhi_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mT_mu1_onRhoOmega_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onRhoOmega_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onRhoOmega_beforeObjectcuts", pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu1_onRhoOmega_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onRhoOmega_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onRhoOmega_beforeWcuts",      pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu1_onRhoOmega_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu2_onRhoOmega_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("mT_mu3_onRhoOmega_afterLoosecuts",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
		
	
	
	
	
	vector<TString> names;
	names.push_back("muons");
	names.push_back("triplet");
	names.push_back("ht");
	names.push_back("lowmassres");
	names.push_back("");          
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";	
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("pT3body"+thisname,      pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("met_et"+thisname,       pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("met_mt_et3mu"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("met_dphi3mu"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
		
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";	
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("mettrk_et"+thisname,       pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mettrk_mt_et3mu"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mettrk_dphi3mu"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mettrk_dphimet"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
		
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=3;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("ht_pt"+thisname,      pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("ht_dR3body"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("ht_mT"+thisname,      pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("ht_mT_track"+thisname,pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("ht_dphimet"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("ht_dphimet_track"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("m3body"+thisname,            pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("m3body_lin"+thisname,        pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("m3body_lin_zoom"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("m3body_sigregion"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("mOS"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mOS1"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mOS2"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mSS"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "mOS2_vs_mOS1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "mOS2_vs_mOS1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "pT3body_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "pT3body_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "mSS_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "m3body_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "m3body_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=3;
		divy=1;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("dROS1"+thisname, pads[increment(pcounter)], channels, histos1);
		drawPadAll("dROS2"+thisname, pads[increment(pcounter)], channels, histos1);
		drawPadAll("dRSS"+thisname,  pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=3;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "dROS1_vs_mOS1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "dROS1_vs_mOS1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "dROS2_vs_mOS2"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "dROS2_vs_mOS2"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "dRSS_vs_mSS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "dRSS_vs_mSS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	
	
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;		
		drawPadSingle("Wtaunu_3mu", "3rdTrk_pt_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "3rdTrk_pt_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "3rdTrk_PtOverPtOS_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "3rdTrk_PtOverPtOS_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	// for(unsigned int i=0 ; i<names.size() ; ++i)
	// {
	// 	TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
	// 	divx=2;
	// 	divy=2;
	// 	makeCnv(divx,divy,false);
	// 	pcounter = -1;		
	// 	drawPadSingle("Wtaunu_3mu", "3rdTrk_chi2dof_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
	// 	drawPadSingle("Data",       "3rdTrk_chi2dof_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
	// 	drawPadSingle("Wtaunu_3mu", "3rdTrk_pval_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
	// 	drawPadSingle("Data",       "3rdTrk_pval_vs_mOS"+thisname, pads[increment(pcounter)], channels, histos2, "col");
	// 	closeCnv(pdffilename);
	// }
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{	
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=3;
		makeCnv(divx,divy,false);
		pcounter = -1;		
		drawPadSingle("Wtaunu_3mu", "mSSsq_vs_mOSsq"+thisname,  pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "mSSsq_vs_mOSsq"+thisname,  pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "mSSsq_vs_mOS1sq"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "mSSsq_vs_mOS1sq"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "mSSsq_vs_mOS2sq"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "mSSsq_vs_mOS2sq"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}

	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=1;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("jet_n"+thisname,       pads[increment(pcounter)], channels, histos1, /*leg,*/ 0.);
		drawPadAll("jet_b"+thisname,       pads[increment(pcounter)], channels, histos1, /*leg,*/ 0.);
		closeCnv(pdffilename);
	}

	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=3;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadAll("jet_pt1_calib"+thisname,     pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("jet_pt2_calib"+thisname,     pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("jet_ptJ1J2_calib"+thisname,  pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("jet_E1_calib"+thisname,      pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("jet_E2_calib"+thisname,      pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("jet_EJ1J2_calib"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}

	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,true);
		pcounter = -1;
		drawPadAll("dPhi3bodyJ1"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("dPhiMETJ1"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("dPhiJ1J2"+thisname,    pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("mJ1J2"+thisname,       pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}

	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "dPhi3bodyJet1_vs_pTjet1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "dPhi3bodyJet1_vs_pTjet1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "dPhiJet1Jet2_vs_sumpTjet12"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "dPhiJet1Jet2_vs_sumpTjet12"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "MV1_vs_pTjet"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "MV1_vs_pTjet"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "MV1_vs_mjet"+thisname,  pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "MV1_vs_mjet"+thisname,  pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}	
	
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,true);
		pcounter = -1;
		drawPadAll("JVF1"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("JVF2"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("JVF3"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		drawPadAll("JVF4"+thisname,   pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=1;
		divy=1;
		makeCnv(divx,divy,true);
		pcounter = -1;
		drawPadAll("JVFall"+thisname, pads[increment(pcounter)], channels, histos1/*, leg*/);
		closeCnv(pdffilename);
	}
	for(unsigned int i=0 ; i<names.size() ; ++i)
	{
		TString thisname = (names[i]!="") ? "_after_"+names[i] : "";
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Wtaunu_3mu", "JVF1_vs_pTjet1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "JVF1_vs_pTjet1"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Wtaunu_3mu", "JVF2_vs_pTjet2"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data",       "JVF2_vs_pTjet2"+thisname, pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
	}
	
	
	divx=2;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("dRmin",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("dRmax",   pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("isolation000",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation001",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation002",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation003",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation004",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation005",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("isolation006",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation007",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation008",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation009",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation010",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation012",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("isolation014",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation016",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation018",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation020",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation022",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation024",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("isolation026",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation028",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("isolation030",       pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVAscore_all",            pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("",                        pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("MVAscore_below_tauMass",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("MVAscore_above_tauMass",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("MVAscore_left_sideband",  pads[increment(pcounter)], channels, histos1/*, leg*/);
	drawPadAll("MVAscore_right_sideband", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "m3body_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3body_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3body_vs_MVAscore",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "pT3body_vs_MVAscore",  pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "mOS_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mOS_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mSS_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu", "mOS1_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mOS1_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mOS2_vs_MVAscore",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	vector<TString> mvaevo;
	mvaevo.push_back("neg_veryloose");
	mvaevo.push_back("neg_loose");
	mvaevo.push_back("neg_medium");
	mvaevo.push_back("neg_tight");
	mvaevo.push_back("zero");
	mvaevo.push_back("pos_loose");
	mvaevo.push_back("pos_medium");
	mvaevo.push_back("pos_tight");
	mvaevo.push_back("pos_verytight");
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
		
		divx=3;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_mOS2_vs_mOS1"+suff,    pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data", "MVAevo_mSS_vs_mOS1"+suff,     pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data", "MVAevo_mSS_vs_mOS2"+suff,     pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data", "MVAevo_m3body_vs_mOS1"+suff,  pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data", "MVAevo_m3body_vs_mOS2"+suff,  pads[increment(pcounter)], channels, histos2, "col");
		drawPadSingle("Data", "MVAevo_m3body_vs_mSS"+suff,   pads[increment(pcounter)], channels, histos2, "col");
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,             pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                              pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_m3body_sigregion"+suff,  pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_pT3body"+suff,           pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,  pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_mSS"+suff,    pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_mOS1"+suff,   pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_mOS2"+suff,   pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,        pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                         pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_cal_et"+suff,   pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_track_et"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,        pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                         pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_cal_mT"+suff,   pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_track_mT"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,             pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_dphiCalTrk"+suff,    pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_cal_dPhi3mu"+suff,   pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_met_track_dPhi3mu"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,          pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_dr3body"+suff, pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_n"+suff,          pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_pt"+suff,      pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,            pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                             pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_mT_metcal"+suff, pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_mT_mettrk"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,             pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                              pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_dphimetcal"+suff, pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_jet_ht_dphimettrk"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	for(unsigned int i=0 ; i<mvaevo.size() ; ++i)
	{
		TString suff = "_"+mvaevo[i];
			
		divx=2;
		divy=2;
		makeCnv(divx,divy,false);
		pcounter = -1;
		drawPadSingle("Data", "MVAevo_score"+suff,                                  pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "",                                                   pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_m3body_sigregion_onPhi_kaon_hypothesis"+suff, pads[increment(pcounter)], channels, histos1);
		drawPadSingle("Data", "MVAevo_m3body_sigregion_onPhi_pion_hypothesis"+suff, pads[increment(pcounter)], channels, histos1);
		closeCnv(pdffilename);
		/////////////////////////////////////////////////////////
	}
	
	


	// blank page
	_INF(visibility,"");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos1/*, leg*/);
	closeCnv(pdffilename+")");
	/////////////////////////////////////////////////////////
}

#endif