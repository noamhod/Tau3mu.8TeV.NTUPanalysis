//////////////////////////////////////////////////////////
//// root -b -l -q cutflow.C++\(\"muons\",\"cuts\"\) /////
//////////////////////////////////////////////////////////

#include "count.h"

ofstream* file;
TString skim1_dat_fname = "";
TString skim1_sig_fname = "";
TString skim1_bkg_fname = "";
TString skim2_fname = "";
TString skim3_fname = "";
TString outname     = "";
string texfilename  = "";

TMapTSP2TH1 histos1;

TString _s(double x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return (TString)str;
}
TString _s(double x, int prcn)
{
	stringstream strm;
	string str;
	strm << setprecision(prcn) << fixed << x;
	strm >> str;
	return (TString)str;
}
void add(TString text)
{
	(*file) << text << endl;
}
void fig(TString name, TString scale, TString label="", TString labelscale="", TString labelcoordinates="")
{
	add("\\begin{tikzpicture}");
	add("\t\\node{\\includegraphics[scale="+scale+"]{"+name+"}};");
	if(label!="") add("\t\\draw("+labelcoordinates+") node[scale="+labelscale+"]{"+label+"};");
	add("\\end{tikzpicture}");
}
void texheader()
{
	_INF(1,"Starting TeX file");
	
	file = new ofstream();
	file->open(texfilename.c_str());
	
	add("\\documentclass[xcolor=dvipsnames,table,final]{beamer}");
	
	add("\\newcommand{\\Sec}[3]{\\section{\\texorpdfstring{#1}{#2}}\\label{#3}}");
	add("\\newcommand{\\Subsec}[3]{\\subsection{\\texorpdfstring{#1}{#2}}\\label{#3}}");
	add("\\def\\vrb#1{\\texttt{\\detokenize{#1}}}");
	
	add("\\RequirePackage{pgfcore}");
	// add("\\pgfplotsset{compat=1.3}");

	add("\\usepackage{lmodern}");
	add("\\usepackage{graphicx}");
	add("\\usepackage{tikz}");
	add("\\usepackage{hyperref}");
	add("\\usepackage[nomessages]{fp}");
	add("\\usepackage{pgfplots}");
	add("\\usepackage{nicefrac,xcolor,verbatim,amsmath,bm,comment,url,epsfig,bm,wrapfig,Tabbing,wasysym,amssymb,axodraw,graphics,graphicx,appendix,multirow,epstopdf}");
	add("\\usepackage[abs]{overpic}");
	add("\\usepackage[english]{babel}");
	add("\\usepackage[latin1]{inputenc}");
	add("\\usepackage{times}");
	add("\\usepackage{etex}");
	add("\\usepackage{booktabs}");
	add("\\usepackage[toc]{multitoc}");
	
	add("\\usetikzlibrary{mindmap}");
	add("\\usetikzlibrary{shapes,decorations,arrows,patterns,automata}");
	add("\\usetikzlibrary{decorations.pathmorphing}");
	add("\\usetikzlibrary{calc,decorations.pathreplacing,decorations.markings}");
	add("\\usetikzlibrary{decorations.shapes}");
	add("\\usetikzlibrary{decorations.text}");
	add("\\usetikzlibrary{positioning}");
	add("\\usetikzlibrary{scopes}");
	add("\\usetikzlibrary{3d,calc}");
	add("\\usepackage{tikz-3dplot}");
	add("\\usetikzlibrary{shadings}");
	
	add("\\usefonttheme{professionalfonts}");
	
	add("\\mode<presentation>");
	add("{");
	add("\t\\usetheme{Montpellier}");
	add("\t\\usecolortheme{beaver}");
	add("\t\\setbeamercolor{title}{fg=Sepia}");
	add("\t\\setbeamercolor{frametitle}{fg=Sepia}");
	add("\t\\setbeamercolor{structure}{fg=Sepia}");
	add("\t\\setbeamercolor{itemize item}{fg=darkred!80!black}");
	add("\t\\setbeamercolor{block title}{fg=cictext, bg=cicbluebg}");
	add("\t\\setbeamerfont{block title}{size=\\large,series=\\bf}");
	add("\t\\setbeamercolor{block body}{parent=normal text,use=block title,bg=block title.bg!10!bg}");
	add("\t\\mode<beamer>{\\setbeamertemplate{blocks}[rounded][shadow=true]}");
	add("\t\\setbeamertemplate{section in toc}[sections numbered]");
	add("\t\\setbeamertemplate{section in toc}[circle]");
	add("\t\\setbeamercolor{title}{bg=white}");
	add("\t\\setbeamercolor{frametitle}{bg=white}");
	add("\t\\setbeamertemplate{bibliography item}[text]");
	add("\t\\author{Noam\\\\{\\scriptsize hod@cern.ch}}");
	add("\t\\title{$\\tau\\longrightarrow 3{\\rm body}$}");
	add("\t\\subtitle{Input for decisions}");
	// add("\t\\institute{\\pgfuseimage{NIKHEF-logo-large}}");
	add("}");
	
	add("\\definecolor{cicbluebg}{HTML}{660000}\%{B7E4FA}");
	add("\\definecolor{cicheadgrey}{HTML}{330000}\%{E0E7E2}");
	add("\\definecolor{cicline}{HTML}{000000}\%{003874}");
	add("\\definecolor{cictext}{HTML}{FFFFFF}\%{6E6E6F}");
	add("\\definecolor{cictitle}{HTML}{FFFFFF}\%{004080}");
	
	add("\\expandafter\\def\\expandafter\\insertshorttitle\\expandafter{\%");
	add("\t\\insertshorttitle\\hfill\%");
	add("\t\\insertframenumber\\,/\\,\\inserttotalframenumber}");
	add("\\setbeamertemplate{navigation symbols}{}");
	add("\\DeclareRobustCommand{\\optimalscorecut}{\\ensuremath{0.958}\\xspace}");
	
	add("\\begin{document}");
}
void table(TH1* hsig, TH1* hdat, TH1* hbkg, TString scale, bool isSkim2)
{
	vector<TString> cuts;
	vector<double>  sig;
	vector<double>  dat;
	vector<double>  bkg;
	vector<double>  eff_sig;
	vector<double>  rej_dat;
	vector<double>  rej_bkg;
	
	_INF(1,"isSkim2="<<isSkim2);
	
	
	
	for(Int_t b=1 ; b<=hsig->GetNbinsX() ; b++)
	{
		TString cutname = hsig->GetXaxis()->GetBinLabel(b);
		TString binname = cutname;
		// cutname.ReplaceAll("nPassing_","");
		TString texcutname = counters_label["nPassing_"+cutname];
		
		/////////////////////////////////////////////////////////////////////////
		_DBG(1,"b="<<b<<" binname="<<binname<<" texcutname="<<texcutname); ////////
		cuts.push_back(texcutname); /////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////
	}
	if(cuts.size()!=(unsigned int)hsig->GetNbinsX()) _ERR(1,"number of cuts do not match");
	
	_INF(1,"hsig->GetNbinsX()="<<hsig->GetNbinsX()<<", cuts.size()="<<cuts.size());
	
	for(unsigned int b=1 ; b<=cuts.size() ; b++)
	{	
		Double_t bPrev = -1.;
		Double_t bCurr = -1.;
		
		bPrev = hsig->GetBinContent(b-1); // if b=1, it'll give the underflow bin content
		bCurr = hsig->GetBinContent(b);
		if(b==1) eff_sig.push_back(-1);
		else     eff_sig.push_back((bCurr/bPrev)*100);
		sig.push_back(bCurr);
		
		bPrev = hdat->GetBinContent(b-1); // if b=1, it'll give the underflow bin content
		bCurr = hdat->GetBinContent(b);
		if(b==1) rej_dat.push_back(-1);
		else     rej_dat.push_back((1-bCurr/bPrev)*100);
		dat.push_back(bCurr);
		
		bPrev = hbkg->GetBinContent(b-1); // if b=1, it'll give the underflow bin content
		bCurr = hbkg->GetBinContent(b);
		if(b==1) rej_bkg.push_back(-1);
		else     rej_bkg.push_back((1-bCurr/bPrev)*100);
		bkg.push_back(bCurr);
	}
	
	add("\\hspace*{+1cm}");
	add("\\begin{tikzpicture}");
	add("\\node[scale="+scale+"]{");
	add("\\begin{tabular}[t]{l|rr|rr|rr}");
	add("\t\\hline\\hline");
	add("\tCut	& $W\\to\\nu\\tau\\to 3\\mu$ & Efficiency & $bb\\to\\mu 4\\mu 4$ & Rejection ($1-\\epsilon$) & Data & Rejection \\\\[3pt]");
	add("\t\\hline");		
	for(unsigned int i=0 ; i<cuts.size() ; i++)
	{
		bool addLine = (cuts[i].EndsWith("|LINE"));
		if(addLine) cuts[i].ReplaceAll("|LINE","");
		TString space = (i!=cuts.size()-1) ? "[3pt]" : "";
		TString effs  = (i<1) ? "--" : _s(eff_sig[i],2)+"\\%"; if(isnan(eff_sig[i])) effs = "--";
		TString rejd  = (i<1) ? "--" : _s(rej_dat[i],2)+"\\%"; if(isnan(rej_dat[i])) rejd = "--";
		TString rejb  = (i<1) ? "--" : _s(rej_bkg[i],2)+"\\%"; if(isnan(rej_bkg[i])) rejb = "--";
		add("\t"+cuts[i]+" &"+_s(sig[i])+" &"+effs+" &"+_s(bkg[i])+" &"+rejb+" &"+_s(dat[i],0)+" &"+rejd+" \\\\"+space);
		if(addLine) add("\t\\hline");
	}
	add("\t\\hline\\hline");
	add("\\end{tabular}};");
	add("\\end{tikzpicture}");
}
void texbody(TString master)
{
	_INF(1,"Filling TeX file");
	
	string master_channel = (master=="cuts") ? "Cut-based" : "MVA-based";
	
	add("\\begin{frame}[noframenumbering]");
	add("\\title{\\Huge{$\\tau\\longrightarrow 3{\\rm body}$ cutflow}}");
	add("\\subtitle{\\Large{"+master_channel+"}}");
	add("\\date{\\today}");
	add("\\titlepage");
	add("\\end{frame}");

	// add("\\Sec{Skim1 cutflow}{Skim1 cutflow}{sec:skim1cutflow}");
	// add("\\begin{frame}{}");
	// add("\\vspace*{-0.1cm}");
	// add("\\begin{columns}");
	// add("\\column{12cm}");
	// add("\\hspace*{+1cm}");
	// fig(outname+".skim1.eps", "0.5");
	// add("\\end{columns}");
	// add("\\end{frame}");
		
	add("\\Sec{Individual absolute cutflow}{Individual absolute cutflow}{sec:individualabsolutecutflow}");
	add("\\begin{frame}{}");
	add("\\vspace*{-0.1cm}");
	add("\\begin{columns}");
	add("\\column{12cm}");
	add("\\hspace*{+1cm}");
	fig(outname+".all.eps", "0.5");
	add("\\end{columns}");
	add("\\end{frame}");
		
	add("\\Sec{Weighted cutflow for 20~fb$^{-1}$}{Weighted cutflow for 20~fb$^{-1}$}{sec:weightedcutflow}");
	add("\\begin{frame}{}");
	add("\\hspace*{-1cm}");
	add("\\begin{columns}");
	add("\\column{12cm}");
	fig(outname+".weighted.eps", "0.6");
	add("\\end{columns}");
	add("\\end{frame}");
		
	add("\\Sec{Normalized cutflow}{Normalized cutflow}{sec:normalizedcutflow}");
	add("\\begin{frame}{}");
	add("\\hspace*{-1cm}");
	add("\\begin{columns}");
	add("\\column{12cm}");
	fig(outname+".normalized.eps", "0.6");
	add("\\end{columns}");
	add("\\end{frame}");
	
	// double nalldat   = 726383357;
	// double nallsig   = 99000;
	// double nskim1dat = histos1["skim2_Data_cutflow_absolute"]->GetBinContent(1);
	// double nskim1sig = histos1["skim2_Wtaunu_3mu_cutflow_absolute"]->GetBinContent(1);
	// add("\\Sec{Skim1}{Skim1}{sec:skim1}");
	// add("\\begin{frame}{}");
	// add("\\begin{columns}");
	// add("\\column{12cm}");
	// add("\\begin{itemize}");
	// add("\t\\item In Skim1 I make the vertex fit on AODs while requireing:");
	// add("\t\\begin{itemize}");
	// add("\t\t\\item GRL");
	// add("\t\t\\item Trigger");
	// add("\t\t\\item At least one muon triplet with all 3 muons having $p_T>2$~GeV and $|\\eta|<2.7$");
	// add("\t\t\\item Muons can be combined, segment tagged or calorimeter tagged");
	// add("\t\t\\item Allowed combinations...");
	// add("\t\t\\item Require $m_{3{\\rm body}}<5$~GeV but do NOT constrain to $|\\Sigma Q_{3{\\rm body}}|=1$");
	// add("\t\t\\item At least one reconstructed triplet vertex (successful fit)");
	// add("\t\\end{itemize}");
	// add("\t\\item Events for 20.2~fb$^{-1}$ (AMI): "+_s(nalldat,0));
	// add("\t\\item Events remaining after skim1 (on AOD): "+_s(nskim1dat,0));
	// add("\t\\item Skim1 rejection for data is: "+_s((1-nskim1dat/nalldat)*100)+"\\%");
	// add("\t\\item Skim1 efficiency for signal is: "+_s((nskim1sig/nallsig)*100)+"\\%");
	// add("\t\\item The output is b-physics formatted ntuples (not D3PDs)");
	// add("\\end{itemize}");
	// add("\\end{columns}");
	// add("\\end{frame}");
	// 
	// add("\\Sec{Skim2}{Skim2}{sec:skim2}");
	// add("\\begin{frame}{}");
	// // add("\\vspace*{-0.3cm}");
	// // add("\\hspace*{+1.5cm}");
	// add("\\begin{columns}");
	// add("\\column{12cm}");
	// add("\\begin{itemize}");
	// add("\t\\item Skim2 input is the b-physics ntuples obtained while running skim1 on AODs");
	// add("\t\\item In this skim level I just reduce the sizes without changing the format");
	// add("\t\\item This is done while enforcing a $W$ enriched sample");
	// add("\t\\item All $3{\\rm body}$ quantities below are taken from the vertex fit");
	// add("\t\\item $E^{\\rm{miss}}_T$ here is RefFinal");
	// add("\t\\item Remember that the acceptance is $\\sim 50\\%$. This is why there's a big difference between the Trigger and the $m_{3{\\rm body}}$ requirements. For $m_{3{\\rm body}}$, have to have 3 reconstructed muons");
	// add("\\end{itemize}");
	// table(histos1["skim2_Wtaunu_3mu_cutflow_absolute"],histos1["skim2_Data_cutflow_absolute"],histos1["skim2_bb_mu4mu4_cutflow_absolute"],"0.4",true);
	// add("\\end{columns}");
	// add("\\end{frame}");	
	
	add("\\Sec{Analysis level}{Analysis level}{sec:analysislevel}");
	add("\\begin{frame}{}");
	add("\\vspace*{-0.3cm}");
	add("\\hspace*{+1.5cm}");
	add("\\begin{columns}");
	add("\\column{12cm}");
	table(histos1["skim3_Wtaunu_3mu_cutflow_absolute"],histos1["skim3_Data_cutflow_absolute"],histos1["skim3_bb_mu4mu4_cutflow_absolute"],"0.4",false);
	add("\\end{columns}");
	add("\\end{frame}");
}
void texfooter()
{
	_INF(1,"Closing TeX file");
	
	add("\\end{document}");
	file->close();
	cout << "Cutflow tex file was written to: " << texfilename << endl;
}


void addHist(TString name, TString skim, TFile* f)
{
	TString newname = skim+"_"+name;
	histos1.insert(make_pair(newname,(TH1*)f->Get(name)->Clone(newname)));
	if(name.Contains("absolute") || name.Contains("weighted")) histos1[newname]->SetMinimum(0.5);
}

void addHist(TString name, TH1* h)
{
	histos1.insert(make_pair(name,h));
}


TH1* getSkim1CutFlow(TString path)
{
	void* dir = gSystem->OpenDirectory(path);
	if(!dir) _FAT("Directory: "<<path<<" could not be opened");
	
	TString sample = gSystem->GetDirEntry(dir);
	TString fullpath = path+"/"+sample+"/";
	dir = gSystem->OpenDirectory(fullpath);
	if(!dir) _FAT("Directory: "<<fullpath<<" could not be opened");
	
	Int_t nFiles = 0;
	TString pattern = ".root";

	char* ent;
	vector<TString> vfiles;
	while((ent=(char*)gSystem->GetDirEntry(dir)))
	{
		TString fn = Form("%s%s",fullpath.Data(),ent);
		if(fn.Contains(pattern))
		{
			// FileStat_t st;
			// if(!gSystem->GetPathInfo(fn,st) && R_ISREG(st.fMode)) chain->Add(fn);
			vfiles.push_back(fn);
			nFiles++;
			
		}
	}
	if(!nFiles) _FAT("Directory: "<<fullpath<<" contains no files with "<<pattern<<" pattern in their name");
	////////////////////////
	//// get the histos ////
	////////////////////////
	TH1*   h  = NULL;
	cout << "------------------------------------------------------------------" << endl;
	cout << "Path:    " << path   << endl;
	cout << "Sample:  " << sample << endl;
	cout << "N files: " << nFiles << endl;
	for(Int_t f=0 ; f<(Int_t)vfiles.size() ; f++)
	{
		cout << "  [" << f+1 << "]" << vfiles[f] << endl;
		TFile* tf = TFile::Open(vfiles[f]);
		if(f==0)
		{
			h = (TH1*)((TH1*)tf->Get("skim1_cutflow")->Clone(sample+"_skim1_cutflow"));
			h->Reset();
		}
		h->Add((TH1*)tf->Get("skim1_cutflow"));
	}
	cout << "------------------------------------------------------------------" << endl;
	gSystem->FreeDirectory(dir);
	return h;
}


void cutflow(TString master = "muons", TString method = "cuts")
{
	// skim2_fname = "/afs/cern.ch/user/h/hod/data/bout/../figs/fig.all."+master+".root";
	// gSystem->Exec("cp -f "+skim2_fname+" .");
	// TFile* fInskim2 = new TFile(skim2_fname,"READ");
	
	skim3_fname = "histos."+master+"."+method+".root"; // histos.muons.cuts.root
	TFile* fInSkim3 = new TFile(skim3_fname,"READ");
	
	outname = "cutflow."+master+"."+method;
	texfilename = (string)outname+".tex";
	
	bool isMVA = (method=="mva");
	initCounters(isMVA,false,true);
	
	// 
	// TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/mrg.bphys/";
	// vector<TString> paths;
	// // paths.push_back(basepath+"data/periodA/");
	// // paths.push_back(basepath+"data/periodB/");
	// paths.push_back(basepath+"data/periodC/");
	// // paths.push_back(basepath+"data/periodD/");
	// paths.push_back(basepath+"data/periodE/");
	// paths.push_back(basepath+"data/periodG/");
	// paths.push_back(basepath+"data/periodH/");
	// paths.push_back(basepath+"data/periodI/");
	// // paths.push_back(basepath+"data/periodJ/");
	// paths.push_back(basepath+"data/periodL/");
	// // paths.push_back(basepath+"data/periodM/");
	// TH1* hSikm1_data = NULL;
	// for(unsigned int p=0 ; p<paths.size() ; p++)
	// {
	// 	if(p==0) hSikm1_data = (TH1*)getSkim1CutFlow(paths[p])->Clone("Data_skim1_cutflow");
	// 	else     hSikm1_data->Add((TH1*)getSkim1CutFlow(paths[p])->Clone(paths[p]+"_skim1_cutflow"));
	// }
	// TH1* hSikm1_sig = (TH1*)getSkim1CutFlow(basepath+"mc/Wtaunu_3mu/")->Clone("Wtaunu_3mu_skim1_cutflow");
	// TH1* hSikm1_bkg = (TH1*)getSkim1CutFlow(basepath+"mc/bb_mu4mu4/")->Clone("bb_mu4mu4_skim1_cutflow");
	
	
	
	
	
	
	// addHist("Data_skim1_cutflow",       hSikm1_data);
	// addHist("Wtaunu_3mu_skim1_cutflow", hSikm1_sig);
	// addHist("bb_mu4mu4_skim1_cutflow",  hSikm1_bkg);
	
	// addHist("Data_cutflow_absolute",            "skim2", fInskim2);
	// addHist("bb_mu4mu4_cutflow_absolute",       "skim2", fInskim2);
	// addHist("bbTotau10_3mu_cutflow_absolute",   "skim2", fInskim2);
	// addHist("Wtaunu_3mu_cutflow_absolute",      "skim2", fInskim2);
	// addHist("Data_cutflow_weighted",            "skim2", fInskim2);
	// addHist("bb_mu4mu4_cutflow_weighted",       "skim2", fInskim2);
	// addHist("bbTotau10_3mu_cutflow_weighted",   "skim2", fInskim2);
	// addHist("Wtaunu_3mu_cutflow_weighted",      "skim2", fInskim2);
	// addHist("Data_cutflow_normalized",          "skim2", fInskim2);
	// addHist("bb_mu4mu4_cutflow_normalized",     "skim2", fInskim2);
	// addHist("bbTotau10_3mu_cutflow_normalized", "skim2", fInskim2);
	// addHist("Wtaunu_3mu_cutflow_normalized",    "skim2", fInskim2);
	
	addHist("Data_cutflow_absolute",            "skim3", fInSkim3);
	addHist("bb_mu4mu4_cutflow_absolute",       "skim3", fInSkim3);
	addHist("bbTotau10_3mu_cutflow_absolute",   "skim3", fInSkim3);
	addHist("Wtaunu_3mu_cutflow_absolute",      "skim3", fInSkim3);
	addHist("Data_cutflow_weighted",            "skim3", fInSkim3);
	addHist("bb_mu4mu4_cutflow_weighted",       "skim3", fInSkim3);
	addHist("bbTotau10_3mu_cutflow_weighted",   "skim3", fInSkim3);
	addHist("Wtaunu_3mu_cutflow_weighted",      "skim3", fInSkim3);
	addHist("Data_cutflow_normalized",          "skim3", fInSkim3);
	addHist("bb_mu4mu4_cutflow_normalized",     "skim3", fInSkim3);
	addHist("bbTotau10_3mu_cutflow_normalized", "skim3", fInSkim3);
	addHist("Wtaunu_3mu_cutflow_normalized",    "skim3", fInSkim3);
	
	
	
	TCanvas* cnv = NULL;
	
	// cnv = new TCanvas("cutflow_skim1", "", 600,400);
	// cnv->SetLogy();
	// histos1["Data_skim1_cutflow"]->Draw("hist text");
	// histos1["bb_mu4mu4_skim1_cutflow"]->Draw("hist same");
	// histos1["Wtaunu_3mu_skim1_cutflow"]->Draw("hist text same");
	// cnv->SaveAs(outname+".cutflow_skim1.eps");
	// cnv->SaveAs(outname+".cutflow.pdf");
	// 
	// delete cnv;
	cnv = new TCanvas("cutflowall", "", 1200,1000);
	vector<TVirtualPad*> pads;
	cnv->Divide(2,2);
	pads.push_back(cnv->cd(1));
	pads.push_back(cnv->cd(2));
	pads.push_back(cnv->cd(3));
	pads.push_back(cnv->cd(4));
	for(unsigned int i=0 ; i<pads.size() ; i++) { pads[i]->SetGridx(); pads[i]->SetLogy(); }
	pads[0]->cd();
	histos1["skim3_Data_cutflow_absolute"]->Draw("hist text");
	pads[1]->cd();
	histos1["skim3_bb_mu4mu4_cutflow_absolute"]->Draw("hist text");
	pads[2]->cd();
	histos1["skim3_Wtaunu_3mu_cutflow_absolute"]->Draw("hist text");
	pads[3]->cd();
	histos1["skim3_bbTotau10_3mu_cutflow_absolute"]->Draw("hist text");
	cnv->SaveAs(outname+".all.eps");
	cnv->SaveAs(outname+".pdf(");

	delete cnv;
	cnv = new TCanvas("cutflownweighted", "", 600,400);
	cnv->SetLogy();
	histos1["skim3_Data_cutflow_weighted"]->Draw("hist text");
	histos1["skim3_bb_mu4mu4_cutflow_weighted"]->Draw("hist same");
	histos1["skim3_bbTotau10_3mu_cutflow_weighted"]->Draw("hist same");
	histos1["skim3_Wtaunu_3mu_cutflow_weighted"]->Draw("hist text same");
	cnv->SaveAs(outname+".weighted.eps");
	cnv->SaveAs(outname+".pdf");

	delete cnv;
	cnv = new TCanvas("cutflownorm", "", 600,400);
	cnv->SetLogy();
	histos1["skim3_Data_cutflow_normalized"]->Draw("hist");
	histos1["skim3_bb_mu4mu4_cutflow_normalized"]->Draw("hist same");
	histos1["skim3_bbTotau10_3mu_cutflow_normalized"]->Draw("hist same");
	histos1["skim3_Wtaunu_3mu_cutflow_normalized"]->Draw("hist same");
	cnv->SaveAs(outname+".normalized.eps");
	cnv->SaveAs(outname+".pdf)");

	TString foname = outname+".root";
	TFile* fOut = new TFile(foname, "RECREATE");
	for(TMapTSP2TH1::iterator it=histos1.begin() ; it!=histos1.end() ; ++it) it->second->Write();
	fOut->Write();
	fOut->Close();
	
	texheader();
	texbody(master);
	texfooter();
	
	// gSystem->Exec("cp -f "+texfilename+"  tex/cutflow."+master+".tex");
	// gSystem->Exec("cd tex/; make clean; make "+master+"; cd ../");
	// cout << "See tex/cutflow."+master+".pdf" << endl;
}