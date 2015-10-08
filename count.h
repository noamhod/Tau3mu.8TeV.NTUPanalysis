#ifndef COUNT_H
#define COUNT_H


#include "std.h"
#include "type.h"

TMapTSui counters;
TMapTSTS counters_label;
TMapTSf  counters_weighted;
TMapTSui counters_evtvisited;
TMapiTS  counters_ordered;

int increment(int& counter)
{
	counter += 1;
	return counter;
}
string str(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm << fixed << x;
	strm >> str;
	return str;
}
void addCounter(TString name, int order, TString label="")
{
	counters.insert(make_pair(name,0));
	counters_label.insert(make_pair(name,label));
	counters_weighted.insert(make_pair(name,0.));
	counters_evtvisited.insert(make_pair(name,0));
	counters_ordered.insert(make_pair(order,name));
}
void initCounters(bool doMVA, bool doMVAout, bool isBlinded)
{	
	//////////////////////////////////////////
	//// cutflow counters (have order>=0) ////
	//////////////////////////////////////////
	
	/*
	if(cutname=="3body_sidebands")  texcutname = "Sidebands: $1500<\\mthreebody<2100$~MeV";
	if(cutname=="3body_blinded")    texcutname = "Blinded: $1700<\\mthreebody<1860$~MeV|LINE";
	if(cutname=="evt_onetriplet")    texcutname = "Exactly one candidate (cut-based only)";
	
	if(cutname=="MVA_optimal")       texcutname = ;
	*/
	
	
	
	//// event preselection
	addCounter("nPassing_evt_all",           0, "All - after skim2");
	addCounter("nPassing_evt_trigger",       1, "Trigger - after skim2");
	addCounter("nPassing_evt_oneCandidate",  2, "Trigger - after skim2");
	addCounter("nPassing_evt_goodtriplets",  3, "Category is $3\\mu$ or $2\\mu+1\\tpa$|LINE");
	
	//// single muon
	addCounter("nPassing_mu_mcp",        20, "Std. MCP ID-hits cut");
	addCounter("nPassing_mu_hits",       21, "TRT+MS hits");
	addCounter("nPassing_mu_pt",         22, "$p_T^{1,2,3}>5.5,3.5,2.5$~GeV");
	addCounter("nPassing_mu_eta",        23, "$|\\eta|<2.5$");
	addCounter("nPassing_mu_sctangsig",  24, "Muon scattering angle sig.$<5$");
	addCounter("nPassing_mu_pbalsig",    25, "Muon momentum balance sig.$<3$");
	addCounter("nPassing_mu_trkquality", 26, "\\tpa track-fit $\\pvalue>0.01$|LINE");
	
	//// triplet cuts
//	addCounter("nPassing_3body_trksprobLoose",     50, "3body track fit probability poduct $\\mathcal{P}>10^{-9}$ (loose)");
	//////addCounter("nPassing_3body_pvalueLoose", 51, "3body vertex fit $\\pvalue>10^{-6}$ (loose)");
//  addCounter("nPassing_3body_vtxcleanLoose",     52, "3body vertex geometry (loose)");
	//////addCounter("nPassing_3body_m",           53, "$\\mthreebody<4$~GeV");
//	addCounter("nPassing_3body_ptLoose",           54, "$\\pTthreebody>10$~GeV (loose)");
	addCounter("nPassing_3body_charge",            55, "$|\\Sigma Q_{3{\\rm body}}|=1$");
//	addCounter("nPassing_3body_m2body",            56, "$m_{2{\\rm body}}>220$~MeV|LINE");
	//////addCounter("nPassing_3body_dalitz",      57, "Peripheral Dalitz cut ($\\mSS^{2}$ vs $\\mOS^{2}$)|LINE");

	//// W-oriented cuts
//	addCounter("nPassing_3body_isolationLoose",   70, "3body track isolation (loose)");
//	addCounter("nPassing_met_metLoose",           71, "$\\etmis>10$~GeV (loose)");
//	addCounter("nPassing_met_mTLoose",            72, "$\\mT>20$~GeV (loose)|LINE");
	//////addCounter("nPassing_ht_htLoose",       73);
	//////addCounter("nPassing_ht_mhTLoose",      74, "$\\mhT>20$~GeV (loose)|LINE");
	
	//// cuts mode ?
	// if(!doMVA) addCounter("nPassing_jets_bjetveto_std", 100);
	// if(!doMVA) addCounter("nPassing_jets_coljetveto",   101, "Collinear jet veto (tight)");
	if(!doMVA) addCounter("nPassing_3body_pt",          102, "$p_T^{3{\\rm body}}>25$~GeV (tight)");
	if(!doMVA) addCounter("nPassing_3body_pvalue",      103, "3body vertex fit $\\pvalue>0.1$ (tight)");
	if(!doMVA) addCounter("nPassing_3body_vtxclean",    104, "3body vertex geometry wrt PV (tight)");
	if(!doMVA) addCounter("nPassing_3body_isolation",   105, "3body track isolation (tight)");
	if(!doMVA) addCounter("nPassing_met_MET",           106, "Calo \\& Track $\\etmis>15$~GeV (tight)");
	if(!doMVA) addCounter("nPassing_met_dphi3bodyMET",  107, "Calo \\& Track $\\dphithreebodymet>2$");
	if(!doMVA) addCounter("nPassing_met_mT",            108, "Calo \\& Track $\\mT>60$~GeV (tight)");
	if(!doMVA) addCounter("nPassing_met_dphiCalTrk",    109, "$\\dphimets<1.5$");
	if(!doMVA) addCounter("nPassing_ht_ht",             110, "$\\HT>25$~GeV");
	if(!doMVA) addCounter("nPassing_ht_dphihtMET",      111, "Calo \\& Track $\\Delta\\phi(\\HT,\\etmis)>2$");
	if(!doMVA) addCounter("nPassing_ht_mhT",            112, "Calo \\& Track $\\mhT>60$~GeV (tight)");
	if(!doMVA) addCounter("nPassing_ht_dr3body",        113, "$\\dRHThreeBody<1.5$|LINE");
	
	if(doMVA && !doMVAout) addCounter("nPassing_optMVA",  300, "MVA score $>\\optimalscorecut$|LINE");
	
	//// Everything eventually need to have only one candidate at the end of the selsction
	if(!doMVAout)              addCounter("nPassing_3body_sidebands",   400, "Sidebands: \\mthreebody in \\SidebandsRegion~MeV");
	if(!doMVAout && isBlinded) addCounter("nPassing_3body_blinded",     401, "Blinded: \\mthreebody in \\BlindedRegion~MeV");
	// if(!doMVAout)              addCounter("nPassing_evt_rhoomegaphi",   402, "$\\rhoomega,\\phi$ removal: $\\left|\\mOS-m_{\\rho/\\omega/\\phi}\\right|<30$~MeV");
	// if(!doMVAout)              addCounter("nPassing_evt_onecandidate",  403, "Exactly one candidate per event|LINE");

	/////////////////////////////////////
	/// other counters (have order<0) ///
	/////////////////////////////////////
}
bool isCounter(TString name)
{
	return (counters.find(name)!=counters.end());
}
void clearCounters()
{
	counters.clear();
	counters_weighted.clear();
	counters_evtvisited.clear();
	counters_ordered.clear();
}
void resetCounterFlags()
{
	for(TMapTSui::iterator it=counters_evtvisited.begin() ; it!=counters_evtvisited.end() ; it++)
	{
		it->second = 0;
	}
}
void incrementCounter(TString name, float weight=1.)
{	
	if(!isCounter(name)) return;
	
	if(!counters_evtvisited[name])
	{	
		counters[name]++;
		counters_weighted[name]+=weight;
		counters_evtvisited[name] = 1;
	}
}
int getCounter(TString name, bool doWeighted=false)
{
	if(!isCounter(name)) return -9999;
	
	if(doWeighted) return counters_weighted[name];
	return counters[name];
}
string printCounters(TString type, TString pname, bool doWeighted=false)
{
	string allLines = "";
	int maxlength = 0;
	int maxcounterlength = 0;
	for(TMapTSui::iterator it = counters.begin() ; it!=counters.end() ; it++)
	{
		TString name = it->first;
		int icounter = counters[name];
		float fcounter = counters_weighted[name];
		TString counter = (doWeighted) ? (TString)str(fcounter,3) : (TString)str(icounter);
		maxlength        = (name.Length()>maxlength)           ? name.Length()    : maxlength;
		maxcounterlength = (counter.Length()>maxcounterlength) ? counter.Length() : maxcounterlength;
	
		// cout << name << " -> " << name.Length() << ", maxlength=" << maxlength << endl;
		// cout << counter << " -> " << counter.Length() << ", maxcounterlength=" << maxcounterlength << endl;
	}
	
	string header = "";
	if     (type=="cutflow" && !doWeighted) header = "================== CUTFLOW "+ pname +" =================";
	else if(type=="cutflow" && doWeighted)  header = "================== WEIGHTED-CUTFLOW "+ pname +" =================";
	else if(type=="objects" && !doWeighted) header = "================== OBJECTS "+pname+" =================";
	else if(type=="objects" && doWeighted)  header = "================== WEIGHTED-OBJECTS "+pname+" =================";
	else _ERR(1,"Unsupported option: "<<type<<" for printCounters(TString,TString,bool) method");
	cout << header << endl;
	allLines += header+"\n";
	
	
	TString previousname = "";
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; it++)
	{
		TString name = it->second;
		int order    = it->first;
		
		int icounter = counters[name];
		float fcounter = counters_weighted[name];
		
		TString thiscounter = (doWeighted) ? (TString)str(fcounter,3) : (TString)str(icounter);
		int thislength = thiscounter.Length();
		int namelength = it->second.Length();
		if(type=="cutflow" && order<0)  continue;
		if(type=="objects" && order>=0) continue;
		
		string spaces = "";
		for(int j=1 ; j<(maxlength-namelength+1) ; j++) spaces+=" ";
		string digits = "";
		for(int j=1 ; j<(maxcounterlength-thislength+1) ; j++) digits+=" ";
		
		// cout << "maxlength-namelength=" << maxlength-namelength << endl;
		// cout << "maxcounterlength-thislength=" << maxcounterlength-thislength << endl;
		
		if(previousname=="") previousname = name; // only for the 1st cut
		float thiscut = (doWeighted) ? fcounter : icounter; // counters[name];
		float prevcut = (doWeighted) ? counters_weighted[previousname] : counters[previousname];
		float reldiff = 100*(prevcut-thiscut)/prevcut;
		
		string sthiscut = (doWeighted) ? str(thiscut,3) : str(thiscut,0);
		string cutline = (string)name+" "+spaces+digits+sthiscut+" [-"+str(reldiff,1)+"\%]";
		cout << cutline << endl;
		allLines += cutline+"\n";
		
		previousname = name;
	}
	string footer = "";
	if(!doWeighted)  footer = "================== CUTFLOW "+pname+" =================";
	else             footer = "=========================== "+pname+" ========================";
	cout << footer << endl;
	allLines += footer+"\n";
	
	return allLines;
}

#endif
