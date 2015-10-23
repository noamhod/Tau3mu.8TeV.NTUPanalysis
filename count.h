#ifndef COUNT_H
#define COUNT_H


#include "std.h"
#include "type.h"

static const bool doLocalCuts = true;

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
	addCounter("nPassing_3body_charge", 50, "$|\\Sigma Q_{3{\\rm body}}|=1|LINE$");
	addCounter("nPassing_3body_mass",   51, "$\\mthreebody<5$~GeV");

	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_TrainingInc", 60, "Training+SB+Blinded: \\TrainingRegion~MeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_m2body",      61, "$\\mOSa>220$~MeV, $\\mOSb>220$~MeV and $\\mSS>220$~MeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_Ptrks",       62, "$\\ptrks>10^{-9}$");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_SLxy",        63, "$-10<\\SLxy<50$");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_Sa0xy",       64, "$\\Sazeroxy<25$");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_pT",          65, "$\\pTthreebody>10$~GeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_metcal",      66, "$10<\\calomet<250$~GeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_mettrk",      67, "$10<\\trackmet<250$~GeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_mtcal",       68, "$\\mTcalo>20$~GeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_mttrk",       69, "$\\mTtrack>20$~GeV");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_iso20",       70, "$\\isoTwenty<0.3$");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_loose_iso30",       71, "$\\isoThirty<1$|LINE");
	
	if(doLocalCuts && doMVA && !doMVAout && isBlinded)  addCounter("nPassing_tight_SB",          81, "SB: \\SidebandLeftBrackets+\\SidebandRightBrackets~MeV");
	if(doLocalCuts && doMVA && !doMVAout && !isBlinded) addCounter("nPassing_tight_SB",          81, "SB: \\SidebandsRegion~MeV");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_m2body",      82, "$\\mOSa>300$~MeV, $\\mOSb>300$~MeV and $\\mSS>300$~MeV");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_Ptrks",       83, "$\\ptrks>8\\times 10^{-9}$");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_SLxy",        84, "$1<\\SLxy<50$");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_mtcal",       85, "$\\mTcalo>45$~GeV");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_mttrk",       86, "$\\mTtrack>45$~GeV");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_pvalue",      87, "$\\pvalue>0.2$");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_dphiHTcal",   88, "$\\dphiHmetcalo>2$");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_dphiHTtrk",   89, "$\\dphiHmettrack>2$");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_rhoomegaphi", 90, "\\rhoomega and $\\phi$ veto");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_Ds",          91, "\\Ds veto|LINE");
	if(doLocalCuts && doMVA && !doMVAout)               addCounter("nPassing_tight_SR",          92, "SR: \\SignalRegion~MeV|LINE");
	                                                   
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_BDT_x0",            200, "BDT score $x>\\xzrovalue$|LINE");
	if(doLocalCuts && doMVA && !doMVAout) addCounter("nPassing_BDT_x1",            201, "BDT score $x>\\xoptvalue$|LINE");
	
	// //// cuts mode ?
	// if(!doMVA) addCounter("nPassing_3body_pt",          102, "$p_T^{3{\\rm body}}>25$~GeV (tight)");
	// if(!doMVA) addCounter("nPassing_3body_pvalue",      103, "3body vertex fit $\\pvalue>0.1$ (tight)");
	// if(!doMVA) addCounter("nPassing_3body_vtxclean",    104, "3body vertex geometry wrt PV (tight)");
	// if(!doMVA) addCounter("nPassing_3body_isolation",   105, "3body track isolation (tight)");
	// if(!doMVA) addCounter("nPassing_met_MET",           106, "Calo \\& Track $\\etmis>15$~GeV (tight)");
	// if(!doMVA) addCounter("nPassing_met_dphi3bodyMET",  107, "Calo \\& Track $\\dphithreebodymet>2$");
	// if(!doMVA) addCounter("nPassing_met_mT",            108, "Calo \\& Track $\\mT>60$~GeV (tight)");
	// if(!doMVA) addCounter("nPassing_met_dphiCalTrk",    109, "$\\dphimets<1.5$");
	// if(!doMVA) addCounter("nPassing_ht_ht",             110, "$\\HT>25$~GeV");
	// if(!doMVA) addCounter("nPassing_ht_dphihtMET",      111, "Calo \\& Track $\\Delta\\phi(\\HT,\\etmis)>2$");
	// if(!doMVA) addCounter("nPassing_ht_mhT",            112, "Calo \\& Track $\\mhT>60$~GeV (tight)");
	// if(!doMVA) addCounter("nPassing_ht_dr3body",        113, "$\\dRHThreeBody<1.5$|LINE");
	


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
