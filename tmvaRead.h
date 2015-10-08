/////////////////////////////////////////////
// usage:
// 1. call initTMVA (on init only)
// 2. call bookMVAvars (on init only)
// 3. call setMVAvars per event
// 4. call setMVAspect per event
// 4. call getMVAscore per event
////////////////////////////////////////////

#include "std.h"
#include "type.h"

// #if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
// #endif

using namespace TMVA;

// #ifdef __CINT__
// 	gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
// #endif

TMVA::Reader* reader; // The Reader object

TMapTSf mva_variables;

TMapiTS mva_var_order;

TMapTSTS mva_var_type;

TMapTSf mva_spectators;
TMapTSTS mva_spe_type;
TMapiTS mva_spe_order;

TString methodName;


void addVariable(TString name, TString type, int& order)
{
	mva_variables.insert(make_pair(name,-999.));
	mva_var_order.insert(make_pair(order,name));
	mva_var_type.insert(make_pair(name,type));
	order++;
}


void addSpectator(TString name, TString type, int& order)
{
	mva_spectators.insert(make_pair(name,-999.));
	mva_spe_order.insert(make_pair(order,name));
	mva_spe_type.insert(make_pair(name,type));
	order++;
}

void clearMVAvars()
{
	mva_variables.clear(); mva_var_order.clear(); mva_var_type.clear();
	
	mva_spectators.clear();
	mva_spe_order.clear();
	mva_spe_type.clear();
}


void bookMVAvars(TString metType)
{
	int order = -1;
	
	//// MVA variables
	order = -1;

	// addVariable("vtx_isolation003", "VF", order);
	// addVariable("vtx_isolation010", "VF", order);
	addVariable("vtx_isolation020", "VF", order);
	addVariable("vtx_pval", "VF", order);
	addVariable("trks_fitprob", "VF", order);
	addVariable("vtx_pt", "VF", order);
	addVariable("ht_pt",      "VF", order);
	addVariable("mets_dphi", "VF", order);
	
	addVariable("met_"+metType+"_et", "F", order);
	addVariable("met_"+metType+"_mT", "VF", order);
	addVariable("mets_dptrelcal", "VF", order);
	// addVariable("mets_dhtrelcal", "VF", order);
	addVariable("met_track_et", "F", order);
	addVariable("met_track_mT", "VF", order);
	addVariable("mets_dptreltrk", "VF", order);
	// addVariable("mets_dhtreltrk", "VF", order);
	
	addVariable("met_"+metType+"_dPhi3mu", "VF", order);
	// addVariable("met_track_dPhi3mu", "VF", order);
	// addVariable("ht_dr3body", "VF", order);
	addVariable("geo_lxySig", "VF", order);
	addVariable("geo_a0xySig", "VF", order);
	addVariable("vtx_pvNtrk", "VI", order);
	
	// addVariable("vtx_code", "VI", order);
	
	// addVariable("ht_mT",      "VF", order);
	// addVariable("ht_dphimet_"+metType, "VF", order);
	// addVariable("ht_mT_mettrk", "VF", order);
	// addVariable("ht_dphimet_track", "VF", order);
    

	// addVariable("vtx_mOS1", "VF", order);
	// addVariable("vtx_mOS2", "VF", order);
	// addVariable("vtx_mSS", "VF", order);
	
	//// MVA spectators
	order = -1;
	// addSpectator("vtx_mSS","VF",order);
}


void initTMVA(TString weightsfilepath, TString method)
{
	TMVA::Tools::Instance(); // This loads the library
	
	// Create the Reader object
	reader = new TMVA::Reader( "!Color:!Silent" );
	
	// Create a set of variables and declare them to the reader
	// the variable names MUST corresponds in name and type to those given in the weight file(s) used
	
	for(TMapiTS::iterator it= mva_var_order.begin() ; it!=mva_var_order.end() ; ++it)
	{
		TString name = it->second;
		reader->AddVariable(name,&mva_variables[name]);
	}

	// Spectator variables declared in the training have to be added to the reader, too
	for(TMapiTS::iterator it= mva_spe_order.begin() ; it!=mva_spe_order.end() ; ++it)
	{
		TString name = it->second;
		reader->AddSpectator(name,&mva_spectators[name]);
	}
	
	// Book MVA method
	methodName = method+" method";
	reader->BookMVA(methodName,weightsfilepath+"tmvaClass_BDTG.weights.xml");
}


void setMVAvars(unsigned int vtx,
				TMapTSi& input_ints, TMapTSf& input_floats, 
				TMapTSP2vi& input_vints, TMapTSP2vf& input_vfloats,
				TMapTSP2vi& input_vints_decorations, TMapTSP2vf& input_vfloats_decorations, float minjetptMeV)
{	
	for(TMapiTS::iterator it= mva_var_order.begin() ; it!=mva_var_order.end() ; ++it)
	{
		TString name = it->second;
		
		//// variables that can be easily calculated
		if(name.Contains("mu_pvaltrkfit := ")) mva_variables[name] = input_vfloats["mu_pvaltrkfit1"]->at(vtx)*input_vfloats["mu_pvaltrkfit2"]->at(vtx)*input_vfloats["mu_pvaltrkfit3"]->at(vtx);
		if(name.Contains("mu_dEdx :="))        mva_variables[name] = input_vfloats["mu_pixeldEdx1"]->at(vtx)+input_vfloats["mu_pixeldEdx2"]->at(vtx)+input_vfloats["mu_pixeldEdx3"]->at(vtx);
		if(name.Contains("vtx_m2body :="))     mva_variables[name] = (input_vfloats["vtx_mOS1"]->at(vtx)+input_vfloats["vtx_mOS2"]->at(vtx))/2.;
		if(name.Contains("jet_n :="))          mva_variables[name] = (input_vfloats["jet_pt1"]->at(vtx)>minjetptMeV)+(input_vfloats["jet_pt2"]->at(vtx)>minjetptMeV)+(input_vfloats["jet_pt3"]->at(vtx)>minjetptMeV)+(input_vfloats["jet_pt4"]->at(vtx)>minjetptMeV);
		if(name.Contains("jet_b :="))          mva_variables[name] = (input_vfloats["jet_MV1w1"]->at(vtx)>0.3511)+(input_vfloats["jet_MV1w2"]->at(vtx)>0.3511)+(input_vfloats["jet_MV1w3"]->at(vtx)>0.3511)+(input_vfloats["jet_MV1w4"]->at(vtx)>0.3511);
		if(name.Contains("vtx_dpfrac113"))     mva_variables[name] = (input_vfloats["mu_pt1"]->at(vtx)-input_vfloats["mu_pt3"]->at(vtx))/input_vfloats["mu_pt3"]->at(vtx);
		if(name.Contains("mu_trkqoverp1 :="))  mva_variables[name] = fabs(input_vfloats["mu_trkqoverp1"]->at(vtx));
		if(name.Contains("mu_trkqoverp2 :="))  mva_variables[name] = fabs(input_vfloats["mu_trkqoverp2"]->at(vtx));
		if(name.Contains("mu_trkqoverp3 :="))  mva_variables[name] = fabs(input_vfloats["mu_trkqoverp3"]->at(vtx));
		if(name.Contains("mu_pvalqoverp1 :=")) mva_variables[name] = fabs(input_vfloats["mu_pvaltrkfit1"]->at(vtx)/input_vfloats["mu_trkqoverp1"]->at(vtx));
		if(name.Contains("mu_pvalqoverp2 :=")) mva_variables[name] = fabs(input_vfloats["mu_pvaltrkfit2"]->at(vtx)/input_vfloats["mu_trkqoverp2"]->at(vtx));
		if(name.Contains("mu_pvalqoverp3 :=")) mva_variables[name] = fabs(input_vfloats["mu_pvaltrkfit3"]->at(vtx)/input_vfloats["mu_trkqoverp3"]->at(vtx));
		
		//// decorations
		if(mva_var_type[name]=="VF")
		{	
			if     (name.BeginsWith("ht_"))    mva_variables[name] = input_vfloats_decorations[name]->at(vtx); //// decorations
			else if(name.BeginsWith("trks_"))  mva_variables[name] = input_vfloats_decorations[name]->at(vtx); //// decorations
			else if(name.BeginsWith("geo_"))   mva_variables[name] = input_vfloats_decorations[name]->at(vtx); //// decorations
			else if(name.BeginsWith("mets_"))  mva_variables[name] = input_vfloats_decorations[name]->at(vtx); //// decorations
			else if(name.BeginsWith("muons_")) mva_variables[name] = input_vfloats_decorations[name]->at(vtx); //// decorations
			else                               mva_variables[name] = input_vfloats[name]->at(vtx);
		}
		if(mva_var_type[name]=="VI")
		{
			if(name.BeginsWith("jets_")) mva_variables[name] = input_vints_decorations[name]->at(vtx); //// decorations
			else                         mva_variables[name] = input_vints[name]->at(vtx);
		}
		
		//// others
		if(mva_var_type[name]=="F")  mva_variables[name] = input_floats[name];
		if(mva_var_type[name]=="I")  mva_variables[name] = input_ints[name];
	}
}


void setMVAspect(unsigned int vtx,
				TMapTSi& input_ints, TMapTSf& input_floats,
				TMapTSP2vi& input_vints, TMapTSP2vf& input_vfloats)
{	
	for(TMapiTS::iterator it= mva_spe_order.begin() ; it!=mva_spe_order.end() ; ++it)
	{
		TString name = it->second;
		if(mva_spe_type[name]=="VF") mva_spectators[name] = input_vfloats[name]->at(vtx);
		if(mva_spe_type[name]=="VI") mva_spectators[name] = input_vints[name]->at(vtx);
		if(mva_spe_type[name]=="F")  mva_spectators[name] = input_floats[name];
		if(mva_spe_type[name]=="I")  mva_spectators[name] = input_ints[name];
	}
}


Double_t getMVAscore()
{
	return reader->EvaluateMVA(methodName);
}
