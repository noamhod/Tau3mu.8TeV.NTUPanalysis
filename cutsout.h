#ifndef CUTSOUT_H
#define CUTSOUT_H

#include "std.h"
#include "type.h"

TDirectory* olddir2 = gDirectory;

TFile* fcutsout;
TTree* tcutsout;
TString fcutsoutname;

float cutsout_weight;
int   cutsout_runnumber;
int   cutsout_eventnumber;
int   cutsout_lbn;
vector<float>* cutsout_m3body;
vector<float>* cutsout_pt3body;
vector<float>* cutsout_mOS1;
vector<float>* cutsout_mOS2;
vector<float>* cutsout_mSS;
vector<float>* cutsout_pval;
vector<float>* cutsout_qop1;
vector<float>* cutsout_qop2;
vector<float>* cutsout_qop3;
vector<float>* cutsout_pval1;
vector<float>* cutsout_pval2;
vector<float>* cutsout_pval3;
vector<float>* cutsout_maxscat;
vector<float>* cutsout_minscat;
vector<float>* cutsout_maxnegb;
vector<float>* cutsout_minnegb;
vector<float>* cutsout_maxpbal;
vector<float>* cutsout_minpbal;
vector<int>*   cutsout_code;


void initCUTSout(TString fname)
{
	//olddir2->cd();
	fcutsoutname = fname;
	fcutsout = new TFile(fname,"RECREATE");
}


void setCUTSout(TString name)
{
	fcutsout->cd();
	tcutsout = new TTree("fltcuts_"+name,"fltcuts_"+name);
	
	if(cutsout_m3body)  delete cutsout_m3body;
	if(cutsout_pt3body) delete cutsout_pt3body;
	if(cutsout_mOS1)    delete cutsout_mOS1;
	if(cutsout_mOS2)    delete cutsout_mOS2;
	if(cutsout_mSS)     delete cutsout_mSS;
	if(cutsout_pval)    delete cutsout_pval;
	if(cutsout_qop1)    delete cutsout_qop1;
	if(cutsout_qop2)    delete cutsout_qop2;
	if(cutsout_qop3)    delete cutsout_qop3;
	if(cutsout_pval1)   delete cutsout_pval1;
	if(cutsout_pval2)   delete cutsout_pval2;
	if(cutsout_pval3)   delete cutsout_pval3;
	if(cutsout_maxscat) delete cutsout_maxscat;
	if(cutsout_minscat) delete cutsout_minscat;
	if(cutsout_maxnegb) delete cutsout_maxnegb;
	if(cutsout_minnegb) delete cutsout_minnegb;
	if(cutsout_maxpbal) delete cutsout_maxpbal;
	if(cutsout_minpbal) delete cutsout_minpbal;
	if(cutsout_code)    delete cutsout_code;
	
	cutsout_m3body          = new vector<float>;
	cutsout_pt3body         = new vector<float>;
	cutsout_mOS1            = new vector<float>;
	cutsout_mOS2            = new vector<float>;
	cutsout_mSS             = new vector<float>;
	cutsout_pval            = new vector<float>;
	cutsout_qop1            = new vector<float>;
	cutsout_qop2            = new vector<float>;
	cutsout_qop3            = new vector<float>;
	cutsout_pval1           = new vector<float>;
	cutsout_pval2           = new vector<float>;
	cutsout_pval3           = new vector<float>;
	cutsout_maxscat         = new vector<float>;
	cutsout_minscat         = new vector<float>;
	cutsout_maxnegb         = new vector<float>;
	cutsout_minnegb         = new vector<float>;
	cutsout_maxpbal         = new vector<float>;
	cutsout_minpbal         = new vector<float>;
	cutsout_code            = new vector<int>;
	
	cutsout_weight      = 0;
	cutsout_runnumber   = 0;
	cutsout_eventnumber = 0;
	cutsout_lbn         = 0;
	
	tcutsout->Branch("weight",          &cutsout_weight);
	tcutsout->Branch("RunNumber",       &cutsout_runnumber);
	tcutsout->Branch("EventNumber",     &cutsout_eventnumber);
	tcutsout->Branch("lbn",             &cutsout_lbn);
	tcutsout->Branch("m3body",          &cutsout_m3body);
	tcutsout->Branch("pt3body",         &cutsout_pt3body);
	tcutsout->Branch("mOS1",            &cutsout_mOS1);
	tcutsout->Branch("mOS2",            &cutsout_mOS2);
	tcutsout->Branch("mSS",             &cutsout_mSS);
	tcutsout->Branch("pval",            &cutsout_pval);
	tcutsout->Branch("qop1",            &cutsout_qop1);
	tcutsout->Branch("qop3",            &cutsout_qop2);
	tcutsout->Branch("qop3",            &cutsout_qop3);
	tcutsout->Branch("pval1",           &cutsout_pval1);
	tcutsout->Branch("pval2",           &cutsout_pval2);
	tcutsout->Branch("pval3",           &cutsout_pval3);
	tcutsout->Branch("maxscat",         &cutsout_maxscat);
	tcutsout->Branch("minscat",         &cutsout_minscat);
	tcutsout->Branch("maxnegb",         &cutsout_maxnegb);
	tcutsout->Branch("minnegb",         &cutsout_minnegb);
	tcutsout->Branch("maxpbal",         &cutsout_maxpbal);
	tcutsout->Branch("minpbal",         &cutsout_minpbal);
	tcutsout->Branch("code",            &cutsout_code);
}

void clearCUTSout()
{
	cutsout_weight      = 0;
	cutsout_runnumber   = 0;
	cutsout_eventnumber = 0;
	cutsout_lbn         = 0;
	
	cutsout_m3body->clear();
	cutsout_pt3body->clear();
	cutsout_mOS1->clear();
	cutsout_mOS2->clear();
	cutsout_mSS->clear();
	cutsout_pval->clear();
	cutsout_qop1->clear();
	cutsout_qop2->clear();
	cutsout_qop3->clear();
	cutsout_pval1->clear();
	cutsout_pval2->clear();
	cutsout_pval3->clear();
	cutsout_maxscat->clear();
	cutsout_minscat->clear();
	cutsout_maxnegb->clear();
	cutsout_minnegb->clear();
	cutsout_maxpbal->clear();
	cutsout_minpbal->clear();
	cutsout_code->clear();
}

void fillCUTSoutVecVars(unsigned int vtx, TMapTSP2vi& vi, TMapTSP2vf& vf, TMapTSP2vi& vid, TMapTSP2vf& vfd)
{
	if(0) cout << "vid.size()=" << vid.size() << endl;

	cutsout_m3body->push_back(vf["vtx_mass"]->at(vtx));
	cutsout_pt3body->push_back(vf["vtx_pt"]->at(vtx));
	cutsout_mOS1->push_back(vf["vtx_mOS1"]->at(vtx));
	cutsout_mOS2->push_back(vf["vtx_mOS2"]->at(vtx));
	cutsout_mSS->push_back(vf["vtx_mSS"]->at(vtx));
	cutsout_pval->push_back(vf["vtx_pval"]->at(vtx));
	cutsout_qop1->push_back(vf["mu_trkqoverp1"]->at(vtx));
	cutsout_qop2->push_back(vf["mu_trkqoverp2"]->at(vtx));
	cutsout_qop3->push_back(vf["mu_trkqoverp3"]->at(vtx));
	cutsout_pval1->push_back(vf["mu_pvaltrkfit1"]->at(vtx));
	cutsout_pval2->push_back(vf["mu_pvaltrkfit2"]->at(vtx));
	cutsout_pval3->push_back(vf["mu_pvaltrkfit3"]->at(vtx));
	cutsout_maxscat->push_back(vfd["muons_maxscatsig"]->at(vtx));
	cutsout_maxscat->push_back(vfd["muons_minscatsig"]->at(vtx));
	cutsout_maxnegb->push_back(vfd["muons_maxnegbsig"]->at(vtx));
	cutsout_maxnegb->push_back(vfd["muons_minnegbsig"]->at(vtx));
	cutsout_maxpbal->push_back(vfd["muons_maxpbalsig"]->at(vtx));
	cutsout_maxpbal->push_back(vfd["muons_minpbalsig"]->at(vtx));
	cutsout_code->push_back(vi["vtx_code"]->at(vtx));
}
void fillCUTSoutTree(TMapTSi& i, TMapTSf& f)
{
	cutsout_weight      = f["wgt_total"];
	cutsout_runnumber   = i["evt_RunNumber"];
	cutsout_eventnumber = i["evt_EventNumber"];
	cutsout_lbn         = i["evt_lbn"];
	
	fcutsout->cd();
	tcutsout->Fill();
}

void finalizeCUTSout()
{
	fcutsout->Write();
	fcutsout->Close();
	delete fcutsout;
	olddir2->cd();
	_INF(1,"Written: "<<fcutsoutname);
}

#endif
