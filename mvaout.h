#ifndef MVAOUT_H
#define MVAOUT_H

#include "std.h"
#include "type.h"

TDirectory* olddir1 = gDirectory;

TFile* fmvaout;
TTree* tmvaout;
TString fmvaoutname;

float mvaout_weight;
int   mvaout_runnumber;
int   mvaout_eventnumber;
int   mvaout_lbn;
vector<float>* mvaout_score;
vector<float>* mvaout_m3body;
vector<float>* mvaout_pt3body;
vector<float>* mvaout_mOS1;
vector<float>* mvaout_mOS2;
vector<float>* mvaout_mSS;
vector<float>* mvaout_pval;
vector<float>* mvaout_qop1;
vector<float>* mvaout_qop2;
vector<float>* mvaout_qop3;
vector<float>* mvaout_pval1;
vector<float>* mvaout_pval2;
vector<float>* mvaout_pval3;
vector<float>* mvaout_maxscat;
vector<float>* mvaout_minscat;
vector<float>* mvaout_maxnegb;
vector<float>* mvaout_minnegb;
vector<float>* mvaout_maxpbal;
vector<float>* mvaout_minpbal;
vector<int>*   mvaout_code;
vector<int>*   mvaout_pass_tight_scat;
vector<int>*   mvaout_pass_tight_pbal;
vector<int>*   mvaout_pass_tight_fitq;

vector<float>* mvaout_mettrk    ;
vector<float>* mvaout_metcal    ;
vector<float>* mvaout_trkspval  ;
vector<float>* mvaout_maxpbalsig;
vector<float>* mvaout_iso003    ;
vector<float>* mvaout_iso004    ;
vector<float>* mvaout_iso005    ;
vector<float>* mvaout_iso006    ;
vector<float>* mvaout_iso007    ;
vector<float>* mvaout_iso008    ;
vector<float>* mvaout_iso009    ;
vector<float>* mvaout_iso010    ;
vector<float>* mvaout_iso012    ;
vector<float>* mvaout_iso014    ;
vector<float>* mvaout_iso016    ;
vector<float>* mvaout_iso018    ;
vector<float>* mvaout_iso020    ;
vector<float>* mvaout_iso022    ;
vector<float>* mvaout_iso024    ;
vector<float>* mvaout_iso026    ;
vector<float>* mvaout_iso028    ;
vector<float>* mvaout_iso030    ;
vector<float>* mvaout_metsdphi  ;
vector<float>* mvaout_metsdptavg;
vector<float>* mvaout_metsdpttrk;
vector<float>* mvaout_metsdptcal;
vector<float>* mvaout_metsdhttrk;
vector<float>* mvaout_metsdhtcal;
vector<float>* mvaout_mtcal     ;
vector<float>* mvaout_dphical   ;
vector<float>* mvaout_mttrk     ;
vector<float>* mvaout_dphitrk   ;
vector<float>* mvaout_lxy       ;
vector<float>* mvaout_dlxy      ;
vector<float>* mvaout_lxysig    ;
vector<int>*   mvaout_pvntrk    ;
vector<float>* mvaout_a0xy      ;
vector<float>* mvaout_da0xy     ;
vector<float>* mvaout_a0xysig   ;
vector<int>*   mvaout_njets     ;
vector<int>*   mvaout_bjets     ;
vector<float>* mvaout_ht        ;
vector<float>* mvaout_mhtcal    ;
vector<float>* mvaout_dphihtcal ;
vector<float>* mvaout_mhttrk    ;
vector<float>* mvaout_dphihttrk ;
vector<float>* mvaout_drht      ;

void initMVAout(TString fname)
{
	//olddir1->cd();
	fmvaoutname = fname;
	fmvaout = new TFile(fname,"RECREATE");
}


void setMVAout(TString name)
{
	fmvaout->cd();
	tmvaout = new TTree("fltmva_"+name,"fltmva_"+name);
	
	if(mvaout_score)   delete mvaout_score;
	if(mvaout_m3body)  delete mvaout_m3body;
	if(mvaout_pt3body) delete mvaout_pt3body;
	if(mvaout_mOS1)    delete mvaout_mOS1;
	if(mvaout_mOS2)    delete mvaout_mOS2;
	if(mvaout_mSS)     delete mvaout_mSS;
	if(mvaout_pval)    delete mvaout_pval;
	if(mvaout_qop1)    delete mvaout_qop1;
	if(mvaout_qop2)    delete mvaout_qop2;
	if(mvaout_qop3)    delete mvaout_qop3;
	if(mvaout_pval1)   delete mvaout_pval1;
	if(mvaout_pval2)   delete mvaout_pval2;
	if(mvaout_pval3)   delete mvaout_pval3;
	if(mvaout_maxscat) delete mvaout_maxscat;
	if(mvaout_minscat) delete mvaout_minscat;
	if(mvaout_maxnegb) delete mvaout_maxnegb;
	if(mvaout_minnegb) delete mvaout_minnegb;
	if(mvaout_maxpbal) delete mvaout_maxpbal;
	if(mvaout_minpbal) delete mvaout_minpbal;
	if(mvaout_code)    delete mvaout_code;
	if(mvaout_pass_tight_scat)    delete mvaout_pass_tight_scat;
	if(mvaout_pass_tight_pbal)    delete mvaout_pass_tight_pbal;
	if(mvaout_pass_tight_fitq)    delete mvaout_pass_tight_fitq;
	
	if(mvaout_mettrk    ) delete mvaout_mettrk    ;
	if(mvaout_metcal    ) delete mvaout_metcal    ;
	if(mvaout_trkspval  ) delete mvaout_trkspval  ;
	if(mvaout_maxpbalsig) delete mvaout_maxpbalsig;
	if(mvaout_iso003    ) delete mvaout_iso003    ;
	if(mvaout_iso004    ) delete mvaout_iso004    ;
	if(mvaout_iso005    ) delete mvaout_iso005    ;
	if(mvaout_iso006    ) delete mvaout_iso006    ;
	if(mvaout_iso007    ) delete mvaout_iso007    ;
	if(mvaout_iso008    ) delete mvaout_iso008    ;
	if(mvaout_iso009    ) delete mvaout_iso009    ;
	if(mvaout_iso010    ) delete mvaout_iso010    ;
	if(mvaout_iso012    ) delete mvaout_iso012    ;
	if(mvaout_iso014    ) delete mvaout_iso014    ;
	if(mvaout_iso016    ) delete mvaout_iso016    ;
	if(mvaout_iso018    ) delete mvaout_iso018    ;
	if(mvaout_iso020    ) delete mvaout_iso020    ;
	if(mvaout_iso022    ) delete mvaout_iso022    ;
	if(mvaout_iso024    ) delete mvaout_iso024    ;
	if(mvaout_iso026    ) delete mvaout_iso026    ;
	if(mvaout_iso028    ) delete mvaout_iso028    ;
	if(mvaout_iso030    ) delete mvaout_iso030    ;
	if(mvaout_metsdphi  ) delete mvaout_metsdphi  ;
	if(mvaout_metsdptavg) delete mvaout_metsdptavg;
	if(mvaout_metsdptcal) delete mvaout_metsdptcal;
	if(mvaout_metsdpttrk) delete mvaout_metsdhttrk;
	if(mvaout_metsdhtcal) delete mvaout_metsdhtcal;
	if(mvaout_metsdhttrk) delete mvaout_metsdpttrk;
	if(mvaout_mtcal     ) delete mvaout_mtcal     ;
	if(mvaout_dphical   ) delete mvaout_dphical   ;
	if(mvaout_mttrk     ) delete mvaout_mttrk     ;
	if(mvaout_dphitrk   ) delete mvaout_dphitrk   ;
	if(mvaout_lxy       ) delete mvaout_lxy       ;
	if(mvaout_dlxy      ) delete mvaout_dlxy      ;
	if(mvaout_lxysig    ) delete mvaout_lxysig    ;
	if(mvaout_pvntrk    ) delete mvaout_pvntrk    ;
	if(mvaout_a0xy      ) delete mvaout_a0xy      ;
	if(mvaout_da0xy     ) delete mvaout_da0xy     ;
	if(mvaout_a0xysig   ) delete mvaout_a0xysig   ;
	if(mvaout_njets     ) delete mvaout_njets     ;
	if(mvaout_bjets     ) delete mvaout_bjets     ;
	if(mvaout_ht        ) delete mvaout_ht        ;
	if(mvaout_mhtcal    ) delete mvaout_mhtcal    ;
	if(mvaout_dphihtcal ) delete mvaout_dphihtcal ;
	if(mvaout_mhttrk    ) delete mvaout_mhttrk    ;
	if(mvaout_dphihttrk ) delete mvaout_dphihttrk ;
	if(mvaout_drht      ) delete mvaout_drht      ;
	
	
	mvaout_score           = new vector<float>;
	mvaout_m3body          = new vector<float>;
	mvaout_pt3body         = new vector<float>;
	mvaout_mOS1            = new vector<float>;
	mvaout_mOS2            = new vector<float>;
	mvaout_mSS             = new vector<float>;
	mvaout_pval            = new vector<float>;
	mvaout_qop1            = new vector<float>;
	mvaout_qop2            = new vector<float>;
	mvaout_qop3            = new vector<float>;
	mvaout_pval1           = new vector<float>;
	mvaout_pval2           = new vector<float>;
	mvaout_pval3           = new vector<float>;
	mvaout_maxscat         = new vector<float>;
	mvaout_minscat         = new vector<float>;
	mvaout_maxnegb         = new vector<float>;
	mvaout_minnegb         = new vector<float>;
	mvaout_maxpbal         = new vector<float>;
	mvaout_minpbal         = new vector<float>;
	mvaout_code            = new vector<int>;
	mvaout_pass_tight_scat = new vector<int>;
	mvaout_pass_tight_pbal = new vector<int>;
	mvaout_pass_tight_fitq = new vector<int>;
	
	mvaout_mettrk    = new vector<float>;
	mvaout_metcal    = new vector<float>;
	mvaout_trkspval  = new vector<float>;
	mvaout_maxpbalsig= new vector<float>;
	mvaout_iso003    = new vector<float>;
	mvaout_iso004    = new vector<float>;
	mvaout_iso005    = new vector<float>;
	mvaout_iso006    = new vector<float>;
	mvaout_iso007    = new vector<float>;
	mvaout_iso008    = new vector<float>;
	mvaout_iso009    = new vector<float>;
	mvaout_iso010    = new vector<float>;
	mvaout_iso012    = new vector<float>;
	mvaout_iso014    = new vector<float>;
	mvaout_iso016    = new vector<float>;
	mvaout_iso018    = new vector<float>;
	mvaout_iso020    = new vector<float>;
	mvaout_iso022    = new vector<float>;
	mvaout_iso024    = new vector<float>;
	mvaout_iso026    = new vector<float>;
	mvaout_iso028    = new vector<float>;
	mvaout_iso030    = new vector<float>;
	mvaout_metsdphi  = new vector<float>;
	mvaout_metsdptavg= new vector<float>;
	mvaout_metsdptcal= new vector<float>;
	mvaout_metsdpttrk= new vector<float>;
	mvaout_metsdhtcal= new vector<float>;
	mvaout_metsdhttrk= new vector<float>;
	mvaout_mtcal     = new vector<float>;
	mvaout_dphical   = new vector<float>;
	mvaout_mttrk     = new vector<float>;
	mvaout_dphitrk   = new vector<float>;
	mvaout_lxy       = new vector<float>;
	mvaout_dlxy      = new vector<float>;
	mvaout_lxysig    = new vector<float>;
	mvaout_pvntrk    = new vector<int>;
	mvaout_a0xy      = new vector<float>;
	mvaout_da0xy     = new vector<float>;
	mvaout_a0xysig   = new vector<float>;
	mvaout_njets     = new vector<int>;
	mvaout_bjets     = new vector<int>;
	mvaout_ht        = new vector<float>;
	mvaout_mhtcal    = new vector<float>;
	mvaout_dphihtcal = new vector<float>;
	mvaout_mhttrk    = new vector<float>;
	mvaout_dphihttrk = new vector<float>;
	mvaout_drht      = new vector<float>;
	
	mvaout_weight      = 0;
	mvaout_runnumber   = 0;
	mvaout_eventnumber = 0;
	mvaout_lbn         = 0;
	
	tmvaout->Branch("weight",          &mvaout_weight);
	tmvaout->Branch("RunNumber",       &mvaout_runnumber);
	tmvaout->Branch("EventNumber",     &mvaout_eventnumber);
	tmvaout->Branch("lbn",             &mvaout_lbn);
	tmvaout->Branch("score",           &mvaout_score);
	tmvaout->Branch("m3body",          &mvaout_m3body);
	tmvaout->Branch("pt3body",         &mvaout_pt3body);
	tmvaout->Branch("mOS1",            &mvaout_mOS1);
	tmvaout->Branch("mOS2",            &mvaout_mOS2);
	tmvaout->Branch("mSS",             &mvaout_mSS);
	tmvaout->Branch("pval",            &mvaout_pval);
	tmvaout->Branch("qop1",            &mvaout_qop1);
	tmvaout->Branch("qop3",            &mvaout_qop2);
	tmvaout->Branch("qop3",            &mvaout_qop3);
	tmvaout->Branch("pval1",           &mvaout_pval1);
	tmvaout->Branch("pval2",           &mvaout_pval2);
	tmvaout->Branch("pval3",           &mvaout_pval3);
	tmvaout->Branch("maxscat",         &mvaout_maxscat);
	tmvaout->Branch("minscat",         &mvaout_minscat);
	tmvaout->Branch("maxnegb",         &mvaout_maxnegb);
	tmvaout->Branch("minnegb",         &mvaout_minnegb);
	tmvaout->Branch("maxpbal",         &mvaout_maxpbal);
	tmvaout->Branch("minpbal",         &mvaout_minpbal);
	tmvaout->Branch("code",            &mvaout_code);
	tmvaout->Branch("pass_tight_scat", &mvaout_pass_tight_scat);
	tmvaout->Branch("pass_tight_pbal", &mvaout_pass_tight_pbal);
	tmvaout->Branch("pass_tight_fitq", &mvaout_pass_tight_fitq);
	
	tmvaout->Branch("mettrk",     &mvaout_mettrk    );
	tmvaout->Branch("metcal",     &mvaout_metcal    );
	tmvaout->Branch("trkspval",   &mvaout_trkspval  );
	tmvaout->Branch("maxpbalsig", &mvaout_maxpbalsig);
	tmvaout->Branch("iso003",     &mvaout_iso003    );
	tmvaout->Branch("iso004",     &mvaout_iso004    );
	tmvaout->Branch("iso005",     &mvaout_iso005    );
	tmvaout->Branch("iso006",     &mvaout_iso006    );
	tmvaout->Branch("iso007",     &mvaout_iso007    );
	tmvaout->Branch("iso008",     &mvaout_iso008    );
	tmvaout->Branch("iso009",     &mvaout_iso009    );
	tmvaout->Branch("iso010",     &mvaout_iso010    );
	tmvaout->Branch("iso012",     &mvaout_iso012    );
	tmvaout->Branch("iso014",     &mvaout_iso014    );
	tmvaout->Branch("iso016",     &mvaout_iso016    );
	tmvaout->Branch("iso018",     &mvaout_iso018    );
	tmvaout->Branch("iso020",     &mvaout_iso020    );
	tmvaout->Branch("iso022",     &mvaout_iso022    );
	tmvaout->Branch("iso024",     &mvaout_iso024    );
	tmvaout->Branch("iso026",     &mvaout_iso026    );
	tmvaout->Branch("iso028",     &mvaout_iso028    );
	tmvaout->Branch("iso030",     &mvaout_iso030    );
	tmvaout->Branch("metsdphi",   &mvaout_metsdphi  );
	tmvaout->Branch("mets_dptrelavg",&mvaout_metsdptavg);
	tmvaout->Branch("mets_dptreltrk",&mvaout_metsdptcal);
	tmvaout->Branch("mets_dptrelcal",&mvaout_metsdpttrk);
	tmvaout->Branch("mets_dhtreltrk",&mvaout_metsdhtcal);
	tmvaout->Branch("mets_dhtrelcal",&mvaout_metsdhttrk);
	tmvaout->Branch("mtcal",      &mvaout_mtcal     );
	tmvaout->Branch("dphical",    &mvaout_dphical   );
	tmvaout->Branch("mttrk",      &mvaout_mttrk     );
	tmvaout->Branch("dphitrk",    &mvaout_dphitrk   );
	tmvaout->Branch("lxy",        &mvaout_lxy       );
	tmvaout->Branch("dlxy",       &mvaout_dlxy      );
	tmvaout->Branch("lxysig",     &mvaout_lxysig    );
	tmvaout->Branch("pvntrk",     &mvaout_pvntrk    );
	tmvaout->Branch("a0xy",       &mvaout_a0xy      );
	tmvaout->Branch("da0xy",      &mvaout_da0xy     );
	tmvaout->Branch("a0xysig",    &mvaout_a0xysig   );
	tmvaout->Branch("njets",      &mvaout_njets     );
	tmvaout->Branch("bjets",      &mvaout_bjets     );
	tmvaout->Branch("ht",         &mvaout_ht        );
	tmvaout->Branch("mhtcal",     &mvaout_mhtcal    );
	tmvaout->Branch("dphihtcal",  &mvaout_dphihtcal );
	tmvaout->Branch("mhttrk",     &mvaout_mhttrk    );
	tmvaout->Branch("dphihttrk",  &mvaout_dphihttrk );
	tmvaout->Branch("drht",       &mvaout_drht      );
}

void clearMVAout()
{
	mvaout_weight      = 0;
	mvaout_runnumber   = 0;
	mvaout_eventnumber = 0;
	mvaout_lbn         = 0;
	
	mvaout_score->clear();
	mvaout_m3body->clear();
	mvaout_pt3body->clear();
	mvaout_mOS1->clear();
	mvaout_mOS2->clear();
	mvaout_mSS->clear();
	mvaout_pval->clear();
	mvaout_qop1->clear();
	mvaout_qop2->clear();
	mvaout_qop3->clear();
	mvaout_pval1->clear();
	mvaout_pval2->clear();
	mvaout_pval3->clear();
	mvaout_maxscat->clear();
	mvaout_minscat->clear();
	mvaout_maxnegb->clear();
	mvaout_minnegb->clear();
	mvaout_maxpbal->clear();
	mvaout_minpbal->clear();
	mvaout_code->clear();
	mvaout_pass_tight_scat->clear();
	mvaout_pass_tight_pbal->clear();
	mvaout_pass_tight_fitq->clear();
	
	mvaout_mettrk    ->clear();
	mvaout_metcal    ->clear();
	mvaout_trkspval  ->clear();
	mvaout_maxpbalsig->clear();
	mvaout_iso003    ->clear();
	mvaout_iso004    ->clear();
	mvaout_iso005    ->clear();
	mvaout_iso006    ->clear();
	mvaout_iso007    ->clear();
	mvaout_iso008    ->clear();
	mvaout_iso009    ->clear();
	mvaout_iso010    ->clear();
	mvaout_iso012    ->clear();
	mvaout_iso014    ->clear();
	mvaout_iso016    ->clear();
	mvaout_iso018    ->clear();
	mvaout_iso020    ->clear();
	mvaout_iso022    ->clear();
	mvaout_iso024    ->clear();
	mvaout_iso026    ->clear();
	mvaout_iso028    ->clear();
	mvaout_iso030    ->clear();
	mvaout_metsdphi  ->clear();
	mvaout_metsdptavg->clear();
	mvaout_metsdptcal->clear();
	mvaout_metsdpttrk->clear();
	mvaout_metsdhtcal->clear();
	mvaout_metsdhttrk->clear();
	mvaout_mtcal     ->clear();
	mvaout_dphical   ->clear();
	mvaout_mttrk     ->clear();
	mvaout_dphitrk   ->clear();
	mvaout_lxy       ->clear();
	mvaout_dlxy      ->clear();
	mvaout_lxysig    ->clear();
	mvaout_pvntrk    ->clear();
	mvaout_a0xy      ->clear();
	mvaout_da0xy     ->clear();
	mvaout_a0xysig   ->clear();
	mvaout_njets     ->clear();
	mvaout_bjets     ->clear();
	mvaout_ht        ->clear();
	mvaout_mhtcal    ->clear();
	mvaout_dphihtcal ->clear();
	mvaout_mhttrk    ->clear();
	mvaout_dphihttrk ->clear();
	mvaout_drht      ->clear();	
}

void fillMVAoutVecVars(unsigned int vtx, float score, TMapTSf& f, TMapTSP2vi& vi, TMapTSP2vf& vf, TMapTSP2vi& vid, TMapTSP2vf& vfd, TString metType)
{
	mvaout_score->push_back(score);
	mvaout_m3body->push_back(vf["vtx_mass"]->at(vtx));
	mvaout_pt3body->push_back(vf["vtx_pt"]->at(vtx));
	mvaout_mOS1->push_back(vf["vtx_mOS1"]->at(vtx));
	mvaout_mOS2->push_back(vf["vtx_mOS2"]->at(vtx));
	mvaout_mSS->push_back(vf["vtx_mSS"]->at(vtx));
	mvaout_pval->push_back(vf["vtx_pval"]->at(vtx));
	mvaout_qop1->push_back(vf["mu_trkqoverp1"]->at(vtx));
	mvaout_qop2->push_back(vf["mu_trkqoverp2"]->at(vtx));
	mvaout_qop3->push_back(vf["mu_trkqoverp3"]->at(vtx));
	mvaout_pval1->push_back(vf["mu_pvaltrkfit1"]->at(vtx));
	mvaout_pval2->push_back(vf["mu_pvaltrkfit2"]->at(vtx));
	mvaout_pval3->push_back(vf["mu_pvaltrkfit3"]->at(vtx));
	mvaout_maxscat->push_back(vfd["muons_maxscatsig"]->at(vtx));
	mvaout_maxscat->push_back(vfd["muons_minscatsig"]->at(vtx));
	mvaout_maxnegb->push_back(vfd["muons_maxnegbsig"]->at(vtx));
	mvaout_maxnegb->push_back(vfd["muons_minnegbsig"]->at(vtx));
	mvaout_maxpbal->push_back(vfd["muons_maxpbalsig"]->at(vtx));
	mvaout_maxpbal->push_back(vfd["muons_minpbalsig"]->at(vtx));
	mvaout_code->push_back(vi["vtx_code"]->at(vtx));
	mvaout_pass_tight_scat->push_back(vid["pass_tight_scat"]->at(vtx));
	mvaout_pass_tight_pbal->push_back(vid["pass_tight_pbal"]->at(vtx));
	mvaout_pass_tight_fitq->push_back(vid["pass_tight_fitq"]->at(vtx));

	mvaout_mettrk    ->push_back(f["met_track_et"]);
	mvaout_metcal    ->push_back(f["met_"+metType+"_et"]);
	mvaout_metsdphi  ->push_back(vfd["mets_dphi"]->at(vtx));

	mvaout_metsdptavg->push_back(vfd["mets_dptrelavg"]->at(vtx));
	mvaout_metsdptcal->push_back(vfd["mets_dptreltrk"]->at(vtx));
	mvaout_metsdpttrk->push_back(vfd["mets_dptrelcal"]->at(vtx));
	mvaout_metsdhtcal->push_back(vfd["mets_dhtreltrk"]->at(vtx));
	mvaout_metsdhttrk->push_back(vfd["mets_dhtrelcal"]->at(vtx));
	
	mvaout_trkspval  ->push_back(vfd["trks_fitprob"]->at(vtx));
	mvaout_maxpbalsig->push_back(vfd["muons_maxpbalsig"]->at(vtx));
	
	mvaout_iso003    ->push_back(vf["vtx_isolation003"]->at(vtx));
	mvaout_iso004    ->push_back(vf["vtx_isolation004"]->at(vtx));
	mvaout_iso005    ->push_back(vf["vtx_isolation005"]->at(vtx));
	mvaout_iso006    ->push_back(vf["vtx_isolation006"]->at(vtx));
	mvaout_iso007    ->push_back(vf["vtx_isolation007"]->at(vtx));
	mvaout_iso008    ->push_back(vf["vtx_isolation008"]->at(vtx));
	mvaout_iso009    ->push_back(vf["vtx_isolation009"]->at(vtx));
	mvaout_iso010    ->push_back(vf["vtx_isolation010"]->at(vtx));
	mvaout_iso012    ->push_back(vf["vtx_isolation012"]->at(vtx));
	mvaout_iso014    ->push_back(vf["vtx_isolation014"]->at(vtx));
	mvaout_iso016    ->push_back(vf["vtx_isolation016"]->at(vtx));
	mvaout_iso018    ->push_back(vf["vtx_isolation018"]->at(vtx));
	mvaout_iso020    ->push_back(vf["vtx_isolation020"]->at(vtx));
	mvaout_iso022    ->push_back(vf["vtx_isolation022"]->at(vtx));
	mvaout_iso024    ->push_back(vf["vtx_isolation024"]->at(vtx));
	mvaout_iso026    ->push_back(vf["vtx_isolation026"]->at(vtx));
	mvaout_iso028    ->push_back(vf["vtx_isolation028"]->at(vtx));
	mvaout_iso030    ->push_back(vf["vtx_isolation030"]->at(vtx));
	
	mvaout_mtcal     ->push_back(vf["met_"+metType+"_mT"]->at(vtx));
	mvaout_dphical   ->push_back(vf["met_"+metType+"_dPhi3mu"]->at(vtx));
	
	mvaout_mttrk     ->push_back(vf["met_track_mT"]->at(vtx));
	mvaout_dphitrk   ->push_back(vf["met_track_dPhi3mu"]->at(vtx));
	
	mvaout_lxy       ->push_back(vf["vtx_lxy"]->at(vtx));
	mvaout_dlxy      ->push_back(vf["vtx_lxyErr"]->at(vtx));
	mvaout_lxysig    ->push_back(vfd["geo_lxySig"]->at(vtx));
	mvaout_pvntrk    ->push_back(vi["vtx_pvNtrk"]->at(vtx));
	mvaout_a0xy      ->push_back(vf["vtx_a0xy"]->at(vtx));
	mvaout_da0xy     ->push_back(vf["vtx_a0xyErr"]->at(vtx));
	mvaout_a0xysig   ->push_back(vfd["geo_a0xySig"]->at(vtx));
	
	mvaout_njets     ->push_back(vid["jets_n"]->at(vtx));
	mvaout_bjets     ->push_back(vid["jets_b"]->at(vtx));
	
	mvaout_ht        ->push_back(vfd["ht_pt"]->at(vtx));
	mvaout_mhtcal    ->push_back(vfd["ht_mT"]->at(vtx));
	mvaout_dphihtcal ->push_back(vfd["ht_dphimet_"+metType]->at(vtx));
	mvaout_mhttrk    ->push_back(vfd["ht_mT_mettrk"]->at(vtx));
	mvaout_dphihttrk ->push_back(vfd["ht_dphimet_track"]->at(vtx));
	mvaout_drht      ->push_back(vfd["ht_dr3body"]->at(vtx));
}
void fillMVAoutTree(TMapTSi& i, TMapTSf& f)
{
	mvaout_weight      = f["wgt_total"];
	mvaout_runnumber   = i["evt_RunNumber"];
	mvaout_eventnumber = i["evt_EventNumber"];
	mvaout_lbn         = i["evt_lbn"];
	
	fmvaout->cd();
	tmvaout->Fill();
}

void finalizeMVAout()
{
	fmvaout->Write();
	fmvaout->Close();
	delete fmvaout;
	olddir1->cd();
	_INF(1,"Written: "<<fmvaoutname);
}

#endif