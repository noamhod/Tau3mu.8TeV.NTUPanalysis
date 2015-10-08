{
	TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.loose.root","READ");
	TTree* tS = (TTree*)f->Get("flatout_Wtaunu_3mu");
	TTree* tD = (TTree*)f->Get("flatout_Data");
	
	TH1* hsMETtrk           = new TH1F("hsMETtrk",";MET;Events",30,0,100000);
	TH1* hsMETsoftNomUncert = new TH1F("hsMETsoftNomUncert",";MET;Uncertainty [\%]",30,0,100000);
	TH1* hsMETsoftNom       = new TH1F("hsMETsoftNom",";MET;Events",30,0,100000);
	TH1* hsMETsoftNom_up    = new TH1F("hsMETsoftNom_up",";MET;Events",30,0,100000);
	TH1* hsMETsoftNom_dn    = new TH1F("hsMETsoftNom_dn",";MET;Events",30,0,100000);
	
	TH1* hdMETtrk           = new TH1F("hdMETtrk",";MET;Events",30,0,100000);
	TH1* hdMETsoftNomUncert = new TH1F("hdMETsoftNomUncert",";MET;Uncertainty [\%]",30,0,100000);
	TH1* hdMETsoftNom       = new TH1F("hdMETsoftNom",";MET;Events",30,0,100000);
	TH1* hdMETsoftNom_up    = new TH1F("hdMETsoftNom_up",";MET;Events",30,0,100000);
	TH1* hdMETsoftNom_dn    = new TH1F("hdMETsoftNom_dn",";MET;Events",30,0,100000);
	
	
	TH1* hsmTtrk           = new TH1F("hsmTtrk",";mT;Events",50,0,150000);
	TH1* hsmTsoftNomUncert = new TH1F("hsmTsoftNomUncert",";mT;Uncertainty [\%]",50,0,150000);
	TH1* hsmTsoftNom       = new TH1F("hsmTsoftNom",";mT;Events",50,0,150000);
	TH1* hsmTsoftNom_up    = new TH1F("hsmTsoftNom_up",";mT;Events",50,0,150000);
	TH1* hsmTsoftNom_dn    = new TH1F("hsmTsoftNom_dn",";mT;Events",50,0,150000);
	
	TH1* hdmTtrk           = new TH1F("hdmTtrk",";mT;Events",50,0,150000);
	TH1* hdmTsoftNomUncert = new TH1F("hdmTsoftNomUncert",";mT;Uncertainty [\%]",50,0,150000);
	TH1* hdmTsoftNom       = new TH1F("hdmTsoftNom",";mT;Events",50,0,150000);
	TH1* hdmTsoftNom_up    = new TH1F("hdmTsoftNom_up",";mT;Events",50,0,150000);
	TH1* hdmTsoftNom_dn    = new TH1F("hdmTsoftNom_dn",";mT;Events",50,0,150000);


	TH1* hsdPhitrk           = new TH1F("hsdPhitrk",";dPhi;Events",20,0,TMath::Pi());
	TH1* hsdPhisoftNomUncert = new TH1F("hsdPhisoftNomUncert",";dPhi;Uncertainty [\%]",20,0,TMath::Pi());
	TH1* hsdPhisoftNom       = new TH1F("hsdPhisoftNom",";dPhi;Events",20,0,TMath::Pi());
	TH1* hsdPhisoftNom_up    = new TH1F("hsdPhisoftNom_up",";dPhi;Events",20,0,TMath::Pi());
	TH1* hsdPhisoftNom_dn    = new TH1F("hsdPhisoftNom_dn",";dPhi;Events",20,0,TMath::Pi());
	
	TH1* hddPhitrk           = new TH1F("hddPhitrk",";dPhi;Events",20,0,TMath::Pi());
	TH1* hddPhisoftNomUncert = new TH1F("hddPhisoftNomUncert",";dPhi;Uncertainty [\%]",20,0,TMath::Pi());
	TH1* hddPhisoftNom       = new TH1F("hddPhisoftNom",";dPhi;Events",20,0,TMath::Pi());
	TH1* hddPhisoftNom_up    = new TH1F("hddPhisoftNom_up",";dPhi;Events",20,0,TMath::Pi());
	TH1* hddPhisoftNom_dn    = new TH1F("hddPhisoftNom_dn",";dPhi;Events",20,0,TMath::Pi());
	
	
	TCut    isLoose    = "pass_loose==1";
	TCut    oneCand    = "@vtx_mass->size()==1";
	TCut    masses2    = "vtx_mOS1>220 && vtx_mOS2>220 && vtx_mSS>220";
	TCut    distances  = "geo_lxySig>-10 && geo_lxySig<50 && geo_a0xySig<25";
	TCut    pvalues    = "trks_fitprob>1e-9";
	TCut    kinematics = "vtx_pt>10000 && (met_muons_et>10000 && met_muons_et<250000) && (met_track_et>10000 && met_track_et<250000)  &&  met_muons_mT>20000 && met_track_mT>20000";
	TCut    iso        = "vtx_isolation020<0.3 && vtx_isolation030<1";
	TCut preTraining   = isLoose+oneCand+masses2+iso+distances+pvalues+kinematics;
	

	tS->Draw("met_track_et>>hsMETtrk",preTraining);
	tS->Draw("met_track_et_softtrk_nom>>hsMETsoftNom",preTraining);
	tS->Draw("met_track_et_softtrk_up>>hsMETsoftNom_up",preTraining);
	tS->Draw("met_track_et_softtrk_dwn>>hsMETsoftNom_dn",preTraining);
	
	tD->Draw("met_track_et>>hdMETtrk",preTraining);
	tD->Draw("met_track_et_softtrk_nom>>hdMETsoftNom",preTraining);
	tD->Draw("met_track_et_softtrk_up>>hdMETsoftNom_up",preTraining);
	tD->Draw("met_track_et_softtrk_dwn>>hdMETsoftNom_dn",preTraining);
	
	
	tS->Draw("met_track_mT>>hsmTtrk",preTraining);
	tS->Draw("met_track_mT_softtrk_nom>>hsmTsoftNom",preTraining);
	tS->Draw("met_track_mT_softtrk_up>>hsmTsoftNom_up",preTraining);
	tS->Draw("met_track_mT_softtrk_dwn>>hsmTsoftNom_dn",preTraining);
	
	tD->Draw("met_track_mT>>hdmTtrk",preTraining);
	tD->Draw("met_track_mT_softtrk_nom>>hdmTsoftNom",preTraining);
	tD->Draw("met_track_mT_softtrk_up>>hdmTsoftNom_up",preTraining);
	tD->Draw("met_track_mT_softtrk_dwn>>hdmTsoftNom_dn",preTraining);


	tS->Draw("met_track_dPhi3mu>>hsdPhitrk",preTraining);
	tS->Draw("met_track_dPhi3mu_softtrk_nom>>hsdPhisoftNom",preTraining);
	tS->Draw("met_track_dPhi3mu_softtrk_up>>hsdPhisoftNom_up",preTraining);
	tS->Draw("met_track_dPhi3mu_softtrk_dwn>>hsdPhisoftNom_dn",preTraining);
	
	tD->Draw("met_track_dPhi3mu>>hddPhitrk",preTraining);
	tD->Draw("met_track_dPhi3mu_softtrk_nom>>hddPhisoftNom",preTraining);
	tD->Draw("met_track_dPhi3mu_softtrk_up>>hddPhisoftNom_up",preTraining);
	tD->Draw("met_track_dPhi3mu_softtrk_dwn>>hddPhisoftNom_dn",preTraining);


	for(Int_t b=1 ; b<hsMETtrk->GetNbinsX() ; ++b)
	{	
		float METsoftNomUp = TMath::Abs(hsMETsoftNom->GetBinContent(b) - hsMETsoftNom_up->GetBinContent(b));
		float METsoftNomDn = TMath::Abs(hsMETsoftNom->GetBinContent(b) - hsMETsoftNom_dn->GetBinContent(b));
		float METsoftNomDif = (METsoftNomUp>METsoftNomDn) ? METsoftNomUp : METsoftNomDn;
		if(hsMETsoftNom->GetBinContent(b)>0) hsMETsoftNomUncert->SetBinContent(b,METsoftNomDif/hsMETsoftNom->GetBinContent(b)*100);	
		hsMETtrk->SetBinError(b,hsMETtrk->GetBinContent(b)*hsMETsoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hsMETtrkOrig = (TH1*)hsMETtrk->Clone("hsMETtrkOrig");
	hsMETtrkOrig->SetLineColor(kBlue);
	hsMETtrk->SetFillColor(kSpring);
	
	for(Int_t b=1 ; b<hdMETtrk->GetNbinsX() ; ++b)
	{	
		float METsoftNomUp = TMath::Abs(hdMETsoftNom->GetBinContent(b) - hdMETsoftNom_up->GetBinContent(b));
		float METsoftNomDn = TMath::Abs(hdMETsoftNom->GetBinContent(b) - hdMETsoftNom_dn->GetBinContent(b));
		float METsoftNomDif = (METsoftNomUp>METsoftNomDn) ? METsoftNomUp : METsoftNomDn;
		if(hdMETsoftNom->GetBinContent(b)>0) hdMETsoftNomUncert->SetBinContent(b,METsoftNomDif/hdMETsoftNom->GetBinContent(b)*100);	
		hdMETtrk->SetBinError(b,hdMETtrk->GetBinContent(b)*hdMETsoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hdMETtrkOrig = (TH1*)hdMETtrk->Clone("hdMETtrkOrig");
	hdMETtrkOrig->SetLineColor(kBlack);
	hdMETtrk->SetFillColor(kYellow);


	for(Int_t b=1 ; b<hsmTtrk->GetNbinsX() ; ++b)
	{	
		float mTsoftNomUp = TMath::Abs(hsmTsoftNom->GetBinContent(b) - hsmTsoftNom_up->GetBinContent(b));
		float mTsoftNomDn = TMath::Abs(hsmTsoftNom->GetBinContent(b) - hsmTsoftNom_dn->GetBinContent(b));
		float mTsoftNomDif = (mTsoftNomUp>mTsoftNomDn) ? mTsoftNomUp : mTsoftNomDn;
		if(hsmTsoftNom->GetBinContent(b)>0) hsmTsoftNomUncert->SetBinContent(b,mTsoftNomDif/hsmTsoftNom->GetBinContent(b)*100);	
		hsmTtrk->SetBinError(b,hsmTtrk->GetBinContent(b)*hsmTsoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hsmTtrkOrig = (TH1*)hsmTtrk->Clone("hsmTtrkOrig");
	hsmTtrkOrig->SetLineColor(kBlue);
	hsmTtrk->SetFillColor(kSpring);
	
	for(Int_t b=1 ; b<hdmTtrk->GetNbinsX() ; ++b)
	{	
		float mTsoftNomUp = TMath::Abs(hdmTsoftNom->GetBinContent(b) - hdmTsoftNom_up->GetBinContent(b));
		float mTsoftNomDn = TMath::Abs(hdmTsoftNom->GetBinContent(b) - hdmTsoftNom_dn->GetBinContent(b));
		float mTsoftNomDif = (mTsoftNomUp>mTsoftNomDn) ? mTsoftNomUp : mTsoftNomDn;
		if(hdmTsoftNom->GetBinContent(b)>0) hdmTsoftNomUncert->SetBinContent(b,mTsoftNomDif/hdmTsoftNom->GetBinContent(b)*100);	
		hdmTtrk->SetBinError(b,hdmTtrk->GetBinContent(b)*hdmTsoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hdmTtrkOrig = (TH1*)hdmTtrk->Clone("hdmTtrkOrig");
	hdmTtrkOrig->SetLineColor(kBlack);
	hdmTtrk->SetFillColor(kSpring);


	for(Int_t b=1 ; b<hsdPhitrk->GetNbinsX() ; ++b)
	{	
		float dPhisoftNomUp = TMath::Abs(hsdPhisoftNom->GetBinContent(b) - hsdPhisoftNom_up->GetBinContent(b));
		float dPhisoftNomDn = TMath::Abs(hsdPhisoftNom->GetBinContent(b) - hsdPhisoftNom_dn->GetBinContent(b));
		float dPhisoftNomDif = (dPhisoftNomUp>dPhisoftNomDn) ? dPhisoftNomUp : dPhisoftNomDn;
		if(hsdPhisoftNom->GetBinContent(b)>0) hsdPhisoftNomUncert->SetBinContent(b,dPhisoftNomDif/hsdPhisoftNom->GetBinContent(b)*100);	
		hsdPhitrk->SetBinError(b,hsdPhitrk->GetBinContent(b)*hsdPhisoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hsdPhitrkOrig = (TH1*)hsdPhitrk->Clone("hsdPhitrkOrig");
	hsdPhitrkOrig->SetLineColor(kBlue);
	hsdPhitrk->SetFillColor(kSpring);
	
	for(Int_t b=1 ; b<hddPhitrk->GetNbinsX() ; ++b)
	{	
		float dPhisoftNomUp = TMath::Abs(hddPhisoftNom->GetBinContent(b) - hddPhisoftNom_up->GetBinContent(b));
		float dPhisoftNomDn = TMath::Abs(hddPhisoftNom->GetBinContent(b) - hddPhisoftNom_dn->GetBinContent(b));
		float dPhisoftNomDif = (dPhisoftNomUp>dPhisoftNomDn) ? dPhisoftNomUp : dPhisoftNomDn;
		if(hddPhisoftNom->GetBinContent(b)>0) hddPhisoftNomUncert->SetBinContent(b,dPhisoftNomDif/hddPhisoftNom->GetBinContent(b)*100);	
		hddPhitrk->SetBinError(b,hddPhitrk->GetBinContent(b)*hddPhisoftNomUncert->GetBinContent(b)/100);
	}
	TH1* hddPhitrkOrig = (TH1*)hddPhitrk->Clone("hddPhitrkOrig");
	hddPhitrkOrig->SetLineColor(kBlack);
	hddPhitrk->SetFillColor(kSpring);
	
	



	TCanvas* cnv = new TCanvas("cnv", "", 1000,1200);
	cnv->Draw();
	cnv->Divide(2,3);
	cnv->cd(1); hsMETtrk->Draw("e2");  hsMETtrkOrig->Draw("hist same");  gPad->RedrawAxis(); gPad->Update();
	cnv->cd(2); hdMETtrk->Draw("e2");  hdMETtrkOrig->Draw("hist same");  gPad->RedrawAxis(); gPad->Update();
	cnv->cd(3); hsmTtrk->Draw("e2");   hsmTtrkOrig->Draw("hist same");   gPad->RedrawAxis(); gPad->Update();
	cnv->cd(4); hdmTtrk->Draw("e2");   hdmTtrkOrig->Draw("hist same");   gPad->RedrawAxis(); gPad->Update();
	cnv->cd(5); hsdPhitrk->Draw("e2"); hsdPhitrkOrig->Draw("hist same"); gPad->RedrawAxis(); gPad->Update();
	cnv->cd(6); hddPhitrk->Draw("e2"); hddPhitrkOrig->Draw("hist same"); gPad->RedrawAxis(); gPad->Update();
	cnv->SaveAs("METtrkUncert.pdf");
	
	
	TFile* fOut = new TFile("METtrkUncerts.root","RECREATE");
	fOut->cd();
	hsMETtrk->Write();
	hdMETtrk->Write();
	hsmTtrk->Write();
	hdmTtrk->Write();
	hsdPhitrk->Write();
	hddPhitrk->Write();
	fOut->Write();
	fOut->Close();
}