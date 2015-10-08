{
	TString s  = "Wtaunu_3mu"; // Data / Wtaunu_3mu
	TString ss = "Wtaunu_3mu"; // periodG / Wtaunu_3mu
	TFile* f = new TFile("fig."+ss+".muons.cuts.analysis.n0.j0.root","READ");
	TH1F* sumet_tool = (TH1F*)f->Get(s+"_mettrk_tool_sumet");       sumet_tool->SetLineColor(kRed);        sumet_tool->SetFillStyle(0);
	TH1F* sumet_ntuple = (TH1F*)f->Get(s+"_mettrk_ntuple_sumet");   sumet_ntuple->SetLineColor(kBlack);    sumet_ntuple->SetFillStyle(0);
	TH1F* sumet_ntuptru = (TH1F*)f->Get(s+"_mettrk_ntuptru_sumet"); sumet_ntuptru->SetLineColor(kOrange+1);    sumet_ntuptru->SetFillStyle(0);
	TH1F* et_tool = (TH1F*)f->Get(s+"_mettrk_tool_et");             et_tool->SetLineColor(kRed);           et_tool->SetFillStyle(0);
	TH1F* et_ntuple = (TH1F*)f->Get(s+"_mettrk_ntuple_et");         et_ntuple->SetLineColor(kBlack);       et_ntuple->SetFillStyle(0);
	TH1F* et_ntuptru = (TH1F*)f->Get(s+"_mettrk_ntuptru_et");         et_ntuptru->SetLineColor(kOrange+1);       et_ntuptru->SetFillStyle(0);
	TH1F* phi_tool = (TH1F*)f->Get(s+"_mettrk_tool_phi");           phi_tool->SetLineColor(kRed);          phi_tool->SetFillStyle(0);
	TH1F* phi_ntuple = (TH1F*)f->Get(s+"_mettrk_ntuple_phi");       phi_ntuple->SetLineColor(kBlack);      phi_ntuple->SetFillStyle(0);
	TH1F* phi_ntuptru = (TH1F*)f->Get(s+"_mettrk_ntuptru_phi");       phi_ntuptru->SetLineColor(kOrange+1);      phi_ntuptru->SetFillStyle(0);
	TH1F* et_diff = (TH1F*)f->Get(s+"_mettrk_toolVSntuple_et_reldiff");   et_diff->SetLineColor(kBlue);    et_diff->SetFillStyle(0);
	TH1F* et_diff_tru = (TH1F*)f->Get(s+"_mettrk_toolVSntuptru_et_reldiff");   et_diff_tru->SetLineColor(kOrange+1);    et_diff_tru->SetFillStyle(0);
	TH1F* phi_diff = (TH1F*)f->Get(s+"_mettrk_toolVSntuple_phi_reldiff"); phi_diff->SetLineColor(kBlue);   phi_diff->SetFillStyle(0);
	TH1F* phi_diff_tru = (TH1F*)f->Get(s+"_mettrk_toolVSntuptru_phi_reldiff"); phi_diff_tru->SetLineColor(kOrange+1);   phi_diff_tru->SetFillStyle(0);



	TLegend* leg = new TLegend(0.55,0.6,0.83,0.83,NULL,"brNDC");
	if(!s.Contains("Data")) leg->AddEntry(sumet_ntuptru,"Truth","l");
	leg->AddEntry(sumet_ntuple,"Ntuple","l");	
	leg->AddEntry(sumet_tool,"Tool","l");	

	TLegend* legdif = new TLegend(0.55,0.6,0.83,0.83,NULL,"brNDC");
	legdif->AddEntry(et_diff,"Tool vs Ntuple","l");
	if(!s.Contains("Data")) legdif->AddEntry(et_diff_tru,"Tool vs Truth","l");


	TCanvas* cnv = new TCanvas("cnv","",1000,1000);
	cnv->Divide(2,2);
	cnv->cd();
	cnv->Draw();

	cnv->cd(1);
	gPad->SetTicks(1,1);
	gPad->SetLogy();
	if(!s.Contains("Data")) sumet_ntuptru->Draw("hist");
	sumet_tool->Draw("hist same");
	sumet_ntuple->Draw("hist same");
	leg->Draw("same");

	cnv->cd(2);
	gPad->SetTicks(1,1);
	if(!s.Contains("Data")) et_ntuptru->Draw("hist");
	et_ntuple->Draw("hist same");
	et_tool->Draw("hist same");
	leg->Draw("same");

	cnv->cd(3);
	gPad->SetTicks(1,1);
	if(!s.Contains("Data")) phi_ntuptru->Draw("hist");
	phi_tool->Draw("hist same");
	phi_ntuple->Draw("hist same");
	leg->Draw("same");

	cnv->SaveAs("mettrk.pdf(");	
	delete cnv;





	cnv = new TCanvas("cnv","",1000,500);
	cnv->Divide(2,1);
	cnv->cd();
	cnv->Draw();

	cnv->cd(1);
	gPad->SetTicks(1,1);
	et_diff->Draw("hist");
	if(!s.Contains("Data")) et_diff_tru->Draw("hist same");
	if(!s.Contains("Data")) legdif->Draw("same");

	cnv->cd(2);
	gPad->SetTicks(1,1);
	phi_diff->Draw("hist");
	if(!s.Contains("Data")) phi_diff_tru->Draw("hist same");

	cnv->SaveAs("mettrk.pdf)");
	delete cnv;
}
