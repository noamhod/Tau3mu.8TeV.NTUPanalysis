{
	TFile* f = new TFile("flatout.periodall+MC.muons.cuts.analysis.n0.j0.loose.root","READ");

	TTree* tS = (TTree*)f->Get("flatout_Wtaunu_3mu");
	TTree* tD = (TTree*)f->Get("flatout_Data");

	TH2F* hS1 = new TH2F("hS1",";m_{OS1}^{2} [MeV^{2}];m_{SS}^{2}  [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*2000.);
	TH2F* hS2 = new TH2F("hS2",";m_{OS1}^{2} [MeV^{2}];m_{OS2}^{2} [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*2000.);
	TH2F* hS3 = new TH2F("hS3",";m_{OS2}^{2} [MeV^{2}];m_{SS}^{2}  [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*2000.);

	TH2F* hD1 = new TH2F("hD1",";m_{OS1}^{2} [MeV^{2}];m_{SS}^{2}  [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*3000.);
	TH2F* hD2 = new TH2F("hD2",";m_{OS1}^{2} [MeV^{2}];m_{OS2}^{2} [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*3000.);
	TH2F* hD3 = new TH2F("hD3",";m_{OS2}^{2} [MeV^{2}];m_{SS}^{2}  [MeV^{2}];Events", 100,0.,3000.*3000., 100,0.,3000.*3000.);

	tS->Draw("vtx_mSS*vtx_mSS:vtx_mOS1*vtx_mOS1>>hS1","pass_loose==1");
	tS->Draw("vtx_mSS*vtx_mSS:vtx_mOS2*vtx_mOS2>>hS3","pass_loose==1");
	tS->Draw("vtx_mOS2*vtx_mOS2:vtx_mOS1*vtx_mOS1>>hS2","pass_loose==1");
	tD->Draw("vtx_mSS*vtx_mSS:vtx_mOS1*vtx_mOS1>>hD1","pass_loose==1");
	tD->Draw("vtx_mSS*vtx_mSS:vtx_mOS2*vtx_mOS2>>hD3","pass_loose==1");
	tD->Draw("vtx_mOS2*vtx_mOS2:vtx_mOS1*vtx_mOS1>>hD2","pass_loose==1");

	TCanvas* cnv1 = new TCanvas("cnv1","",1200,800);
	cnv1->Divide(3,2);

	cnv1->cd(1); hS1->Draw();
	cnv1->cd(2); hS3->Draw();
	cnv1->cd(3); hS2->Draw();

	cnv1->cd(4); hD1->Draw();
	cnv1->cd(5); hD3->Draw();
	cnv1->cd(6); hD2->Draw();

	cnv1->SaveAs("Dalitz.squared.pdf");



	TH2F* h1S1 = new TH2F("h1S1",";m_{OS1} [MeV];m_{SS}  [MeV];Events", 100,0.,4000., 100,0.,4000.);
	TH2F* h1S2 = new TH2F("h1S2",";m_{OS1} [MeV];m_{OS2} [MeV];Events", 100,0.,4000., 100,0.,4000.);
	TH2F* h1S3 = new TH2F("h1S3",";m_{OS2} [MeV];m_{SS}  [MeV];Events", 100,0.,4000., 100,0.,4000.);

	TH2F* h1D1 = new TH2F("h1D1",";m_{OS1} [MeV];m_{SS}  [MeV];Events", 100,0.,4000., 100,0.,4000.);
	TH2F* h1D2 = new TH2F("h1D2",";m_{OS1} [MeV];m_{OS2} [MeV];Events", 100,0.,4000., 100,0.,4000.);
	TH2F* h1D3 = new TH2F("h1D3",";m_{OS2} [MeV];m_{SS}  [MeV];Events", 100,0.,4000., 100,0.,4000.);

	tS->Draw("vtx_mSS:vtx_mOS1>>h1S1","pass_loose==1");
	tS->Draw("vtx_mSS:vtx_mOS2>>h1S3","pass_loose==1");
	tS->Draw("vtx_mOS2:vtx_mOS1>>h1S2","pass_loose==1");
	tD->Draw("vtx_mSS:vtx_mOS1>>h1D1","pass_loose==1");
	tD->Draw("vtx_mSS:vtx_mOS2>>h1D3","pass_loose==1");
	tD->Draw("vtx_mOS2:vtx_mOS1>>h1D2","pass_loose==1");

	TCanvas* cnv2 = new TCanvas("cnv2","",1200,800);
	cnv2->Divide(3,2);

	cnv2->cd(1); h1S1->Draw();
	cnv2->cd(2); h1S3->Draw();
	cnv2->cd(3); h1S2->Draw();

	cnv2->cd(4); h1D1->Draw();
	cnv2->cd(5); h1D3->Draw();
	cnv2->cd(6); h1D2->Draw();

	cnv2->SaveAs("Dalitz.linear.pdf");
	
	
	
	
	TH2F* h2S1 = new TH2F("h2S1",";m_{OS1} [MeV];a0_{xy};Events", 100,0.,4000., 100,0.,1.);
	TH2F* h2S2 = new TH2F("h2S2",";m_{OS2} [MeV]a0_{xy};Events", 100,0.,4000., 100,0.,1.);
	TH2F* h2S3 = new TH2F("h2S3",";m_{SS} [MeV];a0_{xy};Events", 100,0.,4000., 100,0.,1.);

	TH2F* h2D1 = new TH2F("h2D1",";m_{OS1} [MeV];a0_{xy};Events", 100,0.,4000., 100,0.,1.);
	TH2F* h2D2 = new TH2F("h2D2",";m_{OS2} [MeV];a0_{xy};Events", 100,0.,4000., 100,0.,1.);
	TH2F* h2D3 = new TH2F("h2D3",";m_{SS} [MeV];a0_{xy};Events", 100,0.,4000., 100,0.,1.);

	tS->Draw("vtx_pval:vtx_mOS1>>h2S1","pass_loose==1");
	tS->Draw("vtx_pval:vtx_mOS2>>h2S2","pass_loose==1");
	tS->Draw("vtx_pval:vtx_mSS>>h2S3","pass_loose==1");
	tD->Draw("vtx_pval:vtx_mOS1>>h2D1","pass_loose==1");
	tD->Draw("vtx_pval:vtx_mOS2>>h2D2","pass_loose==1");
	tD->Draw("vtx_pval:vtx_mSS>>h2D3","pass_loose==1");

	TCanvas* cnv3 = new TCanvas("cnv3","",1200,800);
	cnv3->Divide(3,2);

	cnv3->cd(1); h2S1->Draw();
	cnv3->cd(2); h2S2->Draw();
	cnv3->cd(3); h2S3->Draw();

	cnv3->cd(4); h2D1->Draw();
	cnv3->cd(5); h2D2->Draw();
	cnv3->cd(6); h2D3->Draw();

	cnv3->SaveAs("Dalitz.others.pdf");
}
