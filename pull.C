{
	TFile* f = new TFile("fig.periodall+MC.muons.cuts.analysis.n0.j0.root","READ");
	TH1F* hX = (TH1F*)f->Get("Wtaunu_3mu_vtx_trumatch_x");
	TH1F* hY = (TH1F*)f->Get("Wtaunu_3mu_vtx_trumatch_y");
	TH1F* hZ = (TH1F*)f->Get("Wtaunu_3mu_vtx_trumatch_z");
	
	TFitResultPtr rX = hX->Fit("gaus","WLMS");
	TFitResultPtr rY = hY->Fit("gaus","WLMS");
	TFitResultPtr rZ = hZ->Fit("gaus","WLMS");
	
	stringstream strm;
	string       str1, str2;
	
	strm << setprecision(3) << fixed << rX->Parameter(2); strm >> str1; strm.clear(); 
	strm << setprecision(4) << fixed << rX->ParError(2);  strm >> str2; strm.clear(); 
	TString sX = "#sigma_{x}="+str1+"#pm"+str2;
	cout << sX << endl;
	strm.clear(); str1.clear(); str2.clear();
	
	strm << setprecision(3) << fixed << rY->Parameter(2); strm >> str1; strm.clear(); 
	strm << setprecision(4) << fixed << rY->ParError(2);  strm >> str2; strm.clear(); 
	TString sY = "#sigma_{y}="+str1+"#pm"+str2;
	cout << sY << endl;
	strm.clear(); str1.clear(); str2.clear();
	
	strm << setprecision(3) << fixed << rZ->Parameter(2); strm >> str1; strm.clear(); 
	strm << setprecision(4) << fixed << rZ->ParError(2);  strm >> str2; strm.clear(); 
	TString sZ = "#sigma_{z}="+str1+"#pm"+str2;
	cout << sZ << endl;
	strm.clear(); str1.clear(); str2.clear();
	
	TPaveText *pX = new TPaveText(0.55,0.75,0.90,0.9,"NDC");
	pX->SetFillStyle(4000); //will be transparent
	pX->SetFillColor(0);
	pX->SetTextFont(42);
	pX->SetBorderSize(0);
	pX->AddText(sX);
	
	TPaveText *pY = new TPaveText(0.55,0.75,0.90,0.9,"NDC");
	pY->SetFillStyle(4000); //will be transparent
	pY->SetFillColor(0);
	pY->SetTextFont(42);
	pY->SetBorderSize(0);
	pY->AddText(sY);
	
	TPaveText *pZ = new TPaveText(0.55,0.75,0.90,0.9,"NDC");
	pZ->SetFillStyle(4000); //will be transparent
	pZ->SetFillColor(0);
	pZ->SetTextFont(42);
	pZ->SetBorderSize(0);
	pZ->AddText(sZ);
	
	
	
	TCanvas* cnv = new TCanvas("cnv","",1500,500);
	cnv->Divide(3,1);
	cnv->Draw();
	cnv->cd(1);
	hX->Draw("p");
	pX->Draw("same");
	cnv->cd(2);
	hY->Draw("p");
	pY->Draw("same");
	cnv->cd(3);
	hZ->Draw("p");
	pZ->Draw("same");
	cnv->RedrawAxis();
	cnv->Update();
	cnv->SaveAs("pull.pdf");
	cnv->SaveAs("pull.eps");
}