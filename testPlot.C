{
	TH1F* h1 = new TH1F("hgaus","hgaus",100,-3,3);
        h1->FillRandom("gaus",10000);

	TPaveText* ptxt = new TPaveText(0.12,0.65,0.37,0.87,"NDC");
	ptxt->SetFillStyle(4000); //will be transparent
	ptxt->SetFillColor(0);
	ptxt->SetTextFont(42);
	ptxt->SetBorderSize(0);
	ptxt->AddText("#bf{#it{ATLAS}} internal");
	ptxt->AddText("#scale[0.55]{#int}Ldt=20.28 fb^{-1}");
	ptxt->AddText("#sqrt{s}=8 TeV");

	TCanvas* c1 = new TCanvas("c1","c1",600,400);
        c1->SetTicks(1,1);
        c1->Draw();
        h1->Draw("p0 e0");
	ptxt->Draw("same");
}
