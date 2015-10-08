{
	TH1F* h = new TH1F("h","Ratios;;Ratio",40,-1,+1);
	Int_t b = 2;
	for(double x=-1.00 ; x<1.00 ; x+=0.05)
	{
		if(x<-0.96) continue;
		cout << "x=" << x << " --> bincenter = " << h->GetBinCenter(b) << endl;
		b++;
	}
}
