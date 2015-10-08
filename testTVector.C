{
	TFile* f = new TFile("BDTfitResults.1450-1690-1870-2290.neg090.pos1000.root","READ");
	TVectorF* xSRmin = (TVectorF*)f->Get("xSRmin");
	TVectorF* xSRmax = (TVectorF*)f->Get("xSRmax");
	TVectorF* nSR00  = (TVectorF*)f->Get("nSR00");
	TVectorF* dnSR00 = (TVectorF*)f->Get("dnSR00");
	TVectorF* chi2SB  = (TVectorF*)f->Get("chi2SB"); 

	Int_t N = xSRmin->GetNoElements();

	cout << "N=" << N << endl;

	for(Int_t i=0 ; i<N ; ++i)
	{
		cout << "i=" << i << ", d=" << ((*xSRmax))[i]-((*xSRmin))[i] << " MeV --> " << "[" << ((*xSRmin))[i] << "," << ((*xSRmax))[i] << "] --> nSR0=" << ((*nSR00))[i] << "+-" << ((*dnSR00))[i] << endl;
	}

	cout << "chi2.0=" << ((*chi2SB))[0] << endl;
	cout << "chi2.1=" << ((*chi2SB))[1] << endl;
	cout << "chi2.2=" << ((*chi2SB))[2] << endl;
	cout << "chi2.3=" << ((*chi2SB))[3] << endl;
	
}
