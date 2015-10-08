{
	TFile* f = new TFile("flatout.EW.muons.cuts.analysis.n0.j0.root","READ");

	TH1F* hMass = new TH1F("m3body",";m_{3body} [MeV];Events",50,0.,4000.);
	TH1F* hpT   = new TH1F("pT3body",";p_{T}^{3body} [MeV];Events",80,0.,100000.);

	
	
}
