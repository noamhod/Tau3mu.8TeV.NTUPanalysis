{
	TString name = "1450-1690-1870-2110.neg090.pos1000";
	cout << "name=" << name << endl;
	
	TString ranges = name;
	ranges.ReplaceAll(".",":");
	ranges.ReplaceAll("-",":");
	
	TPMERegexp regexp(":");
	regexp.Split(ranges);
	regexp.Print("all");
	Double_t mFullMin  = (Double_t)regexp[0].Atof();
	Double_t mBlindMin = (Double_t)regexp[1].Atof();
	Double_t mBlindMax = (Double_t)regexp[2].Atof();
	Double_t mFullMax  = (Double_t)regexp[3].Atof();
	cout << mFullMin << " -> " << mBlindMin << " -> " << mBlindMax << " -> " << mFullMax << endl;

}
