/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// root -b -l -q  compile.C  --chnl=all  --master=muons  --method=cuts  --jetmode=calib --trkmode=calib --metmode=calib  --mettype=muons --blinded=yes  [--ndata=1000]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	TString chnl    = "";
	TString master  = "";
	TString method  = ""; Bool_t bmethod = 0;
	TString jetMode = "";
	TString trkMode = "";
	TString metMode = "";
	TString metType = "";
	TString blinded = "";
	TString nData   = ""; // not mandatory - just for quick tests with fewer events
	Int_t nArgs    = 13;  // + one which is not mandatory
	for(int i=0 ; i<gApplication->Argc() ; i++) printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
	if(gApplication->Argc()<nArgs)
	{
		cout << "gApplication->Argc() = " << gApplication->Argc() << endl;
		cout << "should be: " << nArgs  << " -> exit." << endl;
		exit(-1);
	}
	else
	{
		for(int i=0 ; i<gApplication->Argc() ; i++)
		{
			TString arg = gApplication->Argv(i);
			if(arg.Contains("--chnl="))    chnl    = arg.ReplaceAll("--chnl=","");
			if(arg.Contains("--master="))  master  = arg.ReplaceAll("--master=","");
			if(arg.Contains("--method="))  method  = arg.ReplaceAll("--method=","");
			if(arg.Contains("--jetmode=")) jetMode = arg.ReplaceAll("--jetmode=","");
			if(arg.Contains("--trkmode=")) trkMode = arg.ReplaceAll("--trkmode=","");
			if(arg.Contains("--metmode=")) metMode = arg.ReplaceAll("--metmode=","");
			if(arg.Contains("--mettype=")) metType = arg.ReplaceAll("--mettype=","");
			if(arg.Contains("--blinded=")) blinded = arg.ReplaceAll("--blinded=","");
			if(arg.Contains("--ndata="))   nData   = arg.ReplaceAll("--ndata=",""); // not mandatory - just for quick tests with fewer events
		}
	}
	
	TString args = chnl+"\",\""+master+"\",\""+method+"\",\""+jetMode+"\",\""+trkMode+"\",\""+metMode+"\",\""+metType+"\",\""+blinded;
	if(nData!="") args += +"\",\""+nData;
	TString proc = "NTUPanalysis(\""+args+"\")";
	cout << proc << endl;
	gROOT->ProcessLine(".include .");
	gSystem->Load( "libCintex.so" );
	Cintex::Cintex::Enable();	
	gROOT->ProcessLine(".L Loader.C+");
	gROOT->ProcessLine(".L NTUPanalysis.C++");
	gROOT->ProcessLine(proc);
}
