{
	TFile* fmvaout = new TFile("mvaout.muons.root","READ");
	TTree* tS = (TTree*)fmvaout->Get("fltmva_Wtaunu_3mu");
	

	TString srcuts        = "(m3body>1713 && m3body<1841)";
	
	TString bdtmincuts    = "score>-0.9";
	TString bdtoptcuts    = "score>+0.953";
	TString BDTcuts       = "("+bdtmincuts+" && "+bdtoptcuts+")";
	
	TString ncandcuts     = "(@m3body->size()==1)";
	
	TString iso           = "(iso003<0.1 && iso020<0.3)";
	TString trkspval      = "(trkspval>1.e-9)";
	TString distances     = "(lxysig>-10 && lxysig<50 && a0xysig<25)";
	TString kinematics    = "(pt3body>10000 && (mettrk>10000 && mettrk<250000 && mttrk>20000) && (metcal>10000 && metcal<250000 && mtcal>20000))";
	TString lowmass       = "(mOS1>300. && mOS2>300. && mSS>300.)";
	TString preTraining   = "("+iso+" && "+trkspval+" && "+distances+" && "+kinematics+" && "+lowmass+")";
	
	TString HT            = "(mttrk>60000 && mtcal>60000 && dphihttrk>2 && dphihtcal>2)";
	TString resonances    = "!( ((TMath::Abs(mOS1-1020)<50 || TMath::Abs(mOS2-1020)<50) || (TMath::Abs(mOS1-782)<50 || TMath::Abs(mOS2-782)<50)) && (mettrk<35000 || metcal<35000 || pt3body<35000) )";
	TString postTraining  = HT+" && "+resonances;
	
	TString cuts      = srcuts+" && "+ncandcuts+" && "+BDTcuts+" && "+preTraining+" && "+postTraining;
	cout << "cuts=" << cuts << "\n" << endl;


	Float_t npassedS   = tS->GetEntries(cuts);
	Float_t ninitS     = 99900;
	float AccEffSig    = npassedS/ninitS*100;
	cout << "AccEffSig=" << AccEffSig << "\%" << endl;
}


