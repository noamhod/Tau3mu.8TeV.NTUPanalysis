#ifndef CLONETREE_H
#define CLONETREE_H

#include "type.h"

TMapTSP2TTREE finaltrees;
TDirectory* olddir3 = gDirectory;
TFile* fclonedtrees;

void initFinalTreesFile(TString fname) { fclonedtrees = new TFile(fname,"RECREATE"); }
void initFinalTree(TString name,TTree* t) { finaltrees.insert(make_pair(name,t>CloneTree(0))); }
void finitFinalTrees()
{
	fclonedtrees->cd();
	for(TMapTSP2TTREE::iterator it=finaltrees.begin() ; it!=finaltrees.end() ; it++)
	{
		it->second->AutoSave();
		// it->second->Write("", TObject::kOverwrite);
	}
	fclonedtrees->Write();
	fclonedtrees->Close();
	olddir3->cd();
}

void fillFinalTrees(int N)
{
	fclonedtrees->cd();
	for(TMapTSP2TTREE::iterator it=finaltrees.begin() ; it!=finaltrees.end() ; it++)
	{
		TFile* fout = it->second->GetCurrentFile();
		fout->cd();
		it->second->Fill();
		if(N%10000==0 && N!=0)
		{
			it->second->FlushBaskets();
			// it->second->Write("", TObject::kOverwrite);
		}	
	}
	olddir3->cd();
}

#endif