#include "TF1.h"

#include <iostream>

using namespace std;


Double_t func(Double_t* x, Double_t* par)
{
	return par[0]+par[1]*x[0];
}

void testTF1()
{
	TF1* f = new TF1("func",func,-1,+1);
	f->SetParameter(0,8);
	f->SetParameter(1,9);
	cout << "GetExpFormula: " << f->GetExpFormula("CLING") << endl;
}
