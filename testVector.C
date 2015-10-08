{
	vector<float> v;
	v.push_back(0);
	v.push_back(0);
	v.push_back(0);
	v.push_back(7);
	v.push_back(6);
	v.push_back(5);
	v.push_back(4);
	v.push_back(3);
	v.push_back(2);
	v.push_back(1);
	v.push_back(0);
	v.push_back(0);


	unsigned int i0=0;
	float val0 = -1;
	for(unsigned int i=0 ; i<v.size() ; ++i) { if(v[i]<1.e-3) i0=i; else { val0=v[i]; break; } }
	for(unsigned int j=0 ; j<=i0      ; ++j) { v[j] = val0; }
	
	unsigned int i1=0;
	float val1 = -1;
	for(unsigned int i=(v.size()-1) ; i>=0  ; --i) { if(v[i]<1.e-3) i1=i; else { val1=v[i]; break; } }
	for(unsigned int j=(v.size()-1) ; j>=i1 ; --j) { v[j] = val1; }


	for(unsigned int i=0 ; i<v.size() ; ++i) cout << "v[" << i << "]=" << v[i] << endl;

}
