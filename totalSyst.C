///////////////////////////////////
//// root -b -l -q totalSyst.C ////
///////////////////////////////////
{
	///////////////////////////////
	//// Sidebands=[1450,2290]
	//// Blinded=[1690,1870]
	//// SignalRegion=[1715,1843]
	//// BDT-cutoff=-0.9
	//// BDT-optimal=+0.964
	///////////////////////////////
	float nSB0      = 590;
	float nSR0      = 132.952;
	float f01       = 0.007848;
	///////////////////////////////
	float dnSB0stat = sqrt(nSB0);
	float dnSR0stat = dnSB0stat*(nSR0/nSB0);
	float nSR1      = nSR0*f01;
	float dnSR1stat = dnSR0stat*f01;
	///////////////////////////////
	float dnSR0syst_SBrange   = 19.103;   // 14.368%
	float dnSR0syst_BDTcutoff = 30.285;   // 22.779%
	float df01syst_SBrange    = 0.000493; // 6.280%
	float df01syst_BDTcutoff  = 0.00104;  // 13.25%
	///////////////////////////////
	float dnSR0syst_SBchoice = fabs(nSR0-126.537);   // %
	float df01syst_BDTchoice = fabs(f01-0.005208);   // %
	///////////////////////////////
	float dnSR0syst = sqrt(pow(dnSR0syst_SBrange,2)+pow(dnSR0syst_BDTcutoff,2)+pow(dnSR0syst_SBchoice,2));
	float df01syst  = sqrt(pow(df01syst_SBrange,2)+pow(df01syst_BDTcutoff,2)+pow(df01syst_BDTchoice,2));
	///////////////////////////////
	float dnSR1syst_tot  = nSR1*sqrt(pow(dnSR0syst/nSR0,2)+pow(df01syst/f01,2));
	///////////////////////////////
	cout << "dnSR0syst = "<<100*dnSR0syst_SBrange/nSR0<<"\%(SB range) [+] "<<100*dnSR0syst_BDTcutoff/nSR0<<"\%(BDT cutoff) [+] "<<100*dnSR0syst_SBchoice/nSR0<<"\%(SB fit function choice)"<<endl;
	cout << "df01syst  = "<<100*df01syst_SBrange/f01<<"\%(SB range) [+] "<<100*df01syst_BDTcutoff/f01<<"\%(BDT cutoff) [+] "<<100*df01syst_BDTchoice/f01<<"\%(BDT fit function choice)"<<endl;
	cout << "nSR0 = "<<nSR0<<" +- "<<dnSR0stat<<"(stat. "<<100*dnSR0stat/nSR0<<"\%)"<<" +- "<<dnSR0syst<<"(syst. "<<100*dnSR0syst/nSR0<<"\%)"<<endl;
	cout << "f01  = "<<f01<<" +- "<<df01syst<<"(syst. "<<100*df01syst/f01<<"\%)"<<endl;
	cout << "nSR1 = "<<nSR1<<" +- "<<dnSR1stat<<"(stat. "<<100*dnSR1stat/nSR1<<"\%)"<<" +- "<<dnSR1syst_tot<<"(syst. "<<100*dnSR1syst_tot/nSR1<<"\%)"<<endl;
}