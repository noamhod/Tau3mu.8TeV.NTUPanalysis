# Tau3mu.8TeV.NTUPanalysis
Tau3mu 8TeV NTUPanalysis based on NTUPmaker output

Analysis flow:
-------------
* to get the trunk:
	svn co svn+ssh://svn.cern.ch/reps/atlasphys-exa/Physics/Exotic/Analysis/TauToMuX/ThreeMu/NTUPanalysis/trunk NTUPanalysis
* to tag a version in svn:
	See the latest tag in: https://svnweb.cern.ch/trac/atlasphys-exa/browser/Physics/Exotic/Analysis/TauToMuX
	svn cp svn+ssh://svn.cern.ch/reps/atlasphys-exa/Physics/Exotic/Analysis/TauToMuX/ThreeMu/NTUPanalysis/{trunk,tags/NTUPanalysis-XX-YY-ZZ} -m "tagname"

* need to have flatout.periodall+MC.muons.cuts.analysis.n0.j0.root in this directory
	* can be downloaded from /eos/atlas/user/h/hod/data/flat/
	  like that: eos cp /eos/atlas/user/h/hod/data/flat/flatout.periodall+MC.muons.cuts.analysis.n0.j0.root $HOME/

* root -b -l -q  compile.C	--chnl=all[or a single channel name]\
							--master=muons[or muid]\
							--method=cuts[or mva]\
							--jetmode=calib[uncalib, calib, JESUP, JESDWN, JERUP, JERDWN]\
							--metmode=uncalib[uncalib, calib, JESUP, JESDWN, JERUP, JERDWN]\
							--blinded=yes
	* take a look at:
		* the cutflow
		* figures.cuts.pdf
	* root -b -l -q fit3body.C++\(\"muons\",\"cuts\"\)  # or "muid"
		* take a look at the mass fit plots
	* root -b -l -q cutflow.C++\(\"muons\",\"cuts\"\) # or "muid"
		* take a look at the produced tex file and copy the big cutflow table
	* root -b -l -q objectquality.C++\(\"muons\",\"cuts\"\) # or "muid"
		* take a look at the objqa.pdf file
	
* root -l ./tmvaClass.C\(\"BDTG\"\)
	* produce the MVA plots...
	* root -l getUL.C
	* look for the minimum...
	* root -l UL.root
	* look at the 3 histos above BDTG score=0.9 and chose the best cut
	* update the cut value in fltvtx.C

* root -b -l -q  compile.C  --chnl=all  --master=muons  --method=mva
	* take a look at:
		* the cutflow
		* figures.mva.pdf
	* root -b -l -q fit3body.C++\(\"muons\",\"mva\"\)
		* take a look at the mass fit plots
	* root -b -l -q cutflow.C++\(\"muons\",\"mva\"\) # or "muid"
		* take a look at the produced tex file and copy the big cutflow table (just the line of the BDT cut)

* root -b -l -q  compile.C  --chnl=all  --master=muons  --method=mvaout (produce optimization input file without cutting on the score)
	* copy mvaout.root to the HistFitter directory on lxplus and run the mva optimization there:
	* scp mvaout.muons.root hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37/
	* root -b -l -q readmvaout.C++
		* look on the file: readmvaout.pdf (optimization plot is in the last page)
		* write down the best MVA cut
	* choose the 3body mass fit plot corresponding to the best MVA cut
		* change the expected background events and its error in python/MyUserAnalysis.py (verify all other parameters)
		* choose once the case for a downwards fluctuation in Nobs wrt the expectation
		* only after copying the CLs plot, run the same but for an upwards fluctuation
		* run it with: ./run.sh -l
		* write down the values in the fit result
		* look at the CLs result in file: results/upperlimit_cls_poi_Sig_Asym_CLs_grid_ts3.root.eps

* Proceed to tex build
	* first run: ./getEps.sh  cuts[mva]  muons[muid]
	* make
	
	
	
Short setup for optimization:
	* freeze loose cuts
	* root -b -l -q  compile.C  --chnl=all  --master=muons  --method=cuts    --jetmode=calib  --metmode=calib --mettype=muons --blinded=yes
	* Can train on a range slightly larger than the actual sidebands but without the blinded area:
		-- train: [1300,1690]+[1870,2320]
		-- apply: [1450,1690]+[1870,2290]
		-- note the b3body bin size (chosen as the resolution)
	* chose percentage of events (80%) to use form the signal / background passing events and train the BDT with:
		-- root -l ./tmvaClass.C\(\"1300.\",\"1690.\",\"1870.\",\"2320.\",\"muons\",\"BDTG\"\)
	* look at the training vs test samples
	* root -b -l -q  compile.C  --chnl=all  --master=muons  --method=mvaout  --jetmode=calib  --metmode=calib --mettype=muons --blinded=yes
	* change the min BDT cut and the BDT binnings if necessary (in "const.h")
	* root -l -b -q bdtroofit.C++
	* look on the fit(s) results
	* scp mvaout.muons.root BDTfitResults.root hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37/
	* ssh hod@lxplus.cern.ch
	* cd HistFitter/HistFitter-00-00-37; source setup.sh;
	* root -b -l -q readmvaout.C++\(1500.,1700.,1860.,2100.,20.,0.8000,0.9990,0.001,20.28\)
	* scp hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37/readmvaout.pdf readmvaout.fit.new.pdf
	* look on the optimization and update the optimal BDT cut in const.h
	* root -b -l -q  compile.C  --chnl=all  --master=muons  --method=mva     --jetmode=calib  --metmode=calib --mettype=muons --blinded=yes
	* finalisations:
		-- root -b -l -q cutflow.C++\(\"muons\",\"cuts\"\) # or "muid"
		-- root -b -l -q objectquality.C++\(\"muons\",\"cuts\"\) # or "muid"
		-- root -b -l -q phitest.C++
		-- root -b -l -q sidebandsBias.C++
		-- root -b -l -q sidebandsDef.C++
		-- source runsystsig.sh
		-- source runRangeSyst.sh dry/full
		-- root -b -l -q totalSyst.C
		-- root -b -l -q allFitSystAuto.C++
		-- source copyRangeSyst.sh
	* login lxplus
		-- cd HistFitter/HistFitter-00-00-37/
		-- source setup.sh
		-- root -b -l -q readmvaout.C++\(1450.,1690.,1870.,2290.,50.,0.9000,0.9950,0.005\)
	* back on mac:
		-- scp hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37/Optimization2D.1450-1690-1870-2290.*.pdf
		-- scp hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37/Optimization2D.1450-1690-1870-2290.*.root
		-- root -b -l -q optimization2Dreplot.C++
