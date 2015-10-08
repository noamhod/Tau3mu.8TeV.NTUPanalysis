#!/bin/sh

### trainig:                1250-->2300
### nominal left sideband   1450-->1690
### nominal right sideband  1870-->2110
###   nominal +1bin +2bins +3bins +4bins +5bins +6bins -3bins -2bins -1bin
xmin=(1450    1420  1390   1360   1330   1300   1270   1540   1510   1480)
xmax=(2110    2140  2170   2200   2230   2260   2290   2020   2050   2080)
i=0
while [ $i -lt ${#xmin[*]} ]; do
	scp BDTfitResults.${xmin[$i]}-1690-1870-${xmax[$i]}.neg090.pos1000.root hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
	i=$(($i+1));
done

xcutoff=(-0.89   -0.88   -0.87   -0.86   -0.85   -0.84   -0.83   -0.82   -0.81   -0.80)
scutoff=(neg089  neg088  neg087  neg086  neg085  neg084  neg083  neg082  neg081  neg080)
i=0
while [ $i -lt ${#xcutoff[*]} ]; do
	scp BDTfitResults.1450-1690-1870-2110.${scutoff[$i]}.pos1000.root hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
	i=$(($i+1));
done

scp mvaout.muons.root hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp std.h             hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp type.h            hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp const.h           hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp postBDTcuts.h     hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp allFitSystAuto.h  hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp readmvaout.C      hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp optimization2D.C  hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/
scp MyUserAnalysis.py hod@lxplus.cern.ch:HistFitter/HistFitter-00-00-37_tmp0/python/

