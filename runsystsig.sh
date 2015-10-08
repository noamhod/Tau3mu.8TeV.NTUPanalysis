#!/bin/sh

if [ -z "$1" ] ; then
	echo "ERROR, you should provide the following arguments:"
	echo "source runsystsig.sh Wtaunu_3mu[Data]"
	return
fi

channel=$1
# channel="Wtaunu_3mu"
# channel="Data"

hname=histos.muons.mva

names=(
	"uncalib" "calib"
	"JESUP" "JESDWN"
	"JERUP" "JERDWN" 
)

for name in "${names[@]}"
do
	root -b -l -q  compile.C  --chnl=${channel}  --master=muons  --method=mva  --jetmode=${name}  --metmode=${name}  --mettype=muons;
	/bin/cp -f ${channel}.${hname}.root ${channel}.${hname}.${name}.root
	echo ""
done
root -b -l -q syst.C++\(\"${channel}\"\)
open figures/${channel}.syst.pdf