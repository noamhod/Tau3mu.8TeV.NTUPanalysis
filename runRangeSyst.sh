#!/bin/sh

if [ -z "$1" ] ; then
	echo "ERROR, you should provide the following arguments:"
	echo "source runRangeSyst.sh  dry[full]"
	return
fi


# Below we define the factorial function in bc syntax
relDif="define frel (x,y) {
	return sqrt((100-100*y/x)^2)
}"
absDif="define fabs (x,y) {
	return sqrt((x-y)^2)
}"
fmax () {
	max=-9999999999;
	arr=($@)
	for num in ${arr[@]}; do
		if [ 1 -eq "$(echo "${num} > ${max}" | bc)" ] ;then  
		    max=${num}
		fi
	done
	echo $max
}

### trainig:                1250-->2300
### nominal left sideband   1450-->1690
### nominal right sideband  1870-->2110
###   nominal +1bin +2bins +3bins +4bins +5bins +6bins -3bins -2bins -1bin
xmin=(1450    1420  1390   1360   1330   1300   1270   1540   1510   1480)
xmax=(2110    2140  2170   2200   2230   2260   2290   2020   2050   2080)
nSR0array=()
f01array=()
nSR1array=()
i=0
while [ $i -lt ${#xmin[*]} ]; do
	
	if [[ "$1" == *full* ]] ; then
		root -l -b -q bdtroofit.C++\(${xmin[$i]},${xmax[$i]},-0.9,+1\) > BkgEstimate.${xmin[$i]}-1690-1870-${xmax[$i]}.neg090.pos1000.log; 
		open figures/BkgEstimate.${xmin[$i]}-1690-1870-${xmax[$i]}.neg090.pos1000.pdf
	fi

	j=0
	filename=BkgEstimate.${xmin[$i]}-1690-1870-${xmax[$i]}.neg090.pos1000.txt
	while read line
	do
		if [ $j -eq 14 ] ; then
			x=`echo $line | awk '{print $4}'`;
			nSR0array+=(${x});
			# echo "nSR0=${x}";
		elif [ $j -eq 17 ] ; then
			x=`echo $line | awk '{print $3}'`;
			f01array+=(${x});
			# echo "f01=${x}";
		elif [ $j -eq 19 ] ; then
			x=`echo $line | awk '{print $3}'`;
			nSR1array+=(${x});
			# echo "nSR0=${x}";
		fi
		j=$(($j+1));
	done < <(cat $filename)
	i=$(($i+1));
done




xcutoff=(-0.89   -0.88   -0.87   -0.86   -0.85   -0.84   -0.83   -0.82   -0.81   -0.80)
scutoff=(neg089  neg088  neg087  neg086  neg085  neg084  neg083  neg082  neg081  neg080)
nSR0array1=()
f01array1=()
nSR1array1=()
i=0
while [ $i -lt ${#xcutoff[*]} ]; do
	
	if [[ "$1" == *full* ]] ; then
		root -l -b -q bdtroofit.C++\(1450,2110,${xcutoff[$i]},+1\) > BkgEstimate.1450-1690-1870-2110.${scutoff[$i]}.pos1000.log;
		open figures/BkgEstimate.1450-1690-1870-2110.${scutoff[$i]}.pos1000.pdf
	fi

	j=0
	filename=BkgEstimate.1450-1690-1870-2110.${scutoff[$i]}.pos1000.txt
	while read line
	do
		if [ $j -eq 14 ] ; then
			x=`echo $line | awk '{print $4}'`;
			nSR0array1+=(${x});
			# echo "nSR0=${x}";
		elif [ $j -eq 17 ] ; then
			x=`echo $line | awk '{print $3}'`;
			f01array1+=(${x});
			# echo "f01=${x}";
		elif [ $j -eq 19 ] ; then
			x=`echo $line | awk '{print $3}'`;
			nSR1array1+=(${x});
			# echo "nSR0=${x}";
		fi
		j=$(($j+1));
	done < <(cat $filename)
	i=$(($i+1));
done
echo ""
echo ""




#x1max=(+0.800    +0.900   +0.964)
#sx1max=(pos0800  pos0900  pos0964)
#i=0
#while [ $i -lt ${#x1max[*]} ]; do
#	
#	if [[ "$1" == *full* ]] ; then
#		root -l -b -q bdtroofit.C++\(1450,2110,-0.9,${x1max[$i]}\);
#		open figures/BkgEstimate.1450-1690-1870-2110.neg090.${sx1max[$i]}.pdf
#	fi
#	i=$(($i+1));
#done
#echo ""
#echo ""


#### cut based
#nSRarray=()
#i=0
#while [ $i -lt ${#xmax[*]} ]; do
#	
#	if [[ "$1" == *full* ]] ; then
#		root -l -b -q cutsroofit.C++\(${xmin[$i]},${xmax[$i]}\);
#		open figures/CutBasedFitResults.${xmin[$i]}-1690-1870-${xmax[$i]}.pdf
#	fi
#	
#	j=0
#	filename=CutBasedFitResults.${xmin[$i]}-1690-1870-${xmax[$i]}.txt
#	while read line
#	do
#		if [ $j -eq 12 ] ; then
#			x=`echo $line | awk '{print $3}'`;
#			nSRarray+=(${x});
#			# echo "nSR=${x}";
#		fi
#		j=$(($j+1));
#	done < <(cat $filename)
#	
#	i=$(($i+1));
#done
#echo ""
#echo ""



### put everything in a file
fout=rangesSyst.txt
echo > ${fout}




nSR0relDif=(); nSR0absDif=();
f01relDif=();  f01absDif=(); 
nSR1relDif=(); nSR1absDif=();
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on nSR0 due to SB definition" >> ${fout} 2>&1
# echo ${#nSR0array[@]} #Number of elements in the array
nSR0nom=${nSR0array[0]}
i=0
while [ $i -lt ${#xmin[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${nSR0nom},${nSR0array[$i]})" | bc -l`; nSR0relDif+=(${reldif});
	absdif=`echo "$absDif;scale=3;fabs(${nSR0nom},${nSR0array[$i]})" | bc -l`; nSR0absDif+=(${absdif});
	echo "SB=[${xmin[$i]},${xmax[$i]}]: nSR0nominal=${nSR0nom} --> nSR0array=${nSR0array[$i]}"
	echo "SB=[${xmin[$i]},${xmax[$i]}], BDT-cutoff=-0.9, nSR0nominal=${nSR0nom}:  nSR0=${nSR0array[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${nSR0relDif[@]}`%
maxabsdif=`fmax ${nSR0absDif[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on f01 due to SB definition" >> ${fout} 2>&1
# echo ${#f01array[@]} #Number of elements in the array
f01nom=${f01array[0]}
i=0
while [ $i -lt ${#xmin[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${f01nom},${f01array[$i]})" | bc -l`; nf01relDif+=(${reldif});
	absdif=`echo "$absDif;fabs(${f01nom},${f01array[$i]})" | bc -l`        ; nf01absDif+=(${absdif});
	echo "SB=[${xmin[$i]},${xmax[$i]}], BDT-cutoff=-0.9, f01nominal=${f01nom}:  f01=${f01array[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${nf01relDif[@]}`%
maxabsdif=`fmax ${nf01absDif[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on nSR1 due to SB definition" >> ${fout} 2>&1
# echo ${#nSR1array[@]} #Number of elements in the array
nSR1nom=${nSR1array[0]}
i=0
while [ $i -lt ${#xmin[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${nSR1nom},${nSR1array[$i]})" | bc -l`; nSR1relDif+=(${reldif});
	absdif=`echo "$absDif;scale=3;fabs(${nSR1nom},${nSR1array[$i]})" | bc -l`; nSR1absDif+=(${absdif});
	echo "SB=[${xmin[$i]},${xmax[$i]}], BDT-cutoff=-0.9, nSR1nominal=${nSR1nom}:  nSR1=${nSR1array[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${nSR1relDif[@]}`%
maxabsdif=`fmax ${nSR1absDif[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
echo "" >> ${fout} 2>&1
echo "" >> ${fout} 2>&1



nSR0relDif1=(); nSR0absDif1=();
f01relDif1=();  f01absDif1=(); 
nSR1relDif1=(); nSR1absDif1=();
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on nSR0 due to BDT cutoff" >> ${fout} 2>&1
# echo ${#nSR0array1[@]} #Number of elements in the array
nSR0nom=${nSR0array[0]}
i=0
while [ $i -lt ${#xcutoff[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${nSR0nom},${nSR0array1[$i]})" | bc -l`; nSR0relDif1+=(${reldif});
	absdif=`echo "$absDif;scale=3;fabs(${nSR0nom},${nSR0array1[$i]})" | bc -l`; nSR0absDif1+=(${absdif});
	echo "SB=[1450,2110], BDT-cutoff=${xcutoff[$i]}, nSR0nominal=${nSR0nom}:  nSR0=${nSR0array1[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${nSR0relDif1[@]}`%
maxabsdif=`fmax ${nSR0absDif1[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on f01 due to BDT cutoff" >> ${fout} 2>&1
# echo ${#f01array1[@]} #Number of elements in the array
f01nom=${f01array[0]}
i=0
while [ $i -lt ${#xcutoff[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${f01nom},${f01array1[$i]})" | bc -l`; f01relDif1+=(${reldif});
	absdif=`echo "$absDif;fabs(${f01nom},${f01array1[$i]})" | bc -l`        ; f01absDif1+=(${absdif});
	echo "SB=[1450,2110], BDT-cutoff=${xcutoff[$i]}, f01nominal=${f01nom}:  f01=${f01array1[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${f01relDif1[@]}`%
maxabsdif=`fmax ${f01absDif1[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
echo "----------------------------------------" >> ${fout} 2>&1
echo "Shifts on nSR1 due to BDT cutoff" >> ${fout} 2>&1
# echo ${#nSR1array1[@]} #Number of elements in the array
nSR1nom=${nSR1array[0]}
i=0
while [ $i -lt ${#xcutoff[*]} ]; do
	reldif=`echo "$relDif;scale=3;frel(${nSR1nom},${nSR1array1[$i]})" | bc -l`; nSR1relDif1+=(${reldif});
	absdif=`echo "$absDif;scale=3;fabs(${nSR1nom},${nSR1array1[$i]})" | bc -l`; nSR1absDif1+=(${absdif});
	echo "SB=[1450,2110], BDT-cutoff=${xcutoff[$i]}, nSR1nominal=${nSR1nom}:  nSR1=${nSR1array1[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
	i=$(($i+1));
done
maxreldif=`fmax ${nSR1relDif1[@]}`%
maxabsdif=`fmax ${nSR1absDif1[@]}`
echo "Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
echo "Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1





#nSRrelDif=(); nSRabsDif=();
#echo "----------------------------------------" >> ${fout} 2>&1
#echo "Cut-based: Shifts on nSR due to SB definition" >> ${fout} 2>&1
#nSRnom=${nSRarray[0]}
#i=0
#while [ $i -lt ${#xmin[*]} ]; do
#	reldif=`echo "$relDif;scale=3;frel(${nSRnom},${nSRarray[$i]})" | bc -l`; nSRrelDif+=(${reldif});
#	absdif=`echo "$absDif;scale=3;fabs(${nSRnom},${nSRarray[$i]})" | bc -l`; nSRabsDif+=(${absdif});
#	echo "Cut-based: SB=[${xmin[$i]},${xmax[$i]}], nSRnominal=${nSRnom}:  nSR=${nSRarray[$i]} --> Abs.Error=${absdif}  Rel.Error=${reldif}%" >> ${fout} 2>&1
#	i=$(($i+1));
#done
#maxreldif=`fmax ${nSRrelDif[@]}`%
#maxabsdif=`fmax ${nSRabsDif[@]}`
#echo "Cut-based: Max Rel.Error: ${maxreldif}" >> ${fout} 2>&1
#echo "Cut-based: Max Abs.Error: ${maxabsdif}" >> ${fout} 2>&1
#echo "" >> ${fout} 2>&1
#echo "" >> ${fout} 2>&1






echo "Take a look at ${fout} !"
