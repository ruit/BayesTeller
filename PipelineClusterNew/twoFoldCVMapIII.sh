#!/bin/sh -login
# Tian R. <tianremiATgmail.com>
# Oct 30, 2014
# Nov 5, 2014
# Sample some patients simultaneously from germ /tumor for cross validation.
# Dec 3, always use qsub, always full path, always in parallel
# Dec 12, for panCancer
# Dec 15, 2014
input=$1
maxDist=$2
input2=$3
tempdir=$4
chr=$5


prep(){

path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
inputfile="panCancer.2in1"

cat $inputfile | cut -f3 | sed "s/,/\n/g" | sort | uniq > $tempdir"/"$totalPatientList
cat $inputfile | cut -f1,2 > $inputfile".germ"
cat $inputfile | cut -f1,3 > $inputfile".tumor"
}
#module load mcl


ExtractByChrS ()
	{
	input=$1
	chr=$2
	maxDist=$3
	tempdir=$4
	
	cat $input |sed "s/_/\t/g" | awk -v CHR=$chr '{ if ($1 == CHR) print $0}' > $tempdir"/"$input"_chr"$chr
	
	python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/generateAllPairs.py $tempdir"/"$input"_chr"$chr | awk  -v MAX=$maxDist -v CHR=$chr \
'function abs(x){return ((x < 0.0) ? -x : x)} \
{if (abs($2-$1) <= MAX) print CHR "_"$1 "\t"CHR"_"$2"\t"abs($2-$1)}' > $tempdir"/"$input"_chr"$chr".pair"
	
	rm $tempdir"/"$input"_chr"$chr
	}


#ExtractByChrS $tumorSample X 500 $tempdir
#chr="X"
#mcl $tempdir"/"$input"_chr"$chr".pair" --abc  -o $tempdir"/"$input"_chr"$chr".pair.cluster"


MCLbyChr(){

	module load mcl

	input=$1
	maxDist=$2
	input2=$3
	tempdir=$4
	chr=$5
	
	echo "Chr"$chr
	ExtractByChrS $input $chr $maxDist $tempdir
	mcl $tempdir"/"$input"_chr"$chr".pair" --abc  -o $tempdir"/"$input"_chr"$chr".pair.cluster"
	
	#@Oct 20, 2014 Mon
	grep $chr"_" $input | sed s/$chr"_"//g > $tempdir"/"$input"_chr"$chr
	grep $chr"_" $input2 | sed s/$chr"_"//g > $tempdir"/"$input2"_chr"$chr
	
	#check this for bugs
	cat $tempdir"/"$input"_chr"$chr".pair.cluster" | sed s/$chr"_"//g| awk '{min=$1;max=$1;for (i=1;i<=NF;i++){if($i < min){min=$i}; if ($i > max){max=$i}}; print min"\t"max}' > $tempdir"/"$input"_chr"$chr".pair.cluster.StartEnd"
	
	cat $tempdir"/"$input"_chr"$chr".pair.cluster.StartEnd" | cut -f1 | awk '{print $0"\tA_start"}' >  $tempdir"/"$input"_chr"$chr".pair.cluster.Start" 
	cat $tempdir"/"$input"_chr"$chr".pair.cluster.StartEnd" | cut -f2 | awk '{print $0"\tZ_end"}' >  $tempdir"/"$input"_chr"$chr".pair.cluster.End"

	cat  $tempdir"/"$input"_chr"$chr $tempdir"/"$input"_chr"$chr".pair.cluster.Start" $tempdir"/"$input"_chr"$chr".pair.cluster.End" | sort -k1,1 -n -k2,2 > $tempdir"/"$input"_chr"$chr".raw"
	
	###germline, tumor go together! Nov 5, 2014
	cat  $tempdir"/"$input2"_chr"$chr $tempdir"/"$input"_chr"$chr".pair.cluster.Start" $tempdir"/"$input"_chr"$chr".pair.cluster.End" | sort -k1,1 -n -k2,2 > $tempdir"/"$input2"_chr"$chr".raw"

	rm $tempdir"/"*"_chr"$chr".pair"*
	rm $tempdir"/"*"_chr"$chr
	}


MCLbyChr $input $maxDist $input2 $tempdir $chr

