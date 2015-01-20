#!/bin/sh -login
# Tian R. <tianremiATgmail.com>
# Oct 30, 2014
# Nov 5, 2014
# Sample some patients simultaneously from germ /tumor for cross validation.
# Dec 3, always use qsub, always full path, always in parallel
# Dec 12, for panCancer
# Dec 22, freq, freq, freq
# 1hr


prep(){

path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
inputfile="panCancer.2in1"

tempdir=$path"panCancer_CV_d"$dist"_"$RANDOM"_"`date | cut -d " " -f5`"_"$RANDOM"output"
mkdir $tempdir
totalPatientList="total_"$RANDOM".txt"

cat $inputfile | cut -f3 | sed "s/,/\n/g" | sort | uniq > $tempdir"/"$totalPatientList
cat $inputfile | cut -f1,2 > $inputfile".germ"
cat $inputfile | cut -f1,3 > $inputfile".tumor"

#module load mcl

tumorSample=$inputfile".tumor"
germlineSample=$inputfile".germ"
}


inputfile="panCancer.2in1"
totalPatientList="panCancer.2in1_total_15734.txt"
tempdir="CV2_panCancer_"$RANDOM"output"

germlineSample="panCancer.2in1.germ"
tumorSample="panCancer.2in1.tumor"
mkdir $tempdir

python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/CrossValidationNoClustering.py \
$inputfile $totalPatientList $tempdir



CrossVD(){

	tempdir=$1	
	germlineSample=$2
	tumorSample=$3

	#sorting based on freq mutation in tumor samples, stick to it! Dec 23, 2014
	cat $tempdir"/modelA" | sort -k3 -n -r | gzip > $tempdir"/modelA.sorted.gz"
	rm $tempdir"/modelA"

	cat $tempdir"/modelB" | sort -k3 -n -r | gzip > $tempdir"/modelB.sorted.gz"
	rm $tempdir"/modelB"

	cat $tempdir"/testPatList4A.tab" | while IFS=\t read testPatient
		do
 			#echo $testPatient
			zcat -f $tumorSample | grep -w $testPatient | cut -f1 | gzip >$tempdir"/"$tumorSample"_"$testPatient".snvsa.gz"
			zcat -f $germlineSample | grep -w $testPatient | cut -f1 | gzip >$tempdir"/"$germlineSample"_"$testPatient".snvsa.gz"


		done

	cat $tempdir"/testPatList4B.tab" | while IFS=\t read testPatient
		do
 			#echo $testPatient
			zcat -f $tumorSample | grep -w $testPatient | cut -f1 | gzip >$tempdir"/"$tumorSample"_"$testPatient".snvsb.gz"
			zcat -f $germlineSample | grep -w $testPatient | cut -f1 | gzip >$tempdir"/"$germlineSample"_"$testPatient".snvsb.gz"

		done
		
	}


CrossVD $tempdir $germlineSample $tumorSample


#X_11
#1_222

#A model file
#ChrXXX:111111,111111 Ntumor Ntotal backgroundRate CRscore
 
#Chr17:7578177,7578177,7578190,7578208,7578211,7578212,7578217,7578235,7578253,7578257,7578263,7578280,7578369,7578370,7578388,7578394,7578403,7578406,7578407,7578442,7578445,7578449,7578454,7578455,7578457,7578461,7578503,7578526,7578550	52	220	0	120.63978567313
#Chr12:25398214,25398214,25398218,25398226,25398255,25398262,25398281,25398282,25398284,25398285	75	220	0.00454545	115.906452155921
