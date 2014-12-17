#!/bin/sh -login
#Tian R.
#Dec 15, 2014

inputfile="panCancer.2in1"
tumorSample=$inputfile".tumor"
germlineSample=$inputfile".germ"



dist=$1

tempdir=$path"panCancer_CV_d"$dist"_"$RANDOM"_"`date | sed "s/ /_/g" | sed "s/:/_/g"`"_"$RANDOM"output"
mkdir $tempdir
totalPatientList="total_"$RANDOM".txt"



for chr in X 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22
do
	qsub -j y -cwd -o "/scratch/tianr/queue/Chr"$chr"_"$RANDOM".stdout" -b y "/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/twoFoldCVMapIII.sh $tumorSample $dist $germlineSample $tempdir $chr"	
		
done


