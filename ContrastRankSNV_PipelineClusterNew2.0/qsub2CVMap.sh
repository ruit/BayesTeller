#!/bin/sh -login
#Tian R.
#Dec 15, 2014
#>15hr

inputfile="panCancer.2in1"
tumorSample=$inputfile".tumor"
germlineSample=$inputfile".germ"




totalPatientList="panCancer.2in1_total_15734.txt"
#cat $inputfile | cut -f3 | sed "s/,/\n/g" | sort | uniq > $totalPatientList


runMap (){


	dist=$1

	tempdir=$path"panCancerCV_d"$dist"_"$RANDOM"_"`date | sed "s/ /_/g" | sed "s/:/_/g"`"Map"
	mkdir $tempdir
	cp $totalPatientList $tempdir"/"

	for chr in X 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22
	do
		qsub -j y -cwd -o "/scratch/tianr/queue/Chr"$chr"_"$RANDOM".stdout" -b y "/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/twoFoldCVMapIII.sh $tumorSample $dist $germlineSample $tempdir $chr"	
		
	done
	}	


loop (){
	for i in 1 50 500 2000 10000 100000 1000000
	do
		runMap $i
	done
	}

loop
