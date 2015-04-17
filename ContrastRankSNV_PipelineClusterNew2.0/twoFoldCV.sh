#!/bin/sh -login
# Tian R. <tianremiATgmail.com>
# Oct 30, 2014
# Nov 5, 2014
# Sample some patients simultaneously from germ /tumor for cross validation.


# Dec 3, always use qsub, always full path, always in parallel

cancer=$1
dist=$2
nTest=$3 # number of patients masked
#topN=1000 # default is 1000
sampleSize="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/TCGAsampleN_hash_table.txt"

module load mcl


# path to tempdir and *.filtered
path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
#path="/home/tianr/1Projects/1SNVblocks/RUN_Oct23/old/"


#Dec 8, 2014
#dist, topN are 2 paras
tempdir=$path$cancer"_CV_d"$dist"_"$RANDOM"_"`date | cut -d " " -f5`"_"$RANDOM"output"


total=`cat $sampleSize | awk -v CA=$cancer '{if($1==CA) print $2}'`
NsampleLeft=`expr $total - $nTest`
#echo $NsampleLeft



tumorSample="all"$cancer"Cancertumor.filtered"
germlineSample="all"$cancer"Cancergermline.filtered"

#topN=20
#rank the SNV clusters by ContraskRank score, take the topN as the model





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


RunmclAllChrs () 
	{
	input=$1
	maxDist=$2
	input2=$3
	tempdir=$4

	for chr in X 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22
	do
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
		
		python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/CrossValidation.py $tempdir"/"$input"_chr"$chr".raw" $tempdir"/"$input2"_chr"$chr".raw" $nTest $total $tempdir"/testPatientsList.txt"		


		cat $tempdir"/"$input"_chr"$chr".raw.out"| awk -v CHR=$chr '{print "Chr"CHR":"$0}' >$tempdir"/"$input"_chr"$chr".raw.out2"
		cat $tempdir"/"$input2"_chr"$chr".raw.out"| awk -v CHR=$chr '{print "Chr"CHR":"$0}' >$tempdir"/"$input2"_chr"$chr".raw.out2"

		rm $tempdir"/"*"_chr"$chr".pair"*
		rm $tempdir"/"*"_chr"$chr
	done

	cat $tempdir"/"$input*.raw.out2  > $tempdir"/"$input".snvfreq"
	cat $tempdir"/"$input2*.raw.out2  > $tempdir"/"$input2".snvfreq"

	rm $tempdir"/"*raw*

	}
### Oct 21, 2014 modification focused!!!


ComputePvalue(){

# Tian R.<tianremiATgmail.com>
# Sep 26, 2014
# weigh by num of patients
# filter remaining "common variants"

	tempdir=$1
	cutoff=0.05

#"all"$cancer"Cancertumor.filtered"

	paste $tempdir"/""all"$cancer"Cancergermline.filtered.snvfreq" $tempdir"/""all"$cancer"Cancertumor.filtered.snvfreq" | \
awk -v N=$NsampleLeft -v CUT=$cutoff '{if ($1 == $3 && ($2/N) <= CUT) print $1"\t"$4"\t"N"\t"($2/N)}' >$tempdir"/"$cancer".bloc"

	Rscript /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/bionorm.R $tempdir"/"$cancer".bloc"
	cat $tempdir"/"$cancer".bloc.pval" | sort -k5 -n -r > $tempdir"/"$cancer".bloc.Pval.model"
	rm $tempdir"/"$cancer".bloc" $tempdir"/"$cancer".bloc.pval"
	#echo $cancer "has "$total" patients."
}




CrossVD(){

	modelFile=$1
	tempdir=$2	

	cat $tempdir"/testPatientsList.txt" | while IFS=\t read testPatient
		do
 			#echo $testPatient
			cat $tumorSample | grep -w $testPatient | cut -f1 >$tempdir"/"$tumorSample"_"$testPatient".snvs"
			cat $germlineSample | grep -w $testPatient | cut -f1 >$tempdir"/"$germlineSample"_"$testPatient".snvs"


		done

	}




#-------------------------------------------------------------------------------

run(){

mkdir $tempdir
module load mcl

RunmclAllChrs "all"$cancer"Cancertumor.filtered" $dist "all"$cancer"Cancergermline.filtered" $tempdir

ComputePvalue $tempdir

CrossVD $tempdir"/"$cancer".bloc.Pval.model" $tempdir



}




if [ $# == 3 ]; then
	startTime=$(date +"%T")
	echo "Current time : $startTime"
	echo "Check args:"
	echo "Arg1 is "$1
	echo "Arg2 is "$2
	echo "Arg3 is "$3
	echo "---------------------------------------------------------------------------"
#cancer="colon"
#folder="/projects/common/tcga/coad/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/data/"
### require mcl
	
	run

	startTime=$(date +"%T")
	echo "Current time : $startTime"
	echo "###################################################################"

else
	echo "Error, we need 3 arguments here!"
fi

#X_11
#1_222

#A model file
#ChrXXX:111111,111111 Ntumor Ntotal backgroundRate CRscore
 
#Chr17:7578177,7578177,7578190,7578208,7578211,7578212,7578217,7578235,7578253,7578257,7578263,7578280,7578369,7578370,7578388,7578394,7578403,7578406,7578407,7578442,7578445,7578449,7578454,7578455,7578457,7578461,7578503,7578526,7578550	52	220	0	120.63978567313
#Chr12:25398214,25398214,25398218,25398226,25398255,25398262,25398281,25398282,25398284,25398285	75	220	0.00454545	115.906452155921
