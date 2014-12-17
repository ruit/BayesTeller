#!/bin/sh -login
# Tian R. <tianremiATgmail.com>
# Oct 30, 2014
# Nov 5, 2014
# Sample some patients simultaneously from germ /tumor for cross validation.
# Dec 3, always use qsub, always full path, always in parallel
# Dec 12, for panCancer
# Dec 16, 2014


path=$1
cd $path

input=$2
input2=$3
totalPatientList=$4




crossVD(){

# Tian R.<tianremiATgmail.com>
# Sep 26, 2014
# weigh by num of patients
# filter remaining "common variants"
	input=$1
	input2=$2
	tempdir="cv2_"$RANDOM
	totalPatientList=$3
	cutoff=0.05


	#Dec 16, 2014
	#mkdir tempdir and cat all *raw files	
	mkdir $tempdir 


	for chr in X 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 \
21 22 
	do
		cat $input"_chr"$chr".raw" | awk -v CHR=$chr '{print "Chr"CHR":"$0}' >$tempdir"/"$input"_chr"$chr".raw2"
		cat $input2"_chr"$chr".raw" | awk -v CHR=$chr '{print "Chr"CHR":"$0}' >$tempdir"/"$input2"_chr"$chr".raw2"
	done
	
	cat $tempdir"/"$input"_chr"*".raw2" > $tempdir"/"$input".raw"
	cat $tempdir"/"$input2"_chr"*".raw2" > $tempdir"/"$input2".raw"

	python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/CrossValidationII.py \
$tempdir"/"$input".raw" $tempdir"/"$input2".raw" $totalPatientList $tempdir		


	#remove more ChrX, Chr11 labeling
	cat $tempdir"/"$input".raw.outA"  | sed "s/Chr[0-9X]*://2g" > $tempdir"/"$input".snvfreqA"
	cat $tempdir"/"$input2".raw.outA" | sed "s/Chr[0-9X]*://2g" > $tempdir"/"$input2".snvfreqA"
	
	cat $tempdir"/"$input".raw.outB"  | sed "s/Chr[0-9X]*://2g" > $tempdir"/"$input".snvfreqB"
	cat $tempdir"/"$input2".raw.outB" | sed "s/Chr[0-9X]*://2g" > $tempdir"/"$input2".snvfreqB"

	#rm $tempdir"/"*raw*

	NsampleLeft=`cat $tempdir"/testPatList4B.tab" | wc -l`
#"all"$cancer"Cancertumor.filtered"

	paste $tempdir"/"$input".snvfreqA" $tempdir"/"$input2".snvfreqA" | \
awk -v N=$NsampleLeft -v CUT=$cutoff '{if ($1 == $3 && ($2/N) <= CUT) print $1"\t"$4"\t"N"\t"($2/N)}' >$tempdir"/blocA"

	Rscript /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/bionorm.R $tempdir"/blocA"
	cat $tempdir"/blocA.pval" | sort -k5 -n -r > $tempdir"/blocA.Pval.model"
	rm $tempdir"/blocA" $tempdir"/blocA.pval"


	NsampleLeft=`cat $tempdir"/testPatList4A.tab" | wc -l`
#"all"$cancer"Cancertumor.filtered"

	paste $tempdir"/"$input".snvfreqB" $tempdir"/"$input2".snvfreqB" | \
awk -v N=$NsampleLeft -v CUT=$cutoff '{if ($1 == $3 && ($2/N) <= CUT) print $1"\t"$4"\t"N"\t"($2/N)}' >$tempdir"/blocB"

	Rscript /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/bionorm.R $tempdir"/blocB"
	cat $tempdir"/blocB.pval" | sort -k5 -n -r > $tempdir"/blocB.Pval.model"
	rm $tempdir"/blocB" $tempdir"/blocB.pval"


	cat $tempdir"/testPatList4A.tab" | while IFS=\t read testPatient
		do
 			#echo $testPatient
			cat ../$input2 | grep -w $testPatient | cut -f1 >$tempdir"/"$input2"_"$testPatient".snvsa"
			cat ../$input | grep -w $testPatient | cut -f1 >$tempdir"/"$input"_"$testPatient".snvsa"


		done

	cat $tempdir"/testPatList4B.tab" | while IFS=\t read testPatient
		do
 			#echo $testPatient
			cat ../$input2 | grep -w $testPatient | cut -f1 >$tempdir"/"$input2"_"$testPatient".snvsb"
			cat ../$input | grep -w $testPatient | cut -f1 >$tempdir"/"$input"_"$testPatient".snvsb"

		done

	}






if [ $# == 4 ]; then
	startTime=$(date +"%T")
	echo "Current time : $startTime"
	echo "Check args:"
	echo $1
	echo $2
	echo $3
	echo $4

	crossVD $input $input2 $totalPatientList

	echo "---------------------------------------------------------------------------"

	startTime=$(date +"%T")
	echo "Current time : $startTime"
	echo "###################################################################"

else
	echo "Error, we need 1 arguments here!"
fi

#X_11
#1_222

#A model file
#ChrXXX:111111,111111 Ntumor Ntotal backgroundRate CRscore
 
#Chr17:7578177,7578177,7578190,7578208,7578211,7578212,7578217,7578235,7578253,7578257,7578263,7578280,7578369,7578370,7578388,7578394,7578403,7578406,7578407,7578442,7578445,7578449,7578454,7578455,7578457,7578461,7578503,7578526,7578550	52	220	0	120.63978567313
#Chr12:25398214,25398214,25398218,25398226,25398255,25398262,25398281,25398282,25398284,25398285	75	220	0.00454545	115.906452155921
