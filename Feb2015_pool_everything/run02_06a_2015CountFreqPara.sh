#!/bin/bash
#Feb 6, 2015
#TianR.

#Given a big file, hundreds of MBs, split into smaller parts, run in parallel

inputfile=$1

line=$2 #100000

#panCancer.germ.gz.sorted.org.gz

tempdir=temp_$RANDOM

mkdir $tempdir


split -d -l $line $inputfile "./"$tempdir"/"$inputfile".part_"

for file in `ls "./"$tempdir"/"$inputfile".part_"*`

do
	echo $file

	qsub -j y -cwd -b y "python countFreq3InOne_Feb2_2015.py "./"$tempdir"/"$file  "./"$tempdir"/"$file".out" . []"

done







