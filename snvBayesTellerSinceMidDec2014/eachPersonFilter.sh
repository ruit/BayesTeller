#!/bin/bash
#TianR. Jan22, 2015

pathToIndividualFiles=$1 #
tempdir=$2 #
inFile=$3
eachPersonFilter(){
	pathToIndividualFiles=$1 #
	tempdir=$2 #
	inFile=$3	

	zcat $pathToIndividualFiles$inFile | awk '{if ($4!="0|0" && length($2)==1 && length($3)==1) print $0}' \
| cut -f1 | sort | uniq | gzip >$tempdir$inFile".cleaned.gz"
	zcat $tempdir"1kg_all.fsnvs.rare.gz" $tempdir$inFile".cleaned.gz" \
| sort | uniq -d | gzip >$tempdir$inFile".cleaned.rare.gz" 	
	rm $tempdir$inFile".cleaned.gz"	
	}
eachPersonFilter $pathToIndividualFiles $tempdir $inFile
