#!/bin/bash
# Tian R. tianremiATgmail.com, Jan 22, 2015
# Pick real SNVs, keep <=0.5%, aggregate 1092 genomes

pathToIndividualFiles="/projects/common/1kgenomes/data/individual/"
#HG00096.vcf.gz


snvFreqFile="/projects/common/1kgenomes/data/individual/1kg_all.fsnvs"
tempdir="./Jan22_2015/"


#0.5% MAF cutoff
#cat $snvFreqFile| awk '{if ($3<=0.005) print $1}' | sort | uniq |gzip > $tempdir"1kg_all.fsnvs.rare.gz"


#pathToIndividualFiles=$1 #
#tempdir=$2 #
#inFile=$3
#eachPersonFilter(){
#	pathToIndividualFiles=$1 #
#	tempdir=$2 #
#	inFile=$3	

#	zcat $pathToIndividualFiles$inFile | awk '{if ($4!="0|0" && length($2)==1 && length($3)==1) print $0}' \
#| cut -f1 | sort | uniq | gzip >$tempdir$inFile".cleaned.gz"
#	zcat $tempdir"1kg_all.fsnvs.rare.gz" $tempdir$inFile".cleaned.gz" \
#| sort | uniq -d | gzip >$tempdir$inFile".cleaned.rare.gz" 	
#	rm $tempdir$inFile".cleaned.gz"	
#	}

runAllPersons (){
zcat 1092individuls.txt.gz |while IFS="\t" read inFile

do
	echo $inFile
	#qsub -j y -cwd -o "/scratch/tianr/input/queue/"$inFile".stdout" -b y "./eachPersonFilter.sh $pathToIndividualFiles $tempdir $inFile"
	echo "./eachPersonFilter.sh "$pathToIndividualFiles" "$tempdir" "$inFile >$inFile".q.sh"
	chmod +x $inFile".q.sh"
	qsub -j y -cwd -o ./trash/$inFile".stdout" -b y "./"$inFile".q.sh"
done
}


runAllPersons
