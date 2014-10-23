#Tian R.
# tianremi@gmail.com
# Aug. 20, 2014
# cat 220 colon cancer exome seq vcf file, sort numerically by chrID, by chrom location
# Aug. 21, 2014 Paper, Paper, Paper is everything!!!
#
# Sep 5, 2014, run time 20min for colon cancer

# Sep 15, 2014, Major bugs found in the germline SNVfreq !!!!!!!!!!
# germline SNV block should encompass somatic SNVs!!!
# Sep 23, 2014
# Oct 21, 2014
#---------------------------------------------------------
#the name of cancer type
cancer=$1

#the path of folder where you have *vcf.gz
folder=$2

#how to filter low quality SNVs
filter=$3

#distance for mcl clustering
dist=$4

#SNV freq in 1000 Genomes
### common variants to be removed from cance genome research???
kgenomeVar="/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2"

# Tian R. <tianremi@gmail.com>
# Sep 12, 2014


#for Broad Institute VCF files, but I assume can be used in more general terms

FilterSNVbyMoreParas(){

#filter by quality, by length, by 1 to N, by germline and somatic
	
	filefolder=$1
	cancer=$2
	
	#base quality
	bq=30
	
	#ratio of alt reads
	fa=0.05
	
	#seq depth
	dp=7
	
	mlen=2	

	for vcf in `ls $filefolder | grep vcf.gz`
	do

		
		#10 is normal
		#11 is tumor
		
		echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)!=substr($11, 1, 3)) print $0}' > $vcf".somatic.temp"
		
		cat $vcf".somatic.temp" | awk -v BQ=$bq -v FA=$fa -v DP=$dp -v MLH=$mlen '{split($11, A, ":" ); \
if (A[6]>=BQ && A[4]>=FA && A[3]>=DP && length($4) <= MLH && length($5) <= MLH ) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $vcf".somatic"
		# look at col 11
		cat $vcf".somatic" | cut -f1,2 | sort | uniq > $vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)==substr($11, 1, 3)) print $0}' > $vcf".germline.temp" 
		# look at col 10
		cat $vcf".germline.temp" | awk -v BQ=$bq -v FA=$fa -v DP=$dp -v MLH=$mlen '{split($10, A, ":" ); \
if (A[6]>=BQ && A[4]>=FA && A[3]>=DP && length($4) <= MLH && length($5) <= MLH) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $vcf".germline"
		cat $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".germline.coord"
		
		
		cat $vcf".somatic" $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".tumor.coord"

		rm $vcf".germline.temp"
		rm $vcf".somatic.temp"
	done

#	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
#	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

#	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"

	}





FilterSNVbyPASS(){

#filter by quality, by length, by 1 to N, by germline and somatic
	
	filefolder=$1
	cancer=$2
	mlen=2
	
	for vcf in `ls $filefolder | grep vcf.gz`
	do

		#echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)!=substr($11, 1, 3) && length($4) <= $mlen && length($5) <= $mlen) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".somatic"
		cat $vcf".somatic" | cut -f1,2 | sort | uniq > $vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)==substr($11, 1, 3) && length($4) <= $mlen && length($5) <= $mlen) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".germline" 
		cat $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".germline.coord"
		
		
		cat $vcf".somatic" $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".tumor.coord"

		rm $vcf".germline"
		rm $vcf".somatic"
	done

#	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
#	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

#	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"

	}



#!/bin/bash
#Tian R. <tianremi@gmail.com>
#Aug. 26, 2014
#
#Cat all SNV blocks and each patient data to count in each block this patient is mutated or not
#
#Aug. 27, 2014 Run time, cat all blocks 1hr; for each patient another 4 hrs
#tumor="TCGA-G4-6628-01A-11D-1835-10.vs-10.vcf.gz.tumor.coord"
# tumor file should be sort | uniq first




#CombineAllBlocks "/home/tianr/Projects/SNVblocks/run/block/" "allSNVblocks"
#Prelinary results must be checked up and down using commands, study the composition of a file. Gei rid of bugs!!! Aug26, 2014


LabelPatients(){
	
	patientfolder=$1
	typ=$2
	out=$3

	num=0

	rm $out
	
	for patient in `ls $patientfolder | grep $typ`
	do
		
		num=`expr $num + 1` 	
		#echo $num 
		cat $patientfolder$patient | awk -v ID=$num '{print $0 "\t" "Patient"ID}'> $patient"_"$typ".temp"
		
	done
	
	cat *"_"$typ".temp" | awk '{print $1"_"$2"\t"$3}' > $out
	# Oct 20, Do not sort here, as it will be sorted again

	rm *"_"$typ".temp"
	}




PreprocessVariants(){

mv *.coord $cancer"_run"

echo "#INFO: Labeling patients for tumor."
LabelPatients $cancer"_run/" "tumor" "all"$cancer"Cancertumor"

echo "#INFO: Labeling patients for germline."
LabelPatients $cancer"_run/" "germline" "all"$cancer"Cancergermline"

echo "#INFO: Lableing common variants from 1000 genomes"
cat $kgenomeVar | awk '{print $1"_"$2"\tcomVar"}' > comSNV


cat "all"$cancer"Cancertumor" comSNV| grep -v \< | sed "s/_/\t/g" | sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}' > "all"$cancer"Cancertumor2"
#rm "all"$cancer"Cancertumor"
mv "all"$cancer"Cancertumor2" "all"$cancer"Cancertumor"
python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/listSortedPairs.py "all"$cancer"Cancertumor" | grep -v comVar >"all"$cancer"Cancertumor.filtered"
echo "#INFO: common variants have been filtered out from tumor variants."


cat "all"$cancer"Cancertumor" | cut -f1 | sort | uniq > tumorSNV 
cat "all"$cancer"Cancergermline" | cut -f1 | sort | uniq > germSNV
cat tumorSNV germSNV | sort | uniq -u | awk '{print $1 "\tvoid"}' > somatic
#remove weird <M>, etc.
cat "all"$cancer"Cancergermline" somatic comSNV | grep -v \< | sed "s/_/\t/g" |sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}' > "all"$cancer"Cancergermline2"

rm tumorSNV germSNV  somatic 
mv "all"$cancer"Cancergermline2" "all"$cancer"Cancergermline"
python  /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/listSortedPairs.py "all"$cancer"Cancergermline" | grep -v comVar > "all"$cancer"Cancergermline.filtered"
#rm "all"$cancer"Cancergermline" "all"$cancer"Cancertumor"
echo "#INFO: common variants have been filtered out from germline variants."

rm comSNV

}





ExtractByChrS ()
	{
	input=$1
	chr=$2
	maxDist=$3
	
	cat $input |sed "s/_/\t/g" | awk -v CHR=$chr '{ if ($1 == CHR) print $0}' > $input"_chr"$chr
	
	python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/generateAllPairs.py $input"_chr"$chr | awk  -v MAX=$maxDist -v CHR=$chr \
'function abs(x){return ((x < 0.0) ? -x : x)} \
{if (abs($2-$1) <= MAX) print CHR "_"$1 "\t"CHR"_"$2"\t"abs($2-$1)}' > $input"_chr"$chr".pair"
	
	rm $input"_chr"$chr
	}



RunmclAllChrs () 
	{
	input=$1
	maxDist=$2

	for chr in 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22 X 
	do
		echo "Chr"$chr
		ExtractByChrS $input $chr $maxDist
		mcl $input"_chr"$chr".pair" --abc  -o $input"_chr"$chr".pair.cluster"
		
		#@Oct 20, 2014 Mon
		grep $chr"_" $input | sed s/$chr"_"//g > $input"_chr"$chr
		
		#check this for bugs
		cat $input"_chr"$chr".pair.cluster" | sed s/$chr"_"//g| awk '{min=$1;max=$1;for (i=1;i<=NF;i++){if($i < min){min=$i}; if ($i > max){max=$i}}; print min"\t"max}' > $input"_chr"$chr".pair.cluster.StartEnd"
		
		cat $input"_chr"$chr".pair.cluster.StartEnd" | cut -f1 | awk '{print $0"\tA_start"}' >  $input"_chr"$chr".pair.cluster.Start" 
		cat $input"_chr"$chr".pair.cluster.StartEnd" | cut -f2 | awk '{print $0"\tZ_end"}' >  $input"_chr"$chr".pair.cluster.End"
		cat  $input"_chr"$chr $input"_chr"$chr".pair.cluster.Start" $input"_chr"$chr".pair.cluster.End" | sort -k1,1 -n -k2,2 > $input"_chr"$chr".raw"
		python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/countSNVclusterFreq.py $input"_chr"$chr".raw" 
		
		cat $input"_chr"$chr".raw.out"| awk -v CHR=$chr '{print "Chr"CHR":"$0}' >$input"_chr"$chr".raw.out2"
		rm $input"_chr"$chr".pair"*
		rm $input"_chr"$chr
	done

	cat $input*.raw.out2  > $input".snvfreq"
	rm *raw*

	}
### Oct 21, 2014 modification focused!!!








ComputePvalue(){

# Tian R.<tianremiATgmail.com>
# Sep 26, 2014
# weigh by num of patients
# filter remaining "common variants"


cutoff=0.05



Nbrca=501
Ncolon=220
Nkirc=274
Nkirp=196
Nlgg=468
Nlihc=236
Nluad=625
Nprad=309



for cancer in colon kirc kirp lgg lihc luad prad
do
	if [ $cancer == "colon" ]; then
		total=$Ncolon
	elif [ $cancer == "kirc" ]; then
		total=$Nkirc

	elif [ $cancer == "kirp" ]; then
		total=$Nkirp

	elif [ $cancer == "lgg" ]; then
		total=$Nlgg

	elif [ $cancer == "lihc" ]; then
		total=$Nlihc

	elif [ $cancer == "luad" ]; then
		total=$Nluad

	elif [ $cancer == "prad" ]; then
		total=$Nprad

	fi

#"all"$cancer"Cancertumor.filtered"

	paste "all"$cancer"Cancergermline.filtered.snvfreq" "all"$cancer"Cancertumor.filtered.snvfreq" | \
awk -v N=$total -v CUT=$cutoff '{if ($1 == $3 && ($2/N) <= CUT) print $1"\t"$4"\t"N"\t"($2/N)}' >$cancer".bloc"

	Rscript /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/bionorm.R $cancer".bloc"
	cat $cancer".bloc.pval" | sort -k5 -n -r > $cancer".bloc.Pval"
	rm $cancer".bloc" $cancer".bloc.pval"
	#echo $cancer "has "$total" patients."

done

}


#-------------------------------------------------------------------------------

startTime=$(date +"%T")
echo "Current time : $startTime"
echo "Check args:"
echo "Arg1 is "$1
echo "Arg2 is "$2
echo "Arg3 is "$3
echo "Arg4 is "$4
echo "---------------------------------------------------------------------------"
#cancer="colon"
#folder="/projects/common/tcga/coad/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/data/"
mkdir $cancer"_run"

echo "#INFO: Filtering low quality SNVs from raw vcf files..."
echo "#INFO: Filter option is "$filter


#####space space space!!!!!

if [ "$filter" == "PASS" ];then
	echo "@Filter by PASS"
	FilterSNVbyPASS $folder $cancer
elif [ "$filter" == "other" ];then	
	echo "@Filter by more parameters"
	FilterSNVbyMoreParas $folder $cancer
fi




PreprocessVariants

### require mcl

module load mcl

RunmclAllChrs "all"$cancer"Cancertumor.filtered" $dist
RunmclAllChrs "all"$cancer"Cancergermline.filtered" $dist

rm *vcf.gz.*
rm -r $cancer"_run"

ComputePvalue


startTime=$(date +"%T")
echo "Current time : $startTime"
echo "###################################################################"

