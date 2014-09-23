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
#---------------------------------------------------------
#the name of cancer type
cancer=$1

#the path of folder where you have *vcf.gz
folder=$2

#how to filter low quality SNVs
filter=$3

#SNV freq in 1000 Genomes
kgenomeVar=$4
#/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified

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
	
	
	for vcf in `ls $filefolder | grep vcf.gz`
	do

		
		#10 is normal
		#11 is tumor
		
		echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)!=substr($11, 1, 3)) print $0}' > $vcf".somatic.temp"
		
		cat $vcf".somatic.temp" | awk -v BQ=$bq -v FA=$fa -v DP=$dp '{split($11, A, ":" ); if (A[6]>=BQ && A[4]>=FA && A[3]>=DP) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $vcf".somatic"
		# look at col 11
		cat $vcf".somatic" | cut -f1,2 | sort | uniq > $vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)==substr($11, 1, 3)) print $0}' > $vcf".germline.temp" 
		# look at col 10
		cat $vcf".germline.temp" | awk -v BQ=$bq -v FA=$fa -v DP=$dp '{split($10, A, ":" ); if (A[6]>=BQ && A[4]>=FA && A[3]>=DP) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $vcf".germline"
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
	
	for vcf in `ls $filefolder | grep vcf.gz`
	do

		#echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)!=substr($11, 1, 3)) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".somatic"
		cat $vcf".somatic" | cut -f1,2 | sort | uniq > $vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)==substr($11, 1, 3)) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".germline" 
		cat $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".germline.coord"
		
		
		cat $vcf".somatic" $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".tumor.coord"

		rm $vcf".germline"
		rm $vcf".somatic"
	done

#	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
#	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

#	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"

	}



#SNV blocks, starts from clon cancer data, 220 samples
#annotation free method 

CountSNV (){

	CHR=$1
	inputfile=$2
	###chrom name
	
	cat $inputfile | awk -v CHROMNAME=$CHR '{if ( $1==CHROMNAME ) print $2}' \
| sort | uniq -c | sort -k2 -n  > $inputfile"."$CHR".sum"


	}



something(){
for i in 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22 X
do
	#CountSNV $i $cancer".somatic.stacked.sorted"
	CountSNV $i $cancer".tumor.stacked.sorted"
	CountSNV $i $cancer".germline.stacked.sorted"
done
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
	
	cat *"_"$typ".temp" | sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}' > $out

	rm *"_"$typ".temp"
	}


startTime=$(date +"%T")
echo "Current time : $startTime"

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



Run(){

span=$1
excludeR=$2
repeat=$3


mv *.coord $cancer"_run"

echo "#INFO: Labeling patients for tumor."
LabelPatients $cancer"_run/" "tumor" "all"$cancer"Cancertumor"

echo "#INFO: Labeling patients for germline."
LabelPatients $cancer"_run/" "germline" "all"$cancer"Cancergermline"

echo "#INFO: Lableing common variants from 1000 genomes"
cat $kgenomeVar | awk '{print $1"_"$2"\tcomVar"}' > comSNV


cat "all"$cancer"Cancertumor" comSNV| grep -v \< | sed "s/_/\t/g" | sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}' > "all"$cancer"Cancertumor2"
rm "all"$cancer"Cancertumor"
mv "all"$cancer"Cancertumor2" "all"$cancer"Cancertumor"

cat "all"$cancer"Cancertumor" | cut -f1 | sort | uniq > tumorSNV 
cat "all"$cancer"Cancergermline" | cut -f1 | sort | uniq > germSNV
cat tumorSNV germSNV | sort | uniq -u | awk '{print $1 "\tvoid"}' > somatic
#remove weird <M>, etc.
cat "all"$cancer"Cancergermline" somatic comSNV | grep -v \< | sed "s/_/\t/g" |sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}' > "all"$cancer"Cancergermline2"

rm tumorSNV germSNV  somatic "all"$cancer"Cancergermline"
mv "all"$cancer"Cancergermline2" "all"$cancer"Cancergermline"

echo "#INFO: counting SNV block mutational frequencies for tumor."
python /home/tianr/1Projects/1SNVblocks/pipeline/SNVBlockFreqTCGA.py "all"$cancer"Cancertumor" $span $excludeR $repeat

echo "#INFO: counting SNV block mutational frequencies for germline."
python /home/tianr/1Projects/1SNVblocks/pipeline/SNVBlockFreqTCGA.py "all"$cancer"Cancergermline" $span $excludeR $repeat

rm comSNV

}

Run 500 0 1

#Tian R. <tianremi@gmail.com>
#Sep. 15, 2014 

#germline tumor 
#somatic SNV assign 0 in germline
Summary(){

	#This should be only valid for single SNV, not for SNV blocks!!!!!
	
	cancer=$1
	tumorfile=$2
	germfile=$3

	cat $tumorfile $germfile | cut -f1 | sort | uniq -c | awk '{if($1==1) print $2"\t0"}' >$cancer".somatic.snvfreq"

	cat $germfile $cancer".somatic.snvfreq" | sed "s/_/\t/g" | sort -k1,1 -k2,2 -n > temp1
	cat $tumorfile | sed "s/_/\t/g" | sort -k1,1 -k2,2 -n > temp2
	cut -f1,2 temp1 > temp1.cut1
	cut -f1,2 temp2 > temp2.cut1


	res=`diff temp1.cut1 temp2.cut1`
	if [ $? -eq 0 ]
	then
		
		echo "1st cols of the two files are the same."
		echo "#chr	pos	Nnormal	Ntumor" > $cancer".summary"
		paste temp1 temp2 | cut -f1-3,6	>>$cancer".summary"
	else
		echo "Different!"
	fi 

	rm temp1*
	rm temp2*
	rm *.somatic.snvfreq
}

#Summary luad allluadCancertumor.snvfreq.cn allluadCancergermline.snvfreq.cn

rm *out1
rm *out2
rm -r $cancer"_run"
#Summary

startTime=$(date +"%T")
echo "Current time : $startTime"
echo "#_____________________________________________"

