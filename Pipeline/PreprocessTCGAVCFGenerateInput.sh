#Tian R.
# tianremi@gmail.com
# Aug. 20, 2014
# cat 220 colon cancer exome seq vcf file, sort numerically by chrID, by chrom location
# Aug. 21, 2014 Paper, Paper, Paper is everything!!!
#
# Sep 5, 2014, run time 20min for colon cancer
#---------------------------------------------------------
#the name of cancer type
cancer=$1

#the path of folder where you have *vcf.gz
folder=$2


FilterSNV(){

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

echo "#INFO: filtering out low quatity SNVs from VCF files."
FilterSNV $folder $cancer 

mv *.coord $cancer"_run"

echo "#INFO: Labeling patients for tumor."
LabelPatients $cancer"_run/" "tumor" "all"$cancer"Cancertumor"

echo "#INFO: Labeling patients for germline."
LabelPatients $cancer"_run/" "germline" "all"$cancer"Cancergermline"

echo "#INFO: counting SNV block mutational frequencies for tumor."
python SNVBlockFreqTCGA.py "all"$cancer"Cancertumor" 500

echo "#INFO: counting SNV block mutational frequencies for germline."
python SNVBlockFreqTCGA.py "all"$cancer"Cancergermline" 500


Summary(){
	cat "all"$cancer"Cancertumor.snvfreq" "all"$cancer"Cancergermline.snvfreq" | cut -f1 | sort | uniq -c | awk '{if($1==1) print $2"\t0"}' >$cancer".somatic.snvfreq"

	cat "all"$cancer"Cancergermline.snvfreq" $cancer".somatic.snvfreq" | sort > temp1
	cat "all"$cancer"Cancertumor.snvfreq" | sort > temp2
	cut -f1 temp1 > temp1.cut1
	cut -f1 temp2 > temp2.cut1

	res=`diff temp1.cut1 temp2.cut1`
	if [ $? -eq 0 ]
	then
		
		echo "1st cols are the same."
		echo "#SNVblock	Nnormal	Ntumor" > $cancer".summary"
		cat temp1 temp2 | cut -f1-2,4	>>$cancer".summary"
	else
		echo "Different!"
	fi 

	#rm temp1*
	#rm temp2*
}

rm *out1
rm *out2
rm -r $cancer"_run"
#Summary

startTime=$(date +"%T")
echo "Current time : $startTime"
echo"#_____________________________________________"

