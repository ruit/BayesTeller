#!/bin/bash
#Tian R. <tianremi@gmail.com>
#Aug. 26, 2014
#
#Cat all SNV blocks and each patient data to count in each block this patient is mutated or not
#
#Aug. 27, 2014 Run time, cat all blocks 1hr; for each patient another 4 hrs
#tumor="TCGA-G4-6628-01A-11D-1835-10.vs-10.vcf.gz.tumor.coord"
# tumor file should be sort | uniq first


CountEachPatientForAllSNVblocks(){

	
	allpatient=$1
	SNVblockfile=$2
	out=$3

	cat $allpatient $SNVblockfile | sort -k1,1 -n -k2,2 | awk '{print $0}' | awk '{print $1"_" $2"\t" $3}' > $allpatient".temp.out" 
	
	
	python summarizeSortedPairs2.py $allpatient".temp.out" | grep Patient | awk '{print $2}' | sed "s/,/\t/g" | > $allpatient".temp.out.aggre"

	cat $allpatient".temp.out.aggre" | awk '{for (i=1;i<=NF;i++) {if(substr($i, 1, 7)=="Patient"){j=j+1; patarray[j]=$i} ; \
if (substr($i, 1, 3)!="Pat"){m=m+1; blocarray[m]=$i}}; j=0; m=0; \
for (x=1; x<=length(blocarray); x++){for (y=1;y<=length(patarray);y++){print blocarray[x] "\t" patarray[y] }}}' | sort | uniq > $out

	}


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
		echo $num 
		cat $patientfolder$patient | awk -v ID=$num '{print $0 "\t" "Patient"ID}'> $patient"_"$typ".temp"
		
	done
	
	cat *"_"$typ".temp"  > $out

	#rm *"_"$typ".temp"
	}


#LabelPatients "/home/tianr/Projects/SNVblocks/run/coord/" "tumor" "allClonCancertumor"






