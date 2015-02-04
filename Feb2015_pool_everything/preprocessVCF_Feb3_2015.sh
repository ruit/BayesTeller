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

# Feb 3, 2015 pool all vcf data together, without filtering out commeon variants.

#---------------------------------------------------------
#the name of cancer type
cancer=$1

#the path of folder where you have *vcf.gz
folder=$2

#how to filter low quality SNVs
filter=$3

#distance for mcl clustering
#dist=$4

#SNV freq in 1000 Genomes
### common variants to be removed from cance genome research???
#kgenomeVar="/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2"

# Tian R. <tianremi@gmail.com>
# Sep 12, 2014


#for Broad Institute VCF files, but I assume can be used in more general terms

FilterSNVbyMoreParas(){

#filter by quality, by length, by 1 to N, by germline and somatic
	
	filefolder=$1
	cancer=$2
	outdir=$3	
	#base quality
	bq=30
	
	#ratio of alt reads
	fa=0.05
	
	#seq depth
	dp=7
	
	mlen=1	

	for vcf in `ls $filefolder | grep vcf.gz`
	do

		
		#10 is normal
		#11 is tumor
		
		echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)!=substr($11, 1, 3)) print $0}' \
| awk -v BQ=$bq -v FA=$fa -v DP=$dp -v MLH=$mlen '{split($11, A, ":" ); \
if (A[6]>=BQ && A[4]>=FA && A[3]>=DP && length($4) <= MLH && length($5) <= MLH ) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $outdir"/"$vcf".somatic"
		# look at col 11
		cat $outdir"/"$vcf".somatic" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if (substr($10, 1, 3)==substr($11, 1, 3)) print $0}' \
| awk -v BQ=$bq -v FA=$fa -v DP=$dp -v MLH=$mlen '{split($10, A, ":" ); \
if (A[6]>=BQ && A[4]>=FA && A[3]>=DP && length($4) <= MLH && length($5) <= MLH) print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > $outdir"/"$vcf".germline"
		# look at col 10
		cat $outdir"/"$vcf".germline" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".germline.coord"
		
		
		cat $outdir"/"$vcf".somatic" $outdir"/"$vcf".germline" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".tumor.coord"

	done

#	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
#	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

#	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"

	}





FilterSNVbyPASS(){

#filter by quality, by length, by 1 to N, by germline and somatic
	
	filefolder=$1
	cancer=$2
	outdir=$3

	mlen=1
	
	for vcf in `ls $filefolder | grep vcf.gz`
	do

		#echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)!=substr($11, 1, 3) && length($4) <= $mlen && length($5) <= $mlen) print $1"\t"$2"\t"$4"\t"$5}' > $outdir"/"$vcf".somatic"
		cat $outdir"/"$vcf".somatic" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)==substr($11, 1, 3) && length($4) <= $mlen && length($5) <= $mlen) print $1"\t"$2"\t"$4"\t"$5}' > $outdir"/"$vcf".germline" 
		cat $outdir"/"$vcf".germline" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".germline.coord"
		
		
		cat $outdir"/"$vcf".somatic" $outdir"/"$vcf".germline" | cut -f1,2 | sort | uniq > $outdir"/"$vcf".tumor.coord"

		#rm $outdir"/"$vcf".germline"
		#rm $outdir"/"$vcf".somatic"
	done

#	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
#	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

#	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"

	}




LabelPatients(){
	
	patientfolder=$1
	typ=$2
	out=$3
	

	num=0

	
	for patient in `ls $patientfolder | grep $typ`
	do
		
		num=`expr $num + 1` 	
		#echo $num 
		cat $patientfolder$patient | awk -v ID=$num '{print $0 "\t" "Patient"ID}'> $patientfolder$patient"_"$typ".temp"
		
	done
	
	cat $patientfolder*"_"$typ".temp" | awk '{print $1"_"$2"\t"$3}' | gzip> $out
	# Oct 20, Do not sort here, as it will be sorted again

	rm $patientfolder*"_"$typ".temp"
	}






#-------------------------------------------------------------------------------

startTime=$(date +"%T")
echo "Current time : $startTime"
echo "Check args:"
echo "Arg1 is "$1
echo "Arg2 is "$2
echo "Arg3 is "$3


echo "---------------------------------------------------------------------------"
#cancer="colon"
#folder="/projects/common/tcga/coad/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/data/"
mkdir $cancer"_run"

echo "#INFO: Filtering low quality SNVs from raw vcf files..."
echo "#INFO: Filter option is "$filter


#####space space space!!!!!

if [ "$filter" == "PASS" ];then
	echo "@Filter by PASS"
	FilterSNVbyPASS $folder $cancer $cancer"_run"
elif [ "$filter" == "other" ];then	
	echo "@Filter by more parameters"
	FilterSNVbyMoreParas $folder $cancer $cancer"_run"
fi


echo "#INFO: Labeling patients for tumor."
LabelPatients $cancer"_run/" "tumor" "all_"$cancer"_tumor.loci.gz"

echo "#INFO: Labeling patients for germline."
LabelPatients $cancer"_run/" "germline" "all_"$cancer"_germline.loci.gz"

rm $cancer"_run/"*.coord
gzip $cancer"_run/"*



startTime=$(date +"%T")
echo "Current time : $startTime"
echo "###################################################################"

