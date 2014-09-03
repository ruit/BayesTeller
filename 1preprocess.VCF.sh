#Tian R.
# tianremi@gmail.com
# Aug. 20, 2014
# cat 220 colon cancer exome seq vcf file, sort numerically by chrID, by chrom location
# Aug. 21, 2014 Paper, Paper, Paper is everything!!!

FilterSNV(){

#filter by quality, by length, by 1 to N, by germline and somatic
	
	filefolder=$1
	cancer=$2
	
	for vcf in `ls $filefolder | grep vcf.gz`
	do

		echo $vcf
		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)!=substr($11, 1, 3)) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".somatic"
		cat $vcf".somatic" | cut -f1,2 | sort | uniq > $vcf".somatic.coord"
		

		zcat $filefolder$vcf | grep -v "#" | awk '{if ($7=="PASS" && substr($10, 1, 3)==substr($11, 1, 3)) print $1"\t"$2"\t"$4"\t"$5}' > $vcf".germline" 
		cat $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".germline.coord"
		
		
		cat $vcf".somatic" $vcf".germline" | cut -f1,2 | sort | uniq > $vcf".tumor.coord"

	done

	cat *.somatic.coord | sort -k1,1 -n -k2,2 -n > $cancer".somatic.stacked.sorted"
	cat *.germline.coord | sort -k1,1 -n -k2,2 -n >$cancer".germline.stacked.sorted" 

	cat *.tumor.coord | sort -k1,1 -n -k2,2 -n >$cancer".tumor.stacked.sorted"
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


cancer="colon"
folder="/projects/common/tcga/coad/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/data/"

#FilterSNV $folder $cancer 


for i in 1 2 3 4 5 6 7 8 9 10 \
11 12 13 14 15 16 17 18 19 20 21 22 X
do
	#CountSNV $i $cancer".somatic.stacked.sorted"
	CountSNV $i $cancer".tumor.stacked.sorted"
	CountSNV $i $cancer".germline.stacked.sorted"
done



