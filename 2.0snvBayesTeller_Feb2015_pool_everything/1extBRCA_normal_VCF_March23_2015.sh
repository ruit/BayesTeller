#!/bin/bash -login

#Mar 23, 2014
#Mar 2, 2015
#Normal must be a subset of Tumor SNVs
#Double check and then move on!
#April 13, 2015 Grab every hour, every day.
#normal(germline) should be always a subset of tumor mutations.





dst="/alto/projects/common/tcga/brca/genome.wustl.edu_BRCA.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.3.0/"

dir=`ls $dst | grep ".snv.1.vcf.gz"` 

#908 files


#col 9  format
#GT:DP:DP4:BQ:FA:VAQ:SS:FT
#GT:DP:DP4:BQ:FA:VAQ:SS:FT:GQ:MQ:AD
#GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS:FT
#GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS:FT:DP4

#col 12 is normal

FilterSNV4Tumor(){
	
	#filter by quality, by length, by 1 to N, by germline and somatic
	
	vcf=$1
	outdir=$2	
	#base quality
	bq=30
	#ratio of alt reads
	fa=0.05
	#seq depth
	dp=7
	mlen=1	
	#10 is normal
	#11 is tumor
		
	echo $vcf
	zcat $vcf | grep -v "#" | awk '{if (substr($10, 1, 3) !="./." && substr($10, 1, 3) !="0/0") print $0}' |\
awk -v BQ=$bq '{split($10, A, ":" ); split(A[4], AB, ","); if (AB[1] >=BQ || AB[2] >=BQ) {print $0}}' |\
awk -v DPTHS=$dp '{split($9, FMT, ":" ); split($10, A, ":"); if (FMT[2]=="DP") {xdp=A[2]}; if(FMT[3]=="DP") {xdp=A[3]}; if (xdp >= DPTHS) {print $0} }' |\
awk -v FRAL=$fa '{split($9, FMT, ":" ); split($10, A, ":"); if (FMT[5]=="FA") {xfa=A[5]}; if(FMT[7]=="FA") {xfa=A[7]}; if (xfa >= FRAL) {print $0} }' |\
awk -v LEN=$mlen '{if (length($4)==LEN && length($5)==LEN) print $0}' |gzip > $outdir"/"`basename $vcf`".FLT.gz" 
# March 11, 2015
	}

vcf=$1
outdir=$2

FilterSNV4Tumor $vcf $outdir


#for file in $dir 
#do
#	new=`basename $file`
#	echo $new


#	zcat $dst$new| grep -v \# | grep PASS | gzip >$new".flt2.gz"		
#done





