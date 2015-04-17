# Tian R. <tianremi@gmail.com>
# Sep 22, 2014
# Get the SNV freq in 1000 G and keep those common variants
# So that in the patients genomes, any locus coincide with those common variants should be removed

FilterCommonVar(){
	
	inputfile=$1
	cutoff=$2
	#0.5%
	
	zcat -f $inputfile |cut -f1-8 | grep -v "#" | sed "s/;/\t/g"| awk '{for (i=1;i<=NF;i++){if(substr($i,1,3) == "AF=") print $1"\t"$2"\t"$4"\t"$5"\t"$i}}' \
| sed "s/AF=//g" |awk -v CUTOFF=$cutoff '{if ($5 >= CUTOFF) print $0}'

	}



Run (){
path=$1
out=$2

touch out

for file in ALL.chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz	ALL.chr12.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr14.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr15.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr18.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr19.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr4.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr7.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr8.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ALL.chrMT.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#ALL.chrY.genome_strip_hq.20101123.svs.low_coverage.genotypes.vcf.gz ALL.chrY.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz

do
	echo $path$file
	echo `wc -l $path$file`
	
	FilterCommonVar $path$file 0.005 >> out

done
}




#path="/projects/common/1kgenomes/data/"


#FilterCommonVar ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 0.005

Run "/projects/common/1kgenomes/data/" "/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1000G.0.005.commonVariants"
