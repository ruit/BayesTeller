# Tian R. <tianremi@gmail.com>
# Sep 23, 2014
# Get the SNV freq in 1000 G and keep those common variants
# So that in the patients genomes, any locus coincide with those common variants should be removed

FilterCommonVar(){
	
	inputfile=$1
	cutoff=$2
	#0.5%
	output=$3

	
	zcat -f $inputfile |cut -f1-8 | grep -v "#" | sed "s/;/\t/g"| awk '{for (i=1;i<=NF;i++){if(substr($i,1,3) == "AF=") print $1"\t"$2"\t"$4"\t"$5"\t"$i}}' \
| sed "s/AF=//g" |awk -v CUTOFF=$cutoff '{if ($5 >= CUTOFF && $5 <= (1- CUTOFF) && length($3) <=2 && length($4) <=2) print $0}' > $output

	}


FilterCommonVar $1 $2 $3
