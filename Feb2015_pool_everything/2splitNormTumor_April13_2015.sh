#TianR.
#Apr 13, 2015
#properly filter mutations, split germline mutations from tumor mutation set


inputfile=$1

getGerm() {

	file=$1
	#$10 tumor
	#$12 normal

	zcat $file | awk '{split ($10, A, ":"); split ($12, B, ":"); \
if (A[1] == B[1] && A[1]=="0/1" && B[1] !="./.") \
{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"A[1]"\t"B[1]}}' |gzip > $file".germ.gz"

	}



getSomatic() {

	file=$1
	#$10 tumor
	#$12 normal

	zcat $file | awk '{split ($10, A, ":"); split ($12, B, ":"); \
if (A[1] != B[1] && B[1] !="./.") \
{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"A[1]"\t"B[1]}}' |gzip > $file".somatic.gz"

	}



getGerm $inputfile
getSomatic $inputfile
zcat $inputfile".germ.gz" $inputfile".somatic.gz" | gzip > $inputfile".tumor.gz" 
rm $inputfile".somatic.gz"
