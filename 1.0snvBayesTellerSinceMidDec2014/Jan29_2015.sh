#!/bin/bash
#TianR., Jan29, 2015



splitFile=$1

targetFile="/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/Kgenome1092_rare_aggre_Jan26_2015.gz"


zcat -f $splitFile | while IFS="\t" read var

touch $splitFile".matched.gz"

do 
#parallel --gnu --pipe -j20 
	grep -m 1 -w $var <(zcat $targetFile) | gzip >> $splitFile".matched.gz"
done


 
