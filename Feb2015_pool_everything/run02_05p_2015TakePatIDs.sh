#Feb 5, 2015
#TianR.


getID(){
for file in `ls run_*.sh.o*`

do
	echo $file
	cancer=`cat $file | sed -n 3p | awk '{print $3}'`

	num=`cat $file | grep TCGA |wc -l`
	
	cat $file |grep TCGA > $cancer"_"$num".sampleIDs"

done
}



# sort first by chr then by loci, get rid of <M> and others

mkdir tempFeb5

for input in `ls *.loci.gz`

do
	echo $input

	zcat $input | sed "s/_/\t/g" |grep -v \< \
|sort -k1,1 -k2,2 -n |awk '{print$1"_"$2"\t"$3}' | gzip > tempFeb5/$input".sorted.gz"

	python /home/tianr/1Projects/1SNVblocks/snvBayesTellerSinceMidDec2014/listSortedPairs2.py \
tempFeb5/$input".sorted.gz" "force" | gzip > tempFeb5/$input".sorted.org.gz"


done
