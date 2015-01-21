# Jan 12, EC modified
# Jan 21, 2015 RT


fpr=$1 # for example, 0.5%

#fileNames=$2

cutPoint(){
	inputPerf=$1
	fpr=$2
	acc=`cat $inputPerf | awk -v FPR=$fpr '{if ($1 - FPR >0) {print $1-FPR, $0} else {print FPR-$1, $0} }' | grep -v "^FPR" | sort -gk 1 |awk '{print $2,$3,$4,$5}' | head -n 1`
	echo $file $acc |sed s/.sum//g |sed -e 's/ /\t/g'
	}


#fileArr=(0.000001.sum 0.00001.sum 0.0001.sum 0.001.sum 0.005.sum 0.01.sum 0.05.sum)

#fileArr=(1.sum 5.sum 10.sum 50.sum 100.sum 1000.sum 10000.sum )

echo "NUM FPR TPR FDR AUC" |sed -e 's/ /\t/g'


ls *sum | sed "s/_/\t/g" | sort -k1 -n | sed "s/\t/_/g" >sumList 
#for file in ${fileArr[@]}

while read file;
do
	#Rscript ../perf2all.R $file	
	cutPoint $file".p" $fpr
done < sumList



# should be sorted


