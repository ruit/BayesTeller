#!/bin/sh -login
# Tian R. Dec 8, 2014

#generate ROC

cancer=$1
topN=$2

pkgPath="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"



echo "All tumor samples for test is "`ls "all"$cancer*tumor*.snvs | wc -l`
echo "All normal samples for test is "`ls "all"$cancer*germline*.snvs | wc -l`



for file in `ls "all"$cancer*.snvs`
do
	python $pkgPath"NaiveBayesSNV.py" $cancer".bloc.Pval.model" $file $topN > $file".lod"
	echo $file
done


cat "all"$cancer*"tumor"*".lod" | awk '{print "t\t"$0}' > tumor  

cat "all"$cancer*"germline"*".lod" | awk '{print "g\t"$0}' > germ

cat tumor germ | sort -k2 -n -r > tumor.germ

python $pkgPath"ROC.py" tumor.germ > $cancer"_"$RANDOM".perf"
