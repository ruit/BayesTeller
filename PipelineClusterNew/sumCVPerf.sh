#!/bin/sh -login
#TianR. Dec 8, 2014
#@Dec 9, 2014


#topN=1000

pkgPath="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"
cancer=$1
topN=$2

runNB(){
	cancer=$1
	topN=$2
	pseudo=-0.001
	#this is for NA, where for a given sample no SNV hit in the model


	echo "All tumor samples for test is "`ls "all"$cancer*tumor*.snvs | wc -l`
	echo "All normal samples for test is "`ls "all"$cancer*germline*.snvs | wc -l`



	for file in `ls "all"$cancer*.snvs`
	do
		python $pkgPath"NaiveBayesSNV.py" $cancer".bloc.Pval.model" $file $topN > $file".lod"
		echo $file
	done


	cat "all"$cancer*"tumor"*".lod" | sed "s/NA/$pseudo/g" |awk '{print "t\t"$0}' > tumor  

	cat "all"$cancer*"germline"*".lod" | sed "s/NA/$pseudo/g" |awk '{print "g\t"$0}' > germ
	}


sumTenFold(){
	
	cancer=$1
	topN=$2
	
	for eachFolder in `ls | grep output | grep $cancer`
	do
		distance=`echo $eachFolder | cut -d "_" -f3` 
		cd $eachFolder	
		runNB $cancer $topN
		cat tumor germ >../$cancer"_"$distance"_"$topN"_"$RANDOM".pred"
		
		cd ..
		cat $cancer"_"$distance"_"$topN"_"*".pred" | sort -k2 -n -r >$cancer"_"$distance"_"$topN".pred.sum"
		python $pkgPath"ROC.py" $cancer"_"$distance"_"$topN".pred.sum" > $cancer"_"$distance"_"$topN".pred.sum.avgPerf" 
		
	done
	}


#for cancer in acc blca colon kirc kirp lgg lihc luad paad pcpg prad skcm stad thca

if [ $# -eq 2 ]; then
	runNB $cancer $topN
	sumTenFold $cancer $topN
else
	echo $0" needs 2 args: <cancer><topN>!"
fi


