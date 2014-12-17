#!/bin/sh -login
#TianR. Dec 8, 2014
#@Dec 9, 2014
#@Dec 17, 2014

#topN=1000

pkgPath="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"
topN=$1

runNB(){
	topN=$1
	pseudo=-0.001
	#this is for NA, where for a given sample no SNV hit in the model


	echo "All tumor samples for test is "`ls *tumor*.snvs* | wc -l`
	echo "All normal samples for test is "`ls *germ*.snvs* | wc -l`



	for file in `ls *.snvsa`
	do
		python $pkgPath"NaiveBayesSNV.py" "blocA.Pval.model" $file $topN > $file".lod"

		echo $file
	done


	#the other half	
	for file in `ls *.snvsb`
	do
		python $pkgPath"NaiveBayesSNV.py" "blocB.Pval.model" $file $topN > $file".lod"

		echo $file
	done



	cat *"tumor"*".lod" | sed "s/NA/$pseudo/g" |awk '{print "t\t"$0}' > tumor  

	cat *"germ"*".lod" | sed "s/NA/$pseudo/g" |awk '{print "g\t"$0}' > germ
	}


sumTenFold(){
	
	topN=$1
	
	for eachFolder in `ls | grep "panCancerCV_"`
	do
		distance=`echo $eachFolder | cut -d "_" -f2` 
		cd $eachFolder
		for eachCV in `ls | grep cv2_`
		do	
			cd $eachCV
			runNB $topN
			cat tumor germ >../$distance"_"$topN"_"$RANDOM".pred"
			cd ..
		
		done

		cat $distance"_"$topN"_"*".pred" | sort -k2 -n -r >$distance"_"$topN".pred.sum"
		rm *.pred
		python $pkgPath"ROC.py" $distance"_"$topN".pred.sum" > $distance"_"$topN".pred.sum.avgPerf" 
		
	done
	}


#for cancer in acc blca colon kirc kirp lgg lihc luad paad pcpg prad skcm stad thca

if [ $# -eq 1 ]; then
	#runNB $topN
	sumTenFold $topN
else
	echo $0" needs 1 args: <topN>!"
fi


