#!/bin/sh -login
# Dec 23, 2014
# Jan 13, 2015, All zang and lei jobs
# Jan 14, 2015
# Jan 19, 2015 Grasp every hour, day!


path="/scratch/tianr/input/"
#path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
dir=$1
sortingBy=$2
#sorting of the SNVs (models) by diff methods

array_alpha=(1 5 10 50 100 1000 5000 10000 100000)
suffix=$3

modelSort(){

	modelfile=$1 #".gz"
	patForModeling=$2	
	#total=`cat testPatList4B.tab | wc -l` 
	total=`cat $patForModeling | wc -l`
	sortingBy=$3
		
	zcat -f $modelfile | awk '{print ($3+0.5)/($2+0.5)}' > $modelfile".temp1"$sortingBy # ratio tumor N / germ N 
	zcat -f $modelfile | awk -v TOT=$total '{print $1"\t"$3"\t"TOT"\t"($2/TOT)}' > $modelfile".temp2"$sortingBy
	
		
	Rscript /home/tianr/1Projects/1SNVblocks/snvBayesTellerSinceMidDec2014/binorm2.R $modelfile".temp2"$sortingBy
 
	cat $modelfile".temp2"$sortingBy".pval" | cut -f5 > $modelfile".temp3"$sortingBy # Binormial p
	zcat -f $modelfile | awk '{print ($3-$2)}' > $modelfile".temp4"$sortingBy # tumor N - germ 
	
	zcat -f $modelfile > $modelfile".temp0"$sortingBy

	paste $modelfile".temp0"$sortingBy $modelfile".temp1"$sortingBy \
$modelfile".temp3"$sortingBy $modelfile".temp4"$sortingBy | gzip > $modelfile$sortingBy".all.gz"
	rm $modelfile".temp"*$sortingBy*
	}





moreCutoff(){

	dir=$1
	sortingBy=$2
	array_alpha=$3
	modelfile=$4 #modelA.sorted.gz.diff.all.gz"
	patForModeling=$5 #testPatList4B.tab
	testSuffix=$6 #snvsa.gz
	outSuffix=$7 #lode

	cd $path$dir
	total=`cat $patForModeling | wc -l`  # for model A

	if [ $sortingBy == "ratio" ]; then

		zcat -f $modelfile | sort -k4 -n -r | cut -f1-3 | gzip > $modelfile"_RO.gz"	
			
		for i in ${array_alpha[@]}
		do
			zcat -f $modelfile"_RO.gz" | head -n$i |gzip >$modelfile"_RO""_"$i".gz"	
		done

		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV4.py" $modelfile"_RO""_"$i".gz" $file $total $i "take" > $file"_"$i$outSuffix  			
			done
		done



	elif [ $sortingBy == "binorm" ]; then
		zcat -f $modelfile | sort -k5 -n -r | cut -f1-3 | gzip > $modelfile"_BN.gz"
		for i in ${array_alpha[@]}
		do
			zcat -f $modelfile"_BN.gz" | head -n$i | gzip >$modelfile"_BN""_"$i".gz"	
		done

		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV4.py" $modelfile"_BN""_"$i".gz" \
$file $total $i "take" > $file"_"$i$outSuffix  
			done
		done



	elif [ $sortingBy == "diff" ]; then
		zcat -f $modelfile | sort -k6 -n -r | cut -f1-3| gzip > $modelfile"_DF.gz"
		for i in ${array_alpha[@]}
		do
			zcat -f $modelfile"_DF.gz" | head -n$i |gzip >$modelfile"_DF""_"$i".gz"	
		done
	
		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV4.py" $modelfile"_DF""_"$i".gz" \
$file $total $i "take" > $file"_"$i$outSuffix  
			done
		done


	fi
	
	}


cd $path$dir

modelSort modelA.sorted.gz testPatList4B.tab $sortingBy
modelSort modelB.sorted.gz testPatList4A.tab $sortingBy




moreCutoff $dir $sortingBy $array_alpha  modelA.sorted.gz$sortingBy".all.gz" testPatList4B.tab .snvsa.gz $suffix
moreCutoff $dir $sortingBy $array_alpha  modelB.sorted.gz$sortingBy".all.gz" testPatList4A.tab .snvsb.gz $suffix

#".lodn"

