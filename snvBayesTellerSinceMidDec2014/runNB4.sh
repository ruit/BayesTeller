#!/bin/sh -login
# Dec 23, 2014
# Jan 13, 2015, All zang and lei jobs
# Jan 14, 2015


path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
dir=$1
sortingBy=$2
#sorting of the SNVs (models) by diff methods

#array_alpha=(0.000001) # Jan12, 2015, 10e-3, -6
#array_alpha=(0.05 0.01 0.005 0.001 0.0001 0.00001 0.000001)
array_alpha=(1 5 10 50 100 1000 5000 10000 100000)
suffix=$3

modelSort(){

	modelfile=$1
	patForModeling=$2	
	#total=`cat testPatList4B.tab | wc -l` 
	total=`cat $patForModeling | wc -l`
	sortingBy=$3
		
	cat $modelfile | awk '{print ($3+0.5)/($2+0.5)}' > $modelfile".temp1"$sortingBy # ratio tumor N / germ N 
	cat $modelfile | awk -v TOT=$total '{print $1"\t"$3"\t"TOT"\t"($2/TOT)}' > $modelfile".temp2"$sortingBy
	
		
	Rscript /home/tianr/1Projects/1SNVblocks/snvBayesTellerSinceMidDec2014/binorm2.R $modelfile".temp2"$sortingBy
 
	cat $modelfile".temp2"$sortingBy".pval" | cut -f5 > $modelfile".temp3"$sortingBy # Binormial p
	cat $modelfile | awk '{print ($3-$2)}' > $modelfile".temp4"$sortingBy # tumor N - germ N
	paste $modelfile $modelfile".temp1"$sortingBy $modelfile".temp3"$sortingBy $modelfile".temp4"$sortingBy > $modelfile".out"$sortingBy
	rm $modelfile".temp"*
	}





moreCutoff(){

	dir=$1
	sortingBy=$2
	array_alpha=$3
	modelfile=$4 #modelA.sorted.out
	patForModeling=$5 #testPatList4B.tab
	testSuffix=$6 #snvsa
	outSuffix=$7 #lode

	cd $path$dir
	total=`cat $patForModeling | wc -l`  # for model A

	if [ $sortingBy == "ratio" ]; then
		cat $modelfile | sort -k4 -n -r | cut -f1-3 > $modelfile"_RO"	
		
		for i in ${array_alpha[@]}
		do
			cat $modelfile"_RO" | head -n$i >$modelfile"_RO""_"$i	
		done

		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV3.py" $modelfile"_RO""_"$i $file $total $i "take" > $file"_"$i$outSuffix  			
			done
		done



	elif [ $sortingBy == "binorm" ]; then
		cat $modelfile | sort -k5 -n -r | cut -f1-3 > $modelfile"_BN"
		for i in ${array_alpha[@]}
		do
			cat $modelfile"_BN" | head -n$i >$modelfile"_BN""_"$i	
		done

		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV3.py" $modelfile"_BN""_"$i \
$file $total $i "take" > $file"_"$i$outSuffix  
			done
		done



	elif [ $sortingBy == "diff" ]; then
		cat $modelfile | sort -k6 -n -r | cut -f1-3 > $modelfile"_DF"
		for i in ${array_alpha[@]}
		do
			cat $modelfile"_DF" | head -n$i >$modelfile"_DF""_"$i	
		done
	
		for file in `ls *$testSuffix`
		do
			for i in ${array_alpha[@]} #refiltering alpha, Jan9, 2015
	
			do
				python $path"NBayessingleSNV3.py" $modelfile"_DF""_"$i \
$file $total $i "take" > $file"_"$i$outSuffix  
			done
		done


	fi
	
	}


cd $path$dir

modelSort modelA.sorted testPatList4B.tab $sortingBy
modelSort modelB.sorted testPatList4A.tab $sortingBy

moreCutoff $dir $sortingBy $array_alpha  modelA.sorted.out testPatList4B.tab .snvsa $suffix
moreCutoff $dir $sortingBy $array_alpha  modelB.sorted.out testPatList4A.tab .snvsb $suffix

#".lodn"

