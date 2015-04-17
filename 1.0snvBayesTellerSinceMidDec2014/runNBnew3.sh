#!/bin/sh -login
# Dec 23, 2014

path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
dir=$1

predMassive(){

dir=$1
cd $path$dir
rm *".lod" 
rm *".lodr"
rm *".lode"

total=`cat testPatList4B.tab | wc -l`  # for model A

echo $total



###refiltering based on SNV coverage, EC
for i in 0.05 0.01 0.005 0.001 #r
do
	$path"generateSNVsetalpha.sh" modelA.sorted $i > "modelA.sorted_"$i  
	$path"generateSNVsetalpha.sh" modelB.sorted $i > "modelB.sorted_"$i  
done



#half
for file in `ls *snvsa`
do
	for i in 0.05 0.01 0.005 0.001 #refiltering alpha, Jan9, 2015
	
	do
		python $path"NBayessingleSNV3.py" "modelA.sorted_"$i $file $total `cat "modelA.sorted_"$i | wc -l` "take" > $file"_"$i".lode"  
	done
done




total=`cat testPatList4A.tab | wc -l`  # for model B

echo $total

#another half
for file in `ls *snvsb`
do
	for i in 0.05 0.01 0.005 0.001
	do
		python $path"NBayessingleSNV3.py" "modelB.sorted_"$i $file $total `cat "modelB.sorted_"$i | wc -l` "take" > $file"_"$i".lode"  

	done
done


}



predMassive $dir
