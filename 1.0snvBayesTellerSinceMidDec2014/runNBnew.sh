#!/bin/sh -login
# Dec 23, 2014

path="/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/"
dir=$1

predMassive(){

dir=$1
cd $path$dir
rm *".lod" 
rm *".lodr"

total=`cat testPatList4B.tab | wc -l`  # for model A

echo $total

for file in `ls *snvsa`
do
	for i in 0 1 5 10 50 100 900 950 995
	
	do
		python $path"NBayessingleSNV2.py" modelA.sorted $file $total $i "remove" > $file"_"$i".lodr"  
	done
done




total=`cat testPatList4A.tab | wc -l`  # for model B

echo $total


for file in `ls *snvsb`
do
	for i in 0 1 5 10 50 100 900 995
	do
	
		python $path"NBayessingleSNV2.py" modelB.sorted $file $total $i "remove" > $file"_"$i".lodr"  
	done
done


}



predMassive $dir
