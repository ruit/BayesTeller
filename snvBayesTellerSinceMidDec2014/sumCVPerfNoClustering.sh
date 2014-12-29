#!/bin/sh -login
#TianR. Dec 8, 2014
#@Dec 9, 2014
#@Dec 17, 2014
#@Dec 23, 2014 model should be sorted first!!!
#topN=1000

pkgPath="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"



sumTenFold(){
	dir=$1
	suffix=$2
	pseudo=-0.001
	#this is for NA, where for a given sample no SNV hit in the model

	cd $dir
 	for i in 0 1 5 10 50 100 900 950 995
	do

		cat *"tumor_"*"_"$i$suffix | sed "s/NA/$pseudo/g" |awk '{print "t\t"$0}' > "tumor_"$suffix"."$i  

		cat *"germ_"*"_"$i$suffix | sed "s/NA/$pseudo/g" |awk '{print "g\t"$0}' > "germ_"$suffix"."$i

	
		cat "tumor_"$suffix"."$i "germ_"$suffix"."$i >../$i"_"$RANDOM"_"$suffix".pred"
		rm "tumor_"$suffix"."$i "germ_"$suffix"."$i		
	done
	cd ..
	}	

for dir in `ls  | grep CV2 | grep output`
do
	sumTenFold $dir ".lodr"
done


what(){
for i in 0 1 5 10 50 100 900 950 995
do
	cat $i"_"*".pred" | sort -k2 -n -r >$i".sum"
	python $pkgPath"ROC.py" $i".sum" >$i".avgPerf"
#	rm $i".sum" 
done		
}

what
