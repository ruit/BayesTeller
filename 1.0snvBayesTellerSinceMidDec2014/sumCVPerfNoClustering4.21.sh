#!/bin/sh -login
#TianR. Dec 8, 2014
#@Dec 9, 2014
#@Dec 17, 2014
#@Dec 23, 2014 model should be sorted first!!!
#topN=1000
#@Jan 10, 2015
#Jan 15, 2015


pkgPath="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"
array_alpha=(1 5 10 50 100 1000 5000 10000 100000)
#suffix=".lodn"
suffix=$1

sumTenFold(){
	dir=$1
	suffix=$2
	array_alpha=$3

	pseudo=-0.001
	#this is for NA, where for a given sample no SNV hit in the model

	cd $dir
 	for i in ${array_alpha[@]}
	do

		cat *"tumor_"*"_"$i$suffix | sed "s/NA/$pseudo/g" |awk '{print "t\t"$1}' > "tumor_"$suffix"."$i  

		cat *"germ_"*"_"$i$suffix | sed "s/NA/$pseudo/g" |awk '{print "g\t"$1}' > "germ_"$suffix"."$i

	
		cat "tumor_"$suffix"."$i "germ_"$suffix"."$i >../$i"_"$RANDOM"_"$suffix".pred"
		rm "tumor_"$suffix"."$i "germ_"$suffix"."$i		
	done
	cd ..
	}	


old (){
for dir in `ls  | grep CV2 | grep output`
do
	sumTenFold $dir $suffix $array_alpha
done
}



what(){

array_alpha=$1

for i in ${array_alpha[@]}
do
	cat $i"_"*$suffix".pred" | sort -k2 -n -r >$i"_"$suffix".sum"
	#python $pkgPath"ROC.py" $i"_"$suffix".sum" >$i"_"$suffix".avgPerf"
#	rm $i".sum" 
done		
}


#sumTenFold $dir $suffix $array_alpha
old
what $array_alpha
