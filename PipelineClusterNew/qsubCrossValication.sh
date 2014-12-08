#!/bin/sh -login
#use for looping submit many jobs
#Dec 3, 2014

module load mcl

path="/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/"

#cancer=$1
#dist=$2
#nTest=$3 # number of patients masked


declare -a array_NumSamples
declare -a array_cancerTypes

array_NumSamples=(501 220 274 196 468 236 625 309 91 129 95 186 371 280 429)
array_cancerTypes=(brca colon kirc kirp lgg lihc luad prad acc blca paad pcpg skcm stad thca)

tLen=${#array_NumSamples[@]}



subCVJobs (){

for (( i=0; i<${tLen}; i++));
do
	total=${array_NumSamples[$i]} 
	cancer=${array_cancerTypes[$i]}

	echo $cancer " is " $total


	if [ $cancer == "brca" ] #caution space!
	then
	
		echo "do nothing"
	#elif [ $cancer == "colon" ]
	else
		qsub -j y -cwd -q all.q -o /scratch/tianr/queue/$cancer.stdout -e /scratch/tianr/queue/$cancer.stderr -b y \
"/home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/twoFoldCV.sh $cancer 500 $(($total / 2))"
	fi
done

}



#subCVJobs



for ((j=0; j<10; j++)); #10 time 2fold CV
do
	echo $j
	subCVJobs
#	sleep $((1*60))
done
