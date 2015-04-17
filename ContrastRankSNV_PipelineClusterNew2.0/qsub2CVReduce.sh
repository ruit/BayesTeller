#!/bin/sh -login
#TianR. Dec 17, 2014 
#2hr

route=/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014


runTen2CV (){

path=$1
for k in `seq 1 10`
do

	qsub -j y -cwd -o /scratch/tianr/queue/$k"_"$path".stdout" -b y $route"/twoFoldCVReduceIII.sh $path "panCancer.2in1.germ" "panCancer.2in1.tumor" $route"/"panCancer.2in1_total_15734.txt" 
done
}



for dir in \
panCancerCV_d1000000_4168_Tue_Dec_16_14_15_15_CST_2014Map \
panCancerCV_d100000_19173_Tue_Dec_16_14_15_12_CST_2014Map \
panCancerCV_d10000_26527_Tue_Dec_16_14_15_09_CST_2014Map \
panCancerCV_d1_29621_Tue_Dec_16_14_14_59_CST_2014Map \
panCancerCV_d2000_8312_Tue_Dec_16_14_15_05_CST_2014Map \
panCancerCV_d500_5332_Tue_Dec_16_14_15_03_CST_2014Map \
panCancerCV_d50_30901_Tue_Dec_16_14_15_00_CST_2014Map

do
	runTen2CV $dir
done
