#!/bin/bash
#March 2, 2015
#March 19, 2015

#$1 is the input file
#$2 is the ratio
brcaTest="brcaMar20Test"
snvlist="snv_1.1.ALL"
sampleSize=1954

python bayesTellerFeb27_2015.py $brcaTest"/"$snvlist $1 $brcaTest/germ_freq.A $brcaTest/tumor_freq.A $sampleSize $2



#for file in `ls newFilteredVCFs/*FLT.gz`; do qsub -j y -cwd -b y "./runBayesTeller_Mar19_2015.sh $file 1.5"; done


#merge all files
#paste <(ls newFilteredVCFs/*.gz | cut -d "." -f4) <(cat newFilteredVCFs/*.all |cut -f1-3) <(cat newFilteredVCFs/*.gml |cut -f1-3) <(cat newFilteredVCFs/*.sml |cut -f1-3) >



#Mar 22, 2015
#for rt in 1.2 1.6 2.0; do paste <(ls newFilteredVCFs/*.gz | cut -d "." -f4) <(cat newFilteredVCFs/*_$rt".all" |cut -f1-3) <(cat newFilteredVCFs/*_$rt".gml"|cut -f1-3) <(cat newFilteredVCFs/*_$rt".sml" |cut -f1-3) > resMar22/r$rt"_pred"; done
