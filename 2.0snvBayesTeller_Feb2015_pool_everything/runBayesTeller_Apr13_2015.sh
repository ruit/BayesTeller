#!/bin/bash
#March 2, 2015
#March 19, 2015
#Apr 13, 2015

#$1 is the input file
#$2 is the ratio
brcaTest="brcaMar20Test"
snvlist="snv_1.1.ALL"
sampleSize=1954

python bayesTellerFeb27_2015.py $brcaTest"/"$snvlist $1 $brcaTest/germ_freq.A $brcaTest/tumor_freq.A $sampleSize $2



#Apr 13
#for file in `ls newFilteredVCFs/*FLT.gz.germ.gz`; do qsub -j y -cwd -b y "./runBayesTeller_Apr13_2015.sh $file 1.5"; done
#for file in `ls newFilteredVCFs/*FLT.gz.tumor.gz`; do qsub -j y -cwd -b y "./runBayesTeller_Apr13_2015.sh $file 1.5"; done


#Mar 22, 2015
#for rt in 1.4 1.8; do paste <(ls newFilteredVCFs/*.gz | cut -d "." -f4) <(cat newFilteredVCFs/*_$rt".all" |cut -f1-3) <(cat newFilteredVCFs/*_$rt".gml"|cut -f1-3) <(cat newFilteredVCFs/*_$rt".sml" |cut -f1-3) > resApr13/r$rt"_pred"; done


#April 13, 2015
#for rt in 1.4 1.8; do paste <(ls newFilteredVCFs/tumor/*.gz | cut -d "." -f4) <(cat newFilteredVCFs/tumor/*_$rt".all" |cut -f1-3) <(cat newFilteredVCFs/tumor/*_$rt".gml"|cut -f1-3) <(cat newFilteredVCFs/tumor/*_$rt".sml" |cut -f1-3) > resApr13/r$rt"_tumor.pred"; done
#for rt in 1.4 1.8; do paste <(ls newFilteredVCFs/germ/*.gz | cut -d "." -f4) <(cat newFilteredVCFs/germ/*_$rt".all" |cut -f1-3) <(cat newFilteredVCFs/germ/*_$rt".gml"|cut -f1-3) <(cat newFilteredVCFs/germ/*_$rt".sml" |cut -f1-3) > resApr13/r$rt"_germ.pred"; done
#cat <(cat r1.4_germ.pred | awk '{print "g\t" $0}') <(cat r1.4_tumor.pred | awk '{print "t\t" $0}') | sort -k5 -n -r >1.4.bayes

