# Jan 9, 2015 TianR. 
# Given alpha, select of set of SNVs for cross validation

modelFileOri=$1
alpha=$2

path="/home/tianr/1Projects/1SNVblocks/snvBayesTellerSinceMidDec2014/"

cat $modelFileOri | cut -f2 | sort | uniq -c | sort -k2 -n > $modelFileOri".uniqc"

Rscript $path"ReFilterSNVbyCumulativeCoverage.R" \
$modelFileOri".uniqc" $alpha \
$alpha"_"$modelFileOri"_temp.C"

num=`cat $alpha"_"$modelFileOri"_temp.C"`

rm $alpha"_"$modelFileOri"_temp.C"
rm $modelFileOri".uniqc"

cat $modelFileOri | awk -v N=$num '{ if ($3 >= N) print $0}'
