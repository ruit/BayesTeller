
input=$1

echo $input

        zcat $input | sed "s/_/\t/g" |grep -v \< \
|sort -k1,1 -k2,2 -n |awk '{print$1"_"$2"\t"$3}' | gzip > tempFeb5/$input".sorted.gz"

        python /home/tianr/1Projects/1SNVblocks/snvBayesTellerSinceMidDec2014/listSortedPairs2.py \
tempFeb5/$input".sorted.gz" "force" | gzip > tempFeb5/$input".sorted.org.gz"

