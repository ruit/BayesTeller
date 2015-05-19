#TianR.
#May 14, 2015
#Assign ContrastRank score (cancer associtated gene prioritization method) to uniq set of TFs

cell=$1

#cat cell6312mmc4/$cell/genes-regulate-genes.txt | tr "\t" "\n" | sort | uniq > ./$cell".tf.uniq"
#join <(sort -k1,1 luad_normal_scores.txt) <(sort -k1,1 $cell".tf.uniq") \
#| sort -k5 -n -r >$cell".tf.uniq.score"


#join <(sort -k1,1 luad_normal_scores.txt) <(sort -k1,1 $cell".tf.uniq.diff") \
#| sort -k5 -n -r >$cell".tf.uniq.diff.score"


#tianrMac:3TF_CDS_mutations tianr$ for i in `ls *.tf.uniq.score `; do a=`awk 'BEGIN {n=0} {n=n+$NF} END {print n}' $i`; echo $i $a; done


cat cell6312mmc4/$cell/genes-regulate-genes.txt | awk '{if($1 < $2){print $1"\t"$2} else {print $2"\t"$1}}' \
| sort | uniq |tr "\t" "\n" | sort | uniq -c | sort -k1 -n -r | awk '{print $2"\t"$1}' > $cell".degree.tab"
