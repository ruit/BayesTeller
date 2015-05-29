#May 27, 2015
#May 29, 2015
protein=$1 # should be geneID
uniprot=$2 #should be proID

#get mutation list, sorted by positions
zcat Notsilent.gz | grep -w $protein | grep -v EN | grep -v \? | grep -v \* | grep -v _ | grep -v \> | grep -v del > $protein".temp"
paste $protein".temp" <(cat $protein".temp"| cut -d " " -f2 |sed -r 's/p\.\S{1}([0-9]+)\S{1}/\1/') | sort -k3 -n | awk '{print $2}' | sed "s/p\.//g" >$protein".mut"
rm $protein".temp"


#get protein sequences
zcat human.uniprot.fasta.gz | awk -v UNI=$uniprot '{ if (substr($0,0,1)==">"){split($0, A, "|"); if (A[2]==UNI) {p=1} else {p=0}} else {if (p==1) {print $0}}}'  > $protein".aa"

#run the predictor:
python /share/apps/SNPs-and-GO/SNPs-and-GO.py $protein".aa" -f $protein".mut" -p 2 -s scale.txt > $protein".out" #-p 2 means run PhD SNP

rm $protein".aa.blast" 
rm $protein".aa.hssp"
rm $protein".aa.txt.scale"



#cat hs.HGNC2pro.merged.tab |grep -v "\#" |awk '{print $1"\t"$2}' |while read gene prot;  do qsub -j y -cwd -b y "./runPhDSNP.sh $gene $prot";  done


#check consitency 
#seq=`cat $protein".aa" |grep -v \> | tr -d "\n"`	 
#paste $protein".mut" <(cat $protein".mut" | sed -r 's/\S{1}([0-9]+)\S{1}/\1/')\
#| awk -v SEQ=$seq '{if (substr(SEQ, $2, 1) == substr($1, 1, 1)) {print $1}}' > $protein".mut.checked"
#paste $protein".mut" <(cat $protein".mut" | sed -r 's/\S{1}([0-9]+)\S{1}/\1/')\
#| awk -v SEQ=$seq '{if (substr(SEQ, $2, 1) != substr($1, 1, 1)) {print substr(SEQ, $2, 1), $1}}' >$protein".mut.wrong" 
