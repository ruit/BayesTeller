#!/bin/sh -login
#TianR. Dec 11, 2014
#mark and cat all SNV-patient data to get a pan-cancer data for 14 types

#mkdir panCancer

markByCancer(){

array_cancerTypes=(brca colon kirc kirp lgg lihc luad prad acc blca paad pcpg skcm stad thca)

tLen=${#array_cancerTypes[@]}


#touch panCancer.germ
#touch panCancer.tumor

for (( i=0; i<${tLen}; i++));
do
        cancer=${array_cancerTypes[$i]}

	if [ $cancer == "brca" ]; then
		echo "do nothing"
	else
		echo $cancer
		cat "all"$cancer"Cancergermline.filtered" | sed "s/Patient/Patient@"$cancer"@/g" > "panCancer/all"$cancer"Cancergermline.filtered.marked"
		cat "all"$cancer"Cancertumor.filtered" | sed "s/Patient/Patient@"$cancer"@/g" > "panCancer/all"$cancer"Cancertumor.filtered.marked"
	fi
done
}


catSortOrg () {
cat panCancer/*germline.*marked | sed "s/_/\t/g" |sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}'> panCancer.germ
cat panCancer/*tumor.*marked | sed "s/_/\t/g" |sort -k1,1 -k2,2 -n | awk '{print $1"_"$2"\t"$3}'> panCancer.tumor 


python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/listSortedPairs.py panCancer.germ > panCancer.germ.org
#rm panCancer.germ

python /home/tianr/1Projects/1SNVblocks/NewSinceOct16_2014/listSortedPairs.py panCancer.tumor > panCancer.tumor.org
#rm panCancer.tumor
}

#markByCancer
#catSortOrg


# 14 datasets 3909
# sum(tcga[-1,2])*0.05 
#[1] 195.45

# purify the panCancer 
tmpfile=panCgermN_$RANDOM

cat panCancer.germ.org | cut -f2 | awk -F "," '{print NF}' > $tmpfile

paste panCancer.germ.org $tmpfile panCancer.tumor.org | awk '{if($1==$4 && $3<=200) print $1"\t"$2"\t"$5}' | grep -v Y_ | grep -v MT_ > panCancer.2in1
