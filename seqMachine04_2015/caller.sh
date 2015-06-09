#!/bin/bash

#Rui X. Tian <tianremi@gmail.com>
#June 2, 2015 [April, May2015]
#@June 3,2015
#human genome assembly, by default is hg19 (GRh37)
#bwaIndexes/hg19
#June 8, 2015 samtools version issues for rmdup functionality

########################################################
#4 major steps for variant calling 2015.06.03
#StepI:QC
#StepII:reads alignment
#StepII:mark duplicates
#StepIV:realignment, recalibration, variant calling, filtering
#######################################################


while [[ $# > 0 ]]
do
key="$1"

case $key in
    -f1|--inputfile1)
    input1="$2"
    shift # past argument
    ;;
    -f2|--inputfile2)
    input2="$2"
    shift # past argument
    ;;

    -a|--aligner)
    aligner="$2"
    shift # past argument
    ;;

    -g|--refpath)
    refpath="$2"
    shift # past argument
    ;;

    -c|--caller)
    caller="$2"
    shift #past argument
    ;;

    --fasta)
    fapath="$2"
    shift
    ;;

    *)
            # unknown option
    ;;
esac
shift # past argument or value
done



if [[ $input1 ]] && [[ $input2 ]]; then
    echo "Info: Paired end sequencing data analysis."
    pair="p"
else

    echo "Info: single end sequencing data analysis."
    pair="s"

fi

if [[ $refpath ]]; then
    echo "Info: Path to reference genome is: "${refpath}
else
    echo "Error: Path to reference genome is needed!"

fi

echo $input1
echo $input2
echo $aligner
echo $refpath
echo $caller
echo $fapath


#QA------------------------Step 1 -------------------------------
#QA also make some trimming based on fastqc reports, use trimomatic
#fastqc 1_1.fq -f fastq -o fastqc/


#aligners------------------Step 2 ------------------------------

singleEndBWA () {
    
    INPUT1=$1
    REFPATH=$2
    
    #sudo docker run -v `pwd`:/data -it tianrvarcallers:1.00 
    bwa aln -t 4 -f /data/$INPUT1".sai" -I /data/$REFPATH /data/$INPUT1
    #sudo docker run -v `pwd`:/data -it tianrvarcallers:1.00 
    mkdir /data/output_dir
    bwa samse -r "@RG\tID:<ID>\tLB:<LIBRARY_NAME>\tSM:<SAMPLE_NAME>\tPL:ILLUMINA" -f \
/data/output_dir/$INPUT1".sam" /data/$REFPATH /data/$INPUT1".sai" /data/$INPUT1
   
    samtools view -bS /data/output_dir/$INPUT1".sam" |samtools sort - /data/output_dir/$INPUT1".srt" #sorted 

    #reheadering and remove PCR duplicates
    #samtools view -H /data/output_dir/$INPUT1".srt.bam" |sed "s/SO:unsorted/SO:sorted/g" |samtools reheader - /data/output_dir/$INPUT1".srt.bam"

    #samtools rmdup /data/output_dir/$INPUT1".srt.bam" /data/output_dir/$INPUT1".srt.cleaned.bam"
 
    rm /data/$INPUT1".sai"
    rm /data/output_dir/$INPUT1".sam"
}


pairedEndBWA () {

    INPUT1=$1
    INPUT2=$2
    REFPATH=$3

    echo $INPUT1 
    echo $INPUT2

    bwa aln -t 4 -f /data/$INPUT1".sai" -I /data/$REFPATH /data/$INPUT1
    bwa aln -t 4 -f /data/$INPUT2".sai" -I /data/$REFPATH /data/$INPUT2

    #sudo docker run -v `pwd`:/data -it tianrvarcallers:1.00 
    mkdir /data/output_dir
    bwa sampe -r "@RG\tID:<ID>\tLB:<LIBRARY_NAME>\tSM:<SAMPLE_NAME>\tPL:ILLUMINA" -f \
/data/output_dir/$INPUT1"_"$INPUT2".sam" /data/$REFPATH /data/$INPUT1".sai" /data/$INPUT2".sai" \
/data/$INPUT1 /data/$INPUT2
    samtools view -bS /data/output_dir/$INPUT1"_"$INPUT2".sam" |samtools sort - /data/output_dir/$INPUT1"_"$INPUT2".srt" #sorted
  
    #reheadering and remove PCR duplicates
    samtools view -H /data/output_dir/$INPUT1"_"$INPUT2".srt.bam" |sed "s/SO:unsorted/SO:sorted/g" |samtools reheader - /data/output_dir/$INPUT1"_"$INPUT2".srt.bam" >/data/output_dir/$INPUT1"_"$INPUT2".reh.srt.bam"
    
    #bugs!!!!Here:
    /opt/samtools-0.1.19/samtools rmdup /data/output_dir/$INPUT1"_"$INPUT2".srt.bam" /data/output_dir/$INPUT1"_"$INPUT2".srt.cleaned.bam"

    rm /data/$INPUT1".sai" 
    rm /data/$INPUT2".sai"
    rm /data/output_dir/$INPUT1"_"$INPUT2".sam"
    rm /data/output_dir/$INPUT1"_"$INPUT2".srt.bam"
    rm /data/output_dir/$INPUT1"_"$INPUT2".reh.srt.bam"

}


pairedENDbowtie () {
    INPUT1=$1
    INPUT2=$2    
    pathbowtie2ind=$3

    bowtie2 --end-to-end --very-fast --rg-id \
"@RG\tID:1\tLB:project_name\tSM:1\tPL:ILLUMINA" \
-x /data/$pathbowtie2ind -q -1 /data/$INPUT1 -2 /data/$INPUT2 | \
samtools view - -Sb -o /data/output_dir/$INPUT1"_"$INPUT2".bam" #hg19Ref/hg19-bt2.fa.fai 
#bowtie2_indexes/hg19

}





#mark duplicates---------------------------------------------------

    #java -jar /opt/picard-tools-1.119/MarkDuplicates.jar \
#INPUT=/data/ouput_dir/$INPUT1"_"$INPUT2".srt.bam" \
#OUTPUT=/data/output_dir/$INPUT1"_"$INPUT2".srt.cleaned.bam" \
#METRICS_FILE=/data/output_dir/metrics.txt




#variant callers-------------Step 3-------------------------------

hliCaller () {

    nonredunBAM=$1
    path2RefFA=$2
    
    samtools mpileup -B -u -f $path2RefFA $nonredunBAM| bcftools call - -m -v >$nonredunBAM".vcf"
    #hg19Ref/hg19-bt2.fa 
}


#freeBayes by Erik Garrison

freeBayesCaller () {
    
    nonredunBAM=$1
    path2RefFA=$2
    samtools index $nonredunBAM
    freebayes -f $path2RefFA $nonredunBAM >$nonredunBAM".vcf"
    
}




#Main----------------------------------------------------------------

#paired end:
if [[ ${pair} == "p" ]] && [[ ${aligner} == "bwa" ]]; then
    pairedEndBWA $input1 $input2 $refpath

elif [[ ${pair} == "p" ]] &&  [[ ${aligner} == "bowtie2" ]]; then 
    pairedENDbowtie $input1 $input2 $refpath
   
fi   
#single end
#elif [[ ${pair} =="s" ]] && [[ $refpath ]]; then
#    singleEndBWA $input1 $refpath
    #echo "do nothing"
#fi

if [[ ${caller} == "hli" ]] && [[ ${pair}=="p" ]] && [[ ${fapath} ]]; then
    echo "Variant caller is set to be hli, samtools mpileup plus bcftools"
    hliCaller /data/output_dir/$INPUT1"_"$INPUT2".srt.cleaned.bam" /data/$fapath

elif [[ ${caller} == "freebayes" ]] && [[ ${pair}=="p" ]] && [[ ${fapath} ]]; then
    echo "Variant caller is set to be freebayes"
    freeBayesCaller /data/output_dir/$INPUT1"_"$INPUT2".srt.cleaned.bam" /data/$fapath


fi
