# Tian R. <tianremi@gmail.com>
# Sep 23, 2014
#Feb 4, 2015



#script="/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh"
#Kgenomefile="/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2"
script=./preprocessVCF_Feb3_2015.sh

echo "./preprocessVCF_Feb3_2015.sh colon /home/tianr/1Projects/1SNVblocks/DATA/tcga/colon/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/ PASS" >run_$RANDOM.sh
echo "$script acc /home/tianr/1Projects/1SNVblocks/DATA/tcga/acc/ other" > run_$RANDOM.sh
echo "$script blca /home/tianr/1Projects/1SNVblocks/DATA/tcga/blca/ other"> run_$RANDOM.sh
echo "$script paad /home/tianr/1Projects/1SNVblocks/DATA/tcga/paad/ other"> run_$RANDOM.sh
echo "$script pcpg /home/tianr/1Projects/1SNVblocks/DATA/tcga/pcpg/ other"> run_$RANDOM.sh
echo "$script stad /home/tianr/1Projects/1SNVblocks/DATA/tcga/stad/ other"> run_$RANDOM.sh
echo "$script thca /home/tianr/1Projects/1SNVblocks/DATA/tcga/thca/ other"> run_$RANDOM.sh
echo "$script skcm /home/tianr/1Projects/1SNVblocks/DATA/tcga/skcm/ other"> run_$RANDOM.sh
echo "$script luad /home/tianr/1Projects/1SNVblocks/DATA/tcga/luad/broad.mit.edu_LUAD.IlluminaGA_DNASeq_Cont.Level_2.0.4.0/ other"> run_$RANDOM.sh 
echo "$script prad /home/tianr/1Projects/1SNVblocks/DATA/tcga/prad/broad.mit.edu_PRAD.IlluminaGA_DNASeq_Cont_curated.Level_2.1.4.0/ other"> run_$RANDOM.sh 
echo "$script kirc /home/tianr/1Projects/1SNVblocks/DATA/tcga/kirc/broad.mit.edu_KIRC.IlluminaGA_DNASeq_Cont_automated.Level_2.1.5.0/ other"> run_$RANDOM.sh
echo "$script kirp /home/tianr/1Projects/1SNVblocks/DATA/tcga/kirp/broad.mit.edu_KIRP.IlluminaGA_DNASeq_Cont_curated.Level_2.1.1.0/ other"> run_$RANDOM.sh
echo "$script lgg /home/tianr/1Projects/1SNVblocks/DATA/tcga/lgg/broad.mit.edu_LGG.IlluminaGA_DNASeq_Cont_automated.Level_2.1.4.0/ other"> run_$RANDOM.sh
echo "$script lihc /home/tianr/1Projects/1SNVblocks/DATA/tcga/lihc/broad.mit.edu_LIHC.IlluminaGA_DNASeq_Cont_automated.Level_2.1.0.0/ other"> run_$RANDOM.sh

chmod +x run*.sh

for file in `ls run_*.sh` 
do
	qsub -j y -cwd -b y "./$file"
done


