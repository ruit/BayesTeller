# Tian R. <tianremi@gmail.com>
# Sep 23, 2014

script="/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh"
Kgenomefile="/home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2"


/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh luad /home/tianr/1Projects/1SNVblocks/DATA/tcga/luad/broad.mit.edu_LUAD.IlluminaGA_DNASeq_Cont.Level_2.0.4.0/ other /home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2

/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh prad /home/tianr/1Projects/1SNVblocks/DATA/tcga/prad/broad.mit.edu_PRAD.IlluminaGA_DNASeq_Cont_curated.Level_2.1.4.0/ other /home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2

/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh colon /home/tianr/1Projects/1SNVblocks/DATA/tcga/colon/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/ PASS /home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified2

$script brca /home/tianr/1Projects/1SNVblocks/DATA/tcga/brca/ucsc.edu_BRCA.IlluminaGA_DNASeq_Cont.Level_2.1.1.0/ PASS $Kgenomefile

$script kirc /home/tianr/1Projects/1SNVblocks/DATA/tcga/kirc/broad.mit.edu_KIRC.IlluminaGA_DNASeq_Cont_automated.Level_2.1.5.0/ other $Kgenomefile

$script kirp /home/tianr/1Projects/1SNVblocks/DATA/tcga/kirp/broad.mit.edu_KIRP.IlluminaGA_DNASeq_Cont_curated.Level_2.1.1.0/ other $Kgenomefile

$script lgg /home/tianr/1Projects/1SNVblocks/DATA/tcga/lgg/broad.mit.edu_LGG.IlluminaGA_DNASeq_Cont_automated.Level_2.1.4.0/ other $Kgenomefile

$script lihc /home/tianr/1Projects/1SNVblocks/DATA/tcga/lihc/broad.mit.edu_LIHC.IlluminaGA_DNASeq_Cont_automated.Level_2.1.0.0/ other $Kgenomefile




#/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh colon /home/tianr/1Projects/1SNVblocks/DATA/tcga/colon/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq_Cont.Level_2.1.5.0/ PASS /home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified

#/home/tianr/1Projects/1SNVblocks/pipeline/PreprocessTCGAVCFGenerateInput.sh colon test/ PASS /home/tianr/1Projects/1SNVblocks/DATA/1kgenomes/1Kgenome.SNV.freq.purified

