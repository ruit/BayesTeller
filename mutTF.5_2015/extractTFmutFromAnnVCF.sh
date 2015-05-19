#May 19, 2015
#Tian R.
#extract TF from annotated exonic VCFs for those with stopgain and nonsynonymous mutations (severity can be predicted)

pathTovcf="/alto/projects/common/tcga/luad/broad.mit.edu_LUAD.IlluminaGA_DNASeq_Cont.Level_2.0.4.0/data/"

#tf=$1

tmpdir=$RANDOM"_"$RANDOM

mkdir $tmpdir

for file in `ls $pathTovcf*vcf.avinput.exonic_variant_function.gz`

do

	root=`basename $file`
	zcat $file | grep Somatic > $tmpdir/$root".soma.tmp" #the somatic mutations
		
	for tf in `cat allTF.uniq`
	do
  
		cat $tmpdir/$root".soma.tmp"| grep -w $tf | grep stopgain >$tmpdir/$root".stpg"
		cat $tmpdir/$root".soma.tmp"| grep -w $tf | grep nonsynonymous >$tmpdir/$root".nonsyn"
		cat $tmpdir/$root".stpg" $tmpdir/$root".nonsyn" > $tmpdir/$root"_"$tf".mut" 
		rm $tmpdir/$root".stpg"
		rm $tmpdir/$root".nonsyn"

	done

	cat $tmpdir/$root"_"*".mut" > $tmpdir/$root".mut"
	rm $tmpdir/$root"_"*".mut"
	rm $tmpdir/$root".soma.tmp"
done
		
