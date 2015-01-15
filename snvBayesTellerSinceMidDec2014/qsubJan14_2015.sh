#Dec 23, 2014
#run time 30min


rank="ratio"

for dir in CV2_panCancer_17125output \
CV2_panCancer_17413output \
CV2_panCancer_23157output \
CV2_panCancer_23742output \
CV2_panCancer_24985output \
CV2_panCancer_28120output \
CV2_panCancer_31334output \
CV2_panCancer_66output \
CV2_panCancer_8114output \
CV2_panCancer_8963output


do
	rank="ratio"; suffix=".lodr"
	qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"

	rank="binorm"; suffix=".lodb"
	qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"
	
	rank="diff"; suffix=".lodf"
	qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"

done


