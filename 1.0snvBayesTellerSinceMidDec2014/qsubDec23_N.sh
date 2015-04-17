#Dec 23, 2014
#run time 30min

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
	qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNBnew.sh $dir"
done
