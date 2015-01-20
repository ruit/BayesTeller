#Dec 23, 2014
#run time 30min
#Jan 16, 2015



while read dir;


do echo $dir;

	#rank="ratio"; suffix=".lodr"
	#qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"

	#rank="binorm"; suffix=".lodb"
	#qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"
	
	rank="diff"; suffix=".lodf"
	qsub -j y -cwd -o /scratch/tianr/queue/$dir".stdout" -b y "./runNB4.sh $dir $rank $suffix"

done < $1


