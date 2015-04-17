# Tian R.<tianremiATgmail.com>
# Sep 26, 2014
# weigh by num of patients
# filter remaining "common variants"


cutoff=0.05



Nbrca=501
Ncolon=220
Nkirc=274
Nkirp=196
Nlgg=468
Nlihc=236
Nluad=625
Nprad=309



for cancer in colon kirc kirp lgg lihc luad prad
do
	if [ $cancer == "colon" ]; then
		total=$Ncolon
	elif [ $cancer == "kirc" ]; then
		total=$Nkirc

	elif [ $cancer == "kirp" ]; then
		total=$Nkirp

	elif [ $cancer == "lgg" ]; then
		total=$Nlgg

	elif [ $cancer == "lihc" ]; then
		total=$Nlihc

	elif [ $cancer == "luad" ]; then
		total=$Nluad

	elif [ $cancer == "prad" ]; then
		total=$Nprad

	fi


	paste "all"$cancer"Cancergermline.snvfreq" "all"$cancer"Cancertumor.snvfreq" | \
awk -v N=$total -v CUT=$cutoff '{if ($1 == $3 && ($2/N) <= CUT) print $1"\t"$4"\t"N"\t"($2/N)}' >$cancer".bloc"

	Rscript /home/tianr/1Projects/1SNVblocks/pipeline/bionorm.R $cancer".bloc"

	rm $cancer".bloc"
	#echo $cancer "has "$total" patients."

done


