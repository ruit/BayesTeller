#Tian R.
#March 18, 2015

# run 2CV 10 times
# perf for all and for each cancer type

#python  generateCVdatasets.py setSNV \
#tempDir \
#sampleList \
#germFile \
#tumorFile \
#cutoff \ # ratio


#March 7, 2015, each loop overwrite???

dst="output_March18_2015"
snvList=$1


mkdir $dst

for cutoff in 1.5 2 2.5
do

	for num in 1 2 3 4 5 6 7 8 9 10

	do

		python generateCVdatasetsMar18_2015.py $snvList \
"crossV_"$num \
3909.panCancer.list_Feb19_2015.txt \
panCancer.germ.gz.sorted.org.gz panCancer.tumor.gz.sorted.org.gz \
$cutoff > $dst"/perf_"$cutoff"_"$num".pred"
	done
done
