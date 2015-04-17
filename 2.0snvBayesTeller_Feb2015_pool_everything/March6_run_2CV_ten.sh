#Tian R.
#March 6, 2015

# run 2CV 10 times
# perf for all and for each cancer type

#python  generateCVdatasets.py setSNV \
#tempDir \
#sampleList \
#germFile \
#tumorFile \
#cutoff \ # ratio


#March 7, 2015, each loop overwrite???


mkdir output_March6_2015

for cutoff in 1.5 2
do

	for num in 1 2 3 4 5 6 7 8 9 10

	do

		python generateCVdatasets.py snvList \
"crossV_"$num \
3909.panCancer.list_Feb19_2015.txt \
panCancer.germ.gz.sorted.org.gz panCancer.tumor.gz.sorted.org.gz \
$cutoff > "output_March6_2015/perf_"$cutoff"_"$num".pred"
	done
done
