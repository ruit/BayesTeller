#!/bin/sh -login

for i in 1 2 3 4 5 6 7 8 9 10 

do

	qsub -j y -cwd -o /scratch/tianr/queue/$i"_"$RANDOM".out" \
-b y "/home/tianr/1Projects/1SNVblocks/run_BayesTeller_Oct23_Dec30_2014/twoFoldCVNoClustering.sh"

done
