#!/usr/bin/python

# Tian R. <tianremi@gmail.com>
# June 19, 2014
# June 27, 2014
# Oct 20, 2014
# Jan 26, 2015
<<<<<<< HEAD
# Jan 29, 2015
=======
>>>>>>> 21834cba72bb1a2ad508176d94b82c79d82d3083

import sys, re
import subprocess
import gzip

inputfile=sys.argv[1]
#Mu1	Donor1
<<<<<<< HEAD
check=sys.argv[2]# Jan 29, 2015, whether check the input sorted or not
=======
>>>>>>> 21834cba72bb1a2ad508176d94b82c79d82d3083

# merger col2 sharing the same col1 values
# Aug 26, 2014
# 100K!!!!!!!!!!
# Give up DL, for publications!!!


def Mut2MoreDonors(input):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	'''
	last_mut=""
	current_mut=""
	last_donor=""
	current_donor=""
	f=gzip.open(input, "r")# gzip files, Jan 26, 2015	
	for line in f:
		array=line.strip("\n").split("\t")
		current_mut=array[0]
		current_donor=array[1]
		if last_mut=="":
			last_mut=current_mut
			last_donor=current_donor
		else:
			if current_mut != last_mut:
				print last_mut+"\t"+ last_donor
				last_mut=current_mut
				last_donor=current_donor
			else:
				last_donor=last_donor+","+current_donor

	print last_mut + "\t" +last_donor			
	f.close()



testSorted=subprocess.call('zcat -f '+inputfile+' |sed "s/_/\t/g" |sort -c -k1,1 -k2,2 -n ', shell=True)
#to test whether the aggregated file is sorted first by chr then by locus, Jan 26, 2015


<<<<<<< HEAD

if check=="force":
	Mut2MoreDonors(inputfile)

else:
	if testSorted==0:
	#this file is sorted	
		Mut2MoreDonors(inputfile)

	else:
		print "@ERROR, the input file should be sorted first!"		
=======
if testSorted==0:
	#this file is sorted	
	Mut2MoreDonors(inputfile)

else:
	print "@ERROR, the input file should be sorted first!"		
>>>>>>> 21834cba72bb1a2ad508176d94b82c79d82d3083

#Mut2MoreDonors(inputfile)


#def main ():
#	dict=KeyToListDict(inputfile) 
#	for k, v in sorted(dict.items()):
#		print k+"\t"+v	

#if __name__=="__main__":
#	main()


#Mu1	Donor1
#Mu2	Donor1
#Mu3	Donor2
#Mu2	Donor3
