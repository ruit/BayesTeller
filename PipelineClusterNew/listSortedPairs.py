#!/usr/bin/python

# Tian R. <tianremi@gmail.com>
# June 19, 2014
# June 27, 2014
# Oct 20, 2014

import sys, re
import subprocess


inputfile=sys.argv[1]
#Mu1	Donor1

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
	f=open(input, "r")	
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



testSorted=subprocess.call('cat '+inputfile+' |sed "s/_/\t/g" |sort -c -k1,1 -k2,2 -n ', shell=True)

if testSorted==0:
	#this file is sorted	
	Mut2MoreDonors(inputfile)

else:
	print "@ERROR, the input file should be sorted first!"		

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
