#!/usr/bin/python

# Tian R. <tianremi@gmail.com>
# June 19, 2014
# June 27, 2014
import sys, re

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

Mut2MoreDonors(inputfile)


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
