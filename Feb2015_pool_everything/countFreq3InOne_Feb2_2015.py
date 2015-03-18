#!/bin/python

# Tian R. tianremiATgmail.com
# Feb 2, 2015
import sys


def unique1(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked



def PurifyList(tumorPat, pattern, testList):
	
	import re
	# Feb 2, 2015

	'''
	tumorPat is the input
	newT is the output list
	pattern is the included
	testList is the excluded list
	'''
	
	# patlist might contain "void" (empty) elements
	newT=[]

	for e in tumorPat:
		fuhe=re.search(pattern, e) #"Patient"
		
		if fuhe and e not in testList:
			newT.append(e)
	
	return newT




def CountNonRedundantTrainingFreq(input, outfile, pattern, testList):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#100	P	P22
	#101	void	P4
	#102	P3,P4	P8


	testList=['Patient1','Patient2']

	'''

	import re, gzip
	# June 19, 27, 2014
	

	f=gzip.open(input, "r")	

	out=open(outfile, "w")

	for line in f:

		outStr=""
	
		array=line.strip("\n").split("\t")
		
		for i in range(len(array)):
		
			if i==0:
				chrPos=array[0]
				outStr=chrPos
			else:
				germPat=array[i].split(",")
				
				# patlist might contain "void" (empty) elements
				newG=PurifyList(germPat, pattern, testList)

				outStr=outStr+"\t"+str(len(unique1(newG)))
		
		out.write(outStr+"\n")
		outStr=""

		#print "-------------------------------------------------"	
		#empty the lists
	

	f.close()
	out.close()



	

if len(sys.argv) > 4:

	input=sys.argv[1]
	output=sys.argv[2]
	pattern=sys.argv[3]
	exclusion=sys.argv[4]

	CountNonRedundantTrainingFreq(input, output, pattern, exclusion)	





