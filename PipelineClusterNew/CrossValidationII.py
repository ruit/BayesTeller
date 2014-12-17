#!/usr/bin/python
# Tian R. <tianremiATgmail.com>
# Nov 4, 2014
# Nov 5, 2014
# Dec 12, 2014

# Randomly mask some patients as test data

import random, math
import sys, re


infile1 = sys.argv[1]
infile2 = sys.argv[2]
totalPatListFile = sys.argv[3]
tempdir = sys.argv[4]

# Sep 19, 2014
# Tian R. compare SNVs freq in normal and tumor, build a pridictive model
# Sep 22, 2014. 1000G SNV freq <= 0.5% as the filtering at the beginning!!!

def GenerateTestList (totalPatListFile):
	#totalList is a single column file, containing all patients
	totalL=[]
	for line in open(totalPatListFile, "r"):
		line=line.strip("\n")
		totalL.append(line)

	halfOne=[]
	halfTwo=[]

	while len(halfOne) < round(len(totalL) / 2) :
		singleP=random.choice(totalL)
		if singleP not in halfOne:
			halfOne.append(singleP)

	for anyP in totalL:
		if anyP not in halfOne:
			halfTwo.append(anyP)

	return [halfOne, halfTwo]




def unique1(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked




def SortedKey2GivenStartEndAsList(input, outfile, testList):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#100	A_start
	#100	P2
	#101	void
	#102	P3,P4
	#102	Z_end


	testList=['Patient1','Patient2']

	'''

	import re
	# June 19, 27, 2014
	last_zuo=""
	current_zuo=""
	last_you=""
	current_you=""
	
	f=open(input, "r")	

	out=open(outfile, "w")

	for line in f:
		array=line.strip("\n").split("\t")
		current_zuo=array[0]
		current_you=array[1]
		
		#start with "A_start", Oct 21, 2014		
		if current_you=="A_start":
			last_zuo="what"
			last_you=current_you
		
		if last_you != "":
			#if last_you is not empty
			if current_you != "Z_end":
				if last_zuo != "what":
					if last_zuo=="":
						last_zuo=current_zuo
					else:
						last_zuo=last_zuo+","+current_zuo
					last_you=last_you+","+current_you
				else:
					last_zuo=""
					last_you=last_you+","+current_you
			else:
	
				patlist=last_you.split(",")
		
				# patlist might contain "void" (empty) elements
				newlist=[]
	
				for e in patlist:
					fuhe=re.search("Patient", e)
					if fuhe and e not in testList:
						newlist.append(e)
						#print "Patients counted for modeling are " + e
 
				out.write(last_zuo+"\t"+str(len(unique1(newlist)))+"\n")
				print "Patients counted for modeling are "+ ";".join(unique1(newlist))
				#print "-------------------------------------------------"	
				#empty the lists
				last_zuo=""
				last_you=""
	f.close()
	out.close()


	

def writeTestPatList(testList, outTestPatFile):

	f_out = open(outTestPatFile, "w")
	for p in testList:
		f_out.write(p+"\n")
	f_out.close()


def main():
	pairs=GenerateTestList(totalPatListFile)
	print pairs	
	#complete 2CV!!!	
	SortedKey2GivenStartEndAsList(infile1, infile1+".outA", pairs[0])
	SortedKey2GivenStartEndAsList(infile2, infile2+".outA", pairs[0])
	writeTestPatList(pairs[0], tempdir+"/testPatList4A.tab") ###!!!! double check Dec 12, 2014	
	print "-------------------------------------"	
	SortedKey2GivenStartEndAsList(infile1, infile1+".outB", pairs[1])
	SortedKey2GivenStartEndAsList(infile2, infile2+".outB", pairs[1])
	writeTestPatList(pairs[1], tempdir+"/testPatList4B.tab")

	
	



if __name__=="__main__":
	main()

