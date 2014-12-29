#!/usr/bin/python
# Tian R. <tianremiATgmail.com>
# Nov 4, 2014
# Nov 5, 2014
# Dec 12, 2014
# Dec 22, 2014 @@@@20min for panCancer
# Randomly mask some patients as test data

import random, math
import sys, re


infile1 = sys.argv[1]
totalPatListFile = sys.argv[2]
tempdir = sys.argv[3]

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




def CountNonRedundantTrainingFreq(input, outfile, testList):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#100	P	P22
	#101	void	P4
	#102	P3,P4	P8


	testList=['Patient1','Patient2']

	'''

	import re
	# June 19, 27, 2014
	
	f=open(input, "r")	

	out=open(outfile, "w")

	for line in f:
		array=line.strip("\n").split("\t")
		
		chrPos=array[0]
		germPat=array[1].split(",")
		tumorPat=array[2].split(",")


		# patlist might contain "void" (empty) elements
		newG=[]
		for e in germPat:
			fuhe=re.search("Patient", e)
			if fuhe and e not in testList:
				newG.append(e)

		# patlist might contain "void" (empty) elements
		newT=[]
		for e in tumorPat:
			fuhe=re.search("Patient", e)
			if fuhe and e not in testList:
				newT.append(e)


		out.write(chrPos+"\t"+str(len(unique1(newG)))+"\t"+str(len(unique1(newT)))+"\n")
		#print "-------------------------------------------------"	
		#empty the lists
	

	f.close()
	out.close()


	

def writeTestPatList(testList, outTestPatFile):

	f_out = open(outTestPatFile, "w")
	for p in testList:
		f_out.write(p+"\n")
	f_out.close()


def main():
	pairs=GenerateTestList(totalPatListFile)
	#complete 2CV!!!	
	CountNonRedundantTrainingFreq(infile1, tempdir+"/modelA", pairs[0])
	writeTestPatList(pairs[0], tempdir+"/testPatList4A.tab") ###!!!! double check Dec 12, 2014	
	
	CountNonRedundantTrainingFreq(infile1, tempdir+"/modelB", pairs[1])
	writeTestPatList(pairs[1], tempdir+"/testPatList4B.tab")

	

if __name__=="__main__":
	main()

