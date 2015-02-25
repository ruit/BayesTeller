#!/usr/bin/python

#Feb 20, 2015
#Feb 23, 2015 Use one or two good days to write one major program, do not delay!!!!!!!!!
#Feb 24, 2015
#Tian R. <tianremi@gmail.com>


import sys, os, gzip, re, random

'''
1 A set of SNVs
2 path to SNVs for different individuals
3 splitted lists of random samples
4 A scoring/ ranking method

usage: python generateCVdatasets.py snvTest ./testdir 3909.panCancer.list_Feb19_2015.txt test1.germ.gz test1.tumor.gz

'''

if len(sys.argv)<5:
	raise SystemExit("I need arguments!")



def extractBySNVs(setSNV, aggreData, suffix):
	'''
	Only retrive data based on a list of selected SNVs
	'''
#	import gzip
	
	snvL=[snv.strip("\n") for snv in open(setSNV, "r")]
	f_out=gzip.open(aggreData+suffix+".gz", "wb")
	
	for line in gzip.open(aggreData,"rb"):
		if line.strip("\n").split("\t")[0] in snvL:
			f_out.write(line)
	f_out.close()
	del snvL 
	del line
	
	return 0
		


def extractPatIDs(testList, aggreData, outfilePath, suffix):
	'''
	Only retrieve data based on a list of test patient IDs
	'''
	import gzip
	
	#snvL=[snv.strip("\n") for snv in open(setSNV, "r")]
	patSNV_dict={}
	
	for line in gzip.open(aggreData,"rb"):
		for pat in testList:
			if pat in line.strip("\n").split("\t")[1].split(","):
				if not pat in patSNV_dict:
					patSNV_dict[pat]=[line.strip("\n").split("\t")[0]]
				else:
					patSNV_dict[pat]=patSNV_dict[pat].append(line.strip("\n").split("\t")[0])	
	
	
	for pat in testList:
		f_out=gzip.open(outfilePath+pat+"_"+suffix+".gz", "wb")
		
		if patSNV_dict.get(pat) is None: pass
		else:	
			f_out.write("\n".join(patSNV_dict.get(pat)))
		f_out.close()
	
	del testList
	del patSNV_dict
	del line
	
	return 0


#debug:
#extractPatIDs(['Patient:acc:12'], "panCancer.germ.gz.sorted.org.gz.ext.gz", "./testdir/", ".snvs")



	
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




def CountNonRedundantTrainingFreq(inputfile, outfile, pattern, testList):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#100	P	P22
	#101	void	P4
	#102	P3,P4	P8


	testList=['Patient1','Patient2']

	'''

	#import re, gzip
	# June 19, 27, 2014
	

	f=gzip.open(inputfile, "r")	

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



def writeTestPatList(testList, outTestPatFile):

	f_out = open(outTestPatFile, "w")
	for p in testList:
		f_out.write(p+"\n")
	f_out.close()




try:
	setSNV=sys.argv[1] #the list of selected SNVs
	tempDir=sys.argv[2] #the tempDir to hold individual samples of a collection of SNVs
	sampleList=sys.argv[3] #the total list of samples
	germFile=sys.argv[4] # the aggregated germline file
	tumorFile=sys.argv[5] # the aggregated tumor file
	print sys.argv

except:
	pass



try:
    	os.makedirs(tempDir)
except OSError:
    	if not os.path.isdir(tempDir): raise



#step one, only consider selected SNVs, thereby reducing file sizes

extractBySNVs(setSNV, tumorFile, ".ext")
extractBySNVs(setSNV, germFile, ".ext")

#step two, generate randomly splitted half-half patients list
splittedList=GenerateTestList(sampleList)
#print len(splittedList[0]), len(splittedList[1])

#step three, freq based on train patient list


#use first half for training
CountNonRedundantTrainingFreq(germFile+".ext.gz", tempDir+"/germ_freq.A", ".", splittedList[1]) 

CountNonRedundantTrainingFreq(tumorFile+".ext.gz", tempDir+"/tumor_freq.A", ".", splittedList[1]) 

writeTestPatList(splittedList[1], tempDir+"/testA.tab")



extractPatIDs(splittedList[1], germFile+".ext.gz", tempDir+"/", "germ.snvsa")
extractPatIDs(splittedList[1], tumorFile+".ext.gz", tempDir+"/", "tumor.snvsa")













