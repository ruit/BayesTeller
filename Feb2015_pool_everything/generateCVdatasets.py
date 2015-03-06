#!/usr/bin/python

#Feb 20, 2015
#Feb 23, 2015 Use one or two good days to write one major program, do not delay!!!!!!!!!
#Feb 24, 2015, Feb26 Run over 3909 patients, 2000SNVs, 30min !!!!!
#Tian R. <tianremi@gmail.com>
#March 4, 2015

import sys, os, gzip, re, random
import math

'''
1 A set of SNVs
2 path to SNVs for different individuals
3 splitted lists of random samples
4 A scoring/ ranking method

usage: python generateCVdatasets.py snvTest ./testdir 3909.panCancer.list_Feb19_2015.txt test1.germ.gz test1.tumor.gz

'''

if len(sys.argv)<6:
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
		array=line.strip("\n").split("\t")
		for pat in testList:
			
			if pat in array[1].split(","):
				
				patSNV_dict.setdefault(pat,[]).append(array[0])
				#print pat, array[0]
	
#	print patSNV_dict
	for pat in testList:
		f_out=gzip.open(outfilePath+pat+"_"+suffix+".gz", "wb")
		
		if patSNV_dict.get(pat) is None: pass
		else:	
			f_out.write("\n".join(patSNV_dict.get(pat)))
		f_out.close()

	return patSNV_dict
	
	del testList
	del patSNV_dict
	del line


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
	freqDict={}#Feb 24, 2015

	for line in f:

		outStr=""
		freq=[]
		
		array=line.strip("\n").split("\t")
		
		for i in range(len(array)):
		
			if i==0:
				chrPos=array[0]
				outStr=chrPos
			else:
				germPat=array[i].split(",")
				
				# patlist might contain "void" (empty) elements
				newG=PurifyList(germPat, pattern, testList)
				
				freq.append(len(unique1(newG)))
		
				outStr=outStr+"\t"+str(len(unique1(newG)))
					
		freqDict[chrPos]=freq
		freq=[]

		out.write(outStr+"\n")
		outStr=""

		#print "-------------------------------------------------"	
		#empty the lists
	

	f.close()
	out.close()
	return freqDict



def writeTestPatList(testList, outTestPatFile):

	f_out = open(outTestPatFile, "w")
	for p in testList:
		f_out.write(p+"\n")
	f_out.close()




def bayesTeller(Dict_Bayes, testL, Dict_snv, Dict_freq_germ, Dict_freq_tumor):

	'''
	for a given test list of patients, each patient is represented by a set of SNVs, given model Dict_freq, match each patient with an list of freqs 
	'''

	for pat in testL: #test patients based on a pool of SNVs
		#for "germ" samples
		if Dict_snv.get(pat) is None: Dict_Bayes[pat]=None
		else:
		#print Dict_snv_germ_A.get(pat)
			Dict_Bayes[pat]=[math.log10((0.5+Dict_freq_tumor.get(snv)[0])/(0.5+Dict_freq_germ.get(snv)[0])) \
for snv in Dict_snv.get(pat)] # zero!
	
	return Dict_Bayes




def pred(baydict, label, cutoff):
	'''
	print per patient bayesTeller score, recalibrated
	'''
	keyR=[]

	for pat in baydict:
		if baydict.get(pat) is None:
			print pat+"_"+"0"+"_0"+"\t"+"NA"+"\t"+label
		else:
			logR=baydict.get(pat)
			l=len(logR)
			logR.sort(reverse=True)
			for ratio in logR:
				if ratio>=cutoff:
					keyR.append(ratio)

			#keyR=logR[0:int(topR*l)] # take top 30% hits
			if len(keyR)==0: print pat+"_"+"0"+"_0"+"\t"+"NA"+"\t"+label
			else:
			
				val=reduce(lambda x,y: x+y, keyR)
			#val=float(val/(l*topR))
				print pat+"_"+str(l)+"_"+str(len(keyR))+"\t"+str(val/l)+"\t"+label
				keyR=[]



def getWeights(lista, element):

	return float(len(lista)-lista.index(element))/len(lista)


try:
	setSNV=sys.argv[1] #the list of selected SNVs
	tempDir=sys.argv[2] #the tempDir to hold individual samples of a collection of SNVs
	sampleList=sys.argv[3] #the total list of samples
	germFile=sys.argv[4] # the aggregated germline file
	tumorFile=sys.argv[5] # the aggregated tumor file
	cutoff=float(sys.argv[6]) #only take top 0.3 (30%) for sum (snvlogR)
	cutoff=math.log10(cutoff) # take log10, March 4, 2015
#	print sys.argv

except:
	pass



try:
    	os.makedirs(tempDir)
except OSError:
    	if not os.path.isdir(tempDir): raise

#March 2, 2015
#consider weights for SNVs, top ranked are of higher priorities

snvpool=[snv.strip("\n") for snv in open(setSNV, "r")]





#step one, only consider selected SNVs, thereby reducing file sizes

extractBySNVs(setSNV, tumorFile, ".ext")
extractBySNVs(setSNV, germFile, ".ext")

#step two, generate randomly splitted half-half patients list
splittedList=GenerateTestList(sampleList)
#print len(splittedList[0]), len(splittedList[1])

#step three, freq based on train patient list

#use first half for training, testdat is [1]
Dict_freq_germ_A=CountNonRedundantTrainingFreq(germFile+".ext.gz", tempDir+"/germ_freq.A", ".", splittedList[1]) 
#print Dict_freq_germ_A
Dict_freq_tumor_A=CountNonRedundantTrainingFreq(tumorFile+".ext.gz", tempDir+"/tumor_freq.A", ".", splittedList[1]) 
writeTestPatList(splittedList[1], tempDir+"/testPartA.tab")
Dict_snv_germ_A=extractPatIDs(splittedList[1], germFile+".ext.gz", tempDir+"/", "germ.snvsa")
#print Dict_snv_germ_A
Dict_snv_tumor_A=extractPatIDs(splittedList[1], tumorFile+".ext.gz", tempDir+"/", "tumor.snvsa")
#



Dict_Bayes_germ={}#patient, lists of SNVs-> logRatio
Dict_Bayes_germ=bayesTeller(Dict_Bayes_germ, splittedList[1], Dict_snv_germ_A, Dict_freq_germ_A, Dict_freq_tumor_A)

Dict_Bayes_tumor={}
Dict_Bayes_tumor=bayesTeller(Dict_Bayes_tumor, splittedList[1], Dict_snv_tumor_A, Dict_freq_germ_A, Dict_freq_tumor_A)


#use the second half tor training, testdata is [0]

Dict_freq_germ_B=CountNonRedundantTrainingFreq(germFile+".ext.gz", tempDir+"/germ_freq.B", ".", splittedList[0]) 
Dict_freq_tumor_B=CountNonRedundantTrainingFreq(tumorFile+".ext.gz", tempDir+"/tumor_freq.B", ".", splittedList[0]) 
writeTestPatList(splittedList[0], tempDir+"/testPartB.tab")
Dict_snv_germ_B=extractPatIDs(splittedList[0], germFile+".ext.gz", tempDir+"/", "germ.snvsb")
Dict_snv_tumor_B=extractPatIDs(splittedList[0], tumorFile+".ext.gz", tempDir+"/", "tumor.snvsb")


Dict_Bayes_germ=bayesTeller(Dict_Bayes_germ, splittedList[0], Dict_snv_germ_B, Dict_freq_germ_B, Dict_freq_tumor_B)

Dict_Bayes_tumor=bayesTeller(Dict_Bayes_tumor, splittedList[0], Dict_snv_tumor_B, Dict_freq_germ_B, Dict_freq_tumor_B)



pred(Dict_Bayes_germ, "g", cutoff)
pred(Dict_Bayes_tumor, "t", cutoff)



