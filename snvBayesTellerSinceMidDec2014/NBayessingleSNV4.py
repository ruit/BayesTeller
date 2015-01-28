#!/usr/bin/python
# Nov 3, 2014
# Nov 4, 2014
# @Dec 5, 2014
#X_11
#1_22
#Dec 22, 2014 @Dec 23, 2014
#model file should be sorted based on freq tumor, Dec 23, 2014
#@Jan15, 2015, do not use round(), leave some decimal at least 3

import re, math, sys
from collections import defaultdict
import gzip # Jan 16, 2015


modelFile = sys.argv[1]
patientVariants = sys.argv[2]
total = sys.argv[3]
topN = sys.argv[4]

InExclude= sys.argv[5]
#"take" or "remove" 

def NBmodelInclude (modelFile, patientVariants, total, topN=20, laplace=float(0.5/1092)):
	
	i=0
	snvDict={}
	pNormalDict = defaultdict(lambda: defaultdict(float))
	pTumorDict = defaultdict(lambda: defaultdict(float))
	#print laplace

	for row in gzip.open(modelFile, "r"):
		#print row
		i = i + 1
		#print i
		if i <= topN:
			line = row.strip("\n").split("\t")
			chr, pos = line[0].split("_")
			#posList=positions.split(",")
			if float(line[2])==0.0:
				pTumor=laplace
			else:	
				pTumor=float(line[2])/float(total)
		
			if float(line[1])==0.0:
				pNormal=laplace
			else:
				pNormal=float(line[1])/float(total)
		

			pNormalDict[chr][pos]=pNormal
			pTumorDict[chr][pos]=pTumor

			if chr not in snvDict.keys():
				snvDict[chr]=1
			else:
				snvDict[chr]+=1

		else:
			break
	

	pN=0.0
	pT=0.0
	
	#print snvDict
	#print pNormalDict
	# Read patient/ sample variants
	# 11_2222
	# X_33
	# Nov 4, 2014
	timer=0
	for eachline in gzip.open(patientVariants, "r"):
		chr, pos = eachline.strip("\n").split("_")
		if chr in snvDict.keys():
			if pos in pNormalDict[chr].keys():
				pN=pN+math.log10(pNormalDict[chr][pos])
				pT=pT+math.log10(pTumorDict[chr][pos])
				timer+=1		
		else:
			pass

	if timer == 0:
		return "NA"	#no hit in the model	
	else:		
		return (pT-pN) #the log liklihood ratio, tumor vs normal, Dec 5, 2014
	


#Dec 24, 2014
#do not include top N lines of the model
def NBmodelRemove (modelFile, patientVariants, total, topN=20, laplace=float(0.5/1092)):

	#upLimit=1000#Dec 24, 2014, Jan9, 2015	
	i=0
	snvDict={}
	pNormalDict = defaultdict(lambda: defaultdict(float))
	pTumorDict = defaultdict(lambda: defaultdict(float))
	#print laplace

	for row in gzip.open(modelFile, "r"):
		#print row
		i = i + 1
		#print i
		if i <= topN:
			pass

		elif i > topN and i <= upLimit: #exclude top N lines of the model
			line = row.strip("\n").split("\t")
			chr, pos = line[0].split("_")
			#posList=positions.split(",")
			if float(line[2])==0.0:
				pTumor=laplace
			else:	
				pTumor=float(line[2])/float(total)
		
			if float(line[1])==0.0:
				pNormal=laplace
			else:
				pNormal=float(line[1])/float(total)
		

			pNormalDict[chr][pos]=pNormal
			pTumorDict[chr][pos]=pTumor

			if chr not in snvDict.keys():
				snvDict[chr]=1
			else:
				snvDict[chr]+=1

		elif i > upLimit:
			break
	

	pN=0.0
	pT=0.0
	
	#print snvDict
	#print pNormalDict
	# Read patient/ sample variants
	# 11_2222
	# X_33
	# Nov 4, 2014
	timer=0
	for eachline in gzip.open(patientVariants, "r"):
		chr, pos = eachline.strip("\n").split("_")
		if chr in snvDict.keys():
			if pos in pNormalDict[chr].keys():
				pN=pN+math.log10(pNormalDict[chr][pos])
				pT=pT+math.log10(pTumorDict[chr][pos])
				timer+=1		
		else:
			pass

	if timer == 0:
		return "NA"	#no hit in the model	
	else:		
		return (pT-pN) #the log liklihood ratio, tumor vs normal, Dec 5, 2014
	


def main():

	if InExclude == "take" :
		posterior=NBmodelInclude(modelFile, patientVariants, float(total), int(topN), laplace=float(0.5)/float(total))	
		print str(posterior)	
	elif InExclude == "remove":
		pass
		#posterior=NBmodelRemove(modelFile, patientVariants, float(total), int(topN), laplace=float(0.5)/float(total))	
		#print str(posterior)
	else:
		print "ERROR!"
	


if __name__=="__main__":
	main()

#A model file
#ChrXXX:111111,111111 Ntumor Ntotal backgroundRate CRscore
 
#Chr12:25398214,25398214,25398218,25398226,25398255,25398262,25398281,25398282,25398284,25398285	75	220	0.00454545	115.906452155921

