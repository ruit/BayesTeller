#!/usr/bin/python
# Nov 3, 2014
# Nov 4, 2014
# @Dec 5, 2014
#X_11
#1_222



import re, math, sys
from collections import defaultdict


modelFile = sys.argv[1]
patientVariants = sys.argv[2]
topN = sys.argv[3]



def NBmodel (modelFile, patientVariants, topN=20, laplace=float(0.5/1092)):
	
	i=0
	#pNormalDict={}
	#pTumorDict={}
	snvDict={}
	pNormalDict = defaultdict(lambda: defaultdict(float))
	pTumorDict = defaultdict(lambda: defaultdict(float))
	#print laplace

	for row in open(modelFile, "r"):
		if re.search("^Chr", row):
			#print row
			i = i + 1
			#print i
			if i <= topN:
				line = row.strip("\n").split("\t")
				chr, positions = line[0].split(":")
				posList=positions.split(",")
				
				pTumor=float(line[1])/float(line[2])
			
				if float(line[3])==0.0:
					pNormal=laplace
				else:
					pNormal=float(line[3])
			

				for pos in posList:	
					pNormalDict[chr][pos]=pNormal
					pTumorDict[chr][pos]=pTumor

				if chr not in snvDict.keys():
					snvDict[chr]=[posList]
				else:
					list=snvDict[chr]
					list.append(posList)
					snvDict[chr]=list

			else:
				break
		

	pN=0.0
	pT=0.0
	
	countCluster = defaultdict(lambda: defaultdict(int))

	# Read patient/ sample variants
	# 11_2222
	# X_33
	# Nov 4, 2014

	for eachline in open(patientVariants, "r"):
		num, pos = eachline.strip("\n").split("_")
		chr="Chr"+num
		#print chr+"\t"+pos
		
		if chr in snvDict.keys():
			for j in range(len(snvDict[chr])):
				if pos in snvDict[chr][j]:
					countCluster[chr][j] = countCluster[chr][j] + 1 ### snv in the same chr same cluster should be added only once for the prob.

	timer = 0
	# for all chr and corresponding all clusters match the sample, sum the prob
	#print countCluster	
	for chr in countCluster.keys():
		for k in countCluster[chr].keys():
			pos = snvDict[chr][k][0]
			#print pos
			
			timer = timer + 1		
			# add these prob, do not repeat
			pN=pN+math.log10(pNormalDict[chr][pos])	
			pT=pT+math.log10(pTumorDict[chr][pos])
	if timer == 0:
		return "NA"	
	else:		
		return round(pT-pN) #the log liklihood ratio, tumor vs normal, Dec 5, 2014
	

def main():

	posterior=NBmodel(modelFile, patientVariants, int(topN), laplace=float(0.5/1092))	
	print str(posterior)


if __name__=="__main__":
	main()

#A model file
#ChrXXX:111111,111111 Ntumor Ntotal backgroundRate CRscore
 
#Chr17:7578177,7578177,7578190,7578208,7578211,7578212,7578217,7578235,7578253,7578257,7578263,7578280,7578369,7578370,7578388,7578394,7578403,7578406,7578407,7578442,7578445,7578449,7578454,7578455,7578457,7578461,7578503,7578526,7578550	52	220	0	120.63978567313
#Chr12:25398214,25398214,25398218,25398226,25398255,25398262,25398281,25398282,25398284,25398285	75	220	0.00454545	115.906452155921

