#!/usr/bin/python
# Dec 5, 2014
# TianR.
#t 300.0
#t 192.0
#t 94.0
#g 72.0
#g 66.0

import sys


def uniq(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked



def readTwoCol(inputFile):
	"""
	input file, no header, 2 cols, sorted numerically by col2, -r
	"""	  
	label=[]
	score=[]
	for line in open(inputFile, "r"):
		(x,y)=line.strip("\n").split("\t")	
		#if x !="" and y !="":
		label.append(x)
		score.append(float(y))

	return ([label, score])


def setCutOff(label, score, zheng="t", fu="g"):

	#tuple is not mutable, but list is!!! Dec 8, 2014	
	labelOri=tuple(label)

	zhengList=[]
	fuList=[]

	for instance in labelOri:
		if instance==zheng:
			zhengList.append(instance)
		elif instance==fu:
			fuList.append(instance)	

	sensiList=[]
	One_speciList=[]
	for num in uniq(score):
		for i in range(len(score)):
			if score[i] >= num:
				label[i]=zheng#tumor, positive instance
			else:
				label[i]=fu #germline, negative instance
		
		#compare label and labelOri
	#	print label	
	#	print labelOri
		posiposi=0.0
		negneg=0.0
		sensi=0.0
		One_speci=0.0
		for j in range(len(score)):
			
			if labelOri[j]==zheng and label[j]==zheng:
				#true positive
				posiposi += 1
			#elif labelOri[j]==zheng and label[j]==fu:
				#false negatvie
			#	posineg += 1
			elif (labelOri[j]==fu and label[j]==fu):
				#true neg
				negneg += 1
			#elif labelOri[j]==fu and label[j]==zheng:
				#false positive
			#	negposi += 1
	
		
	#	print posiposi, negneg, len(zhengList), len(fuList)	
			
		sensi=round(posiposi/len(zhengList), 3)
		One_speci=round((1 - (negneg/len(fuList))),3)

		#print sensi
		sensiList.append(sensi)
		One_speciList.append(One_speci)

	return [One_speciList, sensiList]



def main():
	infile=sys.argv[1]
	list=readTwoCol(infile)
	out=setCutOff(label=list[0], score=list[1], zheng="t", fu="g")

	print "one_Speci"+"\t"+"Sensi"	
	for i in range(len(out[0])):
		print str(out[0][i])+"\t"+str(out[1][i])


if __name__=="__main__":
	main()
