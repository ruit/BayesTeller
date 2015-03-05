#!/usr/bin/python

#Tian R. <tianremi@gmail.com>
#Feb 27, 2015 

import sys, gzip, re
import math

'''

'''

if len(sys.argv)<5:
	raise SystemExit("I need arguments!")


try:
	setSNV=sys.argv[1] #the list of selected SNVs
	sampleFile=sys.argv[2] #unknow sample of SNVs, gzipped
	germFreq=sys.argv[3] # freq in panCancer germ trained
	tumorFreq=sys.argv[4] # freq in panCancer tumor trained
	cutoff=float(sys.argv[5]) #take top 0.3 (30%) for sum of the snvlogRatio
	cutoff=math.log10(cutoff)
	#print sys.argv

except:
	pass


#read SNV pool
snvList=[snv.strip("\n") for snv in open(setSNV, "r")]



def getWeights(lista, element):

	return float(len(lista)-lista.index(element))/len(lista)




fh=gzip.open(sampleFile, "rb")
hit=[]

for line in fh:
	if re.search("^#", line): pass
	else:
		array=line.split("\t")
		chr_pos=array[0]+"_"+array[1] # chr_pos
		if chr_pos in snvList:
			hit.append(chr_pos)
		else:
			pass			
			
fh.close()


#read freq in training as dicts
Dict_freq_tumor_A={}
Dict_freq_germ_A={}

for line in open(germFreq, "r"):
	x,y=line.strip("\n").split("\t")
	Dict_freq_germ_A[x]=float(y)

for line in open(tumorFreq, "r"):
	x,y=line.strip("\n").split("\t")
	Dict_freq_tumor_A[x]=float(y)



logR=[]
if len(hit)==0:
	logR=[]
else:
	logR=[math.log10((0.5+Dict_freq_tumor_A.get(snv))/(0.5+Dict_freq_germ_A.get(snv))) \
for snv in hit] # zero!
	


if len(logR)==0:
	print "0"+"\t"+"NA"
else:
	keyR=[]
	l=len(logR)
	logR.sort(reverse=True)
	for ratio in logR:
		if ratio>=cutoff:
			keyR.append(ratio)

	#keyR=logR[0:int(topR*l)] # take top 30% hits
	if len(keyR)==0: print "0"+"\t"+"NA"
	else:
	
		val=reduce(lambda x,y: x+y, keyR)
	#val=float(val/(l*topR))
		print str(l)+"\t"+str(val/l)
		keyR=[]



