#!/usr/bin/python

#Tian R. <tianremi@gmail.com>
#Feb 27, 2015 
#March 18b, 2015

import sys, gzip, re
import math

'''

'''

if len(sys.argv)<6:
	raise SystemExit("I need arguments!")


try:
	setSNV=sys.argv[1] #the list of selected SNVs
	sampleFile=sys.argv[2] #unknow sample of SNVs, gzipped
	germFreq=sys.argv[3] # freq in panCancer germ trained
	tumorFreq=sys.argv[4] # freq in panCancer tumor trained
	sampleSize=float(sys.argv[5]) # 3909 or half, March 18, 2015
	cutoff=float(sys.argv[6]) #take top 0.3 (30%) for sum of the snvlogRatio
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
logR_somaLike=[]
logR_germLike=[]

if len(hit)==0:
	logR=[]
	logR_somaLike=[]
	logR_germLike=[]	
else:

	#March 18, 2015, if  not hits, assign zero
	logR=[math.log10((0.5+Dict_freq_tumor_A.get(snv, 0))/(0.5+Dict_freq_germ_A.get(snv, 0))) \
for snv in hit] # zero!
	
	for snv in hit:
		if Dict_freq_germ_A.get(snv, 0) <= (sampleSize * 0.01): #1% as the cutoff, March 18, 2015	
			logR_somaLike.append(math.log10((0.5+Dict_freq_tumor_A.get(snv, 0))/(0.5+Dict_freq_germ_A.get(snv, 0))))# zero!
		else:
			logR_germLike.append(math.log10((0.5+Dict_freq_tumor_A.get(snv, 0))/(0.5+Dict_freq_germ_A.get(snv, 0))))# zero!





def pred(hit, all_logR, logR, cutoff, sampleFile, suffix):
	keyR=[]
	countedLoci=[]
	l=len(logR)
	#logR.sort(reverse=True)
	for ratio in logR:
		if ratio>=cutoff:
			keyR.append(ratio)
			countedLoci.append(hit[all_logR.index(ratio)]) #March 6, 2015, keep record of counted top SNVs

	#keyR=logR[0:int(topR*l)] # take top 30% hits
	
	fh_out=open(sampleFile+suffix, "w")
	if len(keyR)==0: fh_out.write("0_0"+"\t"+"NA"+"\n")
	else:
	
		val=reduce(lambda x,y: x+y, keyR)
	#val=float(val/(l*topR))
		fh_out.write(str(l)+"\t"+str(len(keyR))+"\t"+str(val/l)+"\t"+";".join(countedLoci)+"\n")
#":"+";".join(countedLoci)+

		keyR=[]
	fh_out.close()



if len(logR)==0:
	fh_out=open(sampleFile+".all", "w")
	fh_out("0"+"\t"+"NA"+"\n")
	fh_out.close()
else:

	pred(hit, logR, logR, cutoff, sampleFile, ".all")
	pred(hit, logR, logR_somaLike, cutoff, sampleFile, ".sml")
	pred(hit, logR, logR_germLike, cutoff, sampleFile, ".gml")
