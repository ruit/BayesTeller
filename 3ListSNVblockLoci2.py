#!/usr/bin/python

# TianR.<tianremi@gmail.com>

# Aug 31, 2014

# input file: colon tumor all possible SNVs, sorted, uniq

# output: Chr11_10000_bloc: 11_10000, 11_10001,...

# a typcial input file is about 5 Mb
# each chrom at least 10 min

import sys

snvdict={}

blockdict={}


def ReadAllPossibleSNVs (infile, span=100):
	'''
	'''
	for line in open (infile, "r"):
		
		chr, pos=line.strip("\n").split("\t")
		
		key=chr+"_"+pos

		snvdict[key]=chr+":"+ pos

#	print snvdict	

	for key, value in snvdict.iteritems():
		
		chr, pos = key.split("_")
		
		blockdict[key]=value
		
		for snvar in snvdict.values():	
			snv=snvar.split(":")
			if snv[0]==chr and int(pos)-span<=int(snv[1])<=int(pos)+span and snv[1]!=value.split(":")[1]:
				blockdict[key]=blockdict[key]+","+snvar
#				print snv

	for key, value in blockdict.iteritems():
		print key+"\t"+value



def FromSortedSNVs2Blocks (infile, span=500):
	'''
	input file should be sorted by chr and by position numerically
	'''

	chrlist=[]
	snvdict={}

	nearlist=[]
	uplist=[]
	
	for line in open (infile, "r"):
		
		line=line.strip("\n")
		chrlist.append(line)
	
	#loop the list to compare each position with left and right neighbours
	for i in range(len(chrlist)):

		chr, pos=chrlist[i].split("\t")
		
		key=chr+"_"+pos
		snvdict[key]=chr+":"+pos
		

		if i==0:
			m=1
			nearlist=chrlist[i+m].split("\t")
			
			while 0 < int(nearlist[1])-int(pos) <= span:
				#print m
				#make sure the chr is the same, be cautious with two chrs boundaries
				if nearlist[0]==chr:
					snvdict[key]=snvdict[key]+","+nearlist[0]+":"+nearlist[1]
					m=m+1
					nearlist=chrlist[i+m].split("\t")
				else:
					break

				
		else:			
			
	
			n=1# scan lines on the top
			uplist=chrlist[i-n].split("\t")
			
			while 0< int(pos)-int(uplist[1]) <=span:
				#print n
				if (uplist[0])==chr:
					snvdict[key]=snvdict[key]+","+uplist[0]+":"+uplist[1]
					n=n+1 	
					if (i-n)>=0:
						uplist=chrlist[i-n].split("\t")	
					else:
						break
				else:
					break
			
			m=1# scan lines blow
			
			if i !=len(chrlist)-1:
				nearlist=chrlist[i+m].split("\t")
			
				while 0 < int(nearlist[1])-int(pos) <= span:
				#print m
					if nearlist[0]==chr:					
						snvdict[key]=snvdict[key]+","+nearlist[0]+":"+nearlist[1]
						m=m+1
						if (i+m) <= (len(chrlist)-1):
							nearlist=chrlist[i+m].split("\t")
						else:
							break
					else:
						break
			else:
				break
			
 

	for key, value in snvdict.iteritems():
 		print key+"\t"+value


def AllPatientsSortedSNVs2Blocks (infile, span=500):
	'''
	input file should be sorted by chr and by position numerically
	An example input file should be like this:
	----------------------------------------------------------------

	X_155386	Patient158,Patient167,Patient179,Patient202,Patient218
	X_173062	Patient175
	X_174296	Patient11,Patient132,Patient7	
	'''

	chrlist=[]
	snvdict={}

	nearlist=[]
	uplist=[]
	
	for line in open (infile, "r"):
		
		line=line.strip("\n")
		chrlist.append(line)
	
	#print chrlist	

	#loop the list to compare each position with left and right neighbours
	for i in range(len(chrlist)):

		chrpos, pat=chrlist[i].split("\t")
		
		chr, pos=chrpos.split("_")
		
		snvdict[chrpos]=pat
		key=chrpos		

		if i==0:
			m=1
			nearlist=chrlist[i+m].split("\t")[0].split("_")
			#print nearlist	
			#print nearlist[1]
			#print nearlist[0]
			while 0 < int(nearlist[1])-int(pos) <= span:
				print m
				#make sure the chr is the same, be cautious with two chrs boundaries
				if nearlist[0]==chr:
					snvdict[key]=snvdict[key]+","+chrlist[i+m].split("\t")[1]
					m=m+1
					nearlist=chrlist[i+m].split("\t")[0].split("_")
				else:
					break

				
		else:			
			
	
			n=1# scan lines on the top
			uplist=chrlist[i-n].split("\t")[0].split("_")
			#print uplist	
			while 0< int(pos)-int(uplist[1]) <=span:
				#print n
				if (uplist[0])==chr:
					snvdict[key]=snvdict[key]+","+chrlist[i-n].split("\t")[1]
					n=n+1 	
					if (i-n)>=0:
						uplist=chrlist[i-n].split("\t")[0].split("_")	
					else:
						break
				else:
					break
			
			m=1# scan lines blow
			
			if i !=len(chrlist)-1:
				nearlist=chrlist[i+m].split("\t")[0].split("_")
			#	print nearlist	
				while 0 < int(nearlist[1])-int(pos) <= span:
				#print m
					if nearlist[0]==chr:					
						snvdict[key]=snvdict[key]+","+chrlist[i+m].split("\t")[1]
						m=m+1
						if (i+m) <= (len(chrlist)-1):
							nearlist=chrlist[i+m].split("\t")[0].split("_")
						else:
							break
					else:
						break
			else:
				break
			
 

	for key, value in snvdict.iteritems():
 		print key+"\t"+value


def main():

# This repeated looping requires a lot of memory and it is slow.
# results can be validated by more than 1 strategies are good.
# Any script should be tested in multiple question settings.
# Sep 1, 2014	
#	ReadAllPossibleSNVs(infile=sys.argv[1], span=int(sys.argv[2]))

#	FromSortedSNVs2Blocks (infile=sys.argv[1], span=int(sys.argv[2]))
	AllPatientsSortedSNVs2Blocks (infile=sys.argv[1], span=int(sys.argv[2]))

if __name__=="__main__":

	main()
