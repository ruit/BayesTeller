#!/usr/bin/python
# Tian R. <tianremiATgmail.com>
# Nov 4, 2014
# Nov 5, 2014


# Randomly mask some patients as test data

import random, math
import sys


infile1 = sys.argv[1]
infile2 = sys.argv[2]
nTest = int (sys.argv[3])
N = int (sys.argv[4])
testPatients = sys.argv[5]


# Sep 19, 2014
# Tian R. compare SNVs freq in normal and tumor, build a pridictive model
# Sep 22, 2014. 1000G SNV freq <= 0.5% as the filtering at the beginning!!!

def GenerateTestList (nTest, N):

	testList = [] # the list of test patients IDs.	
	
	if float(nTest)/float(N) >0.0 and float(nTest)/float(N) <1 :
		#print N
		totalList=range(1, N+1)		
		
		#print totalList
		while len(testList) < nTest:
			num=random.choice(totalList)
			pID="Patient"+str(num)
			if pID not in testList:
				testList.append(pID)
	
	return testList	
	
	


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
					if re.match("Patient", e) and e not in testList:
						newlist.append(e) 
				out.write(last_zuo+"\t"+str(len(unique1(newlist)))+"\n")
	
				#empty the lists
				last_zuo=""
				last_you=""
	f.close()
	out.close()




def main():
	testList = GenerateTestList(nTest, N)

	f_out = open (testPatients, "w")
	
	for p in testList:
		f_out.write(p+"\n")
	f_out.close()
	
	SortedKey2GivenStartEndAsList(infile1, infile1+".out", testList)
	
	SortedKey2GivenStartEndAsList(infile2, infile2+".out", testList)
	
	



if __name__=="__main__":
	main()

