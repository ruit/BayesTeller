#!/usr/bin/python

# TianR.<tianremi@gmail.com>
# Sep 4, 2014
# @Sep5
###################################
#
#
#
#
#
####################################


import sys, re, subprocess
import random as rd
import time

# merger col2 sharing the same col1 values
# Aug 26, 2014 # 300K!!!!!!!!!!
# Give up DL, for publications!!!

def SortedKey2MultipleValuesAsList(input, outfile):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#Mu1	P1
	#Mu1	P2
	#Mu2	P3
	'''
	# June 19, 27, 2014
	last_mut=""
	current_mut=""
	last_donor=""
	current_donor=""
	
	f=open(input, "r")	

	out=open(outfile, "w")

	for line in f:
		array=line.strip("\n").split("\t")
		current_mut=array[0]
		current_donor=array[1]
		if last_mut=="":
			last_mut=current_mut
			last_donor=current_donor
		else:
			if current_mut != last_mut:
				out.write(last_mut+"\t"+ last_donor+"\n")

				last_mut=current_mut
				last_donor=current_donor
			else:
				last_donor=last_donor+","+current_donor

	out.write(last_mut + "\t" +last_donor+"\n")			

	f.close()
	out.close()



def AllPatientsSortedSNVs2Blocks (infile, span, outfile):
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
	
	blockdict={}
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
				#print m
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
			
	
	f_out=open(outfile, "w") 

	for key in sorted(snvdict.keys()):
		
		value=snvdict[key]
 		f_out.write(key+"\t"+value+"\n")
	
	f_out.close()


# results can be validated by more than 1 strategies are good.
# Any script should be tested in multiple question settings.
# Sep 1, 2014	
# For a given SNV block, it is possible that a given patient has been reported multiple times
# These repeated patients should be counted as one


def unique1(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked




def f3(seq):
   # Not order preserving
   dict = {}
   for e in seq:
       dict[e] = 1
   return dict.keys()





def CountNonRedundantPatients(infile, outfile):
	'''
	input file should be something like this:
	1st and 2rd is delimited by tab; then within 2rd column, each patient is seperated by comma
	------------------------------------------
	10_100008701	Patient216
	10_100010909	Patient78,Patient119
	10_100010921	Patient119,Patient78
	10_100011442	Patient12
	somewhere	void (for somatic SNV in germline counting)
	'''

	#Sep 15, 2014	

	patlist=[]
	f_out=open(outfile, "w")

	for line in open(infile, "r"):
		chrpos, patients=line.strip("\n").split("\t")

		patlist=patients.split(",")

		# patlist might contain "void" (empty) elements
	
		newlist=[]
		for e in patlist:	
			if e != "void":
				newlist.append(e)	
			
		f_out.write(chrpos+"\t"+str(len(unique1(newlist)))+"\n")
		
	
	f_out.close()



def CountNonRedundantPatientsTestRemoved(infile, outfile, testpatientList):
	'''
	input file should be something like this:
	1st and 2rd is delimited by tab; then within 2rd column, each patient is seperated by comma
	------------------------------------------
	10_100008701	Patient216
	10_100010909	Patient78,Patient119
	10_100010921	Patient119,Patient78
	10_100011442	Patient12
	somewhere	void (for somatic SNV in germline counting)
	#
	#Sep 17, 2014
	#exclude a list of patients ["Patient1", "Patient2"]
	#excluded patients are to be used as the test data

	'''

	#Sep 15, 2014	

	patlist=[]
	f_out=open(outfile, "w")

	for line in open(infile, "r"):
		chrpos, patients=line.strip("\n").split("\t")

		patlist=patients.split(",")

		# patlist might contain "void" (empty) elements
	
		newlist=[]
		for e in patlist:	
			if e != "void" and e not in testpatientList:
				newlist.append(e)	
			
		f_out.write(chrpos+"\t"+str(len(unique1(newlist)))+"\n")
		
	
	f_out.close()



def main():

	if len(sys.argv)==4:

		start=time.time()

		infile=sys.argv[1]
		outfile1=infile+"_"+str(rd.random())+".out1"
	
		print "@INFO Inputfile is "+infile
		print "@INFO Inputfile should be sorted sequentially by 1st and 2rd columns(numerically)."
	
		testSorted=subprocess.call('cat '+infile+' |sed "s/_/\t/g" |sort -c -k1,1 -k2,2 -n ', shell=True)

		if testSorted==0:
			#this file is sorted	
			SortedKey2MultipleValuesAsList(infile, outfile1)
		
			print "@INFO For a given SNV block, patients show mutations are gathered."		
			AllPatientsSortedSNVs2Blocks (infile=outfile1, span=int(sys.argv[2]), outfile=outfile1+".out2")

			print "@INFO For a given SNV block, number of unique patients are counted."
			print "@INFO The output is "+sys.argv[1]+".snvfreq"


			if int(sys.argv[3]) >0 and int(sys.argv[3]) <=99:

				r=int(sys.argv[3])
				# remove this proportion of patients randomly from both germ and tumor, so that use them for testing the model later
				# Sep 17, 2014
				temppatfile=sys.argv[1]+"_"+str(rd.random())+".temp"
				subprocess.call('cat '+ sys.argv[1]+' | cut -f2 | sort | uniq | grep -v "void" > temppatfile')
				
				testpatlist=[]
				for line in open("temppatfile", "r"):
					line=line.strip("\n")
					testpatlist.append(line)
				
				#shuffle the list of patients				
				newlist=random.shuffle(testpatlist)	
				
				outlist=[]
				for e in newlist:
					if newlist.index(e)<=len(newlist)*r*0.01-1:
						outlist.append(e)
				# take the first r% as the set of patients for testing	

				CountNonRedundantPatientsTestRemoved (infile=outfile1+".out2", outfile=sys.argv[1]+".snvfreq", testpatientList=outlist )
				
	
			elif sys.argv[3] == "0":
				print "Use all the patients to build the model"		
				CountNonRedundantPatients (infile=outfile1+".out2", outfile=sys.argv[1]+".snvfreq")
				end=time.time()
	
				print 'Cong! Total time elapsed: '+str((end-start)/60)+" min." 
				print
			else:
				print "WARNING! The 3rd parameter shoud be within 0-100 integers!"
			#subprocess.call('rm *.out1')
			#subprocess.call('rm *.out2')
		else:
			print "@ERROR The input file should be sorted in Linux first!"
			pass

	else:
		print
		print sys.argv[0]+" needs 3 arguments:"
		print "ARGV[1], input file: a sorted catenated all tumor OR germline SNVs."
		print "ARGV[2], left and right span length in bps, default value is 500."
		print "ARGV[3], proportion of patients to be excluded from building the model"


if __name__=="__main__":
	main()
						
