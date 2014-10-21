#!/usr/bin/python
# Oct 20, 2014
# Tian R. <tianremi@gmail.com>
# count how many patients are mutated in a locus within a cluster of SNVs


import re, sys

def unique1(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked




def SortedKey2GivenStartEndAsList(input, outfile):
	'''
	Input file is of two columns.
	The first col is mut IDs, mut be sorted. No redundancies.
	#100	A_start
	#100	P2
	#101	void
	#102	P3,P4
	#102	Z_end
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
			last_zuo=""
			last_you=current_you
		
		if last_you != "":
			#if last_you is not empty
			if current_you != "Z_end":
				if last_zuo != "":
					last_zuo=last_zuo+","+current_zuo
					last_you=last_you+","+current_you
				else:
					last_zuo=last_zuo+current_zuo
					last_you=last_you+","+current_you
			else:
	
				patlist=last_you.split(",")
		
				# patlist might contain "void" (empty) elements
				newlist=[]
	
				for e in patlist:
					if re.match("Patient", e):
						newlist.append(e) 
				out.write(last_zuo+"\t"+str(len(unique1(newlist)))+"\n")
	
				#empty the lists
				last_zuo=""
				last_you=""
	f.close()
	out.close()



def main():

	SortedKey2GivenStartEndAsList(sys.argv[1], sys.argv[1]+".out")


if __name__=="__main__":

	main()
