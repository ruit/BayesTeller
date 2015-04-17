#!/usr/bin/python

# Tian R. <tianremi@gmail.com>

# For a given SNV block, it is possible that a given patient has been reported multiple times
# These repeated patients should be counted as one

import sys

def unique1(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked



def f3(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()


def CountNonRedundantPatients(infile):
	'''
	input file should be something like this:
	1st and 2rd is delimited by tab; then within 2rd column, each patient is seperated by comma
	------------------------------------------
	10_100008701	Patient216
	10_100010909	Patient78,Patient119
	10_100010921	Patient119,Patient78
	10_100011442	Patient12
	'''

	

	patlist=[]

	for line in open(infile, "r"):
		chrpos, patients=line.strip("\n").split("\t")

		patlist=patients.split(",")
		
		print chrpos+"\t"+str(len(unique1(patlist)))


def main():
	CountNonRedundantPatients (sys.argv[1])


if __name__=="__main__":
	main()
						
