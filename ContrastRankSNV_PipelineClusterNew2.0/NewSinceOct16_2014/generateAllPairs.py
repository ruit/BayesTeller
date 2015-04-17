#!/usr/bin/python

import sys

#17	1638901	3	220	0	3.69516922747676
#L=[1,2,3,4]

def all_pairs(L):
	while L:
		i = L.pop()
		for j in L: print str(i)+"\t"+str(j)


#all_pairs(L)

def main():

	pos=[]
	#count=0
	for line in open(sys.argv[1], "r"):
		line=line.strip("\n").split("\t")
		pos.append(line[1])


	all_pairs(pos)



if __name__=="__main__":
	main()
