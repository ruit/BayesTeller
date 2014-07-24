#/usr/bin/python
#Tian X.R., Feb 4, 2014
#tianremi@gmail.com
#
#Feb.5, 2014
#Feb. 10, 2014
#Feb. 17, 2014
#March 2, 2014 
#-----------------------------------------------------------

#input1#, the Contrast Ranking results list_colorectal_5.txt
#GENE	ID_TUM	GLEN	%MUT_TUM	Z-SCORE_TUM	P-VAL_TUM	ID_NOR	%MUT_NOR	Z-SCORE_NOR	P-VAL_NOR	DIF_%MUT	DIF_Z-SCOR
#KRAS	96	570	  43.636	  127.55	0.00e+00	2	   0.909	    1.53	1.45e-01	0.427273	126.02

#input2#, the gene coordinate in the genome, hg19_refGene_annovar.txt
#661	NM_022787	chr1	+	10003485	10045556	10032131	10042759	5	10003485,10032075,10035649,10041088,10042358,	10003573,10032246,10035833,10041228,10045556,	0	NMNAT1	cmpl	cmpl	-1,0,1,2,1,

#168	NM_017734	chr1	+	100111430	100160097	100111873	100159618	8	100111430,100127874,100133197,100152231,100152485,100152631,100154330,100159574,	100111918,100127955,100133322,100152346,100152519,100152745,100155428,100160097,	0	PALMD	cmpl	cmpl	0,0,0,2,0,1,1,1,


#input3#, the size of each human chrom, hg19
#hg19.chrom.siezes
#chr1	249250621
#chr2	243199373


#input4#
#/Users/tianr/1Rui/Constants/tcga/colon/gene.germline.mut.freq.patients
#ZSCAN20	217
#ZBTB8A	217


#----------------------------------------------------------
#output

#dependent on scipy
import sys, re, operator

	
#declare constants
chrom_size_file="/Users/tianr/1Rui/Constants/hg19.chrom.sizes"
ref="/Users/tianr/1Rui/Constants/hg19_refGene_annovar.txt"

if len(sys.argv) >= 3:
	drivergenelist=sys.argv[1]
	#the flanking range, eg., 100kb
	flanking=int(sys.argv[2])*1000
	mutbg=sys.argv[3]
	#mutumor=sys.argv[4]

else:
	print "------------------------------------------------------------------------------------------------------"
	print
	print "Usage: python nbr.python <driver gene list> <flanking range (Kb)><mutation frequencies in patients file"
	print "driver gene list the output of EC's contrast ranking method for prioritizing cancer associated genes"
	print "flanking range is ## kbs, e.g."
	print
	print "------------------------------------------------------------------------------------------------------"

def ChromSizeDict (hg19_file="hg19.size.file"):

	'''
	given hg19 chrom size file, return a dict
	'''
	chromd={}
	lines=open(hg19_file, "r").readlines()
	for line in lines:
		(chrom, size)=line.strip("\n").split("\t")
		chromd[chrom]=int(size)
	return chromd

#print ChromSizeDict(sys.argv[1])

def GenePosDict (ref="hg19_refGene_annovar.txt"):
	'''
	given gene, return its coordinate in the genome
	'''
	
	genedict={}
	lines=open(ref, "r").readlines()
	for line in lines:
		(num, refseq, chrom, strand, TSS, TTS, CDSstart, CDSend, num_exon, what10, what11, what12, symbol, what14, what15, what16)=line.strip("\n").split("\t")	
		genedict[symbol]=[chrom, TSS, TTS]# what is a symbol corresponds to multiple locations?
	return genedict



def nearby(gene, flanking, chromdict, geneposdict):
	'''
	gene symbol, flanking in bps
	'''
	neargenelist=[]
	gene_pos_list=[]
	value=[]
	(da, xiao)=(0,0)
	if not gene in geneposdict.keys():
		print "Error,genes not in the position dicts"
		return "NULL"
	else:
		
		gene_pos_list=geneposdict[gene]
		#gene symbol, chr, tss, tts
		(chrom, tss, tts) = (gene_pos_list[0], int(gene_pos_list[1]), int(gene_pos_list[2]))
		if tss < tts:
		# strand +
			start = tss
			end = tts
		else:
			start = tts
			end = tss
		up = (start - flanking) 
		down = (end + flanking)
		if up < 0:
			up=0
		if down > int(chromdict[chrom]):
			down= int(chromdict[chrom])
		for key in geneposdict.keys():
			if not geneposdict[key][0]==chrom:
				pass
			else:
				value=geneposdict[key]
				if int(value[1]) < int(value[2]):
					xiao=int(value[1])
					da=int(value[2])
				else:
					xiao=int(value[2])
					da=int(value[1])
				if xiao >= up and da <= down and key!= gene:
					#print key
					neargenelist.append(key)
		return neargenelist							
 
def main():
	
	#the boundary of each chrom
	chromdict=ChromSizeDict(chrom_size_file)
	
	#print chromdict
	geneposdict=GenePosDict(ref)
	
	#gene mut freq in patients
	mutfreqdict=ChromSizeDict(mutbg)	

	array=[]
	query=""
	mutdata=""
	lines=open(drivergenelist, "r").readlines()
	for line in lines:
		if re.match ("^#", line):
			pass
		else:
			array=line.strip("\n").split("\t")
			query=array[0]
			query_linju=nearby(query,flanking,chromdict, geneposdict)
			mutdata=""	
			for nr in range(len(query_linju)):
				if nr != len(query_linju) - 1:
					try:
						mutdata=mutdata+query_linju[nr]+":"+str(mutfreqdict[query_linju[nr]])+","
					except KeyError:
						mutdata=mutdata+query_linju[nr]+":"+"NA"+","
				else:
					try:
						mutdata=mutdata+query_linju[nr]+":"+str(mutfreqdict[query_linju[nr]])
					except KeyError:
						mutdata=mutdata+query_linju[nr]+":"+"NA"
			if query in mutfreqdict.keys():
				print query+":"+str(mutfreqdict[query])+ "\t" +mutdata
			else:
				print query+":NA\t" +mutdata

			
				

if len(sys.argv) < 3:
	pass
else:	
	main()
