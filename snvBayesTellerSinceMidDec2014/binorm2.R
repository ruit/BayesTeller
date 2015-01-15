# Tian R. <tianremiATgmail.com>
# Sep 29, 2014

# calculate 1-accum(k-1,N,p), logP
# Laplacician correction 5e-4

#10_100008701	1	220	0
#10_100010909	2	220	0.00454545
#10_100010921	2	220	0.00454545

# Jan 14, 2015


args<-commandArgs(TRUE)

compbinormP<-function (input=input, p=0.5/2000){
	mat<-read.table(input, header=F)
	head(mat)
	orimat<-mat

	#p<-0.000001
	mat[mat[,4]==0,4] <- p
	
	#compute the K-1 cumulative p

	cumuP_Kminus1<-pbinom((mat[,2]-1), mat[,3], mat[,4], log.p=T, lower.tail=F)/log(10)
	#compute -log10 of P, x>=K
	CRscore<--1*(cumuP_Kminus1)

	outputdf<-data.frame(orimat, CRscore=CRscore)

	write.table(outputdf, file=paste(input,".pval", sep=""), sep="\t", 
quote=F, row.names=F, col.names=F)
	}


input<-args[1]
#p<-args[2]

compbinormP(input)


