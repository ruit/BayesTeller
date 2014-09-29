# Tian R. <tianremiATgmail.com>
# Sep 29, 2014

# calculate 1-accum(k-1,N,p), logP
# Laplacician correction 5e-4

#10_100008701	1	220	0
#10_100010909	2	220	0.00454545
#10_100010921	2	220	0.00454545

args<-commandArgs(TRUE)

pseudo=0.0005

data<-read.table(args[1], header=F)

data[data[,4]==0,4]<-pseudo
	
#compute the K-1 cumulative p

cumuP_Kminus1<-pbinom(
(data[,2]-1), data[,3], data[,4])

#compute -log10 of P, x>=K
CRscore<--log10(1-cumuP_Kminus1)

outputdf<-data.frame(data, CRscore=CRscore)

write.table(outputdf, file=paste(args[1],".pval", sep=""), sep="\t", 
quote=F, row.names=F, col.names=F)
