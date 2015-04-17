# Tian R. <tianremiATgmail.com>
# Sep 29, 2014

# calculate 1-accum(k-1,N,p), logP
# Laplacician correction 5e-4

#10_100008701	1	220	0
#10_100010909	2	220	0.00454545
#10_100010921	2	220	0.00454545



#new inputfile
#10_00000 ptumor pgermline



args<-commandArgs(TRUE)

#Feb 9, 2015


Num=as.integer(args[2])

pseudo=as.numeric(args[3])

#pseudo=0.0005

data<-read.table(args[1], header=F)

oridata<-data

#germline freq
data[data[,3]==0,3]<-pseudo
	
#compute the K-1 cumulative p

cumuP_Kminus1<-pbinom((round(Num*data[,2])-1), Num, data[,3], log.p=T, lower.tail=F)/log(10)

#compute -log10 of P, x>=K
CRscore<--1*(cumuP_Kminus1)

outputdf<-data.frame(oridata, CRscore=CRscore)

write.table(outputdf, file=paste(args[1],".pval", sep=""), sep="\t", 
quote=F, row.names=F, col.names=F)
