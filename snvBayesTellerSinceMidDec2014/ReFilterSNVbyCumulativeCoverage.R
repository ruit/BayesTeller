# Jan 9, 2015 TianR
# Given SNV sort uniq -c file, generate cumulative coverage of SNV, return alpha cutoff 

args=(commandArgs(TRUE))

#requires 3 args
# filename of uniqCfile
# alpha 0.05, 0.01, 0.005, ...
# outputfile name
 
cumulative<-function(uniqCfile="f2"){
    f2<-read.table(uniqCfile, header=F)
    total<-sum(f2[,1])
    what<-c(0, as.vector(f2[-nrow(f2),1]))
    
    for (i in 1:nrow(f2)){
        what[i]<-(total-what[i])
        total<-what[i]
	}
    
    out<-data.frame(numPat=f2[,2], ratio=round(what/what[1], 3))
    	
    return (out)
    }

 

patNrinTumor<-function(uniqCfile="f2", alpha=0.05){
    cumu<-cumulative(uniqCfile)
    N<-head(cumu[cumu$ratio <= alpha,],1)$numPat
    return (N)
    }


num<-patNrinTumor(args[1], args[2])

write.table(num, file=args[3], sep="\t", 
    quote=F, row.names=F, col.names=F)

