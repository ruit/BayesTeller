#TianR.
#Dec 8, 2014


args<-commandArgs(TRUE)

plotROC<-function(cancerType="acc", suffix=".perf"){
	pdf(paste(cancerType, ".pdf", sep=""))
	fileArr=dir(pattern=cancerType)
	fileL=dir(pattern=suffix)
	fileArr<-intersect(fileArr, fileL)

	perf<-read.table(fileArr[1], header=T)

	plot(perf, type="b", lwd=2, col="blue", main=paste(cancerType,"10CV",sep=""))

	for (i in 2:length(fileArr)){perf<-read.table(fileArr[i], header=T)
		points(perf)

		}
	
	dev.off()
	}


plotROC(args[1])

