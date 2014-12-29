#Dec 26, 2014
#TianR.
#dependancy: ROCR


args<-commandArgs(TRUE)

library(ROCR)

calAUC<-function(tenCombinedFile, posiLabel="t", negaLabel="g"){
 	data<-read.table(tenCombinedFile, header=F) # like "10.sum"
 	
	labels<-as.vector(data[,1])
 	labels<-replace(labels, labels=="t", "1")
 	labels<-replace(labels, labels=="g", "0")
 	labels<-factor(labels)

	predictions<-data[,2]

	 pred<-prediction(predictions, labels)
	
	 perf2=performance(pred, "auc")
	 auc=unlist(slot(perf2, "y.values"))
	
	return (auc)
	}

calAUC(args[1])
