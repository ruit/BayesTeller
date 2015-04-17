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


cal_everything<-function(tenCombinedFile, posiLabel="t", negaLabel="g"){
 	data<-read.table(tenCombinedFile, header=F) # like "10.sum"
 	
	labels<-as.vector(data[,1])
 	labels<-replace(labels, labels=="t", "1")
 	labels<-replace(labels, labels=="g", "0")
 	labels<-factor(labels)

	predictions<-data[,2]

	pred<-prediction(predictions, labels)

	#rocperf <- performance(pred, "tpr","fpr")
	#macperf<- performance(pred,"mat",'acc')
	#prperf<- performance(pred,"prec",'rec')
	#negperf<- performance(pred,"spec",'npv')
	
	#area under ROC curve
	aucperf<- performance(pred,"auc")
	auc<-aucperf@y.values[[1]]
	
	#ths=rocperf@alpha.values[[1]]

	#acc
	macperf<- performance(pred,"mat",'acc')
	mac<-macperf@y.values[[1]]
	acc<-macperf@x.values[[1]]
	
	#ROC curve
	rocperf <- performance(pred, "tpr","fpr")
	tpr<-rocperf@y.values[[1]]
	fpr<-rocperf@x.values[[1]]
	ths=rocperf@alpha.values[[1]]

	#pr
	prperf<- performance(pred,"prec",'rec')
	prec<-prperf@y.values[[1]]
	rec<-prperf@x.values[[1]]
	fdr<-1-prec	

	#neg
	negperf<- performance(pred,"spec",'npv')
	spec<-negperf@y.values[[1]]
	npv<-negperf@x.values[[1]]
	
	#fdr=1-prec
	
	#cat(c("AREA: ",auc,"\n"))
	#data<-cbind(ths,acc,fpr,tpr,mac,rec,prec,spec,npv,fdr,rep(auc,length(ths)))
	data<-cbind(fpr,tpr,fdr,rep(auc,length(ths)))
	colnames(data)<-c("FPR", "TPR", "FDR", "AUC")
	#row<-paste(c("TH","Q2","FPR","TPR","C","RECALL","PRECISION","SPEC","NPV","FDR","AUC"),collapse='\t')
	
	#cat(c(row,'\n'))
	#for (i in 3:length(ths)-1){
  	#	row<-paste(data[i,],collapse='\t')
  	#	cat(c(row,'\n'))
	#	}	

	
	#perf2=performance(pred, "auc")
	#auc=unlist(slot(perf2, "y.values"))
	
	return (data)
	}

#calAUC(args[1])
out<-cal_everything(args[1])
write.table(out, file=paste(args[1], ".p", sep=""), quote=F, row.names = F)
