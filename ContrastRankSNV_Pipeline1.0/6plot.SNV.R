#Tian R. <tianremi@gmail.com>
#plot SNV freq over colon cancer 220 patients
#filter germline common variants cutoff 0.5%
#

#

# re cat tumor=germline + somatic do it in linux





Sum1stColbyPos<-function(df1, df2){
	
	#df1: 
	#1 10000
	#3 12000

	#df2:similar

	colnames(df1)<-c("freq", "pos")
	colnames(df2)<-c("freq", "pos")
	
	dfcomb<-rbind(df1, df2)
	
	positions<-unique(sort(dfcomb$pos))
	
	new.mat<-matrix(0, length(positions), 2)
	
	for (i in 1:length(positions)){
		new.mat[i,2]<-positions[i]
		new.mat[i,1]<-sum(dfcomb[dfcomb$pos==positions[i],1])
		}

	return (as.data.frame(new.mat))

	}



#mat1<-matrix(c(1,100,2,150), nrow=2, ncol=2, byrow=T)
#mat2<-matrix(c(1,101,2,150), nrow=2, ncol=2, byrow=T)


PlotSNV<-function (plotName, N=220, p=0.005, upper=0.4){
	
	pdf(plotName)
	
	#N is sample size of the cancer patients cohort
	#p is the cutoff for filtering common variants/ polymorphism

	for (i in c(1:22,"X")){
		gdata<-read.table(paste("colon.germline.stacked.sorted.", paste(i, ".sum", sep=""), sep=""), header=F)
		tdata<-read.table(paste("colon.tumor.stacked.sorted.", paste(i, ".sum", sep=""), sep=""), header=F)

		#tumor<-Sum1stColbyPos(gdata, sdata)
	
		#plot(gdata[gdata[,1]/220<=0.005, 2], gdata[gdata[,1]/220<=0.005, 1]/220, 
		
		#filter out tumor (tdata), where gdata is not rare variant
	
		comVariants<-intersect(tdata[,2], gdata[gdata[,1]/N > p,2])
		
		filter.index<-c()
		for (j in 1:length(comVariants)){
			filter.index<-c(filter.index,which(as.vector(tdata[,2])==comVariants[j]))
			}

		plot(tdata[-filter.index,2], tdata[-filter.index,1]/N, type="l", col="red",main=paste("Chr", i, sep=""), 
xlab="SNV chromosome position", ylab="SNV occurence prob", ylim=c(0, upper))

		points(gdata[gdata[,1]/N<=p,2], gdata[gdata[,1]/N<=p,1]/N, type="p",col="lightblue")
		}	
	dev.off()

	}	


PlotSNV("./tumor_normal.pdf", N=200, p=0.005, upper=0.3)
