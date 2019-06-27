
##' @title Visulization of Various Different Figures
##'        (STEP3 of The Heterotic Model)
##' 
##' @description Visualizing overdominance QTLs by various plots.
##'
##' @param outdir A character. The output directory where generates the folder: "3_Plot".
##' @param PheList_Choose A logic variable. T: Just investigate specified phenotypes; F: Investigate all phenofiles.
##' @param PheList A vector of character. When choose "PheList_Choose=F", the specified phenotype list must be specified.
##' @param covariates_sum A numeric variable. The sum of all covariates.
##' @param RegionMan_chr_whole A logic variable. T: Draw a whole chromosome; F: Draw the specified region.
##' @param RegionMan_chr_region A numeric variable. The length of specified region around the Peak SNP.
##' @param Down_sampling A logic variable. T: Down-sampling points with low logP to speed up the plotting progress; F: Darwing with all loci.
##' @param Down_sampling_logP A numeric variable. The threshold for down-sampling points for rapid rendering of the Manhattan Plots, when Down_sampling is true.
##' @param Down_sampling_distance A numeric variable. The distance of points for equidistant sampling, when Down_sampling is true.
##' @param chrs_sum A numeric variable. The sum of all chromosomes.
##' @param num_nodes A numeric variable. The number of cores used parallelly.
##'
##' @return a folder named "3_Plot" with various plots.
##'
##' @examples 
##' ADDO_Heterotic3_Plot(outdir=outdir, covariates_sum=2, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000)
##'
##' @author Leilei Cui and Bin Yang

ADDO_Heterotic3_Plot <- function(outdir = outdir,
                                 PheList_Choose = F,
                                 PheList = PheList,
                                 covariates_sum = covariates_sum,
                                 RegionMan_chr_whole = F,
                                 RegionMan_chr_region = RegionMan_chr_region,
                                 Down_sampling = F,
                                 Down_sampling_logP = 1,
                                 Down_sampling_distance = 10,
                                 chrs_sum = chrs_sum,
                                 num_nodes = 10){

rawdir = paste(outdir,"/1_PheGen/0_Data",sep="")
Pdir = paste(outdir,"/2_Pvalue",sep="")
setwd(paste(Pdir,"/..",sep=""))
system("mkdir 3_Plot 3_Plot/1_Manhattan 3_Plot/2_QQ_Plot 3_Plot/3_Contour_Plot 3_Plot/4_PeakSNP_RegionMan 3_Plot/5_PeakSNP_BoxPlot")
outdir = paste(outdir,"/3_Plot",sep="")

read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))
read.ped <- function(...) as.data.frame(fread(header=F,colClasses="character",...))

Phe_Clean = read.header(paste(rawdir,"/../2_PheClean/Phe_Clean.txt",sep="")); rownames(Phe_Clean)=Phe_Clean[,1]
Phe_Nor = read.header(paste(rawdir,"/../4_PheNor/Phe_Nor.txt",sep="")); rownames(Phe_Nor)=Phe_Nor[,1]

if(PheList_Choose == T){ PheList = PheList[PheList %in% colnames(Phe_Nor)] }else{ PheList = colnames(Phe_Nor) }

PeakSNPs_all = data.frame()

for(i in (covariates_sum+2):length(PheList)){
try({

	phenotype_name=PheList[i]
	
	tmp_map = read(paste(rawdir,"/Data_clean.bim",sep="")); rownames(tmp_map) = tmp_map[,2]	
	setwd(Pdir); system(paste('tar zxvf Pvalue_',phenotype_name,'.txt.tar.gz',sep=''))
	TEMP_logP_raw = read.header(paste('Pvalue_',phenotype_name,'.txt',sep=''))
	rownames(TEMP_logP_raw) = tmp_map[,2] # Note: Watch out the rownames of TEMP_logP_raw #
	system(paste('rm Pvalue_',phenotype_name,'.txt',sep=''))
	TEMP_logP_raw = subset(TEMP_logP_raw,TEMP_logP_raw[,"chr"]!=0 & TEMP_logP_raw[,"pos"]!=0)
	TEMP_logP_raw = TEMP_logP_raw[complete.cases(TEMP_logP_raw),]
	TEMP_logP = subset(TEMP_logP_raw, TEMP_logP_raw[,"tAB_AA"]*TEMP_logP_raw[,"tAB_BB"]>0) 
	TEMP_logP_all = TEMP_logP_raw 
	
	if(Down_sampling){
		TEMP_logP_low = subset(TEMP_logP_raw, TEMP_logP_raw[,4]<Down_sampling_logP)
		TEMP_logP_high = subset(TEMP_logP_raw, TEMP_logP_raw[,4]>Down_sampling_logP)
		
		uniq_chr_all = unique(TEMP_logP_low[,1]); TEMP_logP_keep = c()
		for(uniq_chr in 1:length(uniq_chr_all)){
			TEMP_logP_tmp = subset(TEMP_logP_low, TEMP_logP_low[,1]==uniq_chr_all[uniq_chr])			
			rownames(TEMP_logP_tmp) = 1:(dim(TEMP_logP_tmp)[1])
			TEMP_logP_tmp_keep = TEMP_logP_tmp[seq(1,dim(TEMP_logP_tmp)[1],Down_sampling_distance),]
			rownames(TEMP_logP_tmp_keep) = paste(TEMP_logP_tmp_keep[,1],TEMP_logP_tmp_keep[,2],sep="_")
			TEMP_logP_keep = rbind(TEMP_logP_keep, TEMP_logP_tmp_keep)
			print(uniq_chr)
		}
		
		TEMP_logP_sampling = rbind(TEMP_logP_keep, TEMP_logP_high)
		TEMP_logP_raw = TEMP_logP_sampling[order(TEMP_logP_sampling[,1],TEMP_logP_sampling[,2]),]
		TEMP_logP = subset(TEMP_logP_raw, TEMP_logP_raw[,"tAB_AA"]*TEMP_logP_raw[,"tAB_BB"]>0)
	}

	NORM_logP_AB_AA = TEMP_logP_raw[,"NORM_logP_AB_AA"]; names(NORM_logP_AB_AA) = rownames(TEMP_logP_raw)	
	NORM_logP_AB_BB = TEMP_logP_raw[,"NORM_logP_AB_BB"]; names(NORM_logP_AB_BB) = rownames(TEMP_logP_raw)
	MVN_logP_Minor = TEMP_logP_raw[,"MVN_logP_Minor"]; names(MVN_logP_Minor) = rownames(TEMP_logP_raw)

	
	if(!file.exists(paste(outdir,'/1_Manhattan/Manhattan_',phenotype_name,'.png',sep=''))){
		try({
			setwd(paste(outdir,"/1_Manhattan",sep=""))
			
			Region_TEMP_logP = cbind(paste("SNP",c(1:dim(TEMP_logP_raw)[1]),sep="_"),TEMP_logP_raw); colnames(Region_TEMP_logP)[1] = "SNP"
			thrgenome = -log10(0.05/nrow(Region_TEMP_logP)); thrsuggest = -log10(1/nrow(Region_TEMP_logP))
			
			plotdat_AB_AA = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_AA)
			plotdat_AB_BB = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_BB) 
			plotdat_AB_AA_BB = cbind(plotdat_AB_AA[,1:3],Pc1df1=MVN_logP_Minor)
			
			y_limit_tmp = c(plotdat_AB_AA[,4],plotdat_AB_BB[,4],plotdat_AB_AA_BB[,4])
			y_limit = max(y_limit_tmp[!is.na(y_limit_tmp)&!is.infinite(y_limit_tmp)])
			
			png(paste("Manhattan_",phenotype_name,".png",sep=""),width=15000,height=10000,res=600,type="cairo")
			layout(matrix(c(1,3,2,3),ncol=2),widths=rep(1,2),heights=c(0.7,1))
			par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=3.5)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_AA),y_limit=y_limit)
			title(paste(phenotype_name,"t(AB-AA)",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_BB),y_limit=y_limit)
			title(paste(phenotype_name,"t(AB-BB)",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_AA_BB),y_limit=y_limit)
			title(paste(phenotype_name,"Minor(|t(AB-AA)|, |t(AB-BB)|)",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			dev.off()
			

			Region_TEMP_logP = cbind(paste("SNP",c(1:dim(TEMP_logP)[1]),sep="_"),TEMP_logP); colnames(Region_TEMP_logP)[1] = "SNP"
			thrgenome = -log10(0.05/nrow(Region_TEMP_logP)); thrsuggest = -log10(1/nrow(Region_TEMP_logP))
			
			plotdat_AB_AA = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_AA[rownames(TEMP_logP)])
			plotdat_AB_BB = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_BB[rownames(TEMP_logP)]) 
			plotdat_AB_AA_BB = cbind(plotdat_AB_AA[,1:3],Pc1df1=MVN_logP_Minor[rownames(TEMP_logP)])
			
			Tvalue_thrgenome = TEMP_logP_raw[rownames(subset(plotdat_AB_AA_BB,plotdat_AB_AA_BB[,"Pc1df1"]>thrgenome)),c("tAB_AA","tAB_BB")]
			Tvalue_thrsuggest = TEMP_logP_raw[rownames(subset(plotdat_AB_AA_BB,plotdat_AB_AA_BB[,"Pc1df1"]>thrsuggest)),c("tAB_AA","tAB_BB")]
			
			y_limit_tmp = c(plotdat_AB_AA[,4],plotdat_AB_BB[,4],plotdat_AB_AA_BB[,4])
			y_limit = max(y_limit_tmp[!is.na(y_limit_tmp)&!is.infinite(y_limit_tmp)])
			
			png(paste("ManhattanFilter_",phenotype_name,".png",sep=""),width=15000,height=10000,res=600,type="cairo")
			layout(matrix(c(1,3,2,3),ncol=2),widths=rep(1,2),heights=c(0.7,1))
			par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=2.8)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_AA),y_limit=y_limit)
			title(paste(phenotype_name,"t(AB-AA) using slected points",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_BB),y_limit=y_limit)
			title(paste(phenotype_name,"t(AB-BB) using slected points",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_AA_BB),y_limit=y_limit)
			title(paste(phenotype_name,"Minor(|t(AB-AA)|, |t(AB-BB)|) using slected points",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			dev.off()
			
			print(paste(phenotype_name,"Manhattan Plot Done","^_^",sep=" "))
		})


		try({
			setwd(paste(outdir,"/2_QQ_Plot",sep=""))
				
			png(paste("QQplot_",phenotype_name,".png",sep=""),width=4600,height=4500,res=200,type="cairo")
			layout(matrix(c(1,4,7,2,5,8,3,6,9),nrow=3),widths=rep(1,3),heights=rep(1,3))
			par(mai=c(0.95,1.1,1,0.4),mgp=c(4.4,1.5,0),cex.main=1)

			Region_TEMP_logP = TEMP_logP_raw

			data_histogram1 = Region_TEMP_logP[,"tAB_AA"]; object_statistic1 = "t(AB_AA)"
			data_histogram1 <- as.numeric(data_histogram1); data_histogram1 <- data_histogram1[!is.na(data_histogram1)]
			h <- hist(data_histogram1,main=paste("Histogram Plot : ",phenotype_name,sep=""),xlab=object_statistic1,ylab="Number",col=c("red"),breaks=170,cex.lab=3.2,cex.axis=3.2,cex.main=3.7)
			xfit = seq(min(data_histogram1),max(data_histogram1),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram1),sd=sd(data_histogram1))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram1)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()

			data_histogram2 = Region_TEMP_logP[,"tAB_BB"]; object_statistic2 = "t(AB_BB)"
			data_histogram2 <- as.numeric(data_histogram2); data_histogram2 <- data_histogram2[!is.na(data_histogram2)]
			h <- hist(data_histogram2,main=paste("Histogram Plot : ",phenotype_name,sep=""),xlab=object_statistic2,ylab="Number",col=c("red"),breaks=170,cex.lab=3.2,cex.axis=3.2,cex.main=3.7)
			xfit = seq(min(data_histogram2),max(data_histogram2),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram2),sd=sd(data_histogram2))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram2)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()

				X_allhap_filter = Region_TEMP_logP
				X_allhap_rmna = X_allhap_filter[complete.cases(X_allhap_filter),]
				data_histogram = apply(X_allhap_rmna[,c("tAB_AA","tAB_BB")],1,function(x){x_abs = abs(x); return(x[x_abs==min(x_abs)])})	
				data_histogram <- as.numeric(data_histogram); data_histogram <- data_histogram[!is.na(data_histogram)]
				
				h <- hist(data_histogram,main=paste("Histogram Plot : ",phenotype_name,sep=""),xlab="Minor_Abs_t",ylab="Number",col=c("red"),breaks=170,cex.lab=3.2,cex.axis=3.2,cex.main=3.7)
				xfit = seq(min(data_histogram),max(data_histogram),length=1000)
				yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
				yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
				lines(xfit,yfit,col="blue",lwd=3)		
				box()

			data_qq1 = Region_TEMP_logP[,"NORM_logP_AB_AA"]; object_statistic1 = "logP(AB_AA)"
			data_qq1 <- as.numeric(data_qq1); data_qq1 <- data_qq1[!is.na(data_qq1)]
			expPvals1 = sort(-log10((1:length(data_qq1))/length(data_qq1))) 
			obsPvals1 = sort(data_qq1)
			plot(expPvals1,obsPvals1,main=paste("QQ Plot : ",phenotype_name,sep=""),xlab=paste("Expected ",object_statistic1,sep=""),ylab=paste("Observed ",object_statistic1,sep=""),cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
			lines(expPvals1,expPvals1,type="l",col="red")

			data_qq2 = Region_TEMP_logP[,"NORM_logP_AB_BB"]; object_statistic2 = "logP(AB_BB)"
			data_qq2 <- as.numeric(data_qq2); data_qq2 <- data_qq2[!is.na(data_qq2)]
			expPvals2 = sort(-log10((1:length(data_qq2))/length(data_qq2))) 
			obsPvals2 = sort(data_qq2)
			plot(expPvals2,obsPvals2,main=paste("QQ Plot : ",phenotype_name,sep=""),xlab=paste("Expected ",object_statistic2,sep=""),ylab=paste("Observed ",object_statistic2,sep=""),cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
			lines(expPvals2,expPvals2,type="l",col="red")

				data_qq = MVN_logP_Minor
				data_qq <- as.numeric(data_qq); data_qq <- data_qq[!is.na(data_qq)]
				expPvals = sort(-log10((1:length(data_qq))/length(data_qq)))
				obsPvals = sort(data_qq)
				plot(expPvals,obsPvals,main=paste("QQ Plot : ",phenotype_name,sep=""),xlab="Expected Minor_logP",ylab="Observed Minor_logP",cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
				lines(expPvals,expPvals,type="l",col="red")	
				
			X_allhap_clean_tmp = Region_TEMP_logP[,c("chr","pos","NORM_logP_AB_AA")]; rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
			Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
			X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
			# X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] ==  Peak_tmp[1,"chr"] & X_allhap_clean_tmp[,2] <= (Peak_tmp[1,"pos"]+10000000) & X_allhap_clean_tmp[,2] >= (Peak_tmp[1,"pos"]-10000000))
			X_allhap_clean = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
			
			data_qq1 = X_allhap_clean[,3]; object_statistic1 = "logP(AB_AA)"
			data_qq1 <- as.numeric(data_qq1); data_qq1 <- data_qq1[!is.na(data_qq1)]
			expPvals1 = sort(-log10((1:length(data_qq1))/length(data_qq1))) 
			obsPvals1 = sort(data_qq1)
			plot(expPvals1,obsPvals1,main=paste("QQ Plot rmchr : ",phenotype_name, sep=""),xlab=paste("Expected ",object_statistic1,sep=""),ylab=paste("Observed ",object_statistic1,sep=""),cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
			lines(expPvals1,expPvals1,type="l",col="red")

			X_allhap_clean_tmp = Region_TEMP_logP[,c("chr","pos","NORM_logP_AB_BB")]; rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
			Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
			X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
			# X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] ==  Peak_tmp[1,"chr"] & X_allhap_clean_tmp[,2] <= (Peak_tmp[1,"pos"]+10000000) & X_allhap_clean_tmp[,2] >= (Peak_tmp[1,"pos"]-10000000))
			X_allhap_clean = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
			
			data_qq2 = X_allhap_clean[,3]; object_statistic2 = "logP(AB_BB)"
			data_qq2 <- as.numeric(data_qq2); data_qq2 <- data_qq2[!is.na(data_qq2)]
			expPvals2 = sort(-log10((1:length(data_qq2))/length(data_qq2))) 
			obsPvals2 = sort(data_qq2)
			plot(expPvals2,obsPvals2,main=paste("QQ Plot rmchr : ",phenotype_name, sep=""),xlab=paste("Expected ",object_statistic2,sep=""),ylab=paste("Observed ",object_statistic2,sep=""),cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
			lines(expPvals2,expPvals2,type="l",col="red")

				X_allhap_filter = cbind(Region_TEMP_logP,MVN_logP_Minor=MVN_logP_Minor)
				X_allhap_clean_tmp = X_allhap_filter[,c("chr","pos","MVN_logP_Minor")]
				rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
				Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
				X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
				# X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] ==  Peak_tmp[1,"chr"] & X_allhap_clean_tmp[,2] <= (Peak_tmp[1,"pos"]+10000000) & X_allhap_clean_tmp[,2] >= (Peak_tmp[1,"pos"]-10000000))
				X_allhap_clean = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
				
				data_qq = X_allhap_clean[,3]; data_qq <- data_qq[!is.na(data_qq)]
				expPvals = sort(-log10((1:length(data_qq))/length(data_qq)))
				obsPvals = sort(data_qq)
				plot(expPvals,obsPvals,main=paste("QQ Plot rmchr: ",phenotype_name, sep=""),xlab="Expected Minor_logP",ylab="Observed Minor_logP",cex.lab=3.2,cex.axis=3.2,cex.main=3.7,cex=1.5)
				lines(expPvals,expPvals,type="l",col="red")

			dev.off()
			
			print(paste(phenotype_name,"QQ Plot Done","^_^",sep=" "))
		})

		
		try({
			setwd(paste(outdir,"/3_Contour_Plot",sep=""))	
			
			data_tvalue_clean = TEMP_logP_raw[,c("tAB_AA","tAB_BB")]
			cor = cor(data_tvalue_clean[,1],data_tvalue_clean[,2]); Sigma = matrix(c(1,cor,cor,1),nrow=2)
			data_tvalue_clean_sim = rmvnorm(nrow(data_tvalue_clean),mean=c(0,0),sigma=Sigma)

			png(paste("ContourPlot_",phenotype_name,".png",sep=""),width=9000,height=9000,res=600,type="cairo")
			layout(matrix(c(1,3,2,4),nrow=2),widths=rep(1,2),heights=rep(1,2))
			par(mai=c(0.95,1,1,0.5),mgp=c(3.6,1.2,0))

			data_histogram1 = data_tvalue_clean[,1]; object_statistic1 = "t(AB-AA)"
			h <- hist(data_histogram1,main=paste("t(AB-AA): ",phenotype_name,sep=""),xlab=object_statistic1,ylab="Number",col=c("red"),breaks=170,cex.lab=2.5,cex.axis=2.5,cex.main=2.8) 
			xfit = seq(min(data_histogram1),max(data_histogram1),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram1),sd=sd(data_histogram1))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram1)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()

			data_histogram2 = data_tvalue_clean[,2]; object_statistic2 = "t(AB-BB)"
			h <- hist(data_histogram2,main=paste("t(AB-BB): ",phenotype_name,sep=""),xlab=object_statistic2,ylab="Number",col=c("red"),breaks=170,cex.lab=2.5,cex.axis=2.5,cex.main=2.8) 
			xfit = seq(min(data_histogram2),max(data_histogram2),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram2),sd=sd(data_histogram2))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram2)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()

			colors <- densCols(data_tvalue_clean_sim)
			plot(data_tvalue_clean_sim,col=colors,pch=20,main="Scatter & Contour Plot (Simulated)",xlab="Simulated: t(AB-AA)",ylab="Simulated: t(AB-BB)",xlim=c(-8,8),ylim=c(-8,8),cex=1,cex.axis=2.5,cex.lab=2,cex.main=2.3)
			legend("topleft",legend=c(paste("cor(X,Y)=",round(cor,3),sep="")),bty="n",cex=2)
			bivn.kde.sim = kde2d(data_tvalue_clean_sim[,1], data_tvalue_clean_sim[,2], n = 100)
			contour(bivn.kde.sim,col="red",lwd=1, labcex = 0.5, nlevels=10, add=TRUE)

			colors <- densCols(data_tvalue_clean)
			plot(data_tvalue_clean,col=colors,pch=20,main="Scatter & Contour Plot (Observed)",xlab="Observed: t(AB-AA)",ylab="Observed: t(AB-BB)",xlim=c(-8,8),ylim=c(-8,8),cex=1,cex.axis=2.5,cex.lab=2,cex.main=2.3)
			points(Tvalue_thrsuggest,col="red",pch=17,cex=1.3)
			points(Tvalue_thrgenome,col="red4",pch=17,cex=1.3)
			legend("topleft",legend=c(paste("cor(t(AB-AA),t(AB-BB))=",round(cor,3),sep="")),bty="n",cex=2)
			bivn.kde = kde2d(data_tvalue_clean[,1], data_tvalue_clean[,2], n = 100)
			contour(bivn.kde,col="red",lwd=1, labcex = 0.5, nlevels=10, add=TRUE)
			dev.off()
		
			print(paste(phenotype_name,"Contour Plot Done","^_^",sep=" "))
		})
		

		try({

			setwd(paste(outdir,"/4_PeakSNP_RegionMan",sep=""))		
			tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,"MVN_logP_Minor"] == max(TEMP_logP[,"MVN_logP_Minor"]))[1,1:2]; ch <- tmp_chr_pos[1,1]
			snp_name = as.character(subset(tmp_map,tmp_map[,1]==tmp_chr_pos[1,1] & tmp_map[,4]==tmp_chr_pos[1,2])[1,2])
			
			regiondat <- subset(TEMP_logP,TEMP_logP[,1]==ch)[,c(1,2,3)]
			regiondat = cbind(snpname=rownames(regiondat),regiondat)
			colnames(regiondat) = c("SNP","chr","pos","logP")
			regionlogp <- regiondat[,4]; maxlogp <- max(regionlogp)[1]; peakSNP <- snp_name
			if(RegionMan_chr_whole==T){
				posh <- regiondat[1,"pos"]; post <- regiondat[nrow(regiondat),"pos"]
				SNPh <- regiondat[1,"SNP"]; SNPt <- regiondat[nrow(regiondat),"SNP"]
			}else{
				peakpos <- regiondat[which(regiondat$logP==maxlogp)[1],"pos"]
				posh <- regiondat[which(regiondat$pos<(peakpos-RegionMan_chr_region)),"pos"]
				post <- regiondat[which(regiondat$pos>(peakpos+RegionMan_chr_region)),"pos"]
				SNPh <- regiondat[which(regiondat$pos==posh[length(posh)]),"SNP"]
				SNPt <- regiondat[which(regiondat$pos==post[1]),"SNP"]
				if(length(posh)==0){posh <- regiondat[1,"pos"]; SNPh <- regiondat[1,"SNP"]}
				if(length(post)==0){post <- regiondat[nrow(regiondat),"pos"]; SNPt <- regiondat[nrow(regiondat),"SNP"] }
			}
			regiondat_used = subset(regiondat, regiondat[,"pos"]>=posh[length(posh)] & regiondat[,"pos"]<=post[1])

			setwd(rawdir); options("scipen"=100, "digits"=4) # NOTE: Avoid the RegionMan_chr_region become into scientific notation for PLINK #
			system(paste("plink --noweb --silent --bfile Data_clean --r2 --ld-snp ",peakSNP,"--ld-window-kb ",RegionMan_chr_region," --ld-window 99999 --ld-window-r2 0","--out RegionMan_tmp"))
			TEMP.ld <- read.header("RegionMan_tmp.ld"); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]
			system("rm RegionMan_tmp.*"); ldinfo = TEMP.ld[regiondat[,1],c("SNP_B","R2")]
			
			setwd(paste(outdir,"/4_PeakSNP_RegionMan",sep=""))		
			plotdat = data.frame(SNP=regiondat$SNP,chr=regiondat$chr,pos=regiondat$pos,logp=regiondat$logP)
			thrgenome = -log10(0.05/dim(TEMP_logP)[1]); thrsuggest = -log10(1/dim(TEMP_logP)[1])
							 
			png(paste("ManhattanRegion_",phenotype_name,".png",sep=""),width=7000,height=3500,res=600,type="cairo")
			par(mai=c(1,1,0.5,0.5))
			ylim_tmp = subset(regiondat,regiondat[,"chr"]==tmp_chr_pos[1,"chr"] & regiondat[,"pos"]==tmp_chr_pos[1,"pos"])[1,"logP"]*1.15
			plotRegion(chrs=ch,alldat=plotdat,traitIdx=1,from=SNPh,to=SNPt,ldinfo=ldinfo,ylim=c(0,ylim_tmp))
			title(paste0(phenotype_name, " (",tmp_chr_pos[1,1],", ",tmp_chr_pos[1,2],")"), cex.main=1.6)
			abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			dev.off()

			print(paste(phenotype_name,"Region Manhattan Plot Done","^_^",sep=" "))
		})
		
		try({
		
			setwd(paste(outdir,"/5_PeakSNP_BoxPlot",sep=""))
			tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,"MVN_logP_Minor"] == max(TEMP_logP[,"MVN_logP_Minor"]))[1,1:2]
			snp_name = as.character(subset(tmp_map,tmp_map[,1]==tmp_chr_pos[1,1] & tmp_map[,4]==tmp_chr_pos[1,2])[1,2])
			
			setwd(rawdir); system(paste("plink --noweb --silent --bfile Data_clean --snp",snp_name,"--recode --out Boxplot_tmp"))
			TEMP.geno <- read.ped("Boxplot_tmp.ped"); system("rm Boxplot_tmp.*")
			TEMP_geno =  as.data.frame(paste(TEMP.geno[,7],TEMP.geno[,8],sep="/")) # PED Format: FamilyID/IndID/FID/MID/Sex/Phe #
			rownames(TEMP_geno) = TEMP.geno[,2]; colnames(TEMP_geno) = snp_name 			
			TEMP_phe = Phe_Clean[,phenotype_name]
			TEMP_gp = data.frame(TEMP_geno,TEMP_phe); TEMP_gp = na.omit(TEMP_gp); TEMP_gp = subset(TEMP_gp, TEMP_gp[,1]!="0/0")
			phe = TEMP_gp[,2]; geno = TEMP_gp[,1]; unigeno = sort(unique(geno))
			PheGeno_clean = as.data.frame(cbind(phe=phe,geno=geno)); PheGeno_clean[,"phe"]=as.numeric(PheGeno_clean[,"phe"])
			
			setwd(paste(outdir,"/5_PeakSNP_BoxPlot",sep=""))
			png(paste("Boxplot_",phenotype_name,".png",sep=""),width=10000,height=9000,res=600,type="cairo")
			par(mai=c(0.9,0.8,0.9,0.5),mgp=c(4.1,1.5,0))
			layout(matrix(c(1,1,1,2,3,4),ncol=2),widths=c(1,0.7),heights=c(1,1,1))
			
			boxplot(PheGeno_clean[,1]~PheGeno_clean[,2],cex.axis=3.5,ylim=c(min(PheGeno_clean[,1],na.rm=T),1.2*max(PheGeno_clean[,1],na.rm=T)),col=rainbow(12))
			y = PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]]
			y = y[!is.na(y)]
			points(jitter(rep(1,length(which(PheGeno_clean[,2]==unigeno[1]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(2,length(which(PheGeno_clean[,2]==unigeno[2]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[2]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(3,length(which(PheGeno_clean[,2]==unigeno[3]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[3]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			text(1,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[1])),sep=""),cex=4)
			text(2,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[2])),sep=""),cex=4)
			text(3,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[3])),sep=""),cex=4)
			text(2,1.18*max(PheGeno_clean[,1],na.rm=T),paste0("(Chr:",tmp_chr_pos[1,1]," Pos:",tmp_chr_pos[1,2],")"),cex=3)
			title(phenotype_name,cex.main=3.3)
			
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[1])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[1],col=c("red"),breaks=150,cex.lab=3,cex.axis=2.5,cex.main=3.3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[2])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[2],col=c("red"),breaks=150,cex.lab=3,cex.axis=2.5,cex.main=3.3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[3])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[3],col=c("red"),breaks=150,cex.lab=3,cex.axis=2.5,cex.main=3.3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			dev.off()
			
			print(paste(phenotype_name,"Genotype Box Plot Done","^_^",sep=" "))
		})

		PeakSNPs = subset(TEMP_logP,TEMP_logP[,"MVN_logP_Minor"] == max(TEMP_logP[,"MVN_logP_Minor"]))[1,] # Choose the first SNP with the Peak P-value #
		PeakSNPs_output = cbind(Trait=phenotype_name,SNP=rownames(PeakSNPs),PeakSNPs)
		PeakSNPs_all = rbind(PeakSNPs_all,PeakSNPs_output)
	
	}else{
		print(paste(phenotype_name," : All the plots have been done! ^_^",sep=""))
	}

})
}

setwd(outdir); write.table(PeakSNPs_all,file="PeakSNPs_Summary.txt",row.names=T,col.names=T,quote=F)

}

