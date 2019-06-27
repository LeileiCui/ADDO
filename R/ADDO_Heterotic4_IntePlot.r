
##' @title Visulization of An Integrated Figure
##'        (STEP4 of The Heterotic Model)
##' 
##' @description Visualizing overdominance QTLs by an integrated plot.
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
##' ADDO_Heterotic4_IntePlot(outdir=outdir, covariates_sum=2, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000)
##'
##' @author Leilei Cui and Bin Yang

ADDO_Heterotic4_IntePlot <- function(outdir = outdir, 
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
system("mkdir 3_Plot 3_Plot/0_All")
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
	
	if(!file.exists(paste(outdir,'/0_All/All_',phenotype_name,'.png',sep=''))){

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
		

		setwd(paste(outdir,"/0_All",sep=""))
		png(paste("All_",phenotype_name,".png",sep=""),width=19500,height=17500,res=600,type="cairo")
		
		dist = 0.001
		m <- rbind(c(0+dist, 0.8-dist, 0.6+dist, 1-dist), c(0.06+dist, 0.4, 0+dist, 0.34),
			c(0.81, 1-dist, 0.8005, 1-dist), c(0.81, 1-dist, 0.6+dist, 0.7995),
			c(0.01+dist, 0.45, 0.34+dist, 0.6), 
			c(0.45+dist, 0.77, 0+dist, 0.6), 
			c(0.77+dist, 1-dist, 0.4+dist, 0.6), c(0.77+dist, 1-dist, 0.2+dist, 0.4), c(0.77+dist, 1-dist, 0+dist, 0.2))
		split.screen(m)
		
		
		try({			

			Region_TEMP_logP = cbind(paste("SNP",c(1:dim(TEMP_logP)[1]),sep="_"),TEMP_logP); colnames(Region_TEMP_logP)[1] = "SNP"
			thrgenome = -log10(0.05/nrow(Region_TEMP_logP)); thrsuggest = -log10(1/nrow(Region_TEMP_logP))
			
			plotdat_AB_AA = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_AA[rownames(TEMP_logP)])
			plotdat_AB_BB = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=NORM_logP_AB_BB[rownames(TEMP_logP)]) 
			plotdat_AB_AA_BB = cbind(plotdat_AB_AA[,1:3],Pc1df1=MVN_logP_Minor[rownames(TEMP_logP)])
			
			Tvalue_thrgenome = TEMP_logP_raw[rownames(subset(plotdat_AB_AA_BB,plotdat_AB_AA_BB[,"Pc1df1"]>thrgenome)),c("tAB_AA","tAB_BB")]
			Tvalue_thrsuggest = TEMP_logP_raw[rownames(subset(plotdat_AB_AA_BB,plotdat_AB_AA_BB[,"Pc1df1"]>thrsuggest)),c("tAB_AA","tAB_BB")]
			
			y_limit_tmp = c(plotdat_AB_AA_BB[,4])
			y_limit = max(y_limit_tmp[!is.na(y_limit_tmp)&!is.infinite(y_limit_tmp)])
					
			screen(1); par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=2.8)
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(plotdat_AB_AA_BB),y_limit=y_limit,cex_points = 1.2,cex_lab = 3.1,cex_axis = 3)
			title(paste(phenotype_name,"Minor(|t(AB-AA))|,|t(AB-BB)|) & t(AB-AA)*t(AB-BB)>0",sep=" : "))
			abline(h=thrgenome,lty=1)
			abline(h=thrsuggest,lty=2)
			
			print(paste(phenotype_name,"Manhattan Plot Done","^_^",sep=" "))
		})


		try({				

			Region_TEMP_logP = TEMP_logP_raw
			
			data_qq = MVN_logP_Minor
			data_qq <- as.numeric(data_qq); data_qq <- data_qq[!is.na(data_qq)]
			expPvals = sort(-log10((1:length(data_qq))/length(data_qq)))
			obsPvals = sort(data_qq)
			screen(3); par(mai=c(0.95,0.95,0.8,0.5),cex.main=1)
			plot(expPvals,obsPvals,main="All QQ Plot",xlab="Expected -logP",ylab="Observed logP",cex.lab=2.1,cex.axis=2,cex.main=1.9,cex=1.5)
			lines(expPvals,expPvals,type="l",col="red")	

			X_allhap_filter = cbind(Region_TEMP_logP,MVN_logP_Minor=MVN_logP_Minor)
			X_allhap_clean_tmp = X_allhap_filter[,c("chr","pos","MVN_logP_Minor")]
			rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
			Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
			X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
			#X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] ==  Peak_tmp[1,"chr"] & X_allhap_clean_tmp[,2] <= (Peak_tmp[1,"pos"]+10000000) & X_allhap_clean_tmp[,2] >= (Peak_tmp[1,"pos"]-10000000))
			X_allhap_clean = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
			
			data_qq = X_allhap_clean[,3]; data_qq <- data_qq[!is.na(data_qq)]
			expPvals = sort(-log10((1:length(data_qq))/length(data_qq)))
			obsPvals = sort(data_qq)
			screen(4); par(mai=c(0.95,0.95,0.8,0.5),cex.main=1)
			plot(expPvals,obsPvals,main="Filtered QQ Plot",xlab="Expected -logP",ylab="Observed logP",cex.lab=2.1,cex.axis=2,cex.main=1.9,cex=1.5)
			lines(expPvals,expPvals,type="l",col="red")
			
			print(paste(phenotype_name,"QQ Plot Done","^_^",sep=" "))
		})

		
		try({				
			data_tvalue_clean = TEMP_logP_raw[complete.cases(TEMP_logP_raw),][,c("tAB_AA","tAB_BB")]
			cor = cor(data_tvalue_clean[,1],data_tvalue_clean[,2]); Sigma = matrix(c(1,cor,cor,1),nrow=2)
			data_tvalue_clean_sim = rmvnorm(nrow(data_tvalue_clean),mean=c(0,0),sigma=Sigma)
			
			screen(2); par(mai=c(1.1,1,1,1))

			colors <- densCols(data_tvalue_clean)
			plot(data_tvalue_clean,col=colors,pch=20,main="Scatter & Contour Plot",xlab="t(AB-AA)",ylab="t(AB-BB)",xlim=c(-8,8),ylim=c(-8,8),cex=2.3,cex.axis=2.2,cex.lab=2.5,cex.main=2.6)
			points(Tvalue_thrsuggest,col="red",pch=17,cex=2.3)
			points(Tvalue_thrgenome,col="red4",pch=17,cex=2.3)
			legend("topleft",legend=c(paste("cor(t(AB-AA),t(AB-BB))=",round(cor,3),sep="")),bty="n",cex=2.7)
			bivn.kde.sim = kde2d(data_tvalue_clean_sim[,1], data_tvalue_clean_sim[,2], n = 100) 
			contour(bivn.kde.sim,col="red",lwd=2.2, labcex = 0.5, nlevels=10, add=TRUE)
		
			print(paste(phenotype_name,"Contour Plot Done","^_^",sep=" "))
		})
		

		try({

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
			
			plotdat = data.frame(SNP=regiondat$SNP,chr=regiondat$chr,pos=regiondat$pos,logp=regiondat$logP)
			thrgenome = -log10(0.05/dim(TEMP_logP)[1]); thrsuggest = -log10(1/dim(TEMP_logP)[1])
			
			screen(5); par(mai=c(1.1,1,0.5,0.5))
			ylim_tmp = subset(regiondat,regiondat[,"chr"]==tmp_chr_pos[1,"chr"] & regiondat[,"pos"]==tmp_chr_pos[1,"pos"])[1,"logP"]*1.15			
			plotRegion(chrs=ch,alldat=plotdat,traitIdx=1,from=SNPh,to=SNPt,ldinfo=ldinfo,ylim=c(0,ylim_tmp),cex_points = 2.7, cex_points_peak = 4, cex_lab = 2.1 ,cex_axis = 2)
			title(paste0(phenotype_name, " (",tmp_chr_pos[1,1],", ",tmp_chr_pos[1,2],")"), cex.main=2.6)
			abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)

			print(paste(phenotype_name,"Region Manhattan Plot Done","^_^",sep=" "))
		})
		
		
		try({
		
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
			
			screen(6); par(mai=c(0.9,0.8,0.9,0.3),mgp=c(4.1,1.5,0))
			boxplot(PheGeno_clean[,1]~PheGeno_clean[,2],cex.axis=2.6,ylim=c(min(PheGeno_clean[,1],na.rm=T),1.15*max(PheGeno_clean[,1],na.rm=T)),col=rainbow(12))
			y = PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]]
			y = y[!is.na(y)]
			points(jitter(rep(1,length(which(PheGeno_clean[,2]==unigeno[1]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(2,length(which(PheGeno_clean[,2]==unigeno[2]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[2]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(3,length(which(PheGeno_clean[,2]==unigeno[3]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[3]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			text(1,1.05*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[1])),sep=""),cex=2.5)
			text(2,1.05*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[2])),sep=""),cex=2.5)
			text(3,1.05*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[3])),sep=""),cex=2.5)
			text(2,1.12*max(PheGeno_clean[,1],na.rm=T),paste("(CHR: ",tmp_chr_pos[1,1],"  POS: ",tmp_chr_pos[1,2],"bp)",sep=""),cex=2.5)
			title("Peak SNP",cex.main=3.3)
			
			screen(7); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[1])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[1],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[1]," Inds",sep=""))
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(8); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[2])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[2],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[2]," Inds",sep=""))
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(9); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[3])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[3],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[3]," Inds",sep=""))
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			print(paste(phenotype_name,"Genotype Box Plot Done","^_^",sep=" "))
		})

		close.screen(all.screens = TRUE)		
		dev.off()	
		print(paste(phenotype_name,"All Plot Done","^_^",sep=" "))
	
	}else{
		print(paste(phenotype_name," : All the plots have been done! ^_^",sep=""))
	}

})
}

}

