
##' @title Visulization of An Integrated Figure
##'        (STEP4 of The Add-Dom Model)
##' 
##' @description Visualizing additive and non-additive QTLs by an integrated plot.
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
##' @param Plot_model A character. The model to be mainly focused in the integrated plot. Please select from choose "AvsAD" or "NvsAD" or "NvsA" or "NvsD".
##' @param chrs_sum A numeric variable. The sum of all chromosomes.
##' @param num_nodes A numeric variable. The number of cores used parallelly.
##'
##' @return a folder named "3_Plot" with various plots.
##'
##' @examples 
##' ADDO_AddDom4_IntePlot(outdir=outdir, covariates_sum=2, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000)
##'
##' @author Leilei Cui and Bin Yang

ADDO_AddDom4_IntePlot <- function(outdir = outdir, 
                                  PheList_Choose = F, 
                                  PheList = PheList, 
                                  covariates_sum = covariates_sum, 
                                  RegionMan_chr_whole = F,
                                  RegionMan_chr_region = RegionMan_chr_region,
                                  Down_sampling = F,
                                  Down_sampling_logP = 1,
                                  Down_sampling_distance = 10,							 
                                  Plot_model = "NvsAD",
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

Phe_Nor = read.table(paste(rawdir,"/../4_PheNor/Phe_Nor.txt",sep=""),header=T)
if(PheList_Choose == T){ PheList = PheList[PheList %in% colnames(Phe_Nor)] }else{ PheList = colnames(Phe_Nor) }

if(Plot_model == "AvsAD"){ Plot_sort = c(1,2,3); Choose_Pvalue_col = 6; Title_used = "A vs A+D"}
if(Plot_model == "NvsAD"){ Plot_sort = c(2,1,3); Choose_Pvalue_col = 5; Title_used = "Null vs A+D"}
if(Plot_model == "NvsA"){ Plot_sort = c(3,1,2); Choose_Pvalue_col = 3; Title_used = "Null vs A"}
if(Plot_model == "NvsD"){ Plot_sort = c(4,2,3); Choose_Pvalue_col = 4; Title_used = "Null vs D"}

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
		TEMP_logP = subset(TEMP_logP_raw,TEMP_logP_raw[,"chr"]!=0 & TEMP_logP_raw[,"pos"]!=0)
		TEMP_logP = TEMP_logP[complete.cases(TEMP_logP),]
		TEMP_logP_all = TEMP_logP
		
		if(Down_sampling){
			TEMP_logP_low = subset(TEMP_logP, TEMP_logP[,5]<Down_sampling_logP)
			TEMP_logP_high = subset(TEMP_logP, TEMP_logP[,5]>Down_sampling_logP)
			
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
			TEMP_logP = TEMP_logP_sampling[order(TEMP_logP_sampling[,1],TEMP_logP_sampling[,2]),] # Revalue TEMP_logP with sampling poings #
		}

		setwd(paste(outdir,"/0_All",sep=""))
		png(paste("All_",phenotype_name,".png",sep=""),width=19500,height=17500,res=600,type="cairo")
		
		dist = 0.001
		m <- rbind(c(0+dist, 0.8-dist, 0.6+dist, 1-dist), c(0.01+dist, 0.45, 0.19+dist, 0.38), c(0.01+dist, 0.45, 0+dist, 0.19), c(0.81, 1-dist, 0.8005, 1-dist), c(0.81, 1-dist, 0.6+dist, 0.7995), c(0.01+dist, 0.45, 0.38+dist, 0.6), c(0.45+dist, 0.77, 0+dist, 0.6), c(0.77+dist, 1-dist, 0.4+dist, 0.6), c(0.77+dist, 1-dist, 0.2+dist, 0.4), c(0.77+dist, 1-dist, 0+dist, 0.2))
		split.screen(m)
		
		try({
			
			print(date()); Region_TEMP_logP = cbind(paste("SNP",c(1:dim(TEMP_logP)[1]),sep="_"),TEMP_logP); colnames(Region_TEMP_logP)[1] = "SNP"
			Region_TEMP_logP = subset(Region_TEMP_logP,Region_TEMP_logP[,2]!=0 | Region_TEMP_logP[,3]!=0) # Remove SNPs with missing chr or pos #
			Region_TEMP_logP = Region_TEMP_logP[complete.cases(Region_TEMP_logP),]
			thrgenome = -log10(0.05/nrow(TEMP_logP_all)); thrsuggest = -log10(1/nrow(TEMP_logP_all))

			logP_AddDomCode_Add = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,7])	
			logP_AddDomCode_Null = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,6])
			logP_AddCode = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,4])
			logP_DomCode = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,5])
			logP_All = c(); logP_All[[1]] = logP_AddDomCode_Add; logP_All[[2]] = logP_AddDomCode_Null; logP_All[[3]] = logP_AddCode; logP_All[[4]] = logP_DomCode
			Title_All = c("A vs A+D", "Null vs A+D", "Null vs A", "Null vs D")
					
			y_limit_tmp = c(logP_AddCode[,4],logP_AddDomCode_Null[,4],logP_AddDomCode_Add[,4])
			y_limit = max(y_limit_tmp[!is.na(y_limit_tmp)&!is.infinite(y_limit_tmp)])

			screen(1); par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=3.4)
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(logP_All[[Plot_sort[1]]]),y_limit=y_limit,cex_points = 1.2,cex_lab = 3.1,cex_axis = 3)
			title(paste(phenotype_name,Title_All[Plot_sort[1]],sep=" : ")); abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			
			screen(2); par(mai=c(1.1,1,0.5,0.5),cex.main=2.6)
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(logP_All[[Plot_sort[2]]]),y_limit=y_limit,cex_points = 0.8,cex_lab = 2.1,cex_axis = 2)
			title(paste(phenotype_name,Title_All[Plot_sort[2]],sep=" : ")); abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			
			screen(3); par(mai=c(1.1,1,0.5,0.5),cex.main=2.6)
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(logP_All[[Plot_sort[3]]]),y_limit=y_limit,cex_points = 0.8,cex_lab = 2.1,cex_axis = 2)
			title(paste(phenotype_name,Title_All[Plot_sort[3]],sep=" : ")); abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			
			print(paste(phenotype_name,"Manhattan Plot Done","^_^",sep=" "))
			
		})
		
		try({
		
			print(date()); screen(4); par(mai=c(0.95,0.95,0.8,0.5),cex.main=1)
			data_qq4 = TEMP_logP_all[,Choose_Pvalue_col]; data_qq4 <- as.numeric(data_qq4); data_qq4 <- data_qq4[!is.na(data_qq4)]
			expPvals4 = sort(-log10((1:length(data_qq4))/length(data_qq4))); obsPvals4 = sort(data_qq4)
			plot(expPvals4,obsPvals4,main=paste0("All QQ Plot (", Title_used,")"),xlab="Expected -logP",ylab="Observed logP",cex.lab=2.1,cex.axis=2,cex.main=1.9,cex=1.5)
			lines(expPvals4,expPvals4,type="l",col="red")
			
			screen(5); par(mai=c(0.95,0.95,0.8,0.5),cex.main=1)
			X_allhap_clean_tmp = TEMP_logP_all[,c(1,2,Choose_Pvalue_col)]; rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
			Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
			X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
			#X_allhap_rm = subset(X_allhap_rm, X_allhap_rm[,2] <= (Peak_tmp[1,"pos"]+RegionMan_chr_region) & X_allhap_rm[,2] >= (Peak_tmp[1,"pos"]-RegionMan_chr_region))
			X_allhap_clean4 = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
			
			data_qq4 = X_allhap_clean4[,3]; data_qq4 <- as.numeric(data_qq4); data_qq4 <- data_qq4[!is.na(data_qq4)]
			expPvals4 = sort(-log10((1:length(data_qq4))/length(data_qq4))); obsPvals4 = sort(data_qq4)
			plot(expPvals4,obsPvals4,main=paste0("Filtered QQ Plot (", Title_used,")"),xlab="Expected -logP",ylab="Observed logP",cex.lab=2.1,cex.axis=2,cex.main=1.9,cex=1.5)
			lines(expPvals4,expPvals4,type="l",col="red")
			
			print(paste(phenotype_name,"QQ Plot Done","^_^",sep=" "))
						
		})
		
		try({

			print(date()); tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,Choose_Pvalue_col] == max(TEMP_logP[,Choose_Pvalue_col]))[1,1:2]; ch <- tmp_chr_pos[1,1]
			snp_name = as.character(subset(tmp_map,tmp_map[,1]==tmp_chr_pos[1,1] & tmp_map[,4]==tmp_chr_pos[1,2])[1,2])
			
			regiondat <- subset(TEMP_logP,TEMP_logP[,1]==ch)[,c(1,2,Choose_Pvalue_col)]
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
			system("rm RegionMan_tmp.*"); ldinfo = TEMP.ld[regiondat_used[,1],c("SNP_B","R2")]
			
			setwd(paste(outdir,"/0_All",sep=""))
			plotdat = data.frame(SNP=regiondat_used$SNP,chr=regiondat_used$chr,pos=regiondat_used$pos,logp=regiondat_used$logP)
			thrgenome = -log10(0.05/dim(TEMP_logP_all)[1]); thrsuggest = -log10(1/dim(TEMP_logP_all)[1])

			screen(6); par(mai=c(1.1,1,0.5,0.5))
			ylim_tmp = subset(regiondat_used,regiondat_used[,"chr"]==tmp_chr_pos[1,"chr"] & regiondat_used[,"pos"]==tmp_chr_pos[1,"pos"])[1,"logP"]*1.15
			plotRegion(chrs=ch,alldat=plotdat,traitIdx=1,from=SNPh,to=SNPt,ldinfo=ldinfo,ylim=c(0,ylim_tmp),cex_points = 2.7, cex_points_peak = 4, cex_lab = 2.1 ,cex_axis = 2)
			title(paste0("Peak SNP Region (", Title_used,")"),cex.main=2.6)
			abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)

			print(paste(phenotype_name,"Region Manhattan Plot Done","^_^",sep=" "))
		})
			
		try({
		
			print(date()); tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,Choose_Pvalue_col] == max(TEMP_logP[,Choose_Pvalue_col]))[1,1:2]
			snp_name = as.character(subset(tmp_map,tmp_map[,1]==tmp_chr_pos[1,1] & tmp_map[,4]==tmp_chr_pos[1,2])[1,2])

			setwd(rawdir); system(paste("plink --noweb --silent --bfile Data_clean --snp",snp_name,"--recode --out Boxplot_tmp"))
			TEMP.geno <- read.ped("Boxplot_tmp.ped"); system("rm Boxplot_tmp.*")
			TEMP_geno =  as.data.frame(paste(TEMP.geno[,7],TEMP.geno[,8],sep="/"))
			rownames(TEMP_geno) = TEMP.geno[,2]; colnames(TEMP_geno) = snp_name 			
			TEMP_phe = read.table(paste(rawdir,"/../2_PheClean/Phe_Clean.txt",sep=""),header=T)[,phenotype_name]
			TEMP_gp = data.frame(TEMP_geno,TEMP_phe); TEMP_gp = na.omit(TEMP_gp); TEMP_gp = subset(TEMP_gp, TEMP_gp[,1]!="0/0")
			phe = TEMP_gp[,2]; geno = TEMP_gp[,1]; unigeno = sort(unique(geno))
			PheGeno_clean = as.data.frame(cbind(phe=phe,geno=geno)); PheGeno_clean[,"phe"]=as.numeric(PheGeno_clean[,"phe"])
			
			screen(7); par(mai=c(0.9,0.8,0.9,0.3),mgp=c(4.1,1.5,0))
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
			title(paste0("Peak SNP (", Title_used,")"),cex.main=3.3)
			
			screen(8); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[1])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[1],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[1]," Individuals",sep=""))
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(9); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[2])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[2],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[2]," Individuals",sep=""))
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(10); par(mai=c(1,1,0.9,0.5))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[3])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=paste(unigeno[3],sep=""),col=c("red"),breaks=150,cex.lab=2.3,cex.axis=2,cex.main=3,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=paste("Phenotype of ",unigeno[3]," Individuals",sep=""))
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
		print(paste(phenotype_name," Already been done! ^_^",sep=""))
	}
	
})
}

}

