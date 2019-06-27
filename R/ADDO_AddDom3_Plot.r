
##' @title Visulization of Various 4in1 Figures
##'        (STEP3 of The Add-Dom Model)
##' 
##' @description Visualizing additive and non-additive QTLs detected by four different models from ADDO_AddDom2_Pvalue.
##'              (1) 4in1 Manhattan Plot;
##'              (2) 4in1 QQ Plot using all loci with or without those loci located in the chromosome contains Peak SNP;
##'              (3) 4in1 Regional Manhattan Plot of the Peak SNP;
##'              (4) 4in1 Genotype Boxplot of the Peak SNP.
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
##' ADDO_AddDom3_Plot(outdir=outdir, covariates_sum=2, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000)
##'
##' @author Leilei Cui and Bin Yang

ADDO_AddDom3_Plot <- function(outdir = outdir,
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
system("mkdir 3_Plot 3_Plot/1_Manhattan 3_Plot/2_QQ_Plot 3_Plot/2_QQ_Plot/1_Whole 3_Plot/2_QQ_Plot/2_Rmchr 3_Plot/3_PeakSNP_RegionMan 3_Plot/4_PeakSNP_BoxPlot")
outdir = paste(outdir,"/3_Plot",sep="")

read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))
read.ped <- function(...) as.data.frame(fread(header=F,colClasses="character",...))

Phe_Clean = read.header(paste(rawdir,"/../2_PheClean/Phe_Clean.txt",sep="")); rownames(Phe_Clean)=Phe_Clean[,1]
Phe_Nor = read.header(paste(rawdir,"/../4_PheNor/Phe_Nor.txt",sep="")); rownames(Phe_Nor)=Phe_Nor[,1]

if(PheList_Choose == T){ PheList = PheList[PheList %in% colnames(Phe_Nor)] }else{ PheList = colnames(Phe_Nor) }

PeakSNPs_all = list("AvsNull"=c(),"DvsNull"=c(),"ADvsNull"=c(),"ADvsA"=c()); traitlamb_all = c(); traitlamb_all_rm = c()

for(i in (covariates_sum+2):length(PheList)){
try({

	phenotype_name=PheList[i]
		
	tmp_map = read(paste(rawdir,"/Data_clean.bim",sep="")); rownames(tmp_map) = tmp_map[,2]	
	setwd(Pdir); system(paste('tar zxvf Pvalue_',phenotype_name,'.txt.tar.gz',sep=''))
	TEMP_logP_raw = read.header(paste('Pvalue_',phenotype_name,'.txt',sep=''))
	rownames(TEMP_logP_raw) = tmp_map[,2]
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
		TEMP_logP = TEMP_logP_sampling[order(TEMP_logP_sampling[,1],TEMP_logP_sampling[,2]),]
	}
		
	if(!file.exists(paste(outdir,'/1_Manhattan/Manhattan_',phenotype_name,'.png',sep=''))){
	try({
		print(date()); setwd(paste(outdir,"/1_Manhattan",sep=""))
			
		Region_TEMP_logP = cbind(paste("SNP",c(1:dim(TEMP_logP)[1]),sep="_"),TEMP_logP); colnames(Region_TEMP_logP)[1] = "SNP"
		Region_TEMP_logP = subset(Region_TEMP_logP,Region_TEMP_logP[,2]!=0 | Region_TEMP_logP[,3]!=0)
		Region_TEMP_logP = Region_TEMP_logP[complete.cases(Region_TEMP_logP),]
		thrgenome = -log10(0.05/nrow(TEMP_logP_all)); thrsuggest = -log10(1/nrow(TEMP_logP_all))

		logP_AddCode = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,4])
		logP_DomCode = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,5])
		logP_AddDomCode_Null = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,6])
		logP_AddDomCode_Add = data.frame(SNP=Region_TEMP_logP[,1],chr=Region_TEMP_logP[,2], pos=Region_TEMP_logP[,3],Pc1df1=Region_TEMP_logP[,7])		
			
		y_limit_tmp = c(logP_AddCode[,4],logP_DomCode[,4],logP_AddDomCode_Null[,4],logP_AddDomCode_Add[,4])
		y_limit = max(y_limit_tmp[!is.na(y_limit_tmp)&!is.infinite(y_limit_tmp)])

		png(paste("Manhattan_",phenotype_name,".png",sep=""),width=30000,height=12000,res=600,type="cairo")
		layout(cbind(c(1,3),c(2,4)),widths=rep(1,2),heights=rep(1,2))
		par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=3.5)
		
		for(mh_sort in 1:4){
		try({
			mkrPos = plotGWAS(chrs=chrs_sum,traitIdx=1,alldat=na.omit(list(logP_AddCode, logP_DomCode, logP_AddDomCode_Null, logP_AddDomCode_Add)[[mh_sort]]),y_limit=y_limit)
			title(paste(phenotype_name,c("Additive Code VS Null","Dominant Code VS Null","Additive-Dominant Code VS Null","Additive-Dominant Code VS Additive Code")[mh_sort],sep=" : "))
			abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			if(mh_sort ==4){dev.off()}
		})
		}
		
		print(paste(phenotype_name,"Manhattan Plot Done","^_^",sep=" "))
	})
	}else{
		print(paste(phenotype_name," : Manhattan Plot Already done! ^_^",sep=""))
	}
		
	if(!(file.exists(paste(outdir,'/2_QQ_Plot/1_Whole/QQplot_',phenotype_name,'.png',sep='')) & file.exists(paste(outdir,'/2_QQ_Plot/2_Rmchr/QQplot_rmqtl_',phenotype_name,'.png',sep='')))){
	try({

		print(date()); setwd(paste(outdir,"/2_QQ_Plot/1_Whole",sep=""))
				
		png(paste("QQplot_",phenotype_name,".png",sep=""),width=3000,height=3000,res=200,type="cairo")
		layout(matrix(c(1,3,2,4),nrow=2),widths=rep(1,2),heights=rep(1,2))
		par(mai=c(0.95,0.95,1,0.5),cex.main=1)

		traitlamb_all_tmp = phenotype_name
		for(qq_sort in 1:4){
		try({				
			data_qq = TEMP_logP_all[,c(3,4,5,6)[qq_sort]]
			object_statistic = c("-logP: Add_Code VS Null","-logP: Dom_Code VS Full","-logP: AddDome_Code VS Null","-logP: AddDome_Code VS Add_Code")[qq_sort]
			data_qq <- as.numeric(data_qq); data_qq <- data_qq[!is.na(data_qq)]
			expPvals = sort(-log10((1:length(data_qq))/length(data_qq))); obsPvals = sort(data_qq)
			lamb = median(obsPvals)/median(expPvals)
			plot(expPvals,obsPvals,main=paste("QQ Plot : ",phenotype_name,sep=""),xlab=paste("Expected ",object_statistic,sep=""),ylab=paste("Observed ",object_statistic,sep=""),cex.lab=1.5,cex.axis=1.5,cex.main=1.7,cex=1.5)
			lines(expPvals,expPvals,type="l",col="red")
			if(qq_sort ==4){dev.off()}
			
			traitlamb_all_tmp = cbind(traitlamb_all_tmp, lamb)
		})
		}
		
		colnames(traitlamb_all_tmp) = c("Trait","AddCode","DomCode","AddDomCode_Null","AddDomCode_Add")
		traitlamb_all = rbind(traitlamb_all, traitlamb_all_tmp)
		print(paste(phenotype_name,"QQ Plot1 Done","^_^",sep=" "))
		
		setwd(paste(outdir,"/2_QQ_Plot/2_Rmchr",sep=""))

		png(paste("QQplot_rmqtl_",phenotype_name,".png",sep=""),width=3000,height=3000,res=200,type="cairo")
		layout(matrix(c(1,3,2,4),nrow=2),widths=rep(1,2),heights=rep(1,2))
		par(mai=c(0.95,0.95,1,0.5),cex.main=1)

		traitlamb_all_rm_tmp = phenotype_name
		for(qq_rm_sort in 1:4){
		try({
			X_allhap_clean_tmp = TEMP_logP_all[,c(1,2,c(3,4,5,6)[qq_rm_sort])]; rownames(X_allhap_clean_tmp) = paste("SNP",1:nrow(X_allhap_clean_tmp),sep="_")
			Peak_tmp = subset(X_allhap_clean_tmp,X_allhap_clean_tmp[,3]==max(X_allhap_clean_tmp[,3],na.rm=T))
			X_allhap_rm = subset(X_allhap_clean_tmp, X_allhap_clean_tmp[,1] == Peak_tmp[1,"chr"])
			#X_allhap_rm = subset(X_allhap_rm, X_allhap_rm[,2] <= (Peak_tmp[1,"pos"]+RegionMan_chr_region) & X_allhap_rm[,2] >= (Peak_tmp[1,"pos"]-RegionMan_chr_region))
			X_allhap_clean = X_allhap_clean_tmp[!(rownames(X_allhap_clean_tmp) %in% rownames(X_allhap_rm)),]
			
			data_qq = X_allhap_clean[,3]
			object_statistic = c("-logP: Add_Code VS Null","-logP: Dom_Code VS Full","-logP: AddDome_Code VS Null","-logP: AddDome_Code VS Add_Code")[qq_rm_sort]
			data_qq <- as.numeric(data_qq); data_qq <- data_qq[!is.na(data_qq)]
			expPvals = sort(-log10((1:length(data_qq))/length(data_qq))); obsPvals = sort(data_qq)
			lamb = median(obsPvals)/median(expPvals)
			plot(expPvals,obsPvals,main=paste("QQ Plot : ",phenotype_name,sep=""),xlab=paste("Expected ",object_statistic,sep=""),ylab=paste("Observed ",object_statistic,sep=""),cex.lab=1.5,cex.axis=1.5,cex.main=1.7,cex=1.5)
			lines(expPvals,expPvals,type="l",col="red")
			if(qq_rm_sort ==4){dev.off()}
			
			traitlamb_all_rm_tmp = cbind(traitlamb_all_rm_tmp, lamb)
		})
		}
		
		colnames(traitlamb_all_rm_tmp) = c("Trait","AddCode_rm","DomCode_rm","AddDomCode_Null_rm","AddDomCode_Add_rm")
		traitlamb_all_rm = rbind(traitlamb_all_rm, traitlamb_all_rm_tmp)
		print(paste(phenotype_name,"QQ Plot2 Done","^_^",sep=" "))
	})
	}else{
		print(paste(phenotype_name," : QQ Plot1&2 Already done! ^_^",sep=""))
	}
		
	if(!file.exists(paste(outdir,'/3_PeakSNP_RegionMan/ManhattanRegion_',phenotype_name,'.png',sep=''))){
	try({
		print(date()); setwd(paste(outdir,"/3_PeakSNP_RegionMan",sep=""))  
	
		png(paste("ManhattanRegion_",phenotype_name,".png",sep=""),width=14000,height=7000,res=600,type="cairo")
		layout(cbind(c(1,3),c(2,4)),widths=rep(1,2),heights=rep(1,2))
		par(mai=c(1,1,0.5,0.5))
			
		for(rm_sort in 1:4){
		try({

			Choose_Pvalue_col = c(3,4,5,6)[rm_sort]
			Title_used = c(") : Additive Code VS Null",") : Dominant Code VS Null",") : Additive-Dominant Code VS Null",") : Additive-Dominant Code VS Additive Code")[rm_sort]
			tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,Choose_Pvalue_col] == max(TEMP_logP[,Choose_Pvalue_col]))[1,1:2]; ch <- tmp_chr_pos[1,1]
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
			
			setwd(paste(outdir,"/3_PeakSNP_RegionMan",sep=""))
			plotdat = data.frame(SNP=regiondat_used$SNP,chr=regiondat_used$chr,pos=regiondat_used$pos,logp=regiondat_used$logP)
			thrgenome = -log10(0.05/dim(TEMP_logP_all)[1]); thrsuggest = -log10(1/dim(TEMP_logP_all)[1])				
			ylim_tmp = subset(regiondat_used, regiondat_used[,"chr"]==tmp_chr_pos[1,"chr"] & regiondat_used[,"pos"]==tmp_chr_pos[1,"pos"])[1,"logP"]*1.15
			plotRegion(chrs=ch, alldat=plotdat, traitIdx=1, from=SNPh, to=SNPt, ldinfo=ldinfo, ylim=c(0,ylim_tmp))
			title(paste0(phenotype_name, "(",tmp_chr_pos[1,1],", ",tmp_chr_pos[1,2],Title_used), cex.main=1.6)
			abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
			
			if(rm_sort == 4){dev.off()}
		})
		}
		
		print(paste(phenotype_name,"Region Manhattan Plot Done","^_^",sep=" "))
	})
	}else{
		print(paste(phenotype_name," : Region Manhattan Plot Already done! ^_^",sep=""))
	}

		
	if(!file.exists(paste(outdir,'/4_PeakSNP_BoxPlot/Boxplot_',phenotype_name,'.png',sep=''))){
	try({
		print(date()); setwd(paste(outdir,"/4_PeakSNP_BoxPlot",sep=""))
		
		png(paste("Boxplot_",phenotype_name,".png",sep=""),width=21000,height=18900,res=600,type="cairo")
		dist = 0.005
		m1 <- rbind(c(0+dist, 0.35-dist, 0.5+dist, 1-dist), c(0.35+dist, 0.5-dist, 0.5+0.5/3*2+dist, 1-dist), c(0.35+dist, 0.5-dist, 0.5+0.5/3*1+dist, 0.5+0.5/3*2-dist), c(0.35+dist, 0.5-dist, 0.5+dist, 0.5+0.5/3*1-dist))
		m2 <- rbind(c(0.5+dist, 0.85-dist, 0.5+dist, 1-dist), c(0.85+dist, 1-dist, 0.5+0.5/3*2+dist, 1-dist), c(0.85+dist, 1-dist, 0.5+0.5/3*1+dist, 0.5+0.5/3*2-dist), c(0.85+dist, 1-dist, 0.5+dist, 0.5+0.5/3*1-dist))
		m3 <- rbind(c(0+dist, 0.35-dist, 0+dist, 0.5-dist), c(0.35+dist, 0.5-dist, 0.5/3*2+dist, 0.5-dist), c(0.35+dist, 0.5-dist, 0.5/3*1+dist, 0.5/3*2-dist), c(0.35+dist, 0.5-dist, 0+dist, 0.5/3*1-dist))
		m4 <- rbind(c(0.5+dist, 0.85-dist, 0+dist, 0.5-dist), c(0.85+dist, 1-dist, 0.5/3*2+dist, 0.5-dist), c(0.85+dist, 1-dist, 0.5/3*1+dist, 0.5/3*2-dist), c(0.85+dist, 1-dist, 0+dist, 0.5/3*1-dist))
		m_all = rbind(m1,m2,m3,m4); rownames(m_all) = 1:16; split.screen(m_all)
		
		for(bx_sort in 1:4){
		try({
			Choose_Pvalue_col = c(3,4,5,6)[bx_sort]; Title_used = c(": AvsNull",": DvsNull",": ADvsNull",": ADvsA")[bx_sort]
			m_used = m_all[((bx_sort*4)-3):(bx_sort*4),]; screen_num = as.numeric(rownames(m_used))
			tmp_chr_pos = subset(TEMP_logP,TEMP_logP[,Choose_Pvalue_col] == max(TEMP_logP[,Choose_Pvalue_col]))[1,1:2]
			snp_name = as.character(subset(tmp_map,tmp_map[,1]==tmp_chr_pos[1,1] & tmp_map[,4]==tmp_chr_pos[1,2])[1,2])

		
			setwd(rawdir); system(paste("plink --noweb --silent --bfile Data_clean --snp",snp_name,"--recode --out Boxplot_tmp"))
			TEMP.geno <- read.ped("Boxplot_tmp.ped"); system("rm Boxplot_tmp.*")
			TEMP_geno =  as.data.frame(paste(TEMP.geno[,7],TEMP.geno[,8],sep="/"))
			rownames(TEMP_geno) = TEMP.geno[,2]; colnames(TEMP_geno) = snp_name
			TEMP_phe = Phe_Clean[,phenotype_name]
			TEMP_gp = data.frame(TEMP_geno,TEMP_phe); TEMP_gp = na.omit(TEMP_gp); TEMP_gp = subset(TEMP_gp, TEMP_gp[,1]!="0/0")
			phe = TEMP_gp[,2]; geno = TEMP_gp[,1]; unigeno = sort(unique(geno))
			PheGeno_clean = as.data.frame(cbind(phe=phe,geno=geno)); PheGeno_clean[,"phe"]=as.numeric(PheGeno_clean[,"phe"])
			
			screen(screen_num[1]); par(mai=c(0.5,0.9,0.9,0.01),mgp=c(4.1,1.5,0))
			boxplot(PheGeno_clean[,1]~PheGeno_clean[,2],cex.axis=3.5,ylim=c(min(PheGeno_clean[,1],na.rm=T),1.2*max(PheGeno_clean[,1],na.rm=T)),col=rainbow(12))
			y = PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]]; y = y[!is.na(y)]
			points(jitter(rep(1,length(which(PheGeno_clean[,2]==unigeno[1]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[1]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(2,length(which(PheGeno_clean[,2]==unigeno[2]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[2]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			points(jitter(rep(3,length(which(PheGeno_clean[,2]==unigeno[3]))),factor=0.7),PheGeno_clean[,1][PheGeno_clean[,2]==unigeno[3]],pch=19,col=rgb(0,0,1,alpha=0.4),cex=3)
			text(1,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[1])),sep=""),cex=4)
			text(2,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[2])),sep=""),cex=4)
			text(3,1.1*max(PheGeno_clean[,1],na.rm=T),paste("n=",length(which(PheGeno_clean[,2]==unigeno[3])),sep=""),cex=4)
			text(2,1.18*max(PheGeno_clean[,1],na.rm=T),paste0("(Chr:",tmp_chr_pos[1,1]," Pos:",tmp_chr_pos[1,2],")"),cex=3)
			title(paste0(phenotype_name,Title_used),cex.main=3.3)
			
			screen(screen_num[2]); par(mai=c(0.7,0.9,0.5,0.1),mgp=c(2.7,1,0))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[1])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[1],col=c("red"),breaks=150,cex.lab=2,cex.axis=1.7,cex.main=2.6,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(screen_num[3]); par(mai=c(0.7,0.9,0.5,0.1),mgp=c(2.7,1,0))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[2])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[2],col=c("red"),breaks=150,cex.lab=2,cex.axis=1.7,cex.main=2.6,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			screen(screen_num[4]); par(mai=c(0.7,0.9,0.5,0.1),mgp=c(2.7,1,0))
			data_histogram = subset(PheGeno_clean,PheGeno_clean[,2]==unigeno[3])[,1]; data_histogram <- data_histogram[!is.na(data_histogram)]
			h <- hist(data_histogram,main=unigeno[3],col=c("red"),breaks=150,cex.lab=2,cex.axis=1.7,cex.main=2.6,xlim=c(min(PheGeno_clean[,1],na.rm=T),max(PheGeno_clean[,1],na.rm=T)),xlab=phenotype_name)
			xfit = seq(min(data_histogram),max(data_histogram),length=1000)
			yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
			yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
			lines(xfit,yfit,col="blue",lwd=3)	
			box()
			
			if(bx_sort == 4){close.screen(all.screens = TRUE); dev.off()}
			
			PeakSNPs = subset(TEMP_logP,TEMP_logP[,Choose_Pvalue_col] == max(TEMP_logP[,Choose_Pvalue_col]))[1,] # Choose the first SNP with the Peak P-value #
			PeakSNPs_output = cbind(Trait=phenotype_name,SNP=rownames(PeakSNPs),PeakSNPs)
			PeakSNPs_all[[bx_sort]] = rbind(PeakSNPs_all[[bx_sort]],PeakSNPs_output)			
		})
		}
		
		print(paste(phenotype_name,"Genotype Box Plot Done","^_^",sep=" "))
	})	
	}else{
		print(paste(phenotype_name," : Genotype Box Plot Already done! ^_^",sep=""))
	}
		
})
}

setwd(outdir); save(PeakSNPs_all,file="PeakSNPs_Summary.Rdata")
write.table(cbind(as.data.frame(traitlamb_all),as.data.frame(traitlamb_all_rm)[,-1]),file="Traitlambs_Summary.txt",row.names=F,col.names=T,quote=F)

}

