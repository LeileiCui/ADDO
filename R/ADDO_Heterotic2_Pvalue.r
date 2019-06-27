
##' @title Pvalue calculation and Verfication of Overdominance QTLs
##'        (STEP2 of The Heterotic Model)
##' 
##' @description P-value Calculation and Verfication of overdominance (or heterotic) QTLs
##'              (1) Run the whole PLINK file or Run the separated PLINK files ("Run_separated = F");
##'              (2) Indicate Reocde Model (AA:1 0 0/AB:0 1 0/BB:0 0 1) without the "1" column of covariance matrix;
##'              (3) Estimate two T-statistics (t(AB-AA) and t(AB-BB)) to measure the deviation between the effect of heterozygote (AB) and that of two homozygotes (AA and BB);
##'              (4) Generate the P-value based on MVN distribution using the minor(abs(t(AB-AA)),abs(t(AB-BB))) from SNPs with t(AB-AA)*t(AB-BB)>0.
##'
##' @param indir A character. The input directory where contains the input bPLINK or GenABEL data.
##' @param outdir A character. The output directory where generates the folder: "2_Pvalue".
##' @param Input_name A character. The prefixes of the input files.
##' @param Kinship_type A character. The method to generate kinship matrix. Please select from "GenABEL","EMMA","EMMAX","GEMMA", "GCTA", "GCTA_ad", "HOMEBREW_AFW" or "HOMEBREW_AS".
##' @param VarComponent_Method A character. The method to estimate variance components. Please select from "EMMA_a", "GCTA_a" or "GCTA_ad" (When VarComponent_Method is "GCTA_ad", the Kinship_type must be "GCTA_ad").
##' @param PheList_Choose A logic variable. T: Just investigate specified phenotypes; F: Investigate all phenotypes.
##' @param PheList A vector of character. When choose "PheList_Choose=F", the specified phenotype list must be specified.
##' @param Phe_HistogramPlot A logic variable. T: Draw the histogram plots for all phenotypes; F: Aviod the histogram plots.
##' @param Run_separated A logic variable. T: Run the separated genotype files; F: Run the whole genotype file.
##' @param covariates_sum A numeric variable. The sum of all covariates.
##' @param Phe_IndMinimum A numeric variable. Remove phenotypes without enough available individuals.
##' @param GT_IndMinimum A numeric variable. Remove loci with available individuals <GT_IndMinimum for all three genotypes.
##' @param num_nodes A numeric variable. The number of cores used parallelly.
##'
##' @return a folder named "2_Pvalue" with various statistics of each significant SNP for all phenotypes.
##'
##' @examples 
##' ADDO_Heterotic2_Pvalue(indir=indir, outdir=outdir, Input_name="TEST", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", covariates_sum=2)
##'
##' @author Leilei Cui and Bin Yang

ADDO_Heterotic2_Pvalue <- function(indir = indir,
                                   outdir = outdir, 
                                   Input_name = Input_name, 
                                   Kinship_type = Kinship_type,
                                   VarComponent_Method = VarComponent_Method,
                                   PheList_Choose = F, 
                                   PheList = PheList,
                                   Phe_HistogramPlot = F,
                                   Run_separated = F,
                                   covariates_sum = covariates_sum, 
                                   Phe_IndMinimum = 200,
                                   GT_IndMinimum = 10, 
                                   num_nodes = 10){

setwd(outdir); system("mkdir 2_Pvalue")

Plot_Phe_Histogram <- function(x,Phe_data,Phe_names){
try({
	phe = Phe_data[,x]; Phe_type = Phe_names[x]
	data_histogram = as.numeric(phe)
	data_histogram <- data_histogram[!is.na(data_histogram)]
	
	png(paste("HisSca_Plot_",Phe_type,".png",sep=""),width=1500,height=1500,res=200,type="cairo")
	par(mai=c(0.95,0.95,1,0.5))
	h <- hist(data_histogram,main=paste("Histogram Plot : ",Phe_type,sep=""),col=c("red"),breaks=170,cex.lab=1.5,cex.axis=1.5,cex.main=1.7)
	xfit = seq(min(data_histogram),max(data_histogram),length=1000)
	yfit = dnorm(xfit,mean=mean(data_histogram),sd=sd(data_histogram))
	yfit = yfit*diff(h$mids[1:2])*length(data_histogram)
	lines(xfit,yfit,col="blue",lwd=3)	
	box()
	dev.off()
})
}

MVN_Tvalue2Pvalue <- function(i,Tvalues,mean1,mean2,Sigma){
try({	
	x = Tvalues[i]
	area1 = pmvnorm(lower=c(abs(x),abs(x)), upper=c(Inf,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
	area2 = pmvnorm(lower=c(-Inf,-Inf), upper=c(-abs(x),-abs(x)), mean=c(mean1,mean2), sigma=Sigma)[1]
	Area1 = pmvnorm(lower=c(0,0), upper=c(Inf,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
	Area2 = pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0), mean=c(mean1,mean2), sigma=Sigma)[1]
	return( -log10((area1+area2)/(Area1+Area2)) )
})
}

MVN_Tvalue2Pvalue_new <- function(i,Tvalues_new,mean1,mean2,Sigma){
try({	
	x1 = Tvalues_new[i,1]; x2 = Tvalues_new[i,2]
	if(x1>0 && x2>0){
		area1 = pmvnorm(lower=c(x1,x2), upper=c(Inf,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
		Area1 = pmvnorm(lower=c(0,0), upper=c(Inf,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
		return( -log10(area1/Area1) )
	}
	if(x1<0 && x2<0){
		area2 = pmvnorm(lower=c(-Inf,-Inf), upper=c(x1,x2), mean=c(mean1,mean2), sigma=Sigma)[1]
		Area2 = pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0), mean=c(mean1,mean2), sigma=Sigma)[1]
		return( -log10(area2/Area2) )
	}
	if(x1>0 && x2<0){
		area1 = pmvnorm(lower=c(x1,-Inf), upper=c(Inf,x2), mean=c(mean1,mean2), sigma=Sigma)[1]
		Area1 = pmvnorm(lower=c(0,-Inf), upper=c(Inf,0), mean=c(mean1,mean2), sigma=Sigma)[1]
		return( -log10(area1/Area1) )
	}
	if(x1<0 && x2>0){
		area2 = pmvnorm(lower=c(-Inf,x2), upper=c(x1,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
		Area2 = pmvnorm(lower=c(-Inf,0), upper=c(0,Inf), mean=c(mean1,mean2), sigma=Sigma)[1]
		return( -log10(area2/Area2) )
	}
})
}

NORM_Tvalue2Pvalue <- function(i,Tvalues,DF_res){
try({
	x = Tvalues[i]; x_df = DF_res[i]
	return( -log10(2*pt(abs(x),x_df,lower.tail=F,log.p=F)) )
})
}

read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))
read.ped <- function(...) as.data.frame(fread(header=F,colClasses="character",...))

rawdir = paste(outdir,"/1_PheGen/0_Data",sep="")
gkin_order = read(paste(rawdir,"/Data_clean.fam",sep=""))[,2]
	
setwd(rawdir); myfolder_vc = "Data_clean_vc"
if(myfolder_vc %in% dir() == FALSE){
	dir.create(myfolder_vc)
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm --out ",myfolder_vc,"/TEMP",sep=""))
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm-d --out ",myfolder_vc,"/TEMP",sep=""))
	try({system("rm TEMP.hh TEMP.log")})
}
	
if(Kinship_type %in% c("GCTA_ad","HOMEBREW_AFW","HOMEBREW_AS")){
	if(Kinship_type == "GCTA_ad"){ tmp_matrix_name = "GCTA" }
	if(Kinship_type == "HOMEBREW_AFW"){ tmp_matrix_name = "AFW" }
	if(Kinship_type == "HOMEBREW_AS"){ tmp_matrix_name = "AS" }
	gkin_a = gkin = read(paste(rawdir,"/Data_clean_kinship_",tmp_matrix_name,".Additive.txt",sep=""))
	gkin_d = read(paste(rawdir,"/Data_clean_kinship_",tmp_matrix_name,".Dominant.txt",sep=""))
	rownames(gkin) = gkin_order; colnames(gkin) = gkin_order; gkin = as.matrix(gkin)
	rownames(gkin_a) = gkin_order; colnames(gkin_a) = gkin_order; gkin_a = as.matrix(gkin_a)
	rownames(gkin_d) = gkin_order; colnames(gkin_d) = gkin_order; gkin_d = as.matrix(gkin_d)
}else{
	if(Kinship_type == "GenABEL"){gkin = read(paste(rawdir,"/Data_clean_kinship_GenABEL.txt",sep="")); gkin[upper.tri(gkin)] = t(gkin)[upper.tri(gkin)]}
	if(Kinship_type == "EMMA"){gkin = read(paste(rawdir,"/Data_clean_kinship_EMMA.Additive.txt",sep=""))}
	if(Kinship_type == "EMMAX"){gkin = read(paste(rawdir,"/Data_clean_kinship_EMMAX.hIBS.txt",sep=""))}
	if(Kinship_type == "GEMMA"){gkin = read(paste(rawdir,"/Data_clean_kinship_GEMMA.cXX.txt",sep=""))}
	if(Kinship_type == "GCTA"){gkin = read(paste(rawdir,"/Data_clean_kinship_GCTA.Additive.txt",sep=""))}
	rownames(gkin) = gkin_order; colnames(gkin) = gkin_order; gkin = as.matrix(gkin)
}

Phe_Nor = read.header(paste(rawdir,"/../4_PheNor/Phe_Nor.txt",sep="")); rownames(Phe_Nor)=Phe_Nor[,1]

if(PheList_Choose == T){ PheList = PheList[PheList %in% colnames(Phe_Nor)] }else{ PheList = colnames(Phe_Nor) }

Phe_Trans = as.data.frame(matrix(NA,length(gkin_order),length(PheList)))
rownames(Phe_Trans) = gkin_order; colnames(Phe_Trans) = PheList
Phe_Trans[,1:(covariates_sum+1)] = Phe_Nor[rownames(Phe_Trans),1:(covariates_sum+1)]
VC_result = c()

for(args in (covariates_sum+2):length(PheList)){
try({

	phenotype_name=PheList[args]

	if(!file.exists(paste(outdir,'/2_Pvalue/Pvalue_',phenotype_name,'.txt.tar.gz',sep=''))){
		
		keep_inds = rownames(Phe_Nor)[!is.na(Phe_Nor[,phenotype_name])]
		y.used = Phe_Nor[keep_inds,phenotype_name]; names(y.used) = keep_inds
		
		if(VarComponent_Method == "GCTA_ad"){
			kinship.a.used = gkin_a[keep_inds,keep_inds]; kinship.d.used = gkin_d[keep_inds,keep_inds]
			
			setwd(paste0(rawdir,"/",myfolder_vc))
			tmp_id = read("TEMP.grm.id"); tmp_pheno = cbind(tmp_id,NA); rownames(tmp_pheno) = tmp_pheno[,2]; tmp_pheno[names(y.used),3] = y.used
			write.table(tmp_pheno,file=paste("tmp_",phenotype_name,".pheno",sep=""),row.names=F,col.names=F,quote=F) # Prepare "tmp.pheno" # 
			write.table(as.data.frame(c("TEMP","TEMP.d")),file=paste("tmp_",phenotype_name,".txt",sep=""),row.names=F,col.names=F,quote=F) # Prepare "tmp.txt" #
			system(paste("gcta64 --reml --mgrm tmp_",phenotype_name,".txt --pheno tmp_",phenotype_name,".pheno --out tmp_",phenotype_name,sep=""))
			
			tmp_hsq = read.table(paste("tmp_",phenotype_name,".hsq",sep=""),header=T,fill=T)
			rownames(tmp_hsq) = tmp_hsq[,1]
			tmp_result = c("trait" = phenotype_name,"V_Add" = as.numeric(as.character(tmp_hsq["V(G1)","Variance"])),"V_Dom" = as.numeric(as.character(tmp_hsq["V(G2)","Variance"])),"V_Envi" = as.numeric(as.character(tmp_hsq["V(e)","Variance"])), "SE_Add" = as.numeric(as.character(tmp_hsq["V(G1)","SE"])),"SE_Dom" = as.numeric(as.character(tmp_hsq["V(G2)","SE"])),"SE_Envi" = as.numeric(as.character(tmp_hsq["V(e)","SE"]))); system(paste("rm tmp_",phenotype_name,"*",sep=""))
			VC_result = rbind(VC_result,tmp_result)
			V = as.numeric(tmp_result["V_Add"])*kinship.a.used + as.numeric(tmp_result["V_Dom"])*kinship.d.used + as.numeric(tmp_result["V_Envi"])*diag(length(y.used))
		}
		
		if(VarComponent_Method == "GCTA_a"){
			kinship.used = gkin[keep_inds,keep_inds]
			
			setwd(paste0(rawdir,"/",myfolder_vc))
			tmp_id = read("TEMP.grm.id"); tmp_pheno = cbind(tmp_id,NA); rownames(tmp_pheno) = tmp_pheno[,2]; tmp_pheno[names(y.used),3] = y.used
			write.table(tmp_pheno,file=paste("tmp_",phenotype_name,".pheno",sep=""),row.names=F,col.names=F,quote=F)
			write.table(as.data.frame(c("TEMP")),file=paste("tmp_",phenotype_name,".txt",sep=""),row.names=F,col.names=F,quote=F)
			system(paste("gcta64 --reml --mgrm tmp_",phenotype_name,".txt --pheno tmp_",phenotype_name,".pheno --out tmp_",phenotype_name,sep=""))
			
			tmp_hsq = read.table(paste("tmp_",phenotype_name,".hsq",sep=""),header=T,fill=T)
			rownames(tmp_hsq) = tmp_hsq[,1]
			tmp_result = c("trait" = phenotype_name,"V_Gene" = as.numeric(as.character(tmp_hsq["V(G)","Variance"])),"V_Envi" = as.numeric(as.character(tmp_hsq["V(e)","Variance"])), "SE_Gene" = as.numeric(as.character(tmp_hsq["V(G)","SE"])), "SE_Envi" = as.numeric(as.character(tmp_hsq["V(e)","SE"]))); system(paste("rm tmp_",phenotype_name,"*",sep=""))
			VC_result = rbind(VC_result,tmp_result)
			V = as.numeric(tmp_result["V_Gene"])*kinship.used + as.numeric(tmp_result["V_Envi"])*diag(length(y.used))
		}
		
		if(VarComponent_Method == "EMMA_a"){
			kinship.used = gkin[keep_inds,keep_inds]
			emma = try(emma.REMLE(y=y.used, X=as.matrix(rep(1,length(y.used))), K=kinship.used, ngrids=100, llim=-10, ulim=10, esp=1e-10, eig.R=NULL))
			tmp_result = cbind("trait" = phenotype_name,do.call(cbind,emma))
			VC_result = rbind(VC_result,tmp_result)
			V = emma$vg*kinship.used+emma$ve*diag(length(y.used))
		}		
		
		#multiplicator.used = solve(svd(V)$u %*% diag(x=sqrt(svd(V)$d)) %*% t(svd(V)$u)) 
		#multiplicator.used = solve(eigen(V)$vectors %*% diag(x=sqrt(eigen(V)$values)) %*% solve(eigen(V)$vectors)) 
		
		svd=try(svd(V)); if (inherits(svd, "try-error")) {print('error svd');next}
		eigen_basis=svd$u; eigen_values=svd$d 
		teigen_basis=t(eigen_basis); sqrt_lambda_matrix=diag(x=sqrt(eigen_values)) 
		r_mat=eigen_basis%*%sqrt_lambda_matrix%*%teigen_basis
		multiplicator.used=solve(r_mat) 
		#multiplicator.used = eigen_basis %*% diag(1/diag(sqrt_lambda_matrix)) %*% teigen_basis 
		
		Phe_Trans[keep_inds,phenotype_name] = multiplicator.used %*% y.used

		
		GWAS_Chip <- function(i){
		try({
			
			if(length(y.used) >= Phe_IndMinimum){
			
				x = Chip_ped[,(i*2+5):(i*2+6)]; x = x[keep_inds,]
				x[x==1] = "A"; x[x==2] = "B" 

				tmp_ind_geno = as.data.frame(cbind(c(1:nrow(x)),paste(x[,1],x[,2],sep="")))
				rownames(tmp_ind_geno) = rownames(x)				
				tmp_ind_geno_rm00 = subset(tmp_ind_geno,tmp_ind_geno[,2]!="00")
				
				if(length(unique(tmp_ind_geno_rm00[,2])) == 3){
					
					if(min(table(as.data.frame(tmp_ind_geno_rm00[,2]))) >= GT_IndMinimum){
						
						tmp_ind_geno[,2] = as.character(tmp_ind_geno[,2]); tmp_ind_geno = data.frame(tmp_ind_geno,0,0,0,0,0)
						colnames(tmp_ind_geno) = c("Sort","Genotype","Add_Code","Dom_Code","Indicate_Code1","Indicate_Code2","Indicate_Code3")
						tmp_ind_geno[tmp_ind_geno[,2] == "AA",c("Indicate_Code1")] = 1
						tmp_ind_geno[tmp_ind_geno[,2] %in% c("AB","BA"),c("Indicate_Code2")] = 1
						tmp_ind_geno[tmp_ind_geno[,2] == "BB",c("Indicate_Code3")] = 1

						if(nrow(subset(tmp_ind_geno,tmp_ind_geno[,2]=="00")) > 0){
							rmna_ind_geno_num = as.numeric(subset(tmp_ind_geno,tmp_ind_geno[,2]=="00")[,1]) 
							rmna_y.used = y.used[-rmna_ind_geno_num]
							rmna_ind_geno = tmp_ind_geno[-rmna_ind_geno_num,]
							rmna_multiplicator.used = as.matrix(multiplicator.used[-rmna_ind_geno_num,-rmna_ind_geno_num])
						}else{
							rmna_y.used = y.used 
							rmna_ind_geno = tmp_ind_geno
							rmna_multiplicator.used = as.matrix(multiplicator.used)
						}
						
						rmna_yt = rmna_multiplicator.used %*% rmna_y.used
						rmna_covt = rmna_multiplicator.used %*% rep(1,length(rmna_y.used))
						rmna_genot = rmna_multiplicator.used %*% as.matrix(rmna_ind_geno[,3:7])
									
						fit_null = lm(rmna_yt ~ -1 + rmna_covt)
						fit_IndiCode = lm(rmna_yt ~ -1 + rmna_genot[,3:5])
						
						fit_IndiCode_coef = summary(fit_IndiCode)$coefficients
						betaAA = fit_IndiCode_coef[1,"Estimate"]; betaAB = fit_IndiCode_coef[2,"Estimate"]; betaBB = fit_IndiCode_coef[3,"Estimate"]
						
						V_Indi = vcov(fit_IndiCode)
						Z1 = matrix(c(rep(0,(dim(V_Indi)[1]-3)),-1,1,0),1); Z2 = t(Z1); varAB_AA = Z1 %*% V_Indi %*% Z2
						Z3 = matrix(c(rep(0,(dim(V_Indi)[1]-3)),0,1,-1),1); Z4 = t(Z3); varAB_BB = Z3 %*% V_Indi %*% Z4
						
						tAB_AA = (betaAB-betaAA)/sqrt(varAB_AA); tAB_BB = (betaAB-betaBB)/sqrt(varAB_BB)
						logP_Indi = -log10(anova(fit_null,fit_IndiCode)[,'Pr(>F)'][2])
						
						return(c(Chip_map[i,1],Chip_map[i,4],logP_Indi,tAB_AA,tAB_BB,betaAA,betaAB,betaBB,varAB_AA,varAB_BB,fit_IndiCode$df.residual))
						
					}else{
						return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA))
					}
				}else{
					return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA))
				}
			}else{
				return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA))
			}
		})
		}
		

		if(Run_separated == F){			
		
			setwd(rawdir); Chip_map_all = Chip_map = read("Data_clean.bim")
			if(!file.exists("Data_clean_recode12.ped")){system("plink --noweb --silent --bfile Data_clean --recode12 --out Data_clean_recode12")}
			
			print(paste0("Whole-Chrs Loading STA: ",date()))
			Chip_ped = read.big.matrix("Data_clean_recode12.ped", sep = " ", type = "integer")
			Chip_tfam <- read(paste(rawdir,"/Data_clean.fam",sep=""))
			options(bigmemory.allow.dimnames=TRUE); rownames(Chip_ped) = Chip_tfam[,2]
			print(paste0("Whole-Chrs Loading END: ",date()))
			
			print(paste0("Whole-Chrs Computing STA: ",date()))
			logP_Chip = do.call('rbind',mclapply(1:((dim(Chip_ped)[2]-6)/2),GWAS_Chip,mc.cores=num_nodes)); rm(Chip_ped); gc()				
			print(paste0("Whole-Chrs Computing END: ",date()))
		
		}else{
		
			setwd(rawdir); myfolder_separated = "Data_clean_separated"			
			Chip_map_all = read("Data_clean.bim"); Part_num = unique(Chip_map_all[,1])
			if(myfolder_separated %in% dir() == FALSE){
				dir.create(myfolder_separated)
				for(i in 1:length(Part_num)){
					system(paste("plink --noweb --silent --bfile Data_clean --chr ",Part_num[i]," --recode12 --out ",myfolder_separated,"/Part",i,sep=""))						
					print(paste0("Part(Chr)",i," Generated!"))
				}
				system(paste("rm ",myfolder_separated,"/*.log",sep=""))	
			}
			
			logP_Chip = data.frame()
			for(j in 1:length(Part_num)){
				Chip_map <- read(paste(myfolder_separated,"/Part",j,".map",sep=""))
				
				print(paste0("Chr",j," Loading STA: ",date()))
				Chip_ped = read.big.matrix(paste(myfolder_separated,"/Part",j,".ped",sep=""), sep = " ", type = "integer")
				Chip_tfam <- read(paste(rawdir,"/Data_clean.fam",sep=""))
				options(bigmemory.allow.dimnames=TRUE); rownames(Chip_ped) = Chip_tfam[,2]
				print(paste0("Chr",j," Loading END: ",date()))
				
				print(paste0("Chr",j," Computing STA: ",date()))
				logP_Chip_tmp = do.call('rbind',mclapply(1:((dim(Chip_ped)[2]-6)/2),GWAS_Chip,mc.cores=num_nodes)); rm(Chip_ped); gc()
				print(paste0("Chr",j," Computing END: ",date()))
				
				logP_Chip = rbind(logP_Chip,logP_Chip_tmp); rm(logP_Chip_tmp); gc()
				print(paste0("Chr",j," Finished!"))
			}
			
		}
		
		colnames(logP_Chip)=c("chr","pos","-logP(IndiCode)","tAB_AA","tAB_BB","betaAA","betaAB","betaBB","varAB_AA","varAB_BB","df_res")
		rownames(logP_Chip) = Chip_map_all[,2]; TEMP_logP_raw = logP_Chip
		TEMP_logP_raw = subset(TEMP_logP_raw,TEMP_logP_raw[,"chr"]!=0 & TEMP_logP_raw[,"pos"]!=0)
		TEMP_logP_raw = TEMP_logP_raw[complete.cases(TEMP_logP_raw),]
		
		AbsMinorT = apply(TEMP_logP_raw[,c("tAB_AA","tAB_BB")],1,function(x){x_abs = abs(x); return(x[x_abs==min(x_abs)])}); AbsMinorT = as.data.frame(AbsMinorT)		
		MVN_logP_Minor = do.call('rbind', mclapply(1:nrow(AbsMinorT), MVN_Tvalue2Pvalue, Tvalues=AbsMinorT[,1], mean1=mean(TEMP_logP_raw[,"tAB_AA"]), mean2=mean(TEMP_logP_raw[,"tAB_BB"]), Sigma=matrix(c(1,cor(TEMP_logP_raw[,"tAB_AA"],TEMP_logP_raw[,"tAB_BB"]),cor(TEMP_logP_raw[,"tAB_AA"],TEMP_logP_raw[,"tAB_BB"]),1),nrow=2), mc.cores=num_nodes))[,1]; names(MVN_logP_Minor) = rownames(AbsMinorT) # Based on tMinor #
		#MVN_logP_Minor = do.call('rbind', mclapply(1:nrow(TEMP_logP_raw), MVN_Tvalue2Pvalue_new, Tvalues_new=TEMP_logP_raw[,c("tAB_AA","tAB_BB")], mean1=mean(TEMP_logP_raw[,"tAB_AA"]), mean2=mean(TEMP_logP_raw[,"tAB_BB"]), Sigma=matrix(c(1,cor(TEMP_logP_raw[,"tAB_AA"],TEMP_logP_raw[,"tAB_BB"]),cor(TEMP_logP_raw[,"tAB_AA"],TEMP_logP_raw[,"tAB_BB"]),1),nrow=2), mc.cores=num_nodes))[,1]; names(MVN_logP_Minor) = rownames(TEMP_logP_raw) # Based on tAB_AA and tAB_BB #
		
		NORM_logP_AB_AA = do.call('rbind', mclapply(1:nrow(TEMP_logP_raw), NORM_Tvalue2Pvalue, Tvalues=TEMP_logP_raw[,"tAB_AA"], DF_res=TEMP_logP_raw[,"df_res"], mc.cores=num_nodes))[,1]; names(NORM_logP_AB_AA) = rownames(TEMP_logP_raw)	
		NORM_logP_AB_BB = do.call('rbind', mclapply(1:nrow(TEMP_logP_raw), NORM_Tvalue2Pvalue, Tvalues=TEMP_logP_raw[,"tAB_BB"], DF_res=TEMP_logP_raw[,"df_res"], mc.cores=num_nodes))[,1]; names(NORM_logP_AB_BB) = rownames(TEMP_logP_raw)

		logP_Chip = cbind(logP_Chip[,1:2],MVN_logP_Minor=NA,NORM_logP_AB_AA=NA,NORM_logP_AB_BB=NA,logP_Chip[,-c(1:2)])
		logP_Chip[rownames(TEMP_logP_raw),"MVN_logP_Minor"] = MVN_logP_Minor
		logP_Chip[rownames(TEMP_logP_raw),"NORM_logP_AB_AA"] = NORM_logP_AB_AA
		logP_Chip[rownames(TEMP_logP_raw),"NORM_logP_AB_BB"] = NORM_logP_AB_BB
	
		setwd(paste(outdir,"/2_Pvalue",sep=""))
		colnames(logP_Chip)=c("chr","pos","MVN_logP_Minor","NORM_logP_AB_AA","NORM_logP_AB_BB","-logP(IndiCode)","tAB_AA","tAB_BB","betaAA","betaAB","betaBB","varAB_AA","varAB_BB","df_res")
		file_name = paste('Pvalue_',phenotype_name,'.txt',sep='')
		write.table(logP_Chip,file=file_name,row.names=F,col.names=T,quote=F)
		system(paste("tar zcvf ",file_name,".tar.gz"," ",file_name,sep=""))
		system(paste("rm ",file_name,sep="")); rm(logP_Chip, TEMP_logP_raw); gc()
		print(paste(phenotype_name,"Done","^_^",sep=" "))

	}else{
		print(paste(phenotype_name,"Done","^_^",sep=" "))
	}

})
}

#setwd(rawdir); system(paste("rm -r",myfolder_separated,myfolder_vc,sep=" "))
setwd(paste(rawdir,"/../5_PheTrans",sep=""))
write.table(Phe_Trans,file="Phe_Trans.txt",row.names=T,col.names=T,quote=F)
write.table(VC_result,file="VC_result.txt",row.names=F,col.names=T,quote=F)
	
if(Phe_HistogramPlot){ Plot_Phe_nor = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Trans,Phe_names=PheList,mc.cores=num_nodes) }
}
