
##' @title Detection and Classification of QTLs with Various Inheritance Categories
##'        (STEP2 of The Add-Dom Model)
##' 
##' @description Select significant QTLs using a matrix operation strategy and estimate the logP of each QTLs with 4 different models as well as their "tAdd" and "tDom" for inheritance categories classification.
##'              (1) Run the whole PLINK file or Run the separated PLINK files if the genotypes are massive ("Run_separated = F");
##'              (2) Additive Recode Model (AA:0/AB:1/BB:2) vs Null Model;
##'              (3) Dominant Recode Model (AA:0/AB:1/BB:0) vs Null Model;
##'              (4) Add+Dom Recode Model (AA:0 0/AB:1 1/BB:2 0) vs Null Model;
##'              (5) Add+Dom Reocde Model (AA:0 0/AB:1 1/BB:2 0) vs Additive Model.
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
##' @param matrix_acceleration A logic variable. T: Implement the matrix acceleration to select significant loci before mixed model; F: Didn't implement the matrix acceleration.
##' @param logP_threshold A numeric variable. The -logP threshold to select significant loci.
##' @param num_nodes A numeric variable. The number of cores used parallelly.
##'
##' @return a folder named "2_Pvalue" with various statistics of each significant SNP for all phenotypes.
##'
##' @examples 
##' ADDO_AddDom2_Pvalue(indir=indir, outdir=outdir, Input_name="TEST", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", covariates_sum=2)
##'
##' @author Leilei Cui and Bin Yang

ADDO_AddDom2_Pvalue <- function(indir = indir,
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
                                matrix_acceleration = T,
                                logP_threshold = 1,
                                num_nodes = 10){

setwd(outdir); system("mkdir 2_Pvalue")

# Function1 #
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

# Function2 #
read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))
read.ped <- function(...) as.data.frame(fread(header=F,colClasses="character",...))

###############
### STAGE 1 ###
###############

rawdir = paste(outdir,"/1_PheGen/0_Data",sep="")
gkin_order = read(paste(rawdir,"/Data_clean.fam",sep=""))[,2]
	

setwd(rawdir); myfolder_vc = "Data_clean_vc"
if(myfolder_vc %in% dir() == FALSE){
	dir.create(myfolder_vc)
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm --out ",myfolder_vc,"/TEMP",sep=""))
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm-d --out ",myfolder_vc,"/TEMP",sep=""))
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
		
		###############
		### STAGE 2 ###
		###############
		
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
		
		# multiplicator.used = solve(svd(V)$u %*% diag(x=sqrt(svd(V)$d)) %*% t(svd(V)$u)) 
		multiplicator.used = solve(eigen(V)$vectors %*% diag(x=sqrt(eigen(V)$values)) %*% solve(eigen(V)$vectors)) 
		
		Phe_Trans[keep_inds,phenotype_name] = multiplicator.used %*% y.used

		###############
		### STAGE 3 ###
		###############

		GWAS_MatrixOperation <- function(i){
		try({
		
			if(num_nodes == 1){
				analysis_allele = 7:dim(Chip_ped)[2]; analysis_genotype = 1:((dim(Chip_ped)[2]-6)/2)
			}else{
				interval_length = as.integer((ncol(Chip_ped)-6)/num_nodes); if(interval_length %% 2 != 0){interval_length = interval_length-1}
				if(i == 1){ analysis_allele = 7:(interval_length+6); analysis_genotype = 1:((interval_length)/2) }
				if(i != 1 && i!= num_nodes){ analysis_allele = (7+(i-1)*interval_length):(i*interval_length+6); analysis_genotype = ((i-1)*interval_length/2+1):((i*interval_length)/2) }
				if(i == num_nodes){ analysis_allele = (7+(i-1)*interval_length):ncol(Chip_ped); analysis_genotype = ((i-1)*interval_length/2+1):((ncol(Chip_ped)-6)/2) }
			}
			
			ped_allele = Chip_ped[keep_inds,analysis_allele]
			dfFull_anova = length(keep_inds) - 3; dfFull_lm = length(keep_inds) - 2
			
			allele_to_genotype_012 <- function(i){
				x_sum = rowSums(ped_allele[,(2*i-1):(2*i)])
				x_result = x_sum - 2; x_result[x_result < 0] = 0 
				return(x_result)
			}
			
			ped_geno_lma = do.call('cbind',mclapply(1:((dim(ped_allele)[2])/2),allele_to_genotype_012,mc.cores=1)) 
			ped_geno_lmd = ped_geno_lma; ped_geno_lmd[ped_geno_lmd == 2] = 0
			ped_geno_anova1 = as.matrix(data.frame(ped_geno_lma==1))
			ped_geno_anova2 = as.matrix(data.frame(ped_geno_lma==2))
			
			yt = multiplicator.used %*% y.used
			covt = multiplicator.used %*% rep(1,length(y.used))
			genot_anova1 = multiplicator.used %*% ped_geno_anova1
			genot_anova2 = multiplicator.used %*% ped_geno_anova2
			
			rm(ped_geno_lma, ped_geno_lmd, ped_geno_anova1, ped_geno_anova2); gc()
			
			cov_qr = t( qr.Q(qr(covt)) )
			
			phe_ortho = t(yt) - tcrossprod(t(yt),cov_qr) %*% cov_qr; div_phe = sqrt( rowSums(phe_ortho^2) ); phe_ortho = phe_ortho/div_phe
			
			geno_ortho_anova1 = t(genot_anova1) - tcrossprod(t(genot_anova1),cov_qr) %*% cov_qr; 
			div_anova1 = sqrt(rowSums(geno_ortho_anova1^2)); 
			geno_ortho_anova1 = geno_ortho_anova1/div_anova1
			
			geno_ortho_anova2 = t(genot_anova2) - tcrossprod(t(genot_anova2),cov_qr) %*% cov_qr;
			geno_ortho_anova2 = geno_ortho_anova2 - rowSums(geno_ortho_anova2*geno_ortho_anova1) * geno_ortho_anova1; 
			div_anova2 = sqrt(rowSums(geno_ortho_anova2^2)); 
			geno_ortho_anova2 = geno_ortho_anova2/div_anova2
			
			rm(genot_anova1, genot_anova2); gc()
					
			r1 = tcrossprod(phe_ortho, geno_ortho_anova1); r2 = tcrossprod(phe_ortho, geno_ortho_anova2)
			r_anova = r1^2 + r2^2; F = (r_anova/(1-r_anova))*(dfFull_anova/2); 
			pvalue_f = pf(F, 2, dfFull_anova, lower.tail=FALSE)

			rm(geno_ortho_anova1, geno_ortho_anova2); gc()
			
			result_MatrixOperation = t(rbind(Chip_map[analysis_genotype,1], Chip_map[analysis_genotype,4], -log10(pvalue_f)))
			return( subset(result_MatrixOperation, result_MatrixOperation[,3]>logP_threshold) )
		})
		}
		
		###############
		### STAGE 4 ###
		###############

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
						
						tmp_ind_geno[,2] = as.character(tmp_ind_geno[,2]); tmp_ind_geno = data.frame(tmp_ind_geno,0,0)
						colnames(tmp_ind_geno) = c("Sort","Genotype","Add_Code","Dom_Code")
						tmp_ind_geno[tmp_ind_geno[,2] %in% c("AB","BA"),c("Add_Code","Dom_Code")] = 1 # Recode "AB"/"BA" with "11" into "1/1" #
						tmp_ind_geno[tmp_ind_geno[,2] == "BB",c("Add_Code")] = 2 # Recode "BB" into "2/0" #

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
						rmna_genot = rmna_multiplicator.used %*% as.matrix(rmna_ind_geno[,3:4])
									
						fit_null = lm(rmna_yt ~ -1 + rmna_covt)
						fit_AddCode = lm(rmna_yt ~ -1 + rmna_covt + rmna_genot[,1])
						fit_DomCode = lm(rmna_yt ~ -1 + rmna_covt + rmna_genot[,2])
						fit_AddDomCode = lm(rmna_yt ~ -1 + rmna_covt + rmna_genot[,1:2])

						fit_AddDomCode_coef = summary(fit_AddDomCode)$coefficients
						betaAdd = fit_AddDomCode_coef[2,"Estimate"]; betaDom = fit_AddDomCode_coef[3,"Estimate"]
						seAdd = fit_AddDomCode_coef[2,"Std. Error"]; seDom = fit_AddDomCode_coef[3,"Std. Error"]
						tAdd = betaAdd/seAdd; tDom = betaDom/seDom
						
						logP1 = -log10(anova(fit_null,fit_AddCode)[,'Pr(>F)'][2])
						logP2 = -log10(anova(fit_null,fit_DomCode)[,'Pr(>F)'][2])
						logP3_1 = -log10(anova(fit_null,fit_AddDomCode)[,'Pr(>F)'][2])
						logP3_2 = -log10(anova(fit_AddCode,fit_AddDomCode)[,'Pr(>F)'][2])
						return(c(Chip_map[i,1],Chip_map[i,4],logP1,logP2,logP3_1,logP3_2,tAdd,tDom,betaAdd,betaDom,seAdd,seDom))
						
					}else{
						return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)) # Loci with genotype type = 3, but one of them < GT_IndMinimum #
					}
				}else{
					return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)) # Loci with genotype type < 3 #
				}
			}else{
				return(c(Chip_map[i,1],Chip_map[i,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)) # Loci from Phenotype with nona-individuals num < Phe_IndMinimum #
			}
		})
		}
		
		###############
		### STAGE 5 ###
		###############

		if(Run_separated == F){			
		
			setwd(rawdir); Chip_map_all = Chip_map = read("Data_clean.bim")
			if(!file.exists("Data_clean_recode12.ped")){system("plink --noweb --silent --bfile Data_clean --recode12 --out Data_clean_recode12")}
			
			print(paste0("Whole-Chrs Loading STA: ",date()))
			Chip_ped = read.big.matrix("Data_clean_recode12.ped", sep = " ", type = "integer")
			Chip_fam <- read(paste(rawdir,"/Data_clean.fam",sep=""))
			options(bigmemory.allow.dimnames=TRUE); rownames(Chip_ped) = Chip_fam[,2]
			print(paste0("Whole-Chrs Loading END: ",date()))
			
			print(paste0("Whole-Chrs Computing STA: ",date()))
			if(matrix_acceleration == T){ 
				Sig_SNPs = do.call('rbind',mclapply(1:num_nodes,GWAS_MatrixOperation,mc.cores=num_nodes))
				logP_Chip = as.data.frame(matrix(NA,nrow(Chip_map),12)); logP_Chip[,1:2] = Chip_map[,c(1,4)]; rownames(logP_Chip) = paste0("SNP",1:nrow(Chip_map))
				Sig_SNPs_order = rownames(logP_Chip)[paste(logP_Chip[,1],logP_Chip[,2],sep="_") %in% paste(Sig_SNPs[,1],Sig_SNPs[,2],sep="_")]
				logP_Chip_Sig = do.call('rbind',mclapply(as.numeric(gsub("SNP","",Sig_SNPs_order)),GWAS_Chip,mc.cores=num_nodes))
				logP_Chip[Sig_SNPs_order,] = logP_Chip_Sig; rm(Chip_ped); gc()
			}else{
				logP_Chip = do.call('rbind',mclapply(1:((ncol(Chip_ped)-6)/2),GWAS_Chip,mc.cores=num_nodes)); rm(Chip_ped); gc()
			}
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
				Chip_fam <- read(paste(rawdir,"/Data_clean.fam",sep=""))
				options(bigmemory.allow.dimnames=TRUE); rownames(Chip_ped) = Chip_fam[,2]
				print(paste0("Chr",j," Loading END: ",date()))
				
				print(paste0("Chr",j," Computing STA: ",date()))
				if(matrix_acceleration == T){ 
					Sig_SNPs_tmp = do.call('rbind',mclapply(1:num_nodes,GWAS_MatrixOperation,mc.cores=num_nodes))
					logP_Chip_tmp = as.data.frame(matrix(NA,nrow(Chip_map),12)); logP_Chip_tmp[,1:2] = Chip_map[,c(1,4)]; rownames(logP_Chip_tmp) = paste0("SNP",1:nrow(Chip_map))
					Sig_SNPs_order_tmp = rownames(logP_Chip_tmp)[paste(logP_Chip_tmp[,1],logP_Chip_tmp[,2],sep="_") %in% paste(Sig_SNPs_tmp[,1],Sig_SNPs_tmp[,2],sep="_")]
					logP_Chip_Sig_tmp = do.call('rbind',mclapply(as.numeric(gsub("SNP","",Sig_SNPs_order_tmp)),GWAS_Chip,mc.cores=num_nodes))
					logP_Chip_tmp[Sig_SNPs_order_tmp,] = logP_Chip_Sig_tmp; rm(Chip_ped); gc()
				}else{
					logP_Chip_tmp = do.call('rbind',mclapply(1:((ncol(Chip_ped)-6)/2),GWAS_Chip,mc.cores=num_nodes)); rm(Chip_ped); gc()
				}
				print(paste0("Chr",j," Computing END: ",date()))
				
				logP_Chip = rbind(logP_Chip,logP_Chip_tmp); rm(logP_Chip_tmp); gc()
				print(paste0("Chr",j," Finished!"))
			}
			
		}
		
		setwd(paste(outdir,"/2_Pvalue",sep=""))
		colnames(logP_Chip)=c("chr","pos","-logP(AddCode)","-logP(DomCode)","-logP(AddDomCode_Null)","-logP(AddDomCode_Add)","tAdd","tDom","betaAdd","betaDom","seAdd","seDom")
		file_name = paste('Pvalue_',phenotype_name,'.txt',sep='')
		write.table(logP_Chip,file=file_name,row.names=F,col.names=T,quote=F)
		system(paste("tar zcvf ",file_name,".tar.gz"," ",file_name,sep=""))
		system(paste("rm ",file_name,sep="")); rm(logP_Chip); gc()
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
