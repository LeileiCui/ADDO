
##' @title Quality Control of Phenotype and Genotype 
##'        (STEP1 of The Add-Dom Model)
##' 
##' @description Quality control of Phenotype and Genotype (Input file format: PLINK or GenABEL)
##'              (1) Discard phenotypes with <200 individuals or logical variables;
##'              (2) Remove extreme values over threefold sd from the mean;
##'              (3) Remove genotypes with MAF<0.05 or missing rate>0.1;
##'              (4) Normalized phenotypes using "quantile" or "-log2" transforming;
##'              (5) Histogram Plot of the raw, clean, residual, normalized and transformed phenotypes;
##'              (6) Calculate kinship matrix using GenABEL, EMMA, EMMAX, GEMMA, GCTA, HOMEBREW_AFW or HOMEBREW_AS;
##'              (7) Summary the mean, sd and sum of each phenotype.
##'
##' @details NOTE1: PLINK Input Format (1) Genotype File, named "file.bed", "file.bim" & "file.fam" (2) Phenotype File, named "file.phe" (1st column name should be "id"; The covariates columns should be prior than phenotypes; The sex column should coded as female=0 and male=1) (3) Covariates File, named "file.covs" (1st column is phenotype names; 2nd column is corresponding covariates separated by ","). NOTE2: GenABEL Input Format (1) file.ABEL.dat (Just contain one GenABEL type variable named "dat") (2) file.covs (1st column is phenotype name; 2nd column is corresponding covariates and all covariates should be separated by ",") NOTE3: Required Softwares: plink (v.1.90) & gcta64 (or emma/emmax-kin/gemma, only required when specified)
##'
##' @param indir A character. The input directory where contains the input bPLINK or GenABEL data.
##' @param outdir A character. The output directory where generates the folder: "1_PheGen".
##' @param Input_name A character. The prefixes of the input files.
##' @param Input_type A character. The format of input data. Please select from "PLINK" or "GenABEL".
##' @param Kinship_type A character. The method to generate kinship matrix. Please select from "GenABEL","EMMA","EMMAX","GEMMA", "GCTA", "GCTA_ad", "HOMEBREW_AFW" or "HOMEBREW_AS".
##' @param PheList_Choose A logic variable. T: Investigate specified phenotypes; F: Investigate all phenofiles.
##' @param PheList A vector of character. Please specifie a list like c("id","cov1","cov2","phe1","phe2"), when "PheList_Choose=F".
##' @param Phe_ResDone A logic variable. T: The input data has already been residualize, won't correct the covariates effect; F: Correct the covariates effect.
##' @param Phe_NormDone A logic variable. T: The input data has already been normalized, won't implement Log Transforming; F: Implement Log Transforming.
##' @param Normal_method A character. When choose "Phe_NormDone = F", the specified normalized method will be needed, "LOG2" or "QUANTILE".
##' @param covariates_sum A numeric variable. The sum of all covariates.
##' @param covariates_types A vector of character. The type of all covariates. Please select from "n" and "f". "f" stands for factorization.
##' @param Phe_IndMinimum A numeric variable. Remove phenotypes without enough available individuals.
##' @param Phe_Extreme A numeric variable. Phenotype QC2: Remove extreme phenotype values over -/+ Phe_Extreme*sd from mean.
##' @param GT_maf A numeric variable. Genotype QC1: Remove genotypes with MAF<GT_maf.
##' @param GT_missing A numeric variable. Genotype QC2: Remove genotypes with rate>GT_missing.
##' @param num_nodes A numeric variable. The number of cores used parallelly.
##'
##' @return a folder named "1_PheGen" with phenotypes and genotypes after QC.
##'
##' @examples 
##' covariates_types = c("n","f")
##' names(covariates_types) = c("sex","batch")
##' ADDO_AddDom1_QC(indir=indir, outdir=outdir, Input_name="TEST", Kinship_type="GCTA_ad", PheList=c("sex","batch"), covariates_sum=2)
##'
##' @author Leilei Cui and Bin Yang

ADDO_AddDom1_QC <- function(indir = indir,
                            outdir = outdir,
                            Input_name = Input_name,
                            Input_type = "PLINK",
                            Kinship_type = Kinship_type,
                            PheList_Choose = F,
                            PheList = PheList,
                            Phe_ResDone = F,
                            Phe_NormDone = F,
                            Normal_method = "QUANTILE",
                            covariates_sum = covariates_sum,
                            covariates_types = covariates_types,
                            Phe_IndMinimum = 200,
                            Phe_Extreme = 5,
                            GT_maf = 0.05,
                            GT_missing = 0.1,
                            num_nodes = 10){

setwd(outdir); system("mkdir 1_PheGen 1_PheGen/0_Data 1_PheGen/1_PheRaw 1_PheGen/2_PheClean 1_PheGen/3_PheRes 1_PheGen/4_PheNor 1_PheGen/5_PheTrans")

read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))
read.ped <- function(...) as.data.frame(fread(header=F,colClasses="character",...))

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

covs_matrix <- function(data,new_covs,covariates_types){
try({
	x_mat=matrix(1,nrow=dim(data)[1],ncol=1)
	for (c in 1:length(new_covs)) {
		cov=new_covs[c]
		type_cov=covariates_types[cov]
		if (type_cov=='f') {
			if(length(unique(data[,cov])) == 1){
				matrice = rep(1,nrow(data))
			}else{
				old=data[,cov]
				old=as.factor(as.numeric(old))
				levels=levels(old)
				n=length(levels)-1
				matrice=matrix(ncol=n,nrow=length(old),NA)
				colnames(matrice)=levels[1:n]
				for (l in levels[1:n]) {
					matrice[,l]=(old==l)
				}
			}
			x_mat=cbind(x_mat,matrice)
		} else {
			data[,cov]=as.numeric(data[,cov])
			x_mat=cbind(x_mat,data[,cov])
		}
	}
	return(x_mat)
})
}
#batch = model.matrix(~as.factor(data.tmp[,"batch"]))

log2_norm <- function(x){
try({
	x_log2 = x
	for(i in 1:length(x)){
		if(!(is.na(x[i])) && x[i]>0){
			x_log2[i] = log2(x[i])
		}else{
			x_log2[i] = NA
		}
	}
	return(x_log2)
})
}

quantile_norm <- function(x){
try({
	names(x) = paste("Ind",1:length(x),sep="")
	x_clean = x[!is.na(x)]
	x_clean_norm = qnorm(rank(x_clean)/(1+length(x_clean)))
	x[names(x_clean)] = x_clean_norm
	return(x)
})
}


if(Input_type == "PLINK"){
	setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))
	if(!file.exists("Data_clean.bed")){system(paste("plink --noweb --silent --bfile ",indir,"/",Input_name," --maf ",GT_maf," --geno ",GT_missing," --make-bed --out Data_clean",sep=""))}
	if(!file.exists("Data_clean.frq")){system("plink --noweb --silent --bfile Data_clean --freq --out Data_clean")}
	Clean_map = read("Data_clean.bim")
	
	Phe_Raw = read.table(paste(indir,"/",Input_name,".phe",sep=""),header=T)
}

if(Input_type == "GenABEL"){
	load(paste(indir,"/",Input_name,".ABEL.dat",sep="")); dat = dat # Note0: The raw variable name of GenABEL #
	
	setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))
	export.plink(dat,filebasename="Data_raw",phenotypes="all")
	#system("plink --noweb --silent --file Data_raw --transpose --recode --out Data_raw") # Just turn on this command when run on sparse/dense #
	if(!file.exists("Data_clean.bed")){system(paste("plink --noweb --silent --tfile Data_raw --maf ",GT_maf," --geno ",GT_missing," --make-bed --out Data_clean",sep=""))}
	if(!file.exists("Data_clean.frq")){	system("plink --noweb --silent --bfile Data_clean --freq --out Data_clean")}
	Clean_map = read("Data_clean.bim")
	
	Phe_Raw = read.table("Data_raw.phe",header=T)[,-1] # NOTE: When the input files are in GenABEL format, the 1st col of exported phenotype file is "FID" #
	rownames(Phe_Raw) = Phe_Raw[,1]; colnames(Phe_Raw)[1] = "id" # NOTE: When the input files are in GenABEL format, the colname of individual's column is "IID" #
}


if(PheList_Choose == T){ PheList = PheList } else{ PheList = colnames(Phe_Raw) }
	
	PheRaw_judge = rep("clean_phe",ncol(Phe_Raw))
	for(i in 1:ncol(Phe_Raw)){
		if(length(Phe_Raw[,i][unique(Phe_Raw[,i]) %in% c("T","TRUE","F","FALSE")]) > 0 || length(unique(Phe_Raw[,i][!is.na(Phe_Raw[,i])])) < 3){
			PheRaw_judge[i] = "bad_phe"
		}else{
			PheRaw_judge[i] = "clean_phe"
		} 
	}
	names(PheRaw_judge) = colnames(Phe_Raw)
	PheRaw_judge[1:(covariates_sum+1)] = "clean_phe" # In case some covariates are logical or case_control #
	PheList = PheList[!PheList %in% names(PheRaw_judge[PheRaw_judge=="bad_phe"])] 
	Phe_Raw = Phe_Raw[,PheList]
	write.table(Phe_Raw,file="../1_PheRaw/Phe_Raw.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/1_PheRaw",sep=""))
	Plot_Phe_raw = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Raw,Phe_names=PheList,mc.cores=num_nodes)

	Phe_Clean = Phe_Raw[,PheList]; rownames(Phe_Clean) = paste("Ind_",Phe_Clean[,1],sep=""); PheList_rm = c() # In case for the individual names are just number #
	for(i in (covariates_sum+2):length(PheList)){ # Note1: The first 3 cols are "id","sex","batch" #
		Phe_tmp = Phe_Clean[,c("id",PheList[i])]; Phe_tmp_clean = Phe_tmp[!is.na(Phe_tmp[,2]),] # Remove NA inds #
		if(length(unique(Phe_tmp_clean)[,2]) != 2){ # In case for some 0/1 phenotypes #
			Phe_extreme1 = rownames(Phe_tmp_clean[Phe_tmp_clean[,2]>(mean(Phe_tmp_clean[,2],na.rm=T)+sd(Phe_tmp_clean[,2],na.rm=T)*Phe_Extreme),]) 
			Phe_extreme2 = rownames(Phe_tmp_clean[Phe_tmp_clean[,2]<(mean(Phe_tmp_clean[,2],na.rm=T)-sd(Phe_tmp_clean[,2],na.rm=T)*Phe_Extreme),])
			Phe_Clean[c(Phe_extreme1,Phe_extreme2),PheList[i]] = NA # Assign individuals with extreme phenotype to be NA #
		}
		if(length(Phe_Clean[,PheList[i]][!is.na(Phe_Clean[,PheList[i]])]) < Phe_IndMinimum) {PheList_rm = c(PheList_rm,PheList[i])} # Get the phe name with Num_Ind < Phe_IndMinimum #
	}
	rownames(Phe_Clean) = Phe_Clean[,1]

	Phe_Clean = Phe_Clean[,colnames(Phe_Clean)[!colnames(Phe_Clean) %in% PheList_rm]]; PheList = PheList[!PheList %in% PheList_rm]
	write.table(Phe_Clean,file="../0_Data/Data_clean.phe",row.names=T,col.names=T,quote=F)
	write.table(Phe_Clean,file="../2_PheClean/Phe_Clean.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/2_PheClean",sep=""))
	Plot_Phe_clean = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Clean,Phe_names=PheList,mc.cores=num_nodes)

if(Phe_ResDone == T){
	Phe_Res = Phe_Clean
	write.table(Phe_Res,file="../3_PheRes/Phe_Res.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/3_PheRes",sep=""))
	Plot_Phe_res = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Res,Phe_names=PheList,mc.cores=num_nodes)
}else{
	Phe_Res = Phe_Clean; Phe_covs = read.table(paste(indir,"/",Input_name,".covs",sep=""),header=T); rownames(Phe_covs) = Phe_covs[,1]

	for(args in (covariates_sum+2):length(PheList)){		
		phenotype_name = PheList[args]; Phe_covs_tmp = strsplit(Phe_covs[phenotype_name,2],",")[[1]]
			
		keep_inds_tmp = rownames(Phe_Clean)[!is.na(Phe_Clean[,phenotype_name])] # (1) Remove individuals with NA phenotype #
		x_mat.tmp = Phe_Clean[keep_inds_tmp,c("id",Phe_covs_tmp)]
		keep_inds = rownames(x_mat.tmp[complete.cases(x_mat.tmp),]) # (2) Remove individuals with NA covariates #
		y.used = Phe_Clean[keep_inds,phenotype_name]; names(y.used) = keep_inds
		covs.tmp = as.data.frame(Phe_Clean[keep_inds,Phe_covs_tmp]); rownames(covs.tmp) = keep_inds; colnames(covs.tmp) = Phe_covs_tmp
		x_mat.used = covs_matrix(data = covs.tmp, new_covs = Phe_covs_tmp, covariates_types=covariates_types)
		res.used = residuals(lm(y.used ~ x_mat.used[,-1])) # Method1: Using covariates matrix #
		
		# residual_data = cbind(phe=y.used,covs.tmp)
		# for(covs_order in 1:length(Phe_covs_tmp)){ if(covariates_types[Phe_covs_tmp[covs_order]] == "f"){ residual_data[,Phe_covs_tmp[covs_order]] = as.factor(residual_data[,Phe_covs_tmp[covs_order]]) }}
		# res.used = residuals(lm(paste("phe ~ ",paste(Phe_covs_tmp,collapse=" + "),sep=""),data=residual_data)) # Method2: Using factored covariates #
		
		Phe_Res[,phenotype_name] = NA; Phe_Res[names(res.used),phenotype_name] = res.used
	}
	write.table(Phe_Res,file="../3_PheRes/Phe_Res.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/3_PheRes",sep=""))
	Plot_Phe_res = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Res,Phe_names=PheList,mc.cores=num_nodes)	
}

if(Phe_NormDone == T){
	Phe_Nor = Phe_Res
	write.table(Phe_Nor,file="../4_PheNor/Phe_Nor.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/4_PheNor",sep=""))
	Plot_Phe_nor = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Nor,Phe_names=PheList,mc.cores=num_nodes)
}else{	
	# normPvals = sapply(PheList[-c(1:(covariates_sum+1))],function(x){
		# if(length(unique(Phe_Res[,x])[!is.na(unique(Phe_Res[,x]))])>1){ 
			# return(shapiro.test(Phe_Res[,x])$p)
		# }else{
			# return(NA)} # [4-1] Skip phenotypes which just contain one unique value #
	# })
	
	Phe_Nor = Phe_Res
	for(i in (covariates_sum+2):length(PheList)){
	#	if(!is.na(normPvals[i-(covariates_sum+1)])){
	#		if(normPvals[i-(covariates_sum+1)] < 1e-8){ # [4-2] Skip phenotypes with shapiro.test < 1e-8 #
				if(Normal_method == "LOG2"){ Phe_Nor[,PheList[i]] = log2_norm(Phe_Res[,PheList[i]]) }
				if(Normal_method == "QUANTILE"){ Phe_Nor[,PheList[i]] = quantile_norm(Phe_Res[,PheList[i]]) }
	#		}
	#	}
	}
	write.table(Phe_Nor,file="../4_PheNor/Phe_Nor.txt",row.names=T,col.names=T,quote=F)
	setwd(paste(outdir,"/1_PheGen/4_PheNor",sep=""))
	Plot_Phe_nor = mclapply((covariates_sum+2):length(PheList),Plot_Phe_Histogram,Phe_data=Phe_Nor,Phe_names=PheList,mc.cores=num_nodes)
}


try({
	setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))
	MAF_info = read.header("Data_clean.frq")
	Plot_Phe_Histogram(1,as.data.frame(MAF_info[,"MAF"]),paste(Input_name,"_MAF",sep=""))
})


try({
	if(dim(Clean_map)[1] <= 100000){
		setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))
		system("awk '{if($1!=0 && $4!=0) print $2}' Data_clean.bim > LD_tmp.txt")
		system("plink --noweb --silent --bfile Data_clean --extract LD_tmp.txt --recode --out LD_tmp")
		system("plink --noweb --silent --bfile Data_clean --r2 --ld-window 99999 --ld-window-r2 0 --out LD_tmp")
		map_tmp = read("Data_clean.bim")
		dat_tmp = read.header("LD_tmp.ld"); system("rm LD_tmp*")

		distance = abs(dat_tmp[,5] - dat_tmp[,2])
		R2 = dat_tmp[,7]
		Beta.st=c(beta1=0.00001)
		Beta.nonlinear = nls(R2~1/(1+4*beta1*distance),start=Beta.st,control=nls.control(maxiter=500))
		tt=summary(Beta.nonlinear)
		Beta.est = tt$parameters[1]
		Beta.error = tt$parameters[2]

		mkrs = c(0,1e3,2e3,4e3,8e3,12e3,16e3,20e3,40e3,80e3,120e3,160e3,200e3,300e3,460e3,800e3,900e3,1000e3)
		D.est = 1/(1+4*Beta.est*mkrs)

		png(paste("LD_decay_",Input_name,".png",sep=""),width=10000,height=8000,res=600,type="cairo")
		par(mai=c(1.4,1.4,1.6,0.3),mgp=c(4.5,2,0))
		plot(log2(mkrs),D.est,col="black",type="l",ylim=c(0,1),ylab=expression(Predicted~LD~r^2),xaxt="n",xlab="Physical distance(kb)",main=paste("LD Decay: ",Input_name," (",nrow(map_tmp)," SNPs)",sep=""),lwd=3,lty=1,cex.axis=2,cex.lab=2,cex.main=3)
		axis(side=1,at=log2(mkrs),labels=mkrs/1000,cex.axis=2)
		abline(h=0.3,lty=3,col="black",lwd=3)
		dev.off(); rm(map_tmp,dat_tmp)
	}
})

try({
	if(dim(Clean_map)[1] <= 100000){
		setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))
		system("plink --noweb --silent --bfile Data_clean --cluster --matrix --out IBS_tmp")
		mibs_raw = read("IBS_tmp.mibs"); system("rm IBS_tmp*")
		
		plotdat_mibs = mibs_raw[lower.tri(mibs_raw, diag = FALSE)]
		Plot_Phe_Histogram(1,as.data.frame(plotdat_mibs),paste(Input_name,"_IBS_IndPairs",sep=""))
	}
})


setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))

if(Kinship_type == "GenABEL"){
	
	if(!file.exists("Data_clean.tped")){system("plink --noweb --silent --bfile Data_clean --transpose --recode --out Data_clean")}
	convert.snp.tped(tpedfile = "Data_clean.tped", tfamfile = "Data_clean.fam", outfile = "Data_clean.raw")
	dat_clean = load.gwaa.data(phenofile= "Data_clean.phe", genofile= "Data_clean.raw", force =TRUE, makemap = FALSE, sort = TRUE, id = "id") # Note2: The colnames of inds is "id" #
	save(dat_clean,file="Data_clean.ABEL.dat")
	
	system("plink --noweb --silent --bfile Data_clean --indep-pairwise 50 10 0.5 --out tmp_genabel")
	prunedSNP = read("tmp_genabel.prune.in")[,1]
	prunedSNP = as.character(prunedSNP[prunedSNP %in% Clean_map[,2]])
	gkin = ibs(dat_clean[,prunedSNP],weight="freq")
	write.table(gkin,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_GenABEL.txt",sep=""),row.names=F,col.names=F,quote=F)
	system("rm tmp_*")
}

if(Kinship_type == "EMMA"){
	if(!file.exists("tmp_emma.raw")){system("plink --noweb --silent --bfile Data_clean --recodeA --out tmp_emma")}
	emma_snp1 = t(read.header("tmp_emma.raw")); colnames(emma_snp1) = emma_snp1["IID",]
	emma_snp2  = as.matrix(emma_snp1[-c(1:6),]) 
	emma_snp2 = gsub("1","0.5",emma_snp2); emma_snp2 = gsub("2","1",emma_snp2) 
	emma_snp3 = data.matrix(as.data.frame(emma_snp2), rownames.force = NA) 
	emma_kinship_A = gkin = emma.kinship(emma_snp3, method="additive", use = "complete.obs") # Note3: emma_snp3 must be (1) matrix (2) numeric elements (3) AA:0;AB:0.5;BB:1 #
	#emma_kinship_D = emma.kinship(emma_snp3, method="dominant", use = "complete.obs")
	#emma_kinship_R = emma.kinship(emma_snp3, method="recessive", use = "complete.obs")
	write.table(emma_kinship_A,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_EMMA.Additive.txt",sep=""),row.names=F,col.names=F,quote=F)	
	system("rm tmp_*")
}

if(Kinship_type == "EMMAX"){
	if(!file.exists("tmp_emmax.tped")){system("plink --noweb --silent --bfile Data_clean --recode12 --output-missing-genotype 0 --transpose --out tmp_emmax")}
	system("emmax-kin -v -h -s -d 10 tmp_emmax") # *.hIBS, IBS matrix  #
	system("mv tmp_emmax.hIBS.kinf Data_clean_kinship_EMMAX.hIBS.txt")
	#system("emmax-kin -v -h -d 10 tmp_emmax") # *.BN, BN matrix #
	#system("mv tmp_emmax.hBN.kinf Data_clean_kinship_EMMAX.hBN.txt")
	gkin = read("Data_clean_kinship_EMMAX.hIBS.txt")
	system("rm tmp_*")
}

if(Kinship_type == "GEMMA"){
	system("sed -i s/-9/1/g Data_clean.fam")
	system("gemma -bfile Data_clean -gk 1 -o Data_clean_kinship_GEMMA") # *.cXX, centered matrix #
	system("mv output/Data_clean_kinship_GEMMA.cXX.txt ./")
	#system("gemma -bfile Data_clean -gk 2 -o Data_clean_kinship_GEMMA") # *.sXX, standardized matrix #	
	#system("mv output/Data_clean_kinship_GEMMA.sXX.txt ./")
	gkin = read("Data_clean_kinship_GEMMA.cXX.txt")
	system("rm -r output"); system("rm tmp_*")
}

if(Kinship_type == "GCTA"){
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm-gz --out tmp_gcta_add",sep=""))	
	tmp_gcta = read.table("tmp_gcta_add.grm.gz",header=F); tmp_gcta_id = read("tmp_gcta_add.grm.id")
	tmp_gcta_matrix = matrix(0,nrow(tmp_gcta_id),nrow(tmp_gcta_id))

	for(gcta_id in 1:nrow(tmp_gcta_matrix)){
		tmp_gcta_matrix[gcta_id,1:gcta_id] = subset(tmp_gcta,tmp_gcta[,1]==gcta_id)[,4]
	}
	tmp_gcta_matrix[upper.tri(tmp_gcta_matrix)] = t(tmp_gcta_matrix)[upper.tri(tmp_gcta_matrix)]
	gkin = tmp_gcta_matrix
	write.table(tmp_gcta_matrix,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_GCTA.Additive.txt",sep=""),row.names=F,col.names=F,quote=F)
	system("rm tmp_*")
}

if(Kinship_type == "GCTA_ad"){
	
	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm-gz --out tmp_gcta_add",sep=""))
	tmp_gcta = read.table("tmp_gcta_add.grm.gz",header=F); tmp_gcta_id = read("tmp_gcta_add.grm.id")
	tmp_gcta_matrix = matrix(0,nrow(tmp_gcta_id),nrow(tmp_gcta_id))
	for(gcta_id in 1:nrow(tmp_gcta_matrix)){
		tmp_gcta_matrix[gcta_id,1:gcta_id] = subset(tmp_gcta,tmp_gcta[,1]==gcta_id)[,4]
	}
	tmp_gcta_matrix[upper.tri(tmp_gcta_matrix)] = t(tmp_gcta_matrix)[upper.tri(tmp_gcta_matrix)]
	gkin = tmp_gcta_matrix
	write.table(tmp_gcta_matrix,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_GCTA.Additive.txt",sep=""),row.names=F,col.names=F,quote=F)

	system(paste("gcta64 --bfile Data_clean --thread-num ",num_nodes," --autosome --make-grm-d-gz --out tmp_gcta_dom",sep=""))
	tmp_gcta = read.table("tmp_gcta_dom.d.grm.gz",header=F); tmp_gcta_id = read("tmp_gcta_dom.d.grm.id")
	tmp_gcta_matrix = matrix(0,nrow(tmp_gcta_id),nrow(tmp_gcta_id))
	for(gcta_id in 1:nrow(tmp_gcta_matrix)){
		tmp_gcta_matrix[gcta_id,1:gcta_id] = subset(tmp_gcta,tmp_gcta[,1]==gcta_id)[,4]
	}
	tmp_gcta_matrix[upper.tri(tmp_gcta_matrix)] = t(tmp_gcta_matrix)[upper.tri(tmp_gcta_matrix)]
	gkin = tmp_gcta_matrix
	write.table(tmp_gcta_matrix,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_GCTA.Dominant.txt",sep=""),row.names=F,col.names=F,quote=F)
	system("rm tmp_*")
}

if(Kinship_type == "HOMEBREW_AFW"){
	
	system("plink --noweb --silent --bfile Data_clean --recode --out tmp_matrix")
	Chip_map <- read("Data_clean.bim"); Chip_ped <- read.ped("tmp_matrix.ped")
	system("rm tmp_matrix*")
	
	ped_allele_merge <- function(i){
		x = Chip_ped[,(i*2+5):(i*2+6)]
		x_merge = paste(x[,1],x[,2],sep="")
		return(x_merge)
	}
	Chip_ped_merge = do.call('cbind',mclapply(1:((dim(Chip_ped)[2]-6)/2),ped_allele_merge,mc.cores=num_nodes))
	Chip_ped_merge = cbind(Chip_ped[,1:6],Chip_ped_merge)
	colnames(Chip_ped_merge) = c("FID","IID","PID","MID","sex","phe",as.character(Chip_map[,2]))
	Chip_ped_merge_used = Chip_ped_merge[,7:dim(Chip_ped_merge)[2]]
	
	Chip_freq = read.header("Data_clean.frq") # Note: Here we use the minor allele A1 as the reference allele (as well as frequency) #
	rownames(Chip_freq) = Chip_freq[,2]
	
	X_AFW_Calculate <- function(j){
		
		if(Kinship_model == "Additive"){
			g = NA; p = Chip_freq[colnames(X)[j],"MAF"]
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A1"],sep="")){ g = 2 } # Recode A1A2 as 2 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 1 } # Recode A1A2 as 1 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A2"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 0 } # Recode A2A2 as 0 #
			if(is.na(g)){ X_OneInd_tmp = NA } else { x = (g-2*p)/(sqrt(2*p*(1-p))); X_OneInd_tmp = x }
		}
		
		if(Kinship_model == "Dominant"){
			g = NA; p = Chip_freq[colnames(X)[j],"MAF"]
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A1"],sep="")){ g = (4*p)-2 } # Recode A1A2 as 4p-2 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 2*p } # Recode A1A2 as 2p #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A2"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 0 } # Recode A2A2 as 0 #
			if(is.na(g)){ X_OneInd_tmp = NA } else { x = (g-2*p*p)/(2*p*(1-p)); X_OneInd_tmp = x }
		}
		return(X_OneInd_tmp)
	}
	
	Kinship_model = "Additive"; X = matrix(0,dim(Chip_ped_merge_used)[1],dim(Chip_ped_merge_used)[2])
	rownames(X) = paste(Chip_ped_merge[,1],Chip_ped_merge[,2],sep="_"); colnames(X) = colnames(Chip_ped_merge)[-c(1:6)]

	for(i in 1:(dim(X)[1])){ X[i,] = do.call('cbind',mclapply(1:(dim(Chip_ped_merge_used)[2]), X_AFW_Calculate, mc.cores=num_nodes))[1,] }
	SNP_num = dim(X)[2]; K = (X%*%t(X))/SNP_num
	write.table(K,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_AFW.Additive.txt",sep=""),row.names=F,col.names=F,quote=F)

	Kinship_model = "Dominant"; X = matrix(0,dim(Chip_ped_merge_used)[1],dim(Chip_ped_merge_used)[2])
	rownames(X) = paste(Chip_ped_merge[,1],Chip_ped_merge[,2],sep="_"); colnames(X) = colnames(Chip_ped_merge)[-c(1:6)]
	
	for(i in 1:(dim(X)[1])){ X[i,] = do.call('cbind',mclapply(1:(dim(Chip_ped_merge_used)[2]), X_AFW_Calculate, mc.cores=num_nodes))[1,] }
	SNP_num = dim(X)[2]; K = (X%*%t(X))/SNP_num
	write.table(K,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_AFW.Dominant.txt",sep=""),row.names=F,col.names=F,quote=F)	
	
}

if(Kinship_type == "HOMEBREW_AS"){

	system("plink --noweb --silent --bfile Data_clean --recode --out tmp_matrix")
	Chip_map <- read("Data_clean.bim"); Chip_ped <- read.ped("tmp_matrix.ped")
	system("rm tmp_matrix*")
	
	ped_allele_merge <- function(i){
		x = Chip_ped[,(i*2+5):(i*2+6)]
		x_merge = paste(x[,1],x[,2],sep="")
		return(x_merge)
	}
	Chip_ped_merge = do.call('cbind',mclapply(1:((dim(Chip_ped)[2]-6)/2),ped_allele_merge,mc.cores=num_nodes))
	Chip_ped_merge = cbind(Chip_ped[,1:6],Chip_ped_merge)
	colnames(Chip_ped_merge) = c("FID","IID","PID","MID","sex","phe",as.character(Chip_map[,2]))
	Chip_ped_merge_used = Chip_ped_merge[,7:dim(Chip_ped_merge)[2]]
	
	Chip_freq = read.header("Data_clean.frq") # Note: Here we use the minor allele A1 as the reference allele, as well as frequency #
	rownames(Chip_freq) = Chip_freq[,2]
	
	X_AS_Calculate <- function(j){
		
		if(Kinship_model == "Additive"){
			g = NA; p = Chip_freq[colnames(X)[j],"MAF"]
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A1"],sep="")){ g = 2 } # Recode A1A2 as 2 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 1 } # Recode A1A2 as 1 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A2"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 0 } # Recode A2A2 as 0 #
			X_OneInd_tmp = g
		}
		
		if(Kinship_model == "Dominant"){
			g = NA; p = Chip_freq[colnames(X)[j],"MAF"]
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A1"],sep="")){ g = (4*p)-2 } # Recode A1A2 as 4p-2 #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A1"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 2*p } # Recode A1A2 as 2p #
			if(Chip_ped_merge_used[i,j] == paste(Chip_freq[colnames(X)[j],"A2"],Chip_freq[colnames(X)[j],"A2"],sep="")){ g = 0 } # Recode A2A2 as 0 #
			X_OneInd_tmp = g
		}
		return(X_OneInd_tmp)
	}

	Kinship_model = "Additive"; X = matrix(0,dim(Chip_ped_merge_used)[1],dim(Chip_ped_merge_used)[2])
	rownames(X) = paste(Chip_ped_merge[,1],Chip_ped_merge[,2],sep="_"); colnames(X) = colnames(Chip_ped_merge)[-c(1:6)]
	
	for(i in 1:(dim(X)[1])){ X[i,] = scale( do.call('cbind',mclapply(1:(dim(Chip_ped_merge_used)[2]), X_AS_Calculate, mc.cores=num_nodes))[1,] )[,1] }
	SNP_num = dim(X)[2]; K = (X%*%t(X))/SNP_num
	write.table(K,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_AS.Additive.txt",sep=""),row.names=F,col.names=F,quote=F)

	Kinship_model = "Dominant"; X = matrix(0,dim(Chip_ped_merge_used)[1],dim(Chip_ped_merge_used)[2])
	rownames(X) = paste(Chip_ped_merge[,1],Chip_ped_merge[,2],sep="_"); colnames(X) = colnames(Chip_ped_merge)[-c(1:6)]
	
	for(i in 1:(dim(X)[1])){ X[i,] = scale( do.call('cbind',mclapply(1:(dim(Chip_ped_merge_used)[2]), X_AS_Calculate, mc.cores=num_nodes))[1,] )[,1] }
	SNP_num = dim(X)[2]; K = (X%*%t(X))/SNP_num
	write.table(K,file=paste(outdir,"/1_PheGen/0_Data/Data_clean_kinship_AS.Dominant.txt",sep=""),row.names=F,col.names=F,quote=F)
	
}


setwd(paste(outdir,"/1_PheGen/0_Data",sep=""))

trait2analy <- colnames(Phe_Clean)[-c(1:(covariates_sum+1))] # phe_name #
dat2analy <- as.matrix(Phe_Clean[,trait2analy]) # phe_data #
traitSummary <- matrix(0,length(trait2analy),4) # output_results #
colnames(traitSummary) <- c("Mean","SD","N","h2")
rownames(traitSummary) <- trait2analy

traitSummary[,"Mean"] <- apply(dat2analy,2,mean,na.rm=T)
traitSummary[,"SD"]   <- apply(dat2analy,2,sd, na.rm=T)
traitSummary[,"N"]    <- apply(dat2analy,2,function(x){length(which(!is.na(x)))})

setwd(paste(outdir,"/1_PheGen/0_Data",sep="")); try({system("rm *.log *.hh")})
write.table(traitSummary,file=paste(outdir,"/1_PheGen/Data_clean_TraitSum.txt",sep=""),row.names=T,col.names=T,quote=F)

}
