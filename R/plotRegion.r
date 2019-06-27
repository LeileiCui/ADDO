
##' @title A Flexible Function for Manhattan Plot 
##' 
##' @description Plot the GWAS results for a single trait on one particular chromosome, or specific region on one chromosome
##'
##' @param chrs A vector of numbers or characters indicating the chromosomes of markers.
##' @param traitIdx  A numeric variable indicating the index of trait to plot.
##' @param alldat A dataframe with columns containing SNP, chr, pos and association strength (-log10 P value) of markers for one or more traits.
##' @param from A character specify the name of SNP.
##' @param to A character specify the name of SNP.
##' @param main A character specifies the title of the plot.
##' @param ldinfo A data frame generated using --r2 --ld-snp command in PLINK .
##' @param ylim A vector of two numeric variables specify the range of y axis.. 
##' @param cex_points A numeric variable specify the cex of points to plot. 
##' @param cex_points_peak A numeric variable specify the cex of point for the lmost significant marker to plot..
##' @param cex_lab A numeric value specifies the cex of labels in the plot.
##' @param cex_axis A numeric value specifies the cex of axis in the plot.
##' @param cex_main AA numeric value specifies the cex of title in the plot.
##'
##' @return a figure
##'
##' @author Bin Yang and Leilei Cui

plotRegion = function(chrs = chrs,
                      traitIdx = 1,
                      alldat = alldat,
                      from = NULL, 
                      to= NULL,
                      main="",
                      ldinfo=ldinfo,
                      ylim=NULL,
                      cex_points = 2.5, cex_points_peak = 4, cex_lab = 1.3 ,cex_axis = 1.3, cex_main = 1.4){
					

    SNP = alldat[,"SNP"]
    CHR = as.numeric(alldat[,"chr"])
    POS = as.numeric(alldat[,"pos"])
    CHR[CHR == "X"] = chrs[length(chrs)]; CHR[CHR == 23] = chrs[length(chrs)]
    if(max(POS) > 1000){
        POS = POS/1e6
    } 
    chrLen <- tapply(POS,CHR,max)
    chrmk <- numeric()
    for(i in 1:(chrs[length(chrs)]-1)){
        chrmk[i] <- sum(as.numeric(chrLen[1:i]))
    }
    mkrPos = numeric(length(POS))

    for(i in as.numeric(unique(CHR))){
        if(i == 1){
            mkrPos[CHR == i] = POS[CHR == i]
        } else {
            mkrPos[CHR == i] = POS[CHR == i] + chrmk[i-1]
        }
    }

allTraits = colnames(alldat)[!(colnames(alldat) %in% c("SNP","chr","pos"))]
traits2study = allTraits[traitIdx]

    idx = which(alldat[,"chr"] %in% chrs)
    #par(xaxs = "i", yaxs = "i") 

    if(length(chrs) == 1){
        if(!is.null(from)){
            index.pos = which(SNP == from):which(SNP == to)
        } else {
        index.pos = CHR %in% chrs
        }
    idx = index.pos
    positions = POS[index.pos]
    cat("\n","Start to plot the results on a single chromosome...")         

        if(length(traitIdx) == 1){
            cat("\n","Pattern: Plot single trait on single chr")
            pvals2plot = alldat[index.pos,traits2study]
            peakIdx = which(pvals2plot == max(pvals2plot,na.rm = T))[1]
            plotSNP = as.character(SNP[index.pos])
            r2 = numeric(length(plotSNP))
            names(r2) = plotSNP
            r2 = ldinfo[names(r2),2]
            color = character(length(plotSNP))
            color[r2 >= 0.8] = "orange"
            color[peakIdx] = "red"
			color[r2 < 0.8 & r2 >= 0.5] = "yellow"
            color[r2 < 0.5 & r2 >= 0.2] = "green"
            color[r2 < 0.2 & r2 >= 0] = "blue"
			color[r2 < 0] = "black"
            cex = numeric(length(plotSNP))
            cex = rep(cex_points,length(cex))*0.6
            cex[peakIdx] = cex_points_peak*0.6
            cat("\n",positions[peakIdx])
            #par(xpd=T, mar=par()$mar+c(0,0,0,0))
			if(is.null(ylim)){
			    ylim=range(0, 1.4*max(alldat[idx,c(traits2study)],na.rm = TRUE))
			} 
            plot(positions,pvals2plot,lwd = 2, col = color,type = "p", pch = 19,cex = cex,
                 ylab = expression(-log[10](italic(P)~value)),xlab = paste(paste("Chromosome",chrs,sep = ""), "(Mb)",sep = " "),
                 main = main,ylim = ylim, cex.lab=cex_lab,cex.axis=cex_lab,cex.main=cex_main)                
            #text(positions[peakIdx][1], pvals2plot[peakIdx][1],cex=1.3,"*",adj = c(-0.1,-0.6))
            legend("topleft",bty="n",cex=1.5,
                   legend=c("Top significance", expression(r^2~">="~"0.8"),expression("0.5"~"<="~r^2~"<"~"0.8"),
                            expression("0.2"~"<="~r^2~"<"~"0.5"), expression(r^2~"<"~"0.2")),
                   col = c("red","orange","yellow","green","blue"),pch=19)    				   
        }
 }
}

