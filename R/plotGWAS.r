
##' @title A Flexible Function to draw Manhattan Plot 
##' 
##' @description Manhattan Plot for whole genome association analyses for single or multiple traits
##'
##' @param chrs A vector of numbers or characters indicating the chromosomes of markers.
##' @param traitIdx A numeric variable indicating the index of trait to plot.
##' @param plotcolor A vector of characters specify the colors used to plot the signatures on each chromosome.
##' @param alldat A dataframe with columns containing SNP, chr, pos and association strength (-log10 P value) of markers for one or more traits.
##' @param y_limit A numeric variable detail the range of y axis.
##' @param main A character specifies the title of the plot. 
##' @param cex_points A numeric value specifies the cex of points in the plot. 
##' @param cex_lab A numeric value specifies the cex of labels in the plot.
##' @param cex_axis A numeric value specifies the cex of axis in the plot.
##'
##' @return a figure
##'
##' @author Bin Yang and Leilei Cui

plotGWAS = function(chrs = chrs,
                    traitIdx = 1,
                    plotcolor = c("darkgreen"),
                    alldat = alldat,
                    y_limit = "",
                    main = "",
                    cex_points = 1, cex_lab = 3.1 ,cex_axis = 3){
					
    SNP = alldat[,"SNP"]
    CHR = alldat[,"chr"]
    chrs = 1:chrs
    CHR[CHR == "X"] = chrs[length(chrs)]; CHR[CHR == 23] = chrs[length(chrs)]
    CHR = as.numeric(as.character(CHR))
    POS = as.numeric(alldat[,"pos"])
    
    if(max(POS,na.rm=T) > 1000){
        POS = POS/1e6 
    } 
	
    chrLen <- tapply(POS,CHR,max) 
	chrmk <- numeric() 
    for(i in 1:length(chrLen)){
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
    
	color <- c("darkred","darkgreen","cyan","red", "darkblue", "brown",
               "black","orange","darkred","darkgreen","cyan",
			   "red", "darkblue", "brown","black","orange","darkred",
			   "darkgreen","cyan","red","blue","brown","black","orange")

	# color <- rep(c("#28292B","#C68D34","#31C79C","#0161AD"),6) # Colors from the Qian Li Jiang Shan Plot #

    if(length(plotcolor) == 2){      
        color2 <- rep(plotcolor,length(chrs))
    } else {
        color2 <- color
    }
	
    idx = which(alldat[,"chr"] %in% chrs)
	
	if(y_limit == "NA"){
		y_limit = range(0, 1.2*max(alldat[idx,c(traits2study)],na.rm = TRUE))
	}else{
		y_limit = range(0, 1.2*y_limit)
	}

    if(length(chrs) > 1 & length(traitIdx) > 1){
        positions = mkrPos
        plot(positions,alldat[idx,traits2study[1]],lwd = 2, col = color[1],type = "p",pch = 19,
        ylab = "-log10(P value)",xlab = "Genome positions",xlim = c(0,max(positions)),
        main = "",xaxt = "n", cex.lab = 1.5, cex.axis = 1.2,
        ylim = range(2, 1.2*max(alldat[idx,c(traits2study)],na.rm = TRUE)))
        for(phe in 2:length(traits2study)){
            points(positions, alldat[idx,traits2study[phe]], lwd = 2,col = color[phe],pch = 19)
        }
        legend("topleft", legend = traits2study, col=color[1:length(traits2study)], 
               border = FALSE, bty = "n",pch = 19, cex = 1.4)
        abline(v=chrmk,lty=4,lwd=0.0005)
        abline(h=2,lty=3,lwd=0.00005)
        axis(side=1,at=(chrmk-(chrLen/2)) ,labels=c(1:18,"X"))
    }
    
    if(length(chrs) > 1 & length(traitIdx) == 1){ 
        positions = mkrPos
        chridx = which(CHR %in% chrs[1])
        plot(positions[chridx], alldat[chridx, traits2study],lwd = 1, col = color2[1], 
             type = "p", pch = 19,ylab = expression(-log[10](italic(P)~value)),xlab ="Chromosomes",
             xlim = c(0,max(positions,na.rm=T)),main = NULL, xaxt = "n",cex.lab = cex_lab,cex.axis = cex_axis,
             ylim = y_limit,xaxs = "i", yaxs = "i", cex = cex_points)             
        for(chrom in 2:length(chrs)){
             chridx = which(CHR %in% chrs[chrom])
             points(positions[chridx],alldat[chridx, traits2study], lwd = 1, col = color2[chrom], pch = 19, cex = cex_points)                
        }
        axis(side=1,at=(chrmk-(chrLen/2)),labels=unique(CHR),cex.axis = cex_axis)     
    }

	return(mkrPos)

}

