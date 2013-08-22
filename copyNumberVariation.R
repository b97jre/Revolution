# TODO: Add comment
# 
# Author: johanreimegard
###############################################################################


#plot CtrlFREEC results
main <- function(){
  CtrlFreecfiles <- c("Cr1GR1-2-KS3_BWA_primary.bam_ratio.txt","Inter3-1_DNA_BWA_primary.bam_ratio.txt", "Inter4-1_DNA_BWA_primary.bam_ratio.txt", "Inter5-1_DNA_BWA_primary.bam_ratio.txt", "Intra6-3_DNA_BWA_primary.bam_ratio.txt", "Intra7-2_DNA_BWA_primary.bam_ratio.txt", "Intra8-2_DNA_BWA_primary.bam_ratio.txt", "Cr1GR1-2-KS3_Bowtie_primary.bam_ratio.txt","Inter3-1_DNA_Bowtie_primary.bam_ratio.txt", "Inter4-1_DNA_Bowtie_primary.bam_ratio.txt", "Inter5-1_DNA_Bowtie_primary.bam_ratio.txt", "Intra6-3_DNA_Bowtie_primary.bam_ratio.txt", "Intra7-2_DNA_Bowtie_primary.bam_ratio.txt", "Intra8-2_DNA_Bowtie_primary.bam_ratio.txt")
  SNPfiles <-  c("simple.1to8.Cr1GR1-2-KS3.info.Radapted" , "simple.1to8.Cr_39_1.info.Radapted" , "simple.1to8.Inter3-1.info.Radapted" , "simple.1to8.Inter4-1.info.Radapted" , "simple.1to8.Inter5-1.info.Radapted" , "simple.1to8.Intra6-3.info.Radapted" , "simple.1to8.Intra7-2.info.Radapted" , "simple.1to8.Intra8-2.info.Radapted")
  centromericFile <- c("crub_centromere_cutoffs.txt")
  
  for (j in 1:length(SNPfiles)){
    print(paste("Analysis with "," START" ,sep = SNPfiles[j]))
    plotCNVsandSNPdistribution(CtrlFreecfiles[j],SNPfiles[j],centromericFile)
    print(paste("Analysis with "," FINISHED" ,sep = SNPfiles[j]))
  }
  
  
}
plotCNVsandSNPdistribution <-function(CtrlFreecfiles,SNPfile,centromericFile){
  
  #concatenate all files
  plotfile <- data.frame()
  for(i in 1:length(CtrlFreecfiles)){
    ratio <- read.table(paste(CtrlFreecfiles[i]), header=T)
    plotfile <- rbind(plotfile,ratio)
  }
  
  SNPs <- read.table(SNPfile,sep="\t",header=TRUE)
  SNPs.common <- subset(SNPs, Call=="0/0"|Call=="0/1"|Call=="1/1")
  
  centromeric <- read.csv(centromericFile,header=T)
  #one of the entries in this file must be a mistake (one order of magnitude too large)
  centromeric[7,3] <- 17322208
  #uncomment next line if you want to output individual plots for each file
  #pdf(file=paste(files[i],".pdf",sep=""))
  pdfFileName = paste("SNPs_and_CNVs.",".pdf",sep=SNPfile)
  pdf(pdfFileName)
  ploidy=2
  
  
  par(mfrow=c(5,2), bty="l", cex=0.6)
  
  for (i in (1:8)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    temp = subset(SNPs.common,chr==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$Ratio[tt]*ploidy,ylim = c(0,3*ploidy) ,pch = 20,col = colors()[88],xlab = paste ("position, chr",i),ylab="")
      mtext("normalized copy number profile",side=2,line=2,col=1,cex=0.5)
      mtext("SNP sites density",side=4,line=2,col=4,cex=0.5)
      #comment next line for individual plots
      #points(plotfile$Start[tt],plotfile$Ratio[tt]*ploidy,pch = 20,col = colors()[88])
      tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep="")  & plotfile$CopyNumber>ploidy )
      points(plotfile$Start[tt],plotfile$Ratio[tt]*ploidy,pch = 20,col = colors()[136])
      tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep="")  & plotfile$CopyNumber<ploidy )
      points(plotfile$Start[tt],plotfile$Ratio[tt]*ploidy,pch = 20,col = colors()[461])
      tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
      #points(plotfile$Start[tt],plotfile$CopyNumber[tt], pch =20, col = colors()[24])
      cen <- centromeric[centromeric$Chromosome==i,]
      abline(v=cen$Centromere.start,col="blue")
      abline(v=cen$Centromere.end, col="blue")
      # add new plot in same window
      par(new=T)
      #getting distribution of heterozygous calls 
      HeterozygousDensity <- density(temp$position[temp$Call=="0/1"],bw=50000)
      #getting distribution of minor homozygous calls calls 
      HomozygousMinorDensity <- density(temp$position[temp$Call=="1/1"],bw=50000)
      
      #getting distribution of all SNPs 
      AllDensity <-density(temp$position,bw=50000)
      
      ##plotting reads 
      ymax = max(c(HeterozygousDensity$y,HomozygousMinorDensity$y,AllDensity$y))
      plot(HeterozygousDensity,axes=FALSE,,xlab="",ylab="",ylim=c(0,ymax),pch=2,col=4,main="")
      axis(4,col=4,col.lab=4,col.axis = 4)
      lines(HomozygousMinorDensity,pch=2,col=2)
      lines(AllDensity,pch=4)
    }
  }
  par(mar=c(0,0,0,0))
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',type="p")  
  plot_colors <- c(colors()[88], colors()[136],colors()[461])
  
  legend(x = "top",inset = 0,
         legend = c("2 alleles", "1 allele", ">2 alleles"),pch=20, 
         col=plot_colors,  cex=1.5)
  
  plot(1, type = "l", axes=FALSE, xlab="", ylab="")
  plot_colors <- c(4,2, 1)
  legend(x = "top",inset = 0,
         legend = c("Heterozygous", "Minor Homozygous", "All"), 
         col=plot_colors, lwd=2, cex=1.5)
  dev.off()
}

boxScatterPlot <- function(x,y,ylabText="Major readCount",xlabText="Minor readCount",title ="Enhanced Scatterplot"){
  par(fig=c(0,0.8,0,0.8), new=TRUE)
  plot(x, y , xlab=xlabText,
       ylab=ylabText)
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  boxplot(x, horizontal=TRUE, axes=FALSE)
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  boxplot(y, axes=FALSE)
  mtext(title, side=3, outer=TRUE, line=-3)
  
}

boxScatterPlot2 <- function(x,y,ylabText="Major readCount",xlabText="Minor readCount",title ="Enhanced Scatterplot"){
  par(fig=c(0,0.8,0,0.8), new=TRUE)
  plot(x, y , xlab=xlabText,
       ylab=ylabText)
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  boxplot(x, horizontal=TRUE, axes=FALSE)
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  boxplot(y, axes=FALSE)
  mtext(title, side=3, outer=TRUE, line=-3)
  
}


