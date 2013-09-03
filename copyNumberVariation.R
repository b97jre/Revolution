# TODO: Add comment
# 
# Author: johanreimegard
###############################################################################


#plot CtrlFREEC results
main <- function(){
  CtrlFreecfiles <- c("Cr1GR1-2-KS3_BWA_primary.bam_ratio.txt","Inter3-1_DNA_BWA_primary.bam_ratio.txt", "Inter4-1_DNA_BWA_primary.bam_ratio.txt", "Inter5-1_DNA_BWA_primary.bam_ratio.txt", "Intra6-3_DNA_BWA_primary.bam_ratio.txt", "Intra7-2_DNA_BWA_primary.bam_ratio.txt", "Intra8-2_DNA_BWA_primary.bam_ratio.txt", "Cr1GR1-2-KS3_Bowtie_primary.bam_ratio.txt","Inter3-1_DNA_Bowtie_primary.bam_ratio.txt", "Inter4-1_DNA_Bowtie_primary.bam_ratio.txt", "Inter5-1_DNA_Bowtie_primary.bam_ratio.txt", "Intra6-3_DNA_Bowtie_primary.bam_ratio.txt", "Intra7-2_DNA_Bowtie_primary.bam_ratio.txt", "Intra8-2_DNA_Bowtie_primary.bam_ratio.txt")
  SNPfiles <-  c("simple.1to8.Cr1GR1-2-KS3.info.Radapted" , "simple.1to8.Cr_39_1.info.Radapted" , "simple.1to8.Inter3-1.info.Radapted" , "simple.1to8.Inter4-1.info.Radapted" , "simple.1to8.Inter5-1.info.Radapted" , "simple.1to8.Intra6-3.info.Radapted" , "simple.1to8.Intra7-2.info.Radapted" , "simple.1to8.Intra8-2.info.Radapted")
  centromericFile <- c("crub_centromere_cutoffs.txt")
  repeatFile <- c("Crubella_183.fa.out.noheader.bed.bin")
  
  for (j in 1:length(SNPfiles)){
    print(paste("Analysis with "," START" ,sep = SNPfiles[j]))
    plotCNVsandSNPdistribution(CtrlFreecfiles[j],SNPfiles[j],centromericFile)
    print(paste("Analysis with "," FINISHED" ,sep = SNPfiles[j]))
  }
  
  
}

plotCNVsandSNPdistribution <-function(CtrlFreecfiles,SNPfile,centromericFile){
  
  #concatenate all files
  plotfile <- data.frame()
  for(i in 1:length(CtrlFreecfile
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


plotCNVsandSNPandRepeatsDistribution2 <-function(CtrlFreecfile,SNPfile,stepSize,repeatFile){
  
  #concatenate all files
  plotfile <- data.frame()
  ratio <- read.table(paste(CtrlFreecfile), header=T)
  plotfile <- rbind(plotfile,ratio)

  plotfile$Diploid=1
  notDiploid <- which(plotfile$CopyNumber!= 2)
  plotfile$Diploid[notDiploid]= 0
  
  plotfile$homozygot=1
  plotfile$FractionHeterozygozity=0
  
  plotfile$repeatFraction = 0
  plotfile$lowRepeat = 1     
  
  
  plotfile$keptForAnalysis=1
  
  
  
    
  
  SNPs <- read.table(SNPfile,sep="\t",header=TRUE)
  SNPs.common <- subset(SNPs, Call=="0/0"|Call=="0/1"|Call=="1/1")
  
  repeatInfo <- read.table(repeatFile, header=F)

  scaffolds <- unique(plotfile$Chromosome)
  
  
  for (i in (1:length(scaffolds))) {
    tt <- which(plotfile$Chromosome==as.character(scaffolds[[i]]))
    tr <- which(repeatInfo$V1==as.character(scaffolds[[i]]))
    starts <- plotfile$Start[tt]
    breakPoints <- c(0,starts+stepSize)
    
    if (length(tt)>0) {
      temp = subset(SNPs.common,chr==as.character(scaffolds[i]))
      heteroHist <- hist(temp$position[temp$Call=="0/1"],breaks=breakPoints,plot=FALSE)
      majorHist <- hist(temp$position[temp$Call=="0/0"],breaks=breakPoints,plot=FALSE)
      minorHist <- hist(temp$position[temp$Call=="1/1"],breaks=breakPoints,plot=FALSE)
      plotfile$FractionHeterozygozity[tt] = heteroHist$counts/majorHist$counts 
      plotfile$homozygot[tt] = heteroHist$counts/majorHist$counts<0.1 
    
    if(length(tr) == length(tt)){
      plotfile$repeatFraction[tt] = repeatInfo$V4[tr]     
      plotfile$lowRepeat[tt] = repeatInfo$V4[tr]<0.2     
    }
    }
  }

  pdf("CutoffSelectionDistribution.pdf")
  
  par(mfrow=c(2,1), bty="l", cex=0.6)
  
  hist(plotfile$repeatFraction,breaks=40, xlab="Fraction nucleotides annotated as repeats in 50 kb windows", main="Histogram of repeat distribution over chromosomes")
  abline(v=0.2,col="blue")
  
  hist(plotfile$FractionHeterozygozity, breaks=40, xlab="Fraction nucleotides annotated as heterozygous compared to major allele homozygous in 50 kb windows", main="Histogram of hererozygous distribution over chromosomes")
  abline(v=0.1,col="blue")
  dev.off()
  
  
  
  
  

  
  heterozygout <- plotfile[which(plotfile$homozygot==0),]  
  heterozygout$Stop = stepSize+heterozygout$Start 
  heterozygout[3:11]<- list(NULL)
  heterozygout$Name <- paste(heterozygout$Chromosome,heterozygout$Start,sep="_")
  
  HetFileName = paste(SNPfile,".bed",sep=".heterozygousRegions")
  write.table(heterozygout, file = HetFileName, quote = FALSE, sep = "\t",row.names=FALSE,col.names=FALSE)

repeatRegions <- plotfile[which(plotfile$lowRepeat==0),]  
repeatRegions$Stop = stepSize+repeatRegions$Start 
repeatRegions[3:11]<- list(NULL)
repeatRegions$Name <- paste(repeatRegions$Chromosome,repeatRegions$Start,sep="_")

RepeatFileName = paste(repeatFile,".bed",sep=".aggregate")
write.table(repeatRegions, file = RepeatFileName, quote = FALSE, sep = "\t",row.names=FALSE,col.names=FALSE)



  
  for(j in ())
    
  pdf("kept_for_analysis_scaffold_1_to_2 distribution.pdf",paper="a4r")
  ploidy=2
  
  
  par(mfrow=c(4,2), bty="l", cex=0.6)
  
  for (i in (1:2)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$keptForAnalysis[tt],pch = 20,col = colors()[88],xlab = paste ("position, chr",i),ylab="", type='l')
      mtext("Kept for analysis (Total)",side=2,line=2,col=1,cex=0.5)
    }
  }
  for (i in (1:2)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$Diploid[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$CopyNumber[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      
      mtext("Copy number (CTRLfreec)",side=4,line=2,col=colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  
  for (i in (1:2)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$lowRepeat[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$repeatFraction[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      mtext("Fraction repeat(repeatMasker)",side=4,line=2,col= colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  
  for (i in (1:2)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      plot(plotfile$Start[tt],plotfile$homozygot[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$FractionHeterozygozity[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      mtext("Fraction heterozygout SNPs (GATK)",side=4,line=2,col= colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  dev.off()
 
  pdf("kept_for_analysis_scaffold_5_to_8 distribution.pdf")
  ploidy=2
  
  
  par(mfrow=c(4,4), bty="l", cex=0.6)
  
  for (i in (5:8)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$keptForAnalysis[tt],pch = 20,col = colors()[88],xlab = paste ("position, chr",i),ylab="", type='l')
      mtext("Kept for analysis (Total)",side=2,line=2,col=1,cex=0.5)
    }
  }
  for (i in (5:8)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$Diploid[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$CopyNumber[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      
      mtext("Copy number (CTRLfreec)",side=4,line=2,col=colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  
  for (i in (5:8)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      #uncomment next line for individual plots
      plot(plotfile$Start[tt],plotfile$lowRepeat[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$repeatFraction[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      mtext("Fraction repeat(repeatMasker)",side=4,line=2,col= colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  
  for (i in (5:8)) {
    tt <- which(plotfile$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      plot(plotfile$Start[tt],plotfile$homozygot[tt],pch = 20,xlab="",,ylab="", type='l')
      par(new=T)
      plot(plotfile$Start[tt],plotfile$FractionHeterozygozity[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
      axis(4,col=4,col.lab=4,col.axis = colors()[99])
      mtext("Fraction heterozygout SNPs (GATK)",side=4,line=2,col= colors()[99],cex=0.5)
      mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
    }
  }
  dev.off()
  
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



