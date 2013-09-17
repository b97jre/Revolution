# TODO: Still some things to do regarding filters.... to be updated
#
# This file was to identify different reigons that is hard for mapping.
#To run all it should be possible to write 
# main()
# with the correct files and all should be taken care of 
# Author: johanreimegard
###############################################################################


#plot CtrlFREEC results
main <- function(dataDir = "/Users/johanreimegard/Vetenskap/Data/capsella/Copy Number Variation/ctrlfreec_bowtie/ratiofiles"){
  
  
  stepSize=50000
  SNPfile <- c("BWA_genome.raw.vcf.Rfriendly")
  repeatFile <- c("Crubella_183.fa.out.noheader.bed.bin")
  CtrlFreecFileMerged <- c("all_FREEC50k.sorted.merged.bed")
  CtrlFreecFileExample <- c("CtrlFreec.r1GR1-2-KS3_DNA.STAR.ratio.txt")
  
  binData <- loadBinData(CtrlFreecFileMerged,SNPfile,stepSize,repeatFile)
  binData <- setCutoffValues(CtrlFreecFileMerged,SNPfile,stepSize,repeatFile,binData,
                              cumulativeDistributionRepeat=0.7,cumulativeDistributionHetero=0.8)
  
}




ExtractBins <- function(distribution,cumulativeDistribution=0.8){
  EmpiricalDensityFunction <- ecdf(distribution)
  InverseEmpiricalDensityFunction <- getInverseFunction(EmpiricalDensityFunction)
  
  cutoff=InverseEmpiricalDensityFunction(cumulativeDistribution)
  
  return (cutoff$root)
  
}


inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

getInverseFunction <- function(f){
  return (inverse(f, 0.01, 100))
}

loadBinsData <-function(CtrlFreecfile,SNPfile,stepSize,repeatFile){
  
  #concatenate all files
  binData <- data.frame()
  ratio <- read.table(paste(CtrlFreecfile), header=T)
  binData <- rbind(binData,ratio)
  
  binData$Diploid=1
  notDiploid <- which(binData$CopyNumber!= 2)
  binData$Diploid[notDiploid]= 0
  
  binData$homozygot=1
  binData$FractionHeterozygozity=0
  
  binData$repeatFraction = 0
  binData$lowRepeat = 1     
  
  
  binData$keptForAnalysis=1
  
  SNPs <- read.table(SNPfile,sep="\t",header=TRUE)
  SNPs.common <- subset(SNPs, Cr1GR1.2.KS3_Call=="0/0"|Cr1GR1.2.KS3_Call=="0/1"|Cr1GR1.2.KS3_Call=="1/1")
  
  repeatInfo <- read.table(repeatFile, header=F)
  
  scaffolds <- unique(SNPs.common$CHROM)
  
  print("checking scaffolds")
  
  for (i in (1:length(scaffolds))) {
    tt <- which(binData$Chromosome==as.character(scaffolds[[i]]))
    tr <- which(repeatInfo$V1==as.character(scaffolds[[i]]))
    starts <- binData$Start[tt]
    breakPoints <- c(0,starts+stepSize)
    
    if (length(tt)>0) {
      temp = subset(SNPs.common,CHROM==as.character(scaffolds[i]))
      heteroHist <- hist(temp$POS[temp$Cr1GR1.2.KS3_Call=="0/1"],breaks=breakPoints,plot=FALSE)
      majorHist <- hist(temp$POS[temp$Cr1GR1.2.KS3_Call=="0/0"],breaks=breakPoints,plot=FALSE)
      minorHist <- hist(temp$POS[temp$Cr1GR1.2.KS3_Call=="1/1"],breaks=breakPoints,plot=FALSE)
      binData$FractionHeterozygozity[tt] = heteroHist$counts/majorHist$counts 
      
      if(length(tr) == length(tt)){
        binData$repeatFraction[tt] = repeatInfo$V4[tr]     
      }
    }
    cat(scaffolds[[i]],sep = ".")
  }
  print("Done")
  return(binData)
  
  
}

setCutoffValues <- function(SNPfile,stepSize,repeatFile,binData,
                            cumulativeDistributionRepeat=0.7,cumulativeDistributionHetero=0.8){
  
  pdf("CutoffSelectionDistribution.pdf")
  par(mfrow=c(2,1), bty="l", cex=0.6)
  
  EmpiricalDensityFunction <- ecdf(binData$repeatFraction)
  RepeatCutoff= getCutoff(binData$repeatFraction,cumulativeDistributionRepeat)
  plot(EmpiricalDensityFunction,main= "repeat cumulative distribution over chromosomes",xlab="Fraction bp annotated as repeats by repeatmasker")
  abline(h=cumulativeDistributionRepeat, col="red")
  abline(v=RepeatCutoff, col="red")
  
  FractionHeterozygozity <- binData$FractionHeterozygozity
  EmpiricalDensityFunction <- ecdf(FractionHeterozygozity)
  HeterozygozityCutoff= getCutoff(FractionHeterozygozity,cumulativeDistributionHetero)
  plot(EmpiricalDensityFunction,main= "Heterozygozity cumulative distribution over chromosomes",xlab="Fraction heterozygous SNPs")
  abline(h=cumulativeDistributionHetero, col="red")
  abline(v=HeterozygozityCutoff, col="red")
  
  
  
  dev.off()
  
  
  
  binData$homozygot = binData$FractionHeterozygozity<HeterozygozityCutoff 
  binData$lowRepeat = binData$repeatFraction<RepeatCutoff     
  binData$keptForAnalysis=binData$homozygot*binData$lowRepeat*binData$Diploid
  
  
  fractions <-c(CNVs = length(binData$Diploid[binData$Diploid==0])/length(binData$Diploid[]), Heterozygousity=length(binData$Diploid[binData$homozygot==0])/length(binData$Diploid[]),
                RepeatRegions=length(binData$Diploid[binData$lowRepeat==0])/length(binData$Diploid[]),total=length(binData$Diploid[binData$keptForAnalysis == 0])/length(binData$Diploid[]))
  
  heterozygout <- binData[which(binData$homozygot==0),]  
  heterozygout$Stop = stepSize+heterozygout$Start 
  heterozygout[3:11]<- list(NULL)
  heterozygout$Name <- paste(heterozygout$Chromosome,heterozygout$Start,sep="_")
  
  HetFileName = paste(SNPfile,".bed",sep=".heterozygousRegions")
  write.table(heterozygout, file = HetFileName, quote = FALSE, sep = "\t",row.names=FALSE,col.names=FALSE)
  
  repeatRegions <- binData[which(binData$lowRepeat==0),]  
  repeatRegions$Stop = stepSize+repeatRegions$Start 
  repeatRegions[3:11]<- list(NULL)
  repeatRegions$Name <- paste(repeatRegions$Chromosome,repeatRegions$Start,sep="_")
  
  RepeatFileName = paste(SNPfile,".bed",sep=".repeatRegoins")
  write.table(repeatRegions, file = RepeatFileName, quote = FALSE, sep = "\t",row.names=FALSE,col.names=FALSE)
  
  
  for(j in seq(1,7, by=2)){
    
    pdfName <- paste("kept_for_analysis_scaffold",j,(j+1),"distribution.pdf",sep= "_")
    pdf(pdfName)    
    
    par(mfrow=c(4,2), bty="l", cex=0.6)
    for (i in seq(j,j+1)) {
      print(i)
      tt <- which(binData$Chromosome==paste("scaffold_",i,sep=""))
      if (length(tt)>0) {
        par(mar=c(5,4,4,4))
        #uncomment next line for individual plots
        plot(binData$Start[tt],binData$keptForAnalysis[tt],pch = 20,col = colors()[88],xlab = paste ("position, chr",i),ylab="", type='l')
        mtext("Kept for analysis (Total)",side=2,line=2,col=1,cex=0.5)
      }
    }
    for (i in seq(j,j+1)) {
      tt <- which(binData$Chromosome==paste("scaffold_",i,sep=""))
      if (length(tt)>0) {
        par(mar=c(5,4,4,4))
        #uncomment next line for individual plots
        plot(binData$Start[tt],binData$Diploid[tt],pch = 20,xlab="",,ylab="", type='l')
        par(new=T)
        plot(binData$Start[tt],binData$CopyNumber[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
        axis(4,col=4,col.lab=4,col.axis = colors()[99])
        
        mtext("Copy number (CTRLfreec)",side=4,line=2,col=colors()[99],cex=0.5)
        mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
      }
    }
    
    for (i in seq(j,j+1)) {
      tt <- which(binData$Chromosome==paste("scaffold_",i,sep=""))
      if (length(tt)>0) {
        par(mar=c(5,4,4,4))
        #uncomment next line for individual plots
        plot(binData$Start[tt],binData$lowRepeat[tt],pch = 20,xlab="",,ylab="", type='l')
        par(new=T)
        plot(binData$Start[tt],binData$repeatFraction[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
        axis(4,col=4,col.lab=4,col.axis = colors()[99])
        mtext("Fraction repeat(repeatMasker)",side=4,line=2,col= colors()[99],cex=0.5)
        mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
      }
    }
    
    for (i in seq(j,j+1)) {
      tt <- which(binData$Chromosome==paste("scaffold_",i,sep=""))
      if (length(tt)>0) {
        par(mar=c(5,4,4,4))
        plot(binData$Start[tt],binData$homozygot[tt],pch = 20,xlab="",,ylab="", type='l')
        par(new=T)
        plot(binData$Start[tt],binData$FractionHeterozygozity[tt],pch = 20,col = colors()[99],axes=FALSE,xlab = paste ("position, chr",i),ylab="", type='l')
        axis(4,col=4,col.lab=4,col.axis = colors()[99])
        mtext("Fraction heterozygout SNPs (GATK)",side=4,line=2,col= colors()[99],cex=0.5)
        mtext("Kept for analysis",side=2,line=2,col=1,cex=0.5)
      }
    }
    dev.off()
  }
  
  
  
  pdfName <- paste("kept_for_analysis_scaffold",".summary.distribution.pdf",sep= "_")
  pdf(pdfName)    
  
  par(mfrow=c(4,2), bty="l", cex=0.6)
  for (i in 1:8) {
    print(i)
    tt <- which(binData$Chromosome==paste("scaffold_",i,sep=""))
    if (length(tt)>0) {
      par(mar=c(5,4,4,4))
      chromsomeSpecific = binData[tt,]
      
      Total=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$keptForAnalysis==0])/length(tt)
      Repeat=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$lowRepeat==FALSE])/length(tt)
      Heterozygout=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$homozygot==FALSE])/length(tt)
      CNVs=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$Diploid==0])/length(tt)
      barplot(c(Total,Repeat,Heterozygout,CNVs),ylim=c(0,1),names.arg=c("Total","Repeats","Heterozygout","CNVs"),main=paste("scaffold_",i,sep=""))
      
    }
  }
  
  dev.off()
  
  
  pdfName <- paste("kept_for_analysis_scaffold",".summary.distribution.pdf",sep= "_All_")
  pdf(pdfName)    
  
  chromsomeSpecific = binData
  
  Total=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$keptForAnalysis==0])/length(chromsomeSpecific$keptForAnalysis)
  Repeat=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$lowRepeat==FALSE])/length(chromsomeSpecific$keptForAnalysis)
  Heterozygout=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$homozygot==FALSE])/length(chromsomeSpecific$keptForAnalysis)
  CNVs=length(chromsomeSpecific$keptForAnalysis[chromsomeSpecific$Diploid==0])/length(chromsomeSpecific$keptForAnalysis)
  barplot(c(Total,Repeat,Heterozygout,CNVs),ylim=c(0,1),names.arg=c("Total","Repeats","Heterozygout","CNVs"),main=paste("All scaffolds"))
  
  dev.off()
  
  
  return (binData)
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


BarplotSNPsCounts <- function(SNPCountfile = "vcfCount.info"){
  countTable <- read.table(SNPCountfile,header=FALSE,sep=" ")
  countTable$fractionAll=countTable$V1/countTable$V1[1]
  countTable$fractionExon=countTable$V1/countTable$V1[9]
  countTable$names <- gsub("\\.vcf","",substring(countTable$V2 , 18)) 
  countTable$names <- gsub("\\."," & ",countTable$names) 
  countTable$names[1] = "None" 
  
  pdf(paste(SNPCountfile,"pdf",sep="."))
  par(mfrow=c(1,2), bty="l", cex=0.6)
  # Fitting Labels 
  par(las=2) # make label text perpendicular to axis
  par(mar=c(15,4,4,2)) # increase y-axis margin.
  
  barplot(countTable$fractionAll[1:8], main="Fraction SNP after filter",  names.arg=countTable$names[1:8], cex.names=0.5)
  barplot(countTable$fractionExon[9:16], main="Fraction SNP after filter in annotated exons",  names.arg=countTable$names[9:16], cex.names=0.5)
  dev.off()
  
}


