# TODO: Add comment
# 
# Author: johanreimegard

###############################################################################


main <-function(fileName="BWA_genome.raw.all_FREEC50k.repeatRegions.heterozygousRegions.Crubella_183_only_exons_unique.vcf.Sample.vcf.Rfriendly", dataDir="/gulo/proj_nobackup/b2012122/private/ASE",
                annotation="BWA_genome.raw.all_FREEC50k.repeatRegions.heterozygousRegions.Crubella_183_only_exons_unique.vcf.Sample.Crubella_183_gene.VCFannotaion",  
                rounds = 10 ,cutoffNominal=0.005, cutoffAdjusted=0.1){
  
  
  file <- paste(dataDir,fileName,sep="/")
  annotationFile <- paste(dataDir,annotation,sep="/")
  
  #load dataset  
  print("Reading data file")
  sampleData <- read.table(file,header=TRUE,sep ="\t")
  print("Reading annotation file")
  annotationData <- read.table(annotationFile,header=FALSE,sep ="\t")
  sampleData$annotation=annotationData$V3
  
  
  # read in different samples kinds od samples will not go through all though
  DNAInterSamples <- c("Inter3.1","Inter4.1","Inter5.1")
  DNAIntraSamples <- c("Intra6.3","Intra7.2","Intra8.2")
  RNAInterSamples <- c("Inter3_1_1_F","Inter3_1_1_L","Inter4_1_1_F","Inter4_1_1_L","Inter4_1_2_F","Inter4_1_2_L","Inter4_1_3_F","Inter4_1_4_L","Inter5_1_1_F","Inter5_1_1_L")
  RNAIntraSamples <- c("Intra6_3_F","Intra6_3_L","Intra7_2_1_F","Intra7_2_1_L","Intra7_2_2_F","Intra7_2_2_L","Intra7_2_3_F","Intra7_2_3_L","Intra8_2_1_F","Intra8_2_1_L")
  AllSamples <- c(DNAInterSamples,DNAIntraSamples,RNAInterSamples,RNAIntraSamples)
  
  
  ASEinfo <- getASEinfo(sampleData,DNAsamples,RNAsamples) 

  
  
  
  
  getASEinfo <- function(sampleData,DNAsamples,RNAsamples, RNAlow=40, RNAhigh=2000,DNAlow=30,DNAhigh=200){
    
    
    ASEinfo= list()
    cutoffs <- c(RNAlow=RNAlow, RNAhigh=RNAhigh,DNAlow=DNAlow,DNAhigh=DNAhigh)
    
    
    ASEinfo <- c(ASEinfo,DNAsamples)
    ASEinfo <- c(ASEinfo,RNAsamples)
    ASEinfo <- c(ASEinfo, sample_Filtered )
    #Filter for SNPs tthat is present in all samples with expression can be changed with RNAlow and RNA high and so on 
    sample_Filtered <- filterData(sampleData,DNAsamples,RNAsamples, RNAlow=40, RNAhigh=2000,DNAlow=30,DNAhigh=200)
    ASEinfo <- c(ASEinfo, sample_Filtered )

    #get OneSNP per Gene
    sample_Filtered_1SNPperGene <- getOneSNPperGene(sample_Filtered)
    #calculate Pvalues
    sample_Filtered_1SNPperGene_pValues <- getPvalues2(sample_Filtered_1SNPperGene,DNAsamples,RNAsamples)
    ASEinfo <- c(ASEinfo, sample_Filtered_1SNPperGene_pValues )
    
    # checking number of genes 
    heterozygousGenes = unique(sample_Filtered_1SNPperGene_pValues$annotation)
    ASEinfo <- c(ASEinfo, heterozygousGenes)
    
    #Union of genes that has ASE
    InterUnionGenesFlowers <- getASEgenesUnion(sample_Filtered_1SNPperGene_pValues,RNAsamples,cutoffAdjusted=0.005)
    #intersection of genes that has ASE
    InterIntersectGenesFlowers <- getASEgenesIntersect(Inter4_1_SampleDataFlowers,Inter4_1_RNA_Flower,cutoffAdjusted=0.005)
    
    
    
    
  }
    
  
  
  # 
  print("Checking for ASE")
 
  Inter4_1_DNA <- c("Inter4.1")
  Inter4_1_RNA <-c("Inter4_1_1_F","Inter4_1_1_L","Inter4_1_2_F","Inter4_1_2_L","Inter4_1_3_F","Inter4_1_4_L")  
  Inter4_1_RNA_Flower <-c("Inter4_1_1_F","Inter4_1_2_F","Inter4_1_3_F")
  Inter4_1_RNA_Leafs <-c("Inter4_1_1_L","Inter4_1_2_L","Inter4_1_4_L")
  
  
  #Filter for SNPs tthat is present in all samples with expression can be changed with RNAlow and RNA high and so on 
  #calculate Pvalues
  
  # checking number of genes
  nrOfheterozygousGenesExpressedInIntraFlowers = length(unique(Inter4_1_SampleDataFlowers$annotation))
  # get One SNP per gene that is the same for all downstream analysis 
  Inter4_1_SampleDataFlowers1SNPperGene <- getOneSNPperGene(Inter4_1_SampleDataFlowers)
  
  # Check number of genes that complied to filter
  nrOfheterozygousGenesExpressedInIntraFlowers = length(unique(Inter4_1_SampleDataFlowers$annotation))

  #Union of genes that has ASE
  InterUnionGenesFlowers <- getASEgenesUnion(Inter4_1_SampleDataFlowers,Inter4_1_RNA_Flower,cutoffAdjusted=0.005)
  #intersection of genes that has ASE
  InterIntersectGenesFlowers <- getASEgenesIntersect(Inter4_1_SampleDataFlowers,Inter4_1_RNA_Flower,cutoffAdjusted=0.005)
  
  
  Inter4_1_SampleDataLeafs <- getPvalues(sampleData,Inter4_1_DNA,Inter4_1_RNA_Leafs)
  nrOfheterozygousGenesExpressedInIntraLeafs = length(unique(Inter4_1_SampleDataLeafs$annotation))
  InterUnionGenesLeafs <- getASEgenesUnion(Inter4_1_SampleDataLeafs,Inter4_1_RNA_Leafs,cutoffAdjusted=0.005)
  InterIntersectGenesLeafs <- getASEgenesIntersect(Inter4_1_SampleDataLeafs,Inter4_1_RNA_Leafs,cutoffAdjusted=0.005)
  
  
  
  Both <- intersect(unique(Inter4_1_SampleDataFlowers$annotation),unique(Inter4_1_SampleDataLeafs$annotation))
  FlowerUnique <- detDiff(unique(Inter4_1_SampleDataFlowers$annotation),Both)
  nrOfheterozygousGenesExpressedInIntra = length(Both)
  BothASEUnion = intersect(UnionGenesFlowers,UnionGenesLeafs)
  BothASEIntersect = intersect(IntersectGenesFlowers,IntersectGenesLeafs)
  
  
  
  
  
  
  Inter4_1_SampleData <- getPvalues(sampleData,Inter4_1_DNA,Inter4_1_RNA)
  ASEinfo <- rbind(ASEinfo, getASEinfo(Inter4_1_SampleData,Inter4_1_RNA,rounds,cutoffNominal,cutoffAdjusted))
  
  Inter5_1_DNA <- c("Inter5.1")
  Inter5_1_RNA <-c("Inter5_1_1_F","Inter5_1_1_L")  
  Inter5_1_SampleData <- getPvalues(sampleData,Inter5_1_DNA,Inter5_1_RNA)
  ASEinfo <- rbind(ASEinfo, getASEinfo(Inter5_1_SampleData,Inter5_1_RNA,rounds,cutoffNominal,cutoffAdjusted))
  
  Intra6_3_DNA <- c("Intra6.3")
  Intra6_3_RNA <-c("Intra6_3_F","Intra6_3_L")
  Intra6_3_SampleData <- getPvalues(sampleData,Intra6_3_DNA,Intra6_3_RNA)
  ASEinfo <- rbind(ASEinfo, getASEinfo(Intra6_3_SampleData,Intra6_3_RNA,rounds,cutoffNominal,cutoffAdjusted))
  
  Intra7_2_DNA <- c("Intra7.2")
  Intra7_2_RNA_Flower <-c("Intra7_2_1_F","Intra7_2_2_F","Intra7_2_3_F")
  Intra7_2_RNA_Leafs <-c("Intra7_2_1_L","Intra7_2_2_L","Intra7_2_3_L")
  
  Intra7_2_SampleDataFlower1 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Flower[1])
  Intra7_2_SampleDataFlower2 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Flower[2])
  Intra7_2_SampleDataFlower3 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Flower[3])
  
  Intra7_2_SampleDataLeafs1 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Leafs[1])
  Intra7_2_SampleDataLeafs2 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Leafs[2])
  Intra7_2_SampleDataLeafs3 <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Leafs[3])
  
  
  #
  Intra7_2_SampleDataFlowers <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Flower)
  nrOfheterozygousGenesExpressedInIntraFlowers = length(unique(Intra7_2_SampleDataFlowers$annotation))
  IntraUnionGenesFlowers <- getASEgenesUnion(Intra7_2_SampleDataFlowers,Intra7_2_RNA_Flower,cutoffAdjusted=0.005)
  IntraIntersectGenesFlowers <- getASEgenesIntersect(Intra7_2_SampleDataFlowers,Intra7_2_RNA_Flower,cutoffAdjusted=0.005)
  
  
  Intra7_2_SampleDataLeafs <- getPvalues(sampleData,Intra7_2_DNA,Intra7_2_RNA_Leafs)
  nrOfheterozygousGenesExpressedInIntraLeafs = length(unique(Intra7_2_SampleDataLeafs$annotation))
  IntraUnionGenesLeafs <- getASEgenesUnion(Intra7_2_SampleDataLeafs,Intra7_2_RNA_Leafs,cutoffAdjusted=0.005)
  IntraIntersectGenesLeafs <- getASEgenesIntersect(Intra7_2_SampleDataLeafs,Intra7_2_RNA_Leafs,cutoffAdjusted=0.005)
  
  Both <- intersect(unique(Intra7_2_SampleDataFlowers$annotation),unique(Intra7_2_SampleDataLeafs$annotation))
  nrOfheterozygousGenesExpressedInIntra = length(Both)
  BothASEUnion = intersect(UnionGenesFlowers,UnionGenesLeafs)
  BothASEIntersect = intersect(IntersectGenesFlowers,IntersectGenesLeafs)

  InterUnionGenesFlowers
  InterUnionGenesLeafs
  IntraUnionGenesFlowers
  IntraUnionGenesLeafs
  
  genes <- unique(sampleData$annotation)
  InterFlowersGenes = unique(Inter4_1_SampleDataFlowers$annotation)
  InterLeavesGenes = unique(Inter4_1_SampleDataLeafs$annotation)
  IntraFlowersGenes = unique(Intra7_2_SampleDataFlowers$annotation)
  IntraLeafsGenes = unique(Intra7_2_SampleDataLeafs$annotation)
  
  InterFlowersGenesIndex = match(InterFlowersGenes,genes)
  InterLeavesGenesIndex  = match(InterLeavesGenes,genes)
  IntraFlowersGenesIndex  = match(IntraFlowersGenes,genes)
  IntraLeafsGenesIndex  = match(IntraLeafsGenes,genes)

  
  
  
  length(setdiff(InterUnionGenesFlowers,InterUnionGenesLeafs))
  length(setdiff(InterUnionGenesLeafs,InterUnionGenesFlowers))
  length(union(InterUnionGenesFlowers,InterUnionGenesLeafs))
  
  length(setdiff(IntraUnionGenesFlowers,IntraUnionGenesLeafs))
  length(setdiff(IntraUnionGenesLeafs,IntraUnionGenesFlowers))
  length(union(IntraUnionGenesFlowers,IntraUnionGenesLeafs))
  length(intersect(IntraUnionGenesFlowers,IntraUnionGenesLeafs))
  
  
  IntraFlowersASE = IntraUnionGenesFlowers
  IntraLeafsASE = IntraUnionGenesLeafs
  
  
  ASEinfo <- rbind(ASEinfo, getASEinfo(Intra7_2_SampleData,Intra7_2_RNA,rounds,cutoffNominal,cutoffAdjusted))
  
  Intra8_2_DNA <- c("Intra8.2")
  Intra8_2_RNA <-c("Intra8_2_1_F","Intra8_2_1_L")
  Intra8_2_SampleData <- getPvalues(sampleData,Intra8_2_DNA,Intra8_2_RNA)
  ASEinfo <- rbind(ASEinfo, getASEinfo(Intra8_2_SampleData,Intra8_2_RNA,rounds,cutoffNominal,cutoffAdjusted))
  
  #Under development
  #printPvalues(DataSetInterHeteroAbundantBOTHwithPvalue,RNAInterSamples,DNAInterSamples,"test.pdf")	
  return (ASEinfo)
  
}


getNrOfGenes <-function(CtrlFreecfile,SNPs.common,stepSize,repeatFile){
  
  #concatenate all files
  binData <- data.frame()
  ratio <- read.table(paste(CtrlFreecfile), header=T)
  binData <- rbind(binData,ratio)

  binData$nrOfGenesGenes=0
  
  
  
  binData$homozygot=1
  binData$FractionHeterozygozity=0
  binData$totalCounts = 0
  binData$heteroHistCounts = 0
  binData$majorHistCounts = 0 
  binData$minorHistCounts = 0 
  
  
  
  binData$repeatFraction = 0
  binData$lowRepeat = 1     
  
  
  binData$keptForAnalysis=1
  
  SNPs.common 
  
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
      binData$FractionHeterozygozity[tt] = heteroHist$counts/ totalCounts 
      
      
      totalCounts = majorHist$counts+ majorHist$counts+minorHist$counts
      binData$totalCounts[tt] = totalCounts
      binData$heteroHistCounts[tt] = heteroHist$counts 
      binData$majorHistCounts[tt] = majorHist$counts 
      binData$minorHistCounts[tt] = minorHist$counts 
      
      if(length(tr) == length(tt)){
        binData$repeatFraction[tt] = repeatInfo$V4[tr]     
      }
    }
    cat(scaffolds[[i]],sep = ".")
  }
  print("Done")
  return(binData)
  
  
  
}


filterData <- function(sampleData,DNAsamples,RNAsamples, RNAlow=40, RNAhigh=2000,DNAlow=30,DNAhigh=200){
  
  
  # Different filters
  # Only keep heterozygous DNAsamples
  DataSetFiltered <- filterCallData(sampleData, DNAsamples)  
  
  # Only keep heterozygous samples with DNA total count above 30 in all heterozygous samples
  DataSetFilteredAbundant <- filterTotalData(DataSetFiltered, DNAsamples,c(DNAlow,DNAhigh))  
  
  # Only keep heterozygous samples with RNA total count above 40 in all heterozygous samples
  DataSetFilteredAbundantBOTH <- filterTotalData(DataSetFilteredAbundant, RNAsamples,c(RNAlow,RNAhigh)) 
  return (DataSetFilteredAbundantBOTH)
}


getPvalues2 <- function(sampleData,DNAsamples,RNAsamples){
  # Calculate pValues 
  sampleDataPvalue <- Pvaluescalculation(sampleData, RNAsamples,DNAsamples)  
  
  # Correct for multiple testing
  sampleDataPvalue<-as.data.frame(sampleDataPvalue)
  sampleDataPvalueCorrected <- MultipleTestCorrection(sampleDataPvalue,RNAsamples)
  
  return(sampleDataPvalueCorrected)
  
}





printPvalues <- function(Dataset,RNAsamples,DNAsamples, pdfFileName){
  
  pdf(pdfFileName)
  par(mfrow=c(length(RNAsamples),2), bty="l", cex=0.6)
  for(i in 1:length(RNAsamples)){
    for(j in 1:length(DNAsamples[j])){
      RNAsample = strsplit(RNAsamples[i], '_')
      DNAsample = strsplit(DNAsamples[j], '\\.')
      if( RNAsample[[1]][1]==DNAsample[[1]][1]){
        
        name = paste(RNAsamples[i], "pValue", sep="_")
        hist(Dataset[name],breaks=40)		
        name = paste(RNAsamples[i], "fraction", sep="_")
        name2 = paste(DNAsamples[j], "fraction", sep="_")
        boxScatterPlot(Dataset[name],Dataset[name],ylabText="DNA distribution",xlabText="RNA distribtuion",title = paste( "ASE distribution"))	
      }
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


filterForHeterozygousSites  <- function(sampleData , samples, call="0/1"){
  # filter for heterozygous sites in all samples
  for(i in 1:length(samples)){
    name = paste(samples[i], "Call", sep="_")
    sampleData <- sampleData[sampleData[name]==call,]		
  }
  return (testData)
}


filterForTotalCount <- function(Dataset , samples,cutoff=c(40,200)){
  for(i in 1:length(samples)){
    name = paste(samples[i], "Total", sep="_")
    Dataset <- Dataset[Dataset[name]>cutoff[1]	,]		
  }
  return (Dataset)
}



Pvaluescalculation <- function(testData,RNAsamples,DNAsamples){
  for(i in 1:length(RNAsamples)){
    for(j in 1:length(DNAsamples)){
      RNAsample = strsplit(RNAsamples[i], '_')
      DNAsample = strsplit(DNAsamples[j], '\\.')
      if( RNAsample[[1]][1]==DNAsample[[1]][1]){
        print(paste(DNAsamples[j], RNAsamples[i], sep = " comparing "))
        count1RNA=paste(RNAsamples[i], "Count1", sep="_")
        TotalRNA=paste(RNAsamples[i], "Total", sep="_")
        count1DNA=paste(DNAsamples[j], "Count1", sep="_")
        TotalDNA=paste(DNAsamples[j], "Total", sep="_")
        pValuesName= paste(RNAsamples[i], "pValue", sep="_")
        #print(testData[count1RNA]/testData[TotalRNA])
        #print(paste(testData[count1DNA][1], testData[TotalDNA][1], sep = " of total "))
        
        
        testData[pValuesName]<- apply(testData,1, BinTest,  count1RNA, TotalRNA, count1DNA, TotalDNA)
      }
    }
  }
  return (testData) 
  
}


MultipleTestCorrection <- function(testData,RNAsamples) {
  
  for(i in 1:length(RNAsamples)){
    RNAsample = strsplit(RNAsamples[i], '_')
    pValue=paste(RNAsamples[i], "pValue", sep="_")
    pValue2= as.numeric(testData[,which(colnames(testData)==pValue)])
    pValuesCorrectedName=paste(RNAsamples[i], "pValueCorrected", sep="_")
    testData[pValuesCorrectedName]<- p.adjust(pValue2,method="BH")
  }
  return(testData)
}




BinTest <- function(x, count1RNA, TotalRNA,count1DNA,TotalDNA){
  temp <- binom.test(as.numeric(x[count1RNA]),as.numeric(x[TotalRNA]),p=(as.numeric(x[count1DNA])/as.numeric(x[TotalDNA])), alternative="two.sided")
  return (temp$p.value)
}


getASEinfo <- function(Dataset,RNAsamples,rounds=10,cutoffNominal=0.005, cutoffAdjusted=0.1){
  info = data.frame()
  for(i in (1:rounds)){
    DatasetOneSNPperGene <- getOneSNPperGene(Dataset) 
    rowInfo <- getASEinformation(DatasetOneSNPperGene,RNAsamples,cutoffNominal, cutoffAdjusted,i)
    info <- rbind(info,rowInfo)
  }  
  return (info)
}




getASEinformation <- function(DatasetOneSNPperGene,RNAsamples,cutoffNominal=0.005, cutoffAdjusted=0.1,round=1){
  info = data.frame()
  
  for(i in seq(1,length(RNAsamples)-1)){
    
    
    columnName  <- paste(RNAsamples[i],"pValueCorrected",sep="_")
    columnName2  <- paste(RNAsamples[i+1],"pValueCorrected",sep="_")
    
    Flowers <- DatasetOneSNPperGene[which(DatasetOneSNPperGene[columnName] <cutoffAdjusted),]
    
    Leafs <- DatasetOneSNPperGene[which(DatasetOneSNPperGene[columnName2] <cutoffAdjusted),]
    LeavesList(RNAsamples[i])=Leafs$annotation
    nrOfASEAdjustedFlowers <- length (Flowers[[1]])
    nrOfASEAdjustedLeafs <- length (Leafs[[1]])
    GeneIntersect <- intersect(Flowers$annotation,Leafs$annotation)
    GeneFlower <- setdiff(Flowers$annotation,Leafs$annotation)
    GeneLeafs <- setdiff(Leafs$annotation,Flowers$annotation)
    nrOfASEAdjustedFlowersAndLeaves <- length (GeneIntersect)
    nrOfASEAdjustedFlowersOnly <- length (GeneFlower)
    nrOfASEAdjustedLeafsOnly <- length (GeneLeafs)
    
    
    rowinfo <- data.frame("Name"=RNAsamples[i],"NrOfGenes"=nrOfGenes,"Adjusted_Pvalue"=cutoffAdjusted,"Adjusted_count"=nrOfASEAdjusted)
    info <- rbind(info,rowinfo)
    
  }
  return (info)
}

getASEgenesUnion <- function(Dataset,RNAsamples,rounds=10,cutoffNominal=0.005, cutoffAdjusted=0.005){
  DatasetOneSNPperGene <- getOneSNPperGene(Dataset) 
  info=NULL
  for(i in 1:length(RNAsamples)){
    rowInfo <- getASEGenes(DatasetOneSNPperGene,RNAsamples[i],cutoffNominal, cutoffAdjusted,i)
    info <- union(info,rowInfo)
  }
  return (info)
}




getASEgenesIntersect <- function(Dataset,RNAsamples,rounds=10,cutoffNominal=0.005, cutoffAdjusted=0.005){
  DatasetOneSNPperGene <- getOneSNPperGene(Dataset) 
  info=NULL
  for(i in 1:length(RNAsamples)){
    rowInfo <- getASEGenes(DatasetOneSNPperGene,RNAsamples[i],cutoffNominal, cutoffAdjusted,i)
    if(i ==1 ){
      info=rowInfo
    }else{
      info <- intersect(info,rowInfo)
    }
  }
  return (info)
}


getASEGenesAdjusted <- function(DatasetOneSNPperGene,RNAsample, cutoffAdjusted=0.1){
  columnName  <- paste(RNAsample,"pValueCorrected",sep="_")
  dataSample <- DatasetOneSNPperGene[which(DatasetOneSNPperGene[columnName] <cutoffAdjusted),]
  return(Flowers$annotation)

}

getASEGenesNominal <- function(DatasetOneSNPperGene,RNAsample, cutoff=0.1){
  columnName  <- paste(RNAsample,"pValue",sep="_")
  dataSample <- DatasetOneSNPperGene[which(DatasetOneSNPperGene[columnName] <cutoffNominal),]
  return(Flowers$annotation)
  
}


getOneSNPperGene <- function(Dataset,classifier="annotation"){
  Dataset$duplicated <- duplicated.random(as.vector(Dataset[,classifier]))
  DataSetOneSNPperGene <- subset(Dataset, duplicated == FALSE)  
  return (DataSetOneSNPperGene)
}

duplicated.random <- function(x, incomparables = FALSE, ...){
  if ( is.vector(x) )
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else if ( is.matrix(x) )
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else
  {
    stop(paste("duplicated.random() only supports vectors",
               "matrices for now."))
  }
}

  plotVenn <- function(A,B,C,D){
  venn.diagram(
    x = list(A,B,C,D),
    filename = "1D-quadruple_Venn.tiff",
    col = "black",
    lty = "dotted",
    lwd = 4,
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 2.5,
    cat.fontfamily = "serif"
  );
  }
  


