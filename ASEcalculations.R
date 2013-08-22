# TODO: Add comment
# 
# Author: johanreimegard
###############################################################################


main <-function(){
	#load dataset	
	sampleData <- read.table("Unified.output.raw.snps.indels.Final.not_centromere.Crubella_183_only_exons_unique.heterozygous.homozygous.vcf.Sample.vcf.Rfriendly",header=TRUE,sep ="\t")

	# read in different samples
	DNAInterSamples <- c("Inter3.1","Inter4.1","Inter5.1")
	DNAIntraSamples <- c("Intra6.3","Intra7.2","Intra8.2")
	RNAInterSamples <- c("Inter3_1_1_F","Inter3_1_1_L","Inter4_1_1_F","Inter4_1_1_L","Inter4_1_2_F","Inter4_1_2_L","Inter4_1_3_F","Inter4_1_4_L","Inter5_1_1_F","Inter5_1_1_L")
	RNAIntraSamples <- c("Intra6_3_1_F","Intra6_3_1_L","Intra7_2_1_F","Intra7_2_1_L","Intra7_2_2_F","Intra7_2_2_L","Intra7_2_3_F","Intra7_2_3_L","Intra8_2_1_F","Intra8_2_1_L")
	AllSamples <- c(DNAInterSamples,DNAIntraSamples,RNAInterSamples,RNAIntraSamples)
	
	# Different filters
	# Only keep heterozygous samples
	DataSetInterHetero <- filterCallData(sampleData, DNAInterSamples)  
	# Only keep heterozygous samples with DNA total count above 30 in all heterozygous samples
	DataSetInterHeteroAbundant <- filterTotalData(DataSetInterHetero, DNAInterSamples,c(30,200))  
	# Only keep heterozygous samples with DNA total count above 40 in all heterozygous samples
	DataSetInterHeteroAbundantBOTH <- filterTotalData(DataSetInterHeteroAbundant, RNAInterSamples,c(40,2000))  
	
	# Calculate pValues 
	DataSetInterHeteroAbundantBOTHwithPvalue <- Pvaluescalculation(DataSetInterHeteroAbundantBOTH, RNAInterSamples,DNAInterSamples)  
	
	
	#Under development
	printPvalues(DataSetInterHeteroAbundantBOTHwithPvalue,RNAInterSamples,DNAInterSamples,"test.pdf")	
	
	
	
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




filterCallData  <- function(testData , samples, call="0/1"){
	
	# filter for heterozygous sites in all samples
	for(i in 1:length(samples)){
		name = paste(samples[i], "Call", sep="_")
		testData <- testData[testData[name]==call,]		
	}
	return (testData)
}

filterTotalData <- function(Dataset , samples,cutoff=c(40,200)){
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


BinTest <- function(x, count1RNA, TotalRNA,count1DNA,TotalDNA){
	temp <- binom.test(as.numeric(x[count1RNA]),as.numeric(x[TotalRNA]),p=(as.numeric(x[count1DNA])/as.numeric(x[TotalDNA])), alternative="two.sided")
	return (temp$p.value)
}


