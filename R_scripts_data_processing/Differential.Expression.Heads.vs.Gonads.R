library(DESeq2)

#function to run DESeq2 for a given tissue type (one sex at a time)
TissuesDiffExpression <- function(sex, output_name){
  setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.DESeq2')
  filenames <- c(list.files(pattern = paste0(sex, '_heads')), list.files(pattern = paste0(sex, '_gonads')))
  samplenames <- gsub('.{4}$', '', filenames)
  sample.tissue <- factor(c(rep("Heads", 3), rep("Gonads", 3)))
  sample.table <- data.frame(sampleName = samplenames, fileName = filenames, tissue = sample.tissue)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table, design= ~tissue)
  keep.dds = rowMeans(counts(ddsHTSeq)) >= 50
  ddsHTSeq = ddsHTSeq[keep.dds,]
  DESeq.df = DESeq(ddsHTSeq)
  results = results(DESeq.df, contrast=c("tissue", "Heads", "Gonads"))
  results$FlyBaseID = rownames(results)
  results = results[,c(7,1:6)]
  
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs/Between.Tissues/")
  write.table(results, file = paste0('DifferentialGeneExpression.', output_name, ".txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  return(results)
}

TissuesDiffExpression('Female', 'Females.Heads.VS.Gonads.tsv')
TissuesDiffExpression('Male', 'Males.Heads.VS.Gonads.tsv')

# function to assign tissue-bias category to genes in the DESeq2 output
AssignTissueBiasCategory <- function(sex){
  
  #load deseq2 output; remove genes with high standard error in expression
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs/Between.Tissues/")
  sample <- read.table(paste0('DifferentialGeneExpression.', sex , '.Heads.VS.Gonads.txt'), header = T, sep = '\t')
  high.SE.genes <- sample[which(sample$lfcSE > 1.5),]$FlyBaseID
  sample <- sample[!(sample$FlyBaseID %in% high.SE.genes),]
  
  #assign tissue bias category: highly GB genes (-infinity, -2), moderately GB genes (-2, -0.5),
  #unbiased (-0.5, 0.5), moderately HB genes (0.5, 2), highly HB genes (2, infinity)
  high.gonad.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < -2),]$FlyBaseID, 
                                      Tissue.Bias.Category = rep("Highly gonad-biased", length(which(sample$log2FoldChange < -2))),
                                      log2FoldChange = sample[which(sample$log2FoldChange < -2),]$log2FoldChange)
  moderate.gonad.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2),]$FlyBaseID, 
                                          Tissue.Bias.Category = rep("Moderately gonad-biased", 
                                                                     length(which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2))),
                                          log2FoldChange = sample[which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2),]$log2FoldChange)
  unbiased <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5),]$FlyBaseID, 
                               Tissue.Bias.Category = rep("Unbiased", 
                                                          length(which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5))),
                               log2FoldChange = sample[which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5),]$log2FoldChange)
  moderate.head.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5),]$FlyBaseID, 
                                         Tissue.Bias.Category = rep("Moderately head-biased", 
                                                                    length(which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5))),
                                         log2FoldChange = sample[which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5),]$log2FoldChange)
  high.head.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange >= 2),]$FlyBaseID, 
                                     Tissue.Bias.Category = rep("Highly head-biased", length(which(sample$log2FoldChange >= 2))),
                                     log2FoldChange = sample[which(sample$log2FoldChange >= 2),]$log2FoldChange)
  
  tissue.bias.summary <- rbind.data.frame(high.gonad.bias, moderate.gonad.bias, unbiased, moderate.head.bias, 
                                          high.head.bias)
  tissue.bias.summary <- tissue.bias.summary[order(tissue.bias.summary$geneID),]
  write.table(tissue.bias.summary, paste0('Tissue.Biased.Genes.', sex,'.tsv'), quote = F, sep = '\t', col.names = T, row.names = F)
  
  return(tissue.bias.summary)
}

AssignTissueBiasCategory('Females')
AssignTissueBiasCategory('Males')

