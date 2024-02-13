library(DESeq2)

#function to run DESeq2 for a given tissue type
SexesDiffExpression <- function(tissue, output_name){
  setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.DESeq2')
  filenames <- list.files(pattern = tissue)
  samplenames <- gsub('.{4}$', '', filenames)
  sample.sex <- factor(c(rep("Female", 3), rep("Male", 3)))
  sample.table <- data.frame(sampleName = samplenames, fileName = filenames, sex = sample.sex)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table, design= ~sex)
  keep.dds = rowMeans(counts(ddsHTSeq)) >= 50
  ddsHTSeq = ddsHTSeq[keep.dds,]
  DESeq.df = DESeq(ddsHTSeq)
  results = results(DESeq.df, contrast=c("sex", "Male", "Female"))
  results$FlyBaseID = rownames(results)
  results = results[,c(7,1:6)]
  
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs")
  write.table(results, file = paste0('DifferentialGeneExpression.', output_name, ".txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  return(results)
}

# SexesDiffExpression("gonads")
# SexesDiffExpression("heads")
# SexesDiffExpression("whole_bodies_rec")
# SexesDiffExpression("whole_bodies_rep")

# function to assign sex-bias category to genes in the DESeq2 output
AssignSexBiasCategory <- function(tissue){
  
  #load deseq2 output; remove genes with high standard error in expression
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs")
  sample <- read.table(paste0('DifferentialGeneExpression.', tissue , '.tsv'), header = T, sep = '\t')
  high.SE.genes <- sample[which(sample$lfcSE > 1.5),]$FlyBaseID
  sample <- sample[!(sample$FlyBaseID %in% high.SE.genes),]
  
  #assign sex bias category: highly FB genes (-infinity, -2), moderately FB genes (-2, -0.5),
  #unbiased (-0.5, 0.5), moderately MB genes (0.5, 2), highly MB genes (2, infinity)
  high.female.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < -2),]$FlyBaseID, 
                      Sex.Bias.Category = rep("Highly female-biased", length(which(sample$log2FoldChange < -2))),
                      log2FoldChange = sample[which(sample$log2FoldChange < -2),]$log2FoldChange)
  moderate.female.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2),]$FlyBaseID, 
                          Sex.Bias.Category = rep("Moderately female-biased", 
                          length(which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2))),
                          log2FoldChange = sample[which(sample$log2FoldChange < -0.5 & sample$log2FoldChange >= -2),]$log2FoldChange)
  unbiased <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5),]$FlyBaseID, 
              Sex.Bias.Category = rep("Unbiased", 
              length(which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5))),
              log2FoldChange = sample[which(sample$log2FoldChange < 0.5 & sample$log2FoldChange >= -0.5),]$log2FoldChange)
  moderate.male.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5),]$FlyBaseID, 
                        Sex.Bias.Category = rep("Moderately male-biased", 
                        length(which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5))),
                        log2FoldChange = sample[which(sample$log2FoldChange < 2 & sample$log2FoldChange >= 0.5),]$log2FoldChange)
  high.male.bias <- cbind.data.frame(geneID = sample[which(sample$log2FoldChange >= 2),]$FlyBaseID, 
                    Sex.Bias.Category = rep("Highly male-biased", length(which(sample$log2FoldChange >= 2))),
                    log2FoldChange = sample[which(sample$log2FoldChange >= 2),]$log2FoldChange)
  
  sex.bias.summary <- rbind.data.frame(high.female.bias, moderate.female.bias, unbiased, moderate.male.bias, 
                                       high.male.bias)
  sex.bias.summary <- sex.bias.summary[order(sex.bias.summary$geneID),]
  write.table(sex.bias.summary, paste0('Sex.Biased.Genes.', tissue,'.tsv'), quote = F, sep = '\t', col.names = T, row.names = F)
  
  return(sex.bias.summary)
}

AssignSexBiasCategory('heads')
AssignSexBiasCategory('gonads')
AssignSexBiasCategory('reciprocal.whole.bodies')
AssignSexBiasCategory('whole.bodies')

