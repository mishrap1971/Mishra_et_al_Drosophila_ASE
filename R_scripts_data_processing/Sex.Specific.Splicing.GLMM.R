#function to estimate average of total allelic expression for all 3 replicates of a sample 
#estimate for each sex and tissue
Avg.Total.Allelic.Expression <- function(tissue_name){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Sex.Effects.in.AI")
  r1 <- read.table(paste0(tissue_name, '_replicate_1.tsv'), header = T)
  r2 <- read.table(paste0(tissue_name, '_replicate_2.tsv'), header = T)
  r3 <- read.table(paste0(tissue_name, '_replicate_3.tsv'), header = T)
  
  
  #combine read count data from all replicates
  combined <- cbind.data.frame(r1, r2$DGRP.177, r2$SP.159N, r3$DGRP.177, r3$SP.159N)
  colnames(combined) <- c("geneID", "DGRP.177.1", "SP.159N.1", "DGRP.177.2", "SP.159N.2", "DGRP.177.3", 
                          "SP.159N.3")
  
  #estimate average total read count over 6 replicates [remember to divide by 6, instead of doing rowMeans]
  avg.read.counts <- cbind.data.frame(geneID = r1$geneID, avg.read.counts = rowSums(combined[,c(2:7)])/3) 
  
  #estimate logarithmic value
  log.counts <- cbind.data.frame(geneID = avg.read.counts[,1], 
                                 avg.total.allelic.expression = avg.read.counts[,2],
                                 Log.avg.allelic.expression = log(avg.read.counts[,2]))
  
  return(log.counts)
}


#function to make dataframes with geneID, presence of allele-specific splicing, presence of allelic imbalance and sex bias
Allele.Specific.Splicing.Dataframes <- function(tissue_name, sample_name){
  #Set global variables
  FDRThreshold <- 0.05
  mappedReadsThreshold <- 50
  
  #load junctionSeq output, remove genes that are untested (i.e., marked 'NA')
  sample <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Alt.Splicing.Analyses/JunctionSeq.Outputs/",
                              sample_name, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  sample <- sample[!(is.na(sample[,23])),]
  sample <- sample[!(is.na(sample[,24])),]
  
  #Obtain a list of significant and non-significant genes and output into separate txt files
  #first, remove all genes with fewer than 50 reads mapping to exons in either treatment
  sample <- sample[sample[,23] > mappedReadsThreshold & sample[,24] > mappedReadsThreshold,]
  
  #obtain unique geneIDs with significant allele-specific splicing, along with the adjusted p-values
  sig.genes <- unique(sample[sample$geneWisePadj < FDRThreshold,][,c(2,25)])
  nonsig.genes <- unique(sample[sample$geneWisePadj > FDRThreshold,][,c(2,25)])
  
  #make dataframe of geneIDs and whether or not they have significant allele-specific splicing
  sig.genes.df <- cbind.data.frame(geneID = sig.genes$geneID, p.adj = sig.genes$geneWisePadj,
                                   Allele.Specific.Splicing = rep("Present", nrow(sig.genes)))
  nonsig.genes.df <- cbind.data.frame(geneID = nonsig.genes$geneID, p.adj = nonsig.genes$geneWisePadj, 
                                      Allele.Specific.Splicing = rep("Absent", nrow(nonsig.genes)))
  df.AS <- rbind.data.frame(sig.genes.df, nonsig.genes.df)
  df.AS <- df.AS[order(df.AS$geneID),]
  
  
  #load sex bias data (log2FC rati0) for tissue
  df.sexbias <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs/Sex.Biased.Genes.",
                                  tissue_name, ".tsv"), header = T, sep ='\t')
  
  #load data for presence of AI
  df.AI <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/",
                             tissue_name, ".sexAI.df.tsv"), sep = '\t', header = T)[,c(1,3)]
  
  #merge dataframes
  final.df <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df.AS, df.sexbias, df.AI))
  
  return(final.df)
}

ovaries.df <- Allele.Specific.Splicing.Dataframes('Gonads', 'Ovaries')  
testes.df <- Allele.Specific.Splicing.Dataframes('Gonads', 'Testes')  
female.heads.df <- Allele.Specific.Splicing.Dataframes('Heads', 'Female.heads')  
male.heads.df <- Allele.Specific.Splicing.Dataframes('Heads', 'Male.heads')  
female.whole.bodies.df <- Allele.Specific.Splicing.Dataframes('Whole.bodies', 'Female.whole.body')  
male.whole.bodies.df <- Allele.Specific.Splicing.Dataframes('Whole.bodies', 'Male.whole.body')  
female.reciprocal.df <- Allele.Specific.Splicing.Dataframes('Reciprocal.whole.bodies', 'Female.reciprocal.whole.body')  
male.reciprocal.df <- Allele.Specific.Splicing.Dataframes('Reciprocal.whole.bodies', 'Male.reciprocal.whole.body')  

#function to merge male and female dataframes of the same tissue
#resulting dataframes only include genes that have been tested for allele-specific splicing in both sexes
CombineDataframes <- function(female.df, male.df, tissue_name){
  gene.list <- intersect(ovaries.df$geneID, testes.df$geneID)
  
  female.df$Sex <- rep("Female", nrow(female.df))
  male.df$Sex <- rep("Male", nrow(male.df))
  female.df <- female.df[female.df$geneID %in% gene.list,]
  male.df <- male.df[male.df$geneID %in% gene.list,]
  
  df.combined <- rbind.data.frame(female.df, male.df)
  write.table(df.combined,
              paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Alt.Splicing.Analyses/Dataframes.for.Splicing.Analyses/",
                     tissue_name, ".Splicing.df.tsv"), sep = '\t', quote = F, row.names = F, col.names = T)

  return(df.combined)
}

CombineDataframes(ovaries.df, testes.df, 'Gonads')
CombineDataframes(female.heads.df, male.heads.df, 'Heads')
CombineDataframes(female.whole.bodies.df, male.whole.bodies.df, 'Whole.bodies')
CombineDataframes(female.reciprocal.df, male.reciprocal.df, 'Reciprocal.whole.bodies')

##################################################################################################################
#run linear models 
#presence.of.allele.specific.splicing ~  log.total.allelic.expression + sex*(log2FoldChange + I(log2FoldChange^2))
# + presence.of.AI) + (1|geneID)

require(lme4)
require(lmerTest)

RunMixedEffectsModel <- function(tissue_name){
  df <- read.table(paste0('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Alt.Splicing.Analyses/Dataframes.for.Splicing.Analyses/',
                          tissue_name, ".Splicing.df.tsv"), header = T, sep = '\t')
  
  #upload list of genes without evidence of mapping bias; include only genes that occur in this list
  gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv",
                          header = T)
  df <- na.omit(df[df$geneID %in% gene.list$geneID,])
  df$Allele.Specific.Splicing <- as.factor(df$Allele.Specific.Splicing)
  df$Presence.of.AI <- as.factor(df$Presence.of.AI)
  df$Sex <- as.factor(df$Sex)
  
  #run linear model
  model.df <- glmer(Allele.Specific.Splicing ~ scale(Log.avg.allelic.expression) +
                      Sex*(scale(log2FoldChange) + I(scale(log2FoldChange)^2)) + 
                      Presence.of.AI + (1|geneID), family = "binomial", data = df, 
                    control = glmerControl(optimizer='bobyqa', optCtrl = list(maxfun = 200000)))
  model.sum <- summary(model.df)
  
  return(model.sum)
}

########################
#As a glm
RunMixedEffectsModel <- function(tissue_name){
  df <- read.table(paste0('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Alt.Splicing.Analyses/Dataframes.for.Splicing.Analyses/',
                          tissue_name, ".Splicing.df.tsv"), header = T, sep = '\t')
  
  #upload list of genes without evidence of mapping bias; include only genes that occur in this list
  gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv",
                          header = T)
  df <- na.omit(df[df$geneID %in% gene.list$geneID,])
  df$Allele.Specific.Splicing <- as.factor(df$Allele.Specific.Splicing)
  df$Presence.of.AI <- as.factor(df$Presence.of.AI)
  df$Sex <- as.factor(df$Sex)
  
  #run linear model
  model.df <- glm(Allele.Specific.Splicing ~ Log.avg.allelic.expression +
                    Sex*log2FoldChange + Sex*I(log2FoldChange^2) + 
                    Presence.of.AI, family = "binomial", data = df)
  model.sum <- summary(model.df)
  
  return(model.sum)
}

