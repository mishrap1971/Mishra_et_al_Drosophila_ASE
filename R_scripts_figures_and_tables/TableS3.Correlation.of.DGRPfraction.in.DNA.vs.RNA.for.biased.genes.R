#load list of genes found to have evidence of mapping bias (either using genomic data or by analysing parental effects)
gene.list <- read.table('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.with.mapping.bias.or.parental.effects.tsv',header=T)

#load genomic read counts for F1 hybrid males and females
#retain genes found to have mapping bias
m.genomic <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/Genomic.read.counts/DGRP177.SP159N.F1.male.tsv", header=T)
f.genomic <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/Genomic.read.counts/DGRP177.SP159N.F1.female.tsv", header=T)
m.genomic <- m.genomic[m.genomic$geneID %in% gene.list$geneID,]
f.genomic <- f.genomic[f.genomic$geneID %in% gene.list$geneID,]

#load AI dataframes containing DGRP177 expression levels estimated from GLM terms
#retain genes found to have mapping bias
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/"
gonads.df <- read.table(paste0(path, "Gonads.sexAI.df.tsv"), header = T,sep = '\t')
heads.df <- read.table(paste0(path, "Heads.sexAI.df.tsv"), header = T,sep = '\t')
whole.bodies.df <- read.table(paste0(path, "Whole.bodies.sexAI.df.tsv"), header = T,sep = '\t')
reciprocal.df <- read.table(paste0(path, "Reciprocal.whole.bodies.sexAI.df.tsv"), header = T,sep = '\t')

gonads.df <- gonads.df[gonads.df$geneID %in% gene.list$geneID,]
heads.df <- heads.df[heads.df$geneID %in% gene.list$geneID,]
whole.bodies.df <- whole.bodies.df[whole.bodies.df$geneID %in% gene.list$geneID,]
reciprocal.df <- reciprocal.df[reciprocal.df$geneID %in% gene.list$geneID,]

#function to obtain correlation between DGRP177 expression and DGRP177 genomic coverage in both sexes for a given 
#tissue for all genes
CorrelationTestFemalesI <- function(tissue_dataframe){
  f.common.genes <- intersect(f.genomic$geneID,tissue_dataframe$geneID)
  f.genomic.new <- f.genomic[f.genomic$geneID %in% f.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% f.common.genes,]
  cor.female <- cor.test(f.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.female)
}

CorrelationTestMalesI <- function(tissue_dataframe){
  m.common.genes <- intersect(m.genomic$geneID,tissue_dataframe$geneID)
  m.genomic.new <- m.genomic[m.genomic$geneID %in% m.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% m.common.genes,]
  cor.male <- cor.test(m.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.male)
}

CorrelationTestFemalesI(gonads.df)
CorrelationTestFemalesI(heads.df)
CorrelationTestFemalesI(whole.bodies.df)
CorrelationTestFemalesI(reciprocal.df)

CorrelationTestMalesI(gonads.df)
CorrelationTestMalesI(heads.df)
CorrelationTestMalesI(whole.bodies.df)
CorrelationTestMalesI(reciprocal.df)


#function to obtain correlation between DGRP177 expression and DGRP177 genomic coverage in both sexes for a given 
#tissue for genes that would meet the criteria for AI
CorrelationTestFemalesII <- function(tissue_dataframe){
  tissue_dataframe <- tissue_dataframe[which(tissue_dataframe$Presence.of.AI == 1),]
  f.common.genes <- intersect(f.genomic$geneID,tissue_dataframe$geneID)
  f.genomic.new <- f.genomic[f.genomic$geneID %in% f.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% f.common.genes,]
  cor.female <- cor.test(f.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.female)
}

CorrelationTestMalesII <- function(tissue_dataframe){
  tissue_dataframe <- tissue_dataframe[which(tissue_dataframe$Presence.of.AI == 1),]
  m.common.genes <- intersect(m.genomic$geneID,tissue_dataframe$geneID)
  m.genomic.new <- m.genomic[m.genomic$geneID %in% m.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% m.common.genes,]
  cor.male <- cor.test(m.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.male)
}

CorrelationTestFemalesII(gonads.df)
CorrelationTestFemalesII(heads.df)
CorrelationTestFemalesII(whole.bodies.df)
CorrelationTestFemalesII(reciprocal.df)

CorrelationTestMalesII(gonads.df)
CorrelationTestMalesII(heads.df)
CorrelationTestMalesII(whole.bodies.df)
CorrelationTestMalesII(reciprocal.df)


#function to obtain correlation between DGRP177 expression and DGRP177 genomic coverage in both sexes for a given 
#tissue for genes would not meet the criteria for AI
CorrelationTestFemalesIII <- function(tissue_dataframe){
  tissue_dataframe <- tissue_dataframe[which(tissue_dataframe$Presence.of.AI == 0),]
  f.common.genes <- intersect(f.genomic$geneID,tissue_dataframe$geneID)
  f.genomic.new <- f.genomic[f.genomic$geneID %in% f.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% f.common.genes,]
  cor.female <- cor.test(f.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.female)
}

CorrelationTestMalesIII <- function(tissue_dataframe){
  tissue_dataframe <- tissue_dataframe[which(tissue_dataframe$Presence.of.AI == 0),]
  m.common.genes <- intersect(m.genomic$geneID,tissue_dataframe$geneID)
  m.genomic.new <- m.genomic[m.genomic$geneID %in% m.common.genes,]
  tissue.df <- tissue_dataframe[tissue_dataframe$geneID %in% m.common.genes,]
  cor.male <- cor.test(m.genomic.new$Fraction.DGRP, tissue.df$Proportion.of.DGRP177.reads.Females)
  
  return(cor.male)
}

CorrelationTestFemalesIII(gonads.df)
CorrelationTestFemalesIII(heads.df)
CorrelationTestFemalesIII(whole.bodies.df)
CorrelationTestFemalesIII(reciprocal.df)

CorrelationTestMalesIII(gonads.df)
CorrelationTestMalesIII(heads.df)
CorrelationTestMalesIII(whole.bodies.df)
CorrelationTestMalesIII(reciprocal.df)
