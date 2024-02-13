require('lme4')
require('lmerTest')

#################################################################################################
### STATISTICAL TEST FOR ALLELIC IMBALANCE ######################################################
#################################################################################################

#function to convert data to suitable format for binomial/quasibinomial test
convert.data.for.binomial.analysis<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "sex", "sampleID"))
  return(md)
}

#function to run GLM with quasibinomial distribution for one data point (i.e., one gene)
run.quasibinomial.model <-function(dat){
  dat2<-convert.data.for.binomial.analysis(dat)
  options(contrasts = rep("contr.sum", 2))
  binomialLM<-glm(cbind(countA, countB) ~ sex , family = "quasibinomial", data = dat2);
  summary.binomialLM<-summary(binomialLM)
  coef.binomialLM <- summary.binomialLM$coefficients
  ans<- as.numeric(c(coef.binomialLM[1, c(1,4)], coef.binomialLM[2, c(1,4)]))
  return(ans)
}

#function to arrange read count files for all 6 replicates (male + female) and run GLM on each gene
glm.For.AI <- function(female_tissue, male_tissue, output_name){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Sex.Effects.in.AI")
  f1 <- t(read.table(paste0(female_tissue, '_replicate_1.tsv')))[2:3,][,-1]
  f2 <- t(read.table(paste0(female_tissue, '_replicate_2.tsv')))[2:3,][,-1]
  f3 <- t(read.table(paste0(female_tissue, '_replicate_3.tsv')))[2:3,][,-1]
  m1 <- t(read.table(paste0(male_tissue, '_replicate_1.tsv')))[2:3,][,-1]
  m2 <- t(read.table(paste0(male_tissue, '_replicate_2.tsv')))[2:3,][,-1]
  m3 <- t(read.table(paste0(male_tissue, '_replicate_3.tsv')))[2:3,][,-1]
  
  #label sample, sex and allele (A = DGRP-177, B = SP-159N)
  geneID <- t(read.table(paste0(female_tissue, '_replicate_1.tsv')))[1,][-1]
  sampleID <- c("F1", "F1", "F2", "F2", "F3", "F3", "M1", "M1", "M2", "M2", "M3", "M3")
  allele <- rep(c("A", "B"), 6)
  sex <- c(rep("F", 6), rep("M", 6))
  
  #make dataframe with read-count data
  df <- data.frame(rbind(f1, f2, f3, m1, m2, m3))
  df <- cbind(sampleID, sex, allele, df)
  label <- cbind(sampleID, sex, allele)
  colnames(df) <- c("sampleID", "sex", "allele", rep("count", length(geneID)))
  df[,4:(length(geneID)+3)] <- sapply(df[, 4:(length(geneID)+3)], as.numeric)
  
  #run GLM on dataframe
  gene_counts <- c()
  output_df <- c()
  temp_results <- c()
  
  for (i in 4:(length(geneID)+3)){
    gene_counts <- df[,c(1:3,i)]
    temp_results <- run.quasibinomial.model(gene_counts)  
    output_df <- cbind(output_df, temp_results)
  }
  
  #edit output dataframe to add column names
  output_df <- t(output_df)
  output_df <- cbind(geneID, output_df)
  colnames(output_df) <- c("geneID", "Overabundance_of_DGRP177","p-val_(quasibinomial)_1","Sex-difference", "p-val_(quasibinomial)_2")
  
  #save dataframe to file
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Test.for.AI.outputs")
  write.table(output_df, paste0(output_name, ".test.for.AI.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(output_df)
}

glm.For.AI('Ovaries', 'Testes', 'Gonads')
glm.For.AI('Female_heads', 'Male_heads', 'Heads')
glm.For.AI('Female_whole_bodies', 'Male_whole_bodies', 'Whole_body')
glm.For.AI('Female_whole_bodies_reciprocal', 'Male_whole_bodies_reciprocal', 'Reciprocal_whole_body')




