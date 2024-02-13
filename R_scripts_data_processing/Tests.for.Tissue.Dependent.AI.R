require('lme4')
require('lmerTest')

#################################################################################################
### STATISTICAL TEST FOR TISSUE-DEPENDENT ALLELIC IMBALANCE #####################################
###               WITH SINGLE-SEX DATA                      #####################################
#################################################################################################

#function to convert data to suitable format for binomial/quasibinomial test
convert.data.for.binomial.analysis<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "tissue", "sampleID"))
  return(md)
}

#function to run GLM with quasibinomial distribution for one data point (i.e., one gene)
run.quasibinomial.model <-function(dat){
  dat2<-convert.data.for.binomial.analysis(dat)
  options(contrasts = rep("contr.sum", 2))
  binomialLM<-glm(cbind(countA, countB) ~ tissue , family = "quasibinomial", data = dat2);
  summary.binomialLM<-summary(binomialLM)
  coef.binomialLM <- summary.binomialLM$coefficients
  ans<- as.numeric(c(coef.binomialLM[1, c(1,4)], coef.binomialLM[2, c(1,4)]))
  return(ans)
}

#function to arrange read count files for all 6 replicates (heads + gonads) for a given sex and run GLM on each gene
glm.For.AI <- function(head_tissue, gonad_tissue, output_name){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Tissue.Effects.in.Single.Sex")
  h1 <- t(read.table(paste0(head_tissue, '_replicate_1.tsv')))[2:3,][,-1]
  h2 <- t(read.table(paste0(head_tissue, '_replicate_2.tsv')))[2:3,][,-1]
  h3 <- t(read.table(paste0(head_tissue, '_replicate_3.tsv')))[2:3,][,-1]
  g1 <- t(read.table(paste0(gonad_tissue, '_replicate_1.tsv')))[2:3,][,-1]
  g2 <- t(read.table(paste0(gonad_tissue, '_replicate_2.tsv')))[2:3,][,-1]
  g3 <- t(read.table(paste0(gonad_tissue, '_replicate_3.tsv')))[2:3,][,-1]
  
  #label sample, sex and allele (A = DGRP-177, B = SP-159N)
  geneID <- t(read.table(paste0(head_tissue, '_replicate_1.tsv')))[1,][-1]
  sampleID <- c("H1", "H1", "H2", "H2", "H3", "H3", "G1", "G1", "G2", "G2", "G3", "G3")
  allele <- rep(c("A", "B"), 6)
  tissue <- c(rep("H", 6), rep("G", 6))
  
  #make dataframe with read-count data
  df <- data.frame(rbind(h1, h2, h3, g1, g2, g3))
  df <- cbind(sampleID, tissue, allele, df)
  label <- cbind(sampleID, tissue, allele)
  colnames(df) <- c("sampleID", "tissue", "allele", rep("count", length(geneID)))
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
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.single.sex")
  write.table(output_df, paste0(output_name, ".test.for.tissue.effects.in.AI.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(output_df)
}

glm.For.AI('Female_heads', 'Ovaries', 'Females')
glm.For.AI('Male_heads', 'Testes', 'Males')


#################################################################################################
### STATISTICAL TEST FOR TISSUE-DEPENDENT ALLELIC IMBALANCE #####################################
###           WITH DATA FROM BOTH SEXES                     #####################################
#################################################################################################

convert.data.for.binomial.analysis2<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "sex", "sampleID", "tissue"))
  return(md)
}

run.quasibinomial.model2<-function(dat){
  dat2<-convert.data.for.binomial.analysis2(dat)
  
  #manually change contrast coding from the default ('contr.treatments') to 'contr.sum'
  dat2$sex <- as.factor(dat2$sex)
  dat2$tissue <- as.factor(dat2$tissue)
  dat2$allele.x <- as.factor(dat2$allele.x)
  dat2$allele.y <- as.factor(dat2$allele.y)
  contrasts(dat2$sex) = contr.sum(2)  
  contrasts(dat2$tissue) = contr.sum(2)
  
  binomialLM<-glm(cbind(countA, countB) ~ sex*tissue, family = "quasibinomial", data = dat2);
  summary.binomialLM<-summary(binomialLM)
  coef.binomialLM <- summary.binomialLM$coefficients
  
  ans <- as.numeric(c(coef.binomialLM[1, c(1,4)], coef.binomialLM[2, c(1,4)], coef.binomialLM[3, c(1,4)], 
                      coef.binomialLM[4, c(1,4)]))
  return(ans)
}

#load read count data 
setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Tissue.Effects.in.Both.Sexes')

#read counts for heads
f.heads.1 <- t(read.table('Female_heads_replicate_1.tsv', header = T))[2:3,]
f.heads.2 <- t(read.table('Female_heads_replicate_2.tsv', header = T))[2:3,]
f.heads.3 <- t(read.table('Female_heads_replicate_3.tsv', header = T))[2:3,]
m.heads.1 <- t(read.table('Male_heads_replicate_1.tsv', header = T))[2:3,]
m.heads.2 <- t(read.table('Male_heads_replicate_2.tsv', header = T))[2:3,]
m.heads.3 <- t(read.table('Male_heads_replicate_3.tsv', header = T))[2:3,]

#read counts for gonads
ovaries.1 <- t(read.table('Ovaries_replicate_1.tsv', header = T))[2:3,]
ovaries.2 <- t(read.table('Ovaries_replicate_2.tsv', header = T))[2:3,]
ovaries.3 <- t(read.table('Ovaries_replicate_3.tsv', header = T))[2:3,]
testes.1 <- t(read.table('Testes_replicate_1.tsv', header = T))[2:3,]
testes.2 <- t(read.table('Testes_replicate_2.tsv', header = T))[2:3,]
testes.3 <- t(read.table('Testes_replicate_3.tsv', header = T))[2:3,]

#label sample, sex and allele (A = DGRP-177, B = SP-159N)
geneID <- t(read.table('Female_heads_replicate_1.tsv'))[1,][-1]
sampleID <- c("H1", "H1", "H2", "H2", "H3", "H3", "H4", "H4", "H5", "H5", "H6", "H6", "G1", "G1", "G2", "G2", "G3", "G3",
              "G4", "G4", "G5", "G5", "G6", "G6")
allele <- rep(c("A", "B"), 12)
sex <- c(rep("F", 6), rep("M", 6), rep("F", 6), rep("M", 6))
tissue <- c(rep("Heads", 12), rep("Gonads", 12))

#make dataframe with read-count data
df <- data.frame(rbind(f.heads.1, f.heads.2, f.heads.3, m.heads.1, m.heads.2, m.heads.3, 
                       ovaries.1, ovaries.2, ovaries.3, testes.1, testes.2, testes.3))
df <- cbind(sampleID, sex, tissue, allele, df)
label <- cbind(sampleID, sex, tissue, allele)
df[,5:(length(geneID) + 4)] <- sapply(df[,5:(length(geneID) + 4)], as.numeric)
colnames(df) <- c("sampleID", "sex", "tissue", "allele", rep("count", length(geneID)))

#run GLM on each gene in dataframe
gene.counts <- c()
output.df <- c()
temp.results <- c()

for (i in 5:(length(geneID)+4)){
  gene.counts <- df[,c(1:4,i)]
  temp.results <- run.quasibinomial.model2(gene.counts)  
  output.df <- cbind(output.df, temp.results)
}

#edit output dataframe to add column names
output.df <- t(output.df)
output.df <- cbind(geneID, output.df)
colnames(output.df) <- c("geneID", "Overabundance_of_DGRP177", "p-val_(quasibinomial)_1", "Sex.effect", 
                         "p-val_(quasibinomial)_2", "Tissue.effect", "p-val_(quasibinomial)_3", 
                         "Sex.Tissue.interaction.effect", "p-val_(quasibinomial)_4")

#save dataframe
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.combined.sexes/"
write.table(output.df, paste0(path, 'Combined.Sexes.test.for.tissue.effects.in.AI.tsv') , sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


###############################################################################################################
###  genes with significant tissue-dependent AI ##############################################################
###############################################################################################################
#upload list of genes without evidence of mapping bias; include only genes that occur in this list
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv",
                        header = T)

#upload output from single-sex analyses, filter in genes without evidence of mapping bias
path1 <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.single.sex/"
females.df <- read.table(paste0(path1, 'Females.test.for.tissue.effects.in.AI.tsv'), header = T)
males.df <- read.table(paste0(path1, 'Males.test.for.tissue.effects.in.AI.tsv'), header = T)
females.df <- females.df[females.df$geneID %in% gene.list$geneID,]
males.df <- males.df[males.df$geneID %in% gene.list$geneID,]


#output from analysis with data from both sexes
path2 <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.combined.sexes/"
combined.df <- read.table(paste0(path2, 'Combined.Sexes.test.for.tissue.effects.in.AI.tsv'), header = T)
combined.df <- combined.df[combined.df$geneID %in% gene.list$geneID,]

#percentage of genes with significant AI 
#in single-sex analysis, AI is significant if intercept term has p-val < 0.05)
#in combined-sex analysis, AI is significant if intercept has p-val < 0.05)
length(which(females.df$p.val_.quasibinomial._1 < 0.05))/nrow(females.df)
length(which(males.df$p.val_.quasibinomial._1 < 0.05))/nrow(males.df)
length(which(combined.df$p.val_.quasibinomial._1 < 0.05))/nrow(combined.df)


#percentage of genes with significant tissue-dependent AI 
#in single-sex analysis, AI is significant if tissue term has p-val < 0.05)
#in combined-sex analysis, AI is significant if tissue term has p-val < 0.05)
length(which(females.df$p.val_.quasibinomial._2 < 0.05))/nrow(females.df)
length(which(males.df$p.val_.quasibinomial._2 < 0.05))/nrow(males.df)
length(which(combined.df$p.val_.quasibinomial._3 < 0.05))/nrow(combined.df)

#percentage of genes with significant tissue*sex interaction effect on AI 
length(which(combined.df$p.val_.quasibinomial._4 < 0.05))/nrow(combined.df)
