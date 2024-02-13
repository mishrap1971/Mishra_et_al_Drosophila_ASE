require('lme4')
require('lmerTest')

#################################################################################################
### STATISTICAL TEST FOR SEX-REVERSED ALLELIC IMBALANCE #########################################
#################################################################################################

#function to convert data to suitable format for binomial/quasibinomial test
convert.data.for.binomial.analysis<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c("tissue", "sampleID"))
  return(md)
}

#function to run GLM for each gene, separately in gonads and heads
TestTissuesSeparately<-function(dat){
  dat2 <- convert.data.for.binomial.analysis(dat)
  gonad.dat <- dat2[dat2$tissue == "G",]
  head.dat <- dat2[dat2$tissue == "H",]
  
  options(contrasts = rep("contr.sum", 2))
  gonad.binomialLM<-glm(cbind(countA, countB) ~ 1 , family = "quasibinomial", data = gonad.dat);
  
  options(contrasts = rep("contr.sum", 2))
  head.binomialLM<-glm(cbind(countA, countB) ~ 1 , family = "quasibinomial", data = head.dat);
  
  gonad.summary.binomialLM <- summary(gonad.binomialLM)
  gonad.coef.binomialLM<- gonad.summary.binomialLM$coefficients
  head.summary.binomialLM <- summary(head.binomialLM)
  head.coef.binomialLM<- head.summary.binomialLM$coefficients
  
  ans <- as.numeric(c(gonad.coef.binomialLM[1, c(1,4)], head.coef.binomialLM[1, c(1,4)]))
  return(ans)
}

################################################################################
## Testing for tissue-reversed AI in females ###################################
################################################################################
#arrange read count files for all 6 replicates (3 head + 3 gonad replicates) in a dataframe
setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Tissue.Effects.in.Single.Sex")
ovaries.1 <- t(read.table('Ovaries_replicate_1.tsv'))[2:3,][,-1]
ovaries.2 <- t(read.table('Ovaries_replicate_2.tsv'))[2:3,][,-1]
ovaries.3 <- t(read.table('Ovaries_replicate_3.tsv'))[2:3,][,-1]
female.heads.1 <- t(read.table('Female_heads_replicate_1.tsv'))[2:3,][,-1]
female.heads.2 <- t(read.table('Female_heads_replicate_2.tsv'))[2:3,][,-1]
female.heads.3 <- t(read.table('Female_heads_replicate_3.tsv'))[2:3,][,-1]

#label sample, sex, cross and allele (A = DGRP-177, B = SP-159N)
geneID <- t(read.table('Ovaries_replicate_1.tsv'))[1,][-1]
sampleID <- c("G1", "G1", "G2", "G2", "G3", "G3", "H1", "H1", "H2", "H2", "H3", "H3")
allele <- rep(c("A", "B"), 6)
tissue <- c(rep("G", 6), rep("H", 6))

##make dataframe with read-count data
df <- data.frame(rbind(ovaries.1, ovaries.2, ovaries.3, female.heads.1, female.heads.2, female.heads.3))
df <- cbind(sampleID, tissue, allele, df)
label <- cbind(sampleID, tissue, allele)
colnames(df) <- c("sampleID", "tissue", "allele", rep("count", length(geneID)))
df[,4:(length(geneID) + 3)] <- sapply(df[,4:(length(geneID) + 3)], as.numeric)

#run GLM on dataframe
gene.counts <- c()
output.df <- c()
temp.results <- c()

for (i in 4:(length(geneID)+3)){
  gene.counts <- df[,c(1:3,i)]
  temp.results <- TestTissuesSeparately(gene.counts)  
  output.df <- cbind(output.df, temp.results)
}

#edit output dataframe to add column names
output.df <- t(output.df)
output.df <- cbind(geneID, output.df)
colnames(output.df) <- c("geneID", "Overabundance_of_DGRP177.gonads", "p.val.gonads", "Overabundance_of_DGRP177.heads", "p.val.heads")

#save dataframe to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Tissue.Reversed.AI/"
write.table(output.df, paste0(path, "Females.Tissue.Reversed.AI.test.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


################################################################################
## Testing for tissue-reversed AI in males ###################################
################################################################################
#arrange read count files for all 6 replicates (3 head + 3 gonad replicates) in a dataframe
setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Tissue.Effects.in.Single.Sex")
testes.1 <- t(read.table('Testes_replicate_1.tsv'))[2:3,][,-1]
testes.2 <- t(read.table('Testes_replicate_2.tsv'))[2:3,][,-1]
testes.3 <- t(read.table('Testes_replicate_3.tsv'))[2:3,][,-1]
male.heads.1 <- t(read.table('male_heads_replicate_1.tsv'))[2:3,][,-1]
male.heads.2 <- t(read.table('male_heads_replicate_2.tsv'))[2:3,][,-1]
male.heads.3 <- t(read.table('male_heads_replicate_3.tsv'))[2:3,][,-1]

#label sample, sex, cross and allele (A = DGRP-177, B = SP-159N)
geneID <- t(read.table('Testes_replicate_1.tsv'))[1,][-1]
sampleID <- c("G1", "G1", "G2", "G2", "G3", "G3", "H1", "H1", "H2", "H2", "H3", "H3")
allele <- rep(c("A", "B"), 6)
tissue <- c(rep("G", 6), rep("H", 6))

##make dataframe with read-count data
df <- data.frame(rbind(testes.1, testes.2, testes.3, male.heads.1, male.heads.2, male.heads.3))
df <- cbind(sampleID, tissue, allele, df)
label <- cbind(sampleID, tissue, allele)
colnames(df) <- c("sampleID", "tissue", "allele", rep("count", length(geneID)))
df[,4:(length(geneID) + 3)] <- sapply(df[,4:(length(geneID) + 3)], as.numeric)

#run GLM on dataframe
gene.counts <- c()
output.df <- c()
temp.results <- c()

for (i in 4:(length(geneID)+3)){
  gene.counts <- df[,c(1:3,i)]
  temp.results <- TestTissuesSeparately(gene.counts)  
  output.df <- cbind(output.df, temp.results)
}

#edit output dataframe to add column names
output.df <- t(output.df)
output.df <- cbind(geneID, output.df)
colnames(output.df) <- c("geneID", "Overabundance_of_DGRP177.gonads", "p.val.gonads", "Overabundance_of_DGRP177.heads", "p.val.heads")

#save dataframe to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Tissue.Reversed.AI/"
write.table(output.df, paste0(path, "Males.Tissue.Reversed.AI.test.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


###############################################################################################################
#obtain list of genes with significant sex-reversed AI
#a gene has tissue-reversed AI if there's significant AI in BOTH heads and gonads + direction of AI is opposite in 
#direction between the sexes

#load list of genes without evidence of mapping bias
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load GLM output dataframe for sex-reversed AI; filter it to only have genes that pass mapping bias tests
females.TRAI.df <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Tissue.Reversed.AI/Females.Tissue.Reversed.AI.test.tsv", header=T)
males.TRAI.df <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Tissue.Reversed.AI/Males.Tissue.Reversed.AI.test.tsv", header=T)
females.TRAI.df <- females.TRAI.df[females.TRAI.df$geneID %in% gene.list$geneID,]
males.TRAI.df <- males.TRAI.df[males.TRAI.df$geneID %in% gene.list$geneID,]

#list of genes with tissue-reversed AI; save to file
females.TRAI.genes <- females.TRAI.df[which(females.TRAI.df$Overabundance_of_DGRP177.gonads*females.TRAI.df$Overabundance_of_DGRP177.heads < 0 &
                                              females.TRAI.df$p.val.gonads < 0.05 & females.TRAI.df$p.val.heads < 0.05),]
males.TRAI.genes <- males.TRAI.df[which(males.TRAI.df$Overabundance_of_DGRP177.gonads*males.TRAI.df$Overabundance_of_DGRP177.heads < 0 &
                                          males.TRAI.df$p.val.gonads < 0.05 & males.TRAI.df$p.val.heads < 0.05),]

