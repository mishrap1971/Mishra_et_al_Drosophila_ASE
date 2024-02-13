require('lme4')
require('lmerTest')

#################################################################################################
### STATISTICAL TEST FOR SEX-REVERSED ALLELIC IMBALANCE #########################################
#################################################################################################

#function to convert data to suitable format for binomial/quasibinomial test
convert.data.for.binomial.analysis<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "sex", "sampleID", "cross"))
  return(md)
}

#function to run GLM for each gene, separately in males and females
TestSexesSeparately<-function(dat){
  dat2 <- convert.data.for.binomial.analysis(dat)
  maledat <- dat2[dat2$sex == "M",]
  options(contrasts = rep("contr.sum", 2))
  male.binomialLM<-glm(cbind(countA, countB) ~ 1 , family = "quasibinomial", data = maledat)
  
  
  femaledat <- dat2[dat2$sex == "F",]
  options(contrasts = rep("contr.sum", 2))
  female.binomialLM<-glm(cbind(countA, countB) ~ 1 , family = "quasibinomial", data = femaledat)
  
  male.summary.binomialLM <- summary(male.binomialLM)
  male.coef.binomialLM<- male.summary.binomialLM$coefficients
  female.summary.binomialLM <- summary(female.binomialLM)
  female.coef.binomialLM<- female.summary.binomialLM$coefficients
  
  ans <- as.numeric(c(male.coef.binomialLM[1, c(1,4)], female.coef.binomialLM[1, c(1,4)]))
  return(ans)
}


#arrange read count files for all 12 replicates (male + female replicates in reciprocal and main cross) in a dataframe

#load read count data
setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Sex.Reversed.AI")
f.main.1 <- t(read.table('Female_whole_bodies_replicate_1.tsv', header = T))[2:3,]
f.main.2 <- t(read.table('Female_whole_bodies_replicate_2.tsv', header = T))[2:3,]
f.main.3 <- t(read.table('Female_whole_bodies_replicate_3.tsv', header = T))[2:3,]
f.reciprocal.1 <- t(read.table('Female_whole_bodies_reciprocal_replicate_1.tsv', header = T))[2:3,]
f.reciprocal.2 <- t(read.table('Female_whole_bodies_reciprocal_replicate_2.tsv', header = T))[2:3,]
f.reciprocal.3 <- t(read.table('Female_whole_bodies_reciprocal_replicate_3.tsv', header = T))[2:3,]
m.main.1 <- t(read.table('Male_whole_bodies_replicate_1.tsv', header = T))[2:3,]
m.main.2 <- t(read.table('Male_whole_bodies_replicate_2.tsv', header = T))[2:3,]
m.main.3 <- t(read.table('Male_whole_bodies_replicate_3.tsv', header = T))[2:3,]
m.reciprocal.1 <- t(read.table('Male_whole_bodies_reciprocal_replicate_1.tsv', header = T))[2:3,]
m.reciprocal.2 <- t(read.table('Male_whole_bodies_reciprocal_replicate_2.tsv', header = T))[2:3,]
m.reciprocal.3 <- t(read.table('Male_whole_bodies_reciprocal_replicate_3.tsv', header = T))[2:3,]

#label sample, sex, cross and allele (A = DGRP-177, B = SP-159N)
geneID <- t(read.table('Female_whole_bodies_replicate_1.tsv', header = T))[1,]
sampleID <- c("F1", "F1", "F2", "F2", "F3", "F3", "F4", "F4", "F5", "F5", "F6", "F6", "M1", "M1", "M2", "M2", "M3", "M3",
              "M4", "M4", "M5", "M5", "M6", "M6")
allele <- rep(c("A", "B"), 12)
sex <- c(rep("F", 12), rep("M", 12))
cross <- c(rep("main", 6), rep("reciprocal", 6), rep("main", 6), rep("reciprocal", 6))

##make dataframe with read-count data
df <- data.frame(rbind(f.main.1, f.main.2, f.main.3, f.reciprocal.1, f.reciprocal.2, f.reciprocal.3,
                       m.main.1, m.main.2, m.main.3, m.reciprocal.1, m.reciprocal.2, m.reciprocal.3))
df <- cbind(sampleID, sex, cross, allele, df)
label <- cbind(sampleID, sex, cross, allele)
df[,5:(length(geneID) + 4)] <- sapply(df[,5:(length(geneID) + 4)], as.numeric)

#removing gene where model showed convergence warnings
df <- df[,-2461]  
geneID <- geneID[-2457]
colnames(df) <- c("sampleID", "sex", "cross", "allele", rep("count", length(geneID)))

#run GLM on dataframe
gene.counts <- c()
output.df <- c()
temp.results <- c()

for (i in 5:(length(geneID)+4)){
  gene.counts <- df[,c(1:4,i)]
  temp.results <- TestSexesSeparately(gene.counts)  
  output.df <- cbind(output.df, temp.results)
}

#edit output dataframe to add column names
output.df <- t(output.df)
output.df <- cbind(geneID, output.df)
colnames(output.df) <- c("geneID", "Overabundance_of_DGRP177.males", "p.val.males", "Overabundance_of_DGRP177.females", "p.val.females")

#save dataframe to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Sex.Reversed.AI/"
write.table(output.df, paste0(path, "Sex.Reversed.AI.test.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###############################################################################################################
#obtain list of genes with significant sex-reversed AI
#a gene has sex-reversed AI if there's significant AI in BOTH males and females + direction of AI is opposite in direction
#between the sexes

#load list of genes without evidence of mapping bias
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load GLM output dataframe for sex-reversed AI; filter it to only have genes that pass mapping bias tests
SRAI.df <-read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Sex.Reversed.AI/Sex.Reversed.AI.test.tsv", header=T)
SRAI.df <- SRAI.df[SRAI.df$geneID %in% gene.list$geneID,]

#list of genes with sex-reversed AI; save to file
SRAI.genes <- SRAI.df[which(SRAI.df$Overabundance_of_DGRP177.males*SRAI.df$Overabundance_of_DGRP177.females < 0 &
                              SRAI.df$p.val.males < 0.05 & SRAI.df$p.val.females < 0.05),]
write.table(SRAI.genes,
            "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Sex.Reversed.AI/Genes.with.Sex.Reversed.AI.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

