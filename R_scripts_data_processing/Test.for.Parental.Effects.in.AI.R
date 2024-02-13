require('lme4')
require('lmerTest')

#################################################################################################
### STATISTICAL TEST FOR PARENTAL EFFECTS IN       #############################################
### ALLELIC IMBALANCE FOR SEPARATE SEXES           #############################################
#################################################################################################
#function to convert data to suitable format for binomial/quasibinomial test
convert.data.for.binomial.analysis<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "cross", "sampleID"))
  return(md)
}

#function to run GLM with quasibinomial distribution for one data point (i.e., one gene)
run.quasibinomial.model <-function(dat){
  dat2<-convert.data.for.binomial.analysis(dat)
  options(contrasts = rep("contr.sum", 2))
  binomialLM<-glm(cbind(countA, countB) ~ cross , family = "quasibinomial", data = dat2);
  summary.binomialLM<-summary(binomialLM)
  coef.binomialLM <- summary.binomialLM$coefficients
  ans<- as.numeric(c(coef.binomialLM[1, c(1,4)], coef.binomialLM[2, c(1,4)]))
  return(ans)
}

#function to arrange read count files for all 6 replicates (male + female) and run GLM on each gene
glm.for.ParentalEffects.on.AI <- function(sex, output_name){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Parental.Effects.in.AI.in.Single.Sex")
  m1 <- t(read.table(paste0(sex, '_whole_bodies_replicate_1.tsv')))[2:3,][,-1]
  m2 <- t(read.table(paste0(sex, '_whole_bodies_replicate_2.tsv')))[2:3,][,-1]
  m3 <- t(read.table(paste0(sex, '_whole_bodies_replicate_3.tsv')))[2:3,][,-1]
  r1 <- t(read.table(paste0(sex, '_whole_bodies_reciprocal_replicate_1.tsv')))[2:3,][,-1]
  r2 <- t(read.table(paste0(sex, '_whole_bodies_reciprocal_replicate_2.tsv')))[2:3,][,-1]
  r3 <- t(read.table(paste0(sex, '_whole_bodies_reciprocal_replicate_3.tsv')))[2:3,][,-1]
  
  #label sample, sex and allele (A = DGRP-177, B = SP-159N)
  geneID <- t(read.table(paste0(sex, '_whole_bodies_replicate_1.tsv')))[1,][-1]
  sampleID <- c("M1", "M1", "M2", "M2", "M3", "M3", "R1", "R1", "R2", "R2", "R3", "R3")
  allele <- rep(c("A", "B"), 6)
  cross <- c(rep("M", 6), rep("R", 6))
  
  #make dataframe with read-count data
  df <- data.frame(rbind(m1, m2, m3, r1, r2, r3))
  df <- cbind(sampleID, cross, allele, df)
  label <- cbind(sampleID, cross, allele)
  colnames(df) <- c("sampleID", "cross", "allele", rep("count", length(geneID)))
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
  colnames(output_df) <- c("geneID", "Overabundance_of_DGRP177","p-val_(quasibinomial)_1","Cross-difference", "p-val_(quasibinomial)_2")
  
  #save dataframe to file
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Parental.Effects.on.AI.in.separate.sexes")
  write.table(output_df, paste0(sex, "s.test.for.parentalEffects.in.AI.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(output_df)
}

glm.for.ParentalEffects.on.AI("Female")
glm.for.ParentalEffects.on.AI("Male")

#################################################################################################
### STATISTICAL TEST FOR PARENTAL EFFECTS IN       #############################################
### ALLELIC IMBALANCE WITH DATA FROM BOTH SEXES    ############################################
#################################################################################################
require('lme4')
require('lmerTest')

#functions
convert.data.for.binomial.analysis2<-function(data){
  dA<-data[data$allele=="A",]; names(dA)[names(dA) == "count"] <- "countA"
  dB<-data[data$allele=="B",]; names(dB)[names(dB) == "count"] <- "countB"
  md<-merge(dA, dB, by = c( "sex", "sampleID", "cross"))
  return(md)
}

run.quasibinomial.model2<-function(dat){
  dat2<-convert.data.for.binomial.analysis2(dat)
  dat2$sex <- as.factor(dat2$sex)
  dat2$cross <- as.factor(dat2$cross)
  dat2$allele.x <- as.factor(dat2$allele.x)
  dat2$allele.y <- as.factor(dat2$allele.y)
  contrasts(dat2$sex) = contr.sum(2)  
  contrasts(dat2$cross) = contr.sum(2)
  binomialLM<-glm(cbind(countA, countB) ~ sex + cross, family = "quasibinomial", data = dat2);
  summary.binomialLM<-summary(binomialLM)
  coef.binomialLM <- summary.binomialLM$coefficients
  
  ans<- as.numeric(c(coef.binomialLM[1, c(1,4)], coef.binomialLM[2, c(1,4)], coef.binomialLM[3, c(1,4)]))
  return(ans)
}

#load read count data 
setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Parental.Effects.in.AI.in.Both.Sexes')

#females
female.main.1 <- t(read.table('Female_whole_bodies_replicate_1.tsv', header = T))[2:3,]
female.main.2 <- t(read.table('Female_whole_bodies_replicate_2.tsv', header = T))[2:3,]
female.main.3 <- t(read.table('Female_whole_bodies_replicate_3.tsv', header = T))[2:3,]
female.reciprocal.1 <- t(read.table('Female_whole_bodies_reciprocal_replicate_1.tsv', header = T))[2:3,]
female.reciprocal.2 <- t(read.table('Female_whole_bodies_reciprocal_replicate_2.tsv', header = T))[2:3,]
female.reciprocal.3 <- t(read.table('Female_whole_bodies_reciprocal_replicate_3.tsv', header = T))[2:3,]

#males
male.main.1 <- t(read.table('Male_whole_bodies_replicate_1.tsv', header = T))[2:3,]
male.main.2 <- t(read.table('Male_whole_bodies_replicate_2.tsv', header = T))[2:3,]
male.main.3 <- t(read.table('Male_whole_bodies_replicate_3.tsv', header = T))[2:3,]
male.reciprocal.1 <- t(read.table('Male_whole_bodies_reciprocal_replicate_1.tsv', header = T))[2:3,]
male.reciprocal.2 <- t(read.table('Male_whole_bodies_reciprocal_replicate_2.tsv', header = T))[2:3,]
male.reciprocal.3 <- t(read.table('Male_whole_bodies_reciprocal_replicate_3.tsv', header = T))[2:3,]

#label sample, sex and allele (A = DGRP-177, B = SP-159N)
geneID <- t(read.table('Female_whole_bodies_replicate_1.tsv'))[1,][-1]
sampleID <- c("F1", "F1", "F2", "F2", "F3", "F3", "F4", "F4", "F5", "F5", "F6", "F6", "M1", "M1", "M2", "M2", 
              "M3", "M3", "M4", "M4", "M5", "M5", "M6", "M6")
allele <- rep(c("A", "B"), 12)
sex <- c(rep("F", 12), rep("M", 12))
cross <- c(rep("main", 6), rep("reciprocal", 6), rep("main", 6), rep("reciprocal", 6))

#make dataframe with read-count data
df <- data.frame(rbind(female.main.1, female.main.2, female.main.3, female.reciprocal.1, 
                       female.reciprocal.2, female.reciprocal.3, male.main.1, male.main.2, 
                       male.main.3, male.reciprocal.1, male.reciprocal.2, male.reciprocal.3))
df <- cbind(sampleID, sex, cross, allele, df)
label <- cbind(sampleID, sex, cross, allele)
df[,5:(length(geneID) + 4)] <- sapply(df[,5:(length(geneID) + 4)], as.numeric)

#removing rows that showed warning when running GLM through the for loop below (GLM didn't converge for these genes)
df <- df[,-c(308, 1565,2212, 2235, 2439, 2473, 7280)]  
geneID <- geneID[-c(304, 1561, 2208, 2231, 2435, 2469, 7276)]  
colnames(df) <- c("sampleID", "sex", "cross", "allele", rep("count", length(geneID)))

#run GLM on dataframe
gene.counts <- c()
output.df <- c()
temp.results <- c()

for (i in 5:(length(geneID)+4)){
  gene.counts <- df[,c(1:4,i)]
  temp.results <- run.quasibinomial.model(gene.counts)  
  output.df <- cbind(output.df, temp.results)
}

#edit output dataframe to add column names
output.df <- t(output.df)
output.df <- cbind(geneID, output.df)
colnames(output.df) <- c("geneID", "Overabundance_of_DGRP177", "p-val_(quasibinomial)_1", "Sex.effect", "p-val_(quasibinomial)_2",
                         "Cross.effect", "p-val_(quasibinomial)_3")

#save dataframe to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Parental.Effects.on.AI.with.both.sexes/"
write.table(output.df, paste0(path, "Test.for.parentalEffects.in.AI.combined.sexes.tsv"), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


##########################################################################
### GENES WITH SIG. PARENTAL EFFECTS #####################################
##########################################################################
require(dplyr)
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/"

#genes with significant parental effects are likely to have mapping bias - we are going to only analyse
#genes that do not have a sig. parental effect in the single-sex AND combined-sexes tests
#load outputs from GLMs run to test for parental effects on single-sex AND combined-sex data
PE.males <- read.table(paste0(path,
            "Parental.Effects.on.AI.in.separate.sexes/Males.test.for.parentalEffects.in.AI.tsv"), header = T)
PE.females <- read.table(paste0(path,
            "Parental.Effects.on.AI.in.separate.sexes/Females.test.for.parentalEffects.in.AI.tsv"), header = T)
PE.combined <- read.table(paste0(path,
          "Parental.Effects.on.AI.with.both.sexes/Test.for.parentalEffects.in.AI.combined.sexes.tsv"), header = T)


#REMEMBER to remove sex-chromosome genes (both X and Y) from the male 
#Do not remove sex-chromosome genes from combined-sexes analysis; we do not want to exclude X-genes with possible mapping bias
x.genes <- read.table(paste0(path, "X.chromosome.genes.tsv"), header = T)
y.genes <- read.table(paste0(path, "Y.chromosome.genes.tsv"), header = T)
sex.genes <- union(x.genes$geneID, y.genes$geneID)
PE.males <- PE.males[!PE.males$geneID %in% sex.genes,]
#PE.combined <- PE.combined[!PE.combined$geneID %in% sex.genes,]

#function to mark 'PRESENT' against genes that have a significant parental effect in the tests
#presence.of.PE = PRESENT means there is parental effect in the gene
PresenceOfParentalEffectsOnAI <- function(geneID, pval){
  presence.of.PE <- c()
  for (i in 1:length(geneID)){
    if (pval[i] < 0.05){
      presence.of.PE[i] = "Present"
    } else {
      presence.of.PE[i] = "Absent"
    }
  }
  df <- cbind.data.frame(geneID = geneID, presence.of.PE = presence.of.PE)
  return(df)
}

#run above function on male, female and combined test outputs, merge resulting dataframes with full_join
#p-values should correspond to the "cross direction" term
PE.f <- PresenceOfParentalEffectsOnAI(PE.females$geneID, PE.females$p.val_.quasibinomial._2)
PE.m <- PresenceOfParentalEffectsOnAI(PE.males$geneID, PE.males$p.val_.quasibinomial._2)
PE.c <- PresenceOfParentalEffectsOnAI(PE.combined$geneID, PE.combined$p.val_.quasibinomial._3) 

PE.mf <- full_join(PE.m, PE.f, by = "geneID")
PE.all <- full_join(PE.mf, PE.c, by = "geneID")
colnames(PE.all) <- c("geneID", "Presence.of.PE.males", "Presence.of.PE.females", "Presence.of.PE.combined.sexes")


#format the dataframe, save to file
PE.all$Should.gene.be.included <- rep("X", nrow(PE.all))  ##"X" is just a filler string, to be replaced by "TRUE" or "FALSE
PE.all$Should.gene.be.included[PE.all$Presence.of.PE.males == "Present" | 
              PE.all$Presence.of.PE.females == "Present" | PE.all$Presence.of.PE.combined.sexes == "Present"] <- "FALSE"
PE.all$Should.gene.be.included[PE.all$Should.gene.be.included == 'X'] <- "TRUE"

write.table(PE.all, paste0(path, "Parental.Effects.summary.dataframe3.tsv"), sep = '\t', row.names = F,
            col.names = T, quote = F)

###############################################################################################
### Genes to be included in further analysis ##################################################
###############################################################################################
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/"

#genes with mapping bias detected from genomic data
df.genomic <- read.table(paste0(path,
                        "Mapping.Bias.with.F1.DNA.counts/Tests.for.Bias.in.F1.genomic.data.tsv"), header = T)
genes.to.be.included1 <- df.genomic[which(df.genomic$crit1and2 == 1),]
genes.to.be.excluded1 <- df.genomic[which(df.genomic$crit1and2 == 0),]

#genes with mapping bias detected using parental effects analysis 
df.PE <- read.table(paste0(path, 'Parental.Effects.summary.dataframe.tsv'), header = T)
genes.to.be.included2 <- df.PE[which(df.PE$Should.gene.be.included == "TRUE"),]
genes.to.be.excluded2 <- df.PE[which(df.PE$Should.gene.be.included == "FALSE"),]

#final list of genes to be included in further analysis
genes.to.be.included.3 <- union(genes.to.be.included1$geneID, genes.to.be.included2$geneID)
genes.to.be.excluded.3 <- union(genes.to.be.excluded1$geneID, genes.to.be.excluded2$geneID)

genes.to.be.included.FINAL <- genes.to.be.included.3[!genes.to.be.included.3 %in% genes.to.be.excluded.3]

write.table(genes.to.be.included.FINAL, paste0(path, "List.of.genes.for.analysis.tsv"), sep = '\t', row.names = F,
            col.names = T, quote = F)
