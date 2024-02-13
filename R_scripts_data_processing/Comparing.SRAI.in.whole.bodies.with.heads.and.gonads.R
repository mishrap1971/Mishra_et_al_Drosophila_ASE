#load list of genes that have no evidence of mapping bias
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load results of GLM testing for sex-reversed AI for genes found to have sex-reversed AI
SRAI.genes <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Sex.Reversed.AI/Genes.with.Sex.Reversed.AI.tsv", header=T)

#add a column to above dataframe that lists the sex where direction of AI is positive
for (i in 1:nrow(SRAI.genes)){
  if (SRAI.genes[i,2] > 0){
    SRAI.genes$Sex.with.Positive.AI[i] = "M"
  } else {
    SRAI.genes$Sex.with.Positive.AI[i] = "F"
  }
}

#load dataframe obtained from the outputs of the GLM testing for sex effects on AI in heads and gonads
#use filter to only include genes without evidence of mapping bias
heads <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/Heads.sexAI.df.tsv",header=T,sep='\t')
gonads <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/Gonads.sexAI.df.tsv",header=T,sep='\t')
heads <- heads[heads$geneID %in% gene.list$geneID,]
gonads <- gonads[gonads$geneID %in% gene.list$geneID,]

#for both heads and gonads obtain GLM results for genes found to have sex-reversed AI (tested using whole body data)
heads.SRAI <- heads[heads$geneID %in% SRAI.genes$geneID,]
gonads.SRAI <- gonads[gonads$geneID %in% SRAI.genes$geneID,]

#make dataframe with geneID, presence/absence of sex-dependent AI, -(bias in DGRP177 expression in both sexes), and 
#point estimates of sex-reversed AI
#sex-reversed point estimate of AI = bias in DGRP-177 expression in males)*(bias in DGRP-177 expression in females)
#we use negative values of the bias in DGRP177 to be consistent with how SRAI was tested in the whole bodies; AI is
#tested as the 'overabundance in DGRP177', which is deviation of DGRP177 expression from 0.5 (DGRP177 - 0.5) -- bias in 
#DGRP177 expression is defined the other way around, 0.5 - DGRP177
heads.SRAI.df <- cbind.data.frame(geneID = heads.SRAI$geneID,
                                  Presence.of.sex.dependent.AI = heads.SRAI$Presence.of.sex.dependent.AI,
                                  Bias.in.DGRP.expression.male = -heads.SRAI$Bias.in.DGRP.expression.Male,
                                  Bias.in.DGRP.expression.female = -heads.SRAI$Bias.in.DGRP.expression.Female,
                                  SRAI.Point.Estimate = heads.SRAI$Bias.in.DGRP.expression.Male*heads.SRAI$Bias.in.DGRP.expression.Female)
gonads.SRAI.df <- cbind.data.frame(geneID = gonads.SRAI$geneID,
                                   Presence.of.sex.dependent.AI = gonads.SRAI$Presence.of.sex.dependent.AI,
                                   Bias.in.DGRP.expression.male = -gonads.SRAI$Bias.in.DGRP.expression.Male,
                                   Bias.in.DGRP.expression.female = -gonads.SRAI$Bias.in.DGRP.expression.Female,
                                   SRAI.Point.Estimate = gonads.SRAI$Bias.in.DGRP.expression.Male*gonads.SRAI$Bias.in.DGRP.expression.Female)

#add a column to above dataframes for presence/absence of sex-reversed AI (i.e., opposite signs of 
#bias in DGRP177 expression in sexes)
for (i in 1:nrow(heads.SRAI.df)){
  if (heads.SRAI.df[i,3]*heads.SRAI.df[i,4] < 0){
    heads.SRAI.df$Presence.of.sex.reversed.AI[i] = "1"
  } else {
    heads.SRAI.df$Presence.of.sex.reversed.AI[i] = "0"
  }
}

for (i in 1:nrow(gonads.SRAI.df)){
  if (gonads.SRAI.df[i,3]*gonads.SRAI.df[i,4] < 0){
    gonads.SRAI.df$Presence.of.sex.reversed.AI[i] = "1"
  } else {
    gonads.SRAI.df$Presence.of.sex.reversed.AI[i] = "0"
  }
}



#subset above dataframes to only include genes with sex-reversed point estimates of AI
#add a column that lists the sex where direction of AI (i.e., bias in DGRP177 expression) is positive
heads.SRAI.genes.only <- heads.SRAI.df[heads.SRAI.df$Presence.of.sex.reversed.AI == 1,]
gonads.SRAI.genes.only <- gonads.SRAI.df[gonads.SRAI.df$Presence.of.sex.reversed.AI == 1,]

for (i in 1:nrow(heads.SRAI.genes.only)){
  if (heads.SRAI.genes.only[i,3] > 0){
    heads.SRAI.genes.only$Sex.with.Positive.AI[i] = "M"
  } else {
    heads.SRAI.genes.only$Sex.with.Positive.AI[i] = "F"
  }
}

for (i in 1:nrow(gonads.SRAI.genes.only)){
  if (gonads.SRAI.genes.only[i,3] > 0){
    gonads.SRAI.genes.only$Sex.with.Positive.AI[i] = "M"
  } else {
    gonads.SRAI.genes.only$Sex.with.Positive.AI[i] = "F"
  }
}

###################################################################################
#SRAI in genes (from whole body) that also had sex-reversed point estimates of AI in heads or gonads
SRAI.tested.in.heads <- SRAI.genes[SRAI.genes$geneID %in% heads.SRAI.genes.only$geneID,]
SRAI.tested.in.gonads <- SRAI.genes[SRAI.genes$geneID %in% gonads.SRAI.genes.only$geneID,]

heads.SRAI.genes.only$Sex.with.Positive.AI.in.Whole.Body <- SRAI.tested.in.heads$Sex.with.Positive.AI
gonads.SRAI.genes.only$Sex.with.Positive.AI.in.Whole.Body <- SRAI.tested.in.gonads$Sex.with.Positive.AI

#number of genes that can be tested for SRAI in heads and gonads
nrow(heads.SRAI.df)
nrow(gonads.SRAI.df)

#number of genes where the direction of SRAI in heads/gonads is same as SRAI in whole bodies
length(which(heads.SRAI.genes.only$Sex.with.Positive.AI == heads.SRAI.genes.only$Sex.with.Positive.AI.in.Whole.Body))
length(which(gonads.SRAI.genes.only$Sex.with.Positive.AI == gonads.SRAI.genes.only$Sex.with.Positive.AI.in.Whole.Body))
