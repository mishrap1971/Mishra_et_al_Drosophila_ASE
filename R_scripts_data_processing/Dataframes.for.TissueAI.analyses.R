##################################################################################
#### FOR ANALYSIS OF COMBINED DATA FROM MALES & FEMALES  #########################
##################################################################################

#load list of genes without evidence of mapping bias; only the genes in this list will be analysed further
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load output for GLM testing for AI using data from both sexes in heads and gonads 
#use filter to include genes that have no evidence of mapping bias
glm.results <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.combined.sexes/Combined.Sexes.test.for.tissue.effects.in.AI.tsv", header=T)
glm.results <- glm.results[glm.results$geneID %in% gene.list$geneID,]

############################################################################################
#absence/presence of significant AI and significant sex/tissue/interaction effects on AI ##
############################################################################################

#generate 1s for genes with significant AI (intercept term < 0.05) and 0s for genes without
Presence.of.AI <- c()

for (i in 1:nrow(glm.results)){
  if (glm.results[,3][i] < 0.05){
    Presence.of.AI[i] = 1
  } else {
    Presence.of.AI[i] = 0
  }
}

#generate 1s for genes with significant SDAI (sex term < 0.05) and 0s for genes without
Presence.of.SDAI <- c()

for (i in 1:nrow(glm.results)){
  if (glm.results[,5][i] < 0.05){
    Presence.of.SDAI[i] = 1
  } else {
    Presence.of.SDAI[i] = 0
  }
}

#generate 1s for genes with significant tissue-dependent AI (TDAI for short; tissue term < 0.05) and 0s for genes without
Presence.of.TDAI <- c()

for (i in 1:nrow(glm.results)){
  if (glm.results[,7][i] < 0.05){
    Presence.of.TDAI[i] = 1
  } else {
    Presence.of.TDAI[i] = 0
  }
}

#generate 1s for genes with significant sex*tissue interaction effects AI (interaction term < 0.05) and 0s for genes without
Presence.of.interactionEffects <- c()

for (i in 1:nrow(glm.results)){
  if (glm.results[,9][i] < 0.05){
    Presence.of.interactionEffects[i] = 1
  } else {
    Presence.of.interactionEffects[i] = 0
  }
}


####################################################################################################
### Estimating fraction of DGRP-177 expression in heads and gonads for males and females ###########
####################################################################################################

# logit link function: logit(fraction of DGRP177 expression) = intercept + sexEffect*X_sex + 
# tissueEffect*X_tissue + interactionEffect*X_sex.tissue (X_sex = 1 for females, -1 for males,
# X_tissue = 1 for gonads, -1 for heads, X_sex.tissue = X_sex*X_tissue)

#for female gonads
intercept <- 0
sexEffect <- 0
tissueEffect <- 0
interactionEffect <- 0
ovaries.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  tissueEffect <- glm.results$Tissue.effect[i]
  interactionEffect <- glm.results$Sex.Tissue.interaction.effect[i]
  ovaries.DGRP.fractions[i] <- exp(intercept + sexEffect + tissueEffect + interactionEffect)/(1 + exp(intercept + sexEffect + tissueEffect + interactionEffect))
}


#for male gonads
intercept <- 0
sexEffect <- 0
tissueEffect <- 0
interactionEffect <- 0
testes.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  tissueEffect <- glm.results$Tissue.effect[i]
  interactionEffect <- glm.results$Sex.Tissue.interaction.effect[i]
  testes.DGRP.fractions[i] <- exp(intercept - sexEffect + tissueEffect - interactionEffect)/(1 + exp(intercept - sexEffect + tissueEffect - interactionEffect))
}


#for female heads
intercept <- 0
sexEffect <- 0
tissueEffect <- 0
interactionEffect <- 0
fheads.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  tissueEffect <- glm.results$Tissue.effect[i]
  interactionEffect <- glm.results$Sex.Tissue.interaction.effect[i]
  fheads.DGRP.fractions[i] <- exp(intercept + sexEffect - tissueEffect - interactionEffect)/(1 + exp(intercept + sexEffect - tissueEffect - interactionEffect))
}

#for male heads
intercept <- 0
sexEffect <- 0
tissueEffect <- 0
interactionEffect <- 0
mheads.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  tissueEffect <- glm.results$Tissue.effect[i]
  interactionEffect <- glm.results$Sex.Tissue.interaction.effect[i]
  mheads.DGRP.fractions[i] <- exp(intercept - sexEffect - tissueEffect + interactionEffect)/(1 + exp(intercept - sexEffect - tissueEffect + interactionEffect))
}


##################################################################################################
#  make a dataframe containing fraction of DGRP-177 expression for each sex/tissue combination + 
#  presence/absence of AI or tissue/sex/interaction effects on AI

df <- cbind.data.frame(geneID = glm.results$geneID, Presence.of.AI = Presence.of.AI, 
                       Presence.of.sex.dependent.AI = Presence.of.SDAI,
                       Presence.of.tissue.dependent.AI = Presence.of.TDAI,
                       Presence.of.interactionEffects.on.AI = Presence.of.interactionEffects,
                       Proportion.of.DGRP177.reads.Ovaries = ovaries.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Testes = testes.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Female.Heads = fheads.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Male.Heads = mheads.DGRP.fractions)

#save to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Tissue.dependent.AI.both.sexes/"
write.table(df, paste0(path, 'Combined.Sexes.Tissue.Effects.df.tsv'), quote = F, row.names = F, col.names = T, sep = '\t')


###############################################################################################
##################################################################################
#### FOR SEPARTE ANALYSIS OF DATA FROM MALES & FEMALES  #########################
##################################################################################
################################################################################################

#estimate log total allelic expression per gene
#this function outputs log (total allelic expression) for A) all male replicates, B)
#all female replicates & C) all replicates from both sexes
Log.Total.Allelic.Expression <- function(head_tissue, gonad_tissue){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Tissue.Effects.in.Single.Sex/")
  h1 <- read.table(paste0(head_tissue, '_replicate_1.tsv'), header = T)
  h2 <- read.table(paste0(head_tissue, '_replicate_2.tsv'), header = T)
  h3 <- read.table(paste0(head_tissue, '_replicate_3.tsv'), header = T)
  g1 <- read.table(paste0(gonad_tissue, '_replicate_1.tsv'), header = T)
  g2 <- read.table(paste0(gonad_tissue, '_replicate_2.tsv'), header = T)
  g3 <- read.table(paste0(gonad_tissue, '_replicate_3.tsv'), header = T)
  
  
  #combine read count data from all replicates (males + females)
  combined <- cbind.data.frame(h1, h2$DGRP.177, h2$SP.159N, h3$DGRP.177, h3$SP.159N,
                               g1$DGRP.177, g1$SP.159N, g2$DGRP.177, g2$SP.159N, g3$DGRP.177, g3$SP.159N)
  colnames(combined) <- c("geneID", "DGRP.177.1", "SP.159N.1", "DGRP.177.2", "SP.159N.2", "DGRP.177.3", 
                          "SP.159N.3", "DGRP.177.4", "SP.159N.4", "DGRP.177.5", "SP.159N.5", "DGRP.177.6", 
                          "SP.159N.6")
  
  #estimate average total read count over all replicates (6 replicates in all)
  avg.read.counts <- cbind.data.frame(geneID = h1$geneID, avg.read.counts = rowSums(combined[,c(2:13)])/6) 
  
  #estimate logarithmic value of total allelic expression; make dataframe
  all.log.counts <- cbind.data.frame(geneID = avg.read.counts$geneID, 
                                     Log.total.allelic.expression = log(avg.read.counts$avg.read.counts))
  
  return(all.log.counts)
}


##############################################################################################
#absence/presence of significant AI and significant tissue-dependent AI (TDAI) from GLM outputs
presence.of.AI.and.TDAI <- function(sex){
  #load output from GLM testing for AI
  setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Tissue.Effects.in.AI.single.sex/')
  sample <- read.table(paste0(sex, '.test.for.tissue.effects.in.AI.tsv'), header = T)
  
  #generate 1s for genes with significant AI (intercept term < 0.05) and 0s for genes without
  evidence.for.AI <- c()
  
  for (i in 1:nrow(sample)){
    if (sample[,3][i] < 0.05){
      evidence.for.AI[i] = 1
    } else {
      evidence.for.AI[i] = 0
    }
  }
  
  #generate 1s for genes with significant SDAI (sex term < 0.05) and 0s for genes without
  evidence.for.TDAI <- c()
  
  for (i in 1:nrow(sample)){
    if (sample[,5][i] < 0.05){
      evidence.for.TDAI[i] = 1
    } else {
      evidence.for.TDAI[i] = 0
    }
  }
  
  presence.of.AI.TDAI <- cbind.data.frame(geneID = sample$geneID, Presence.of.AI = evidence.for.AI, 
                                          Presence.of.tissue.dependent.AI = evidence.for.TDAI)
  return(presence.of.AI.TDAI)
}


###########################################################################################
#compile outputs above functions to one dataframe for each tissue
#also add information about sex bias (sex bias category + estimates of log2 fold change in expression in head vs gonads)

#log total allele-specific expression
df1.females <- Log.Total.Allelic.Expression('Female_heads', 'Ovaries')
df1.males <- Log.Total.Allelic.Expression('Male_heads', 'Testes')

#presence/absence of AI and tissue-dependent AI
df2.females <- presence.of.AI.and.TDAI('Females')
df2.males <- presence.of.AI.and.TDAI('Males')

#load tissue bias dataframes (specify path to dataframes before that)
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs/Between.Tissues/"
df3.females <- read.table(paste0(path, "Tissue.Biased.Genes.Females.tsv"), header = T, sep = '\t')
df3.males <- read.table(paste0(path, "Tissue.Biased.Genes.Males.tsv"), header = T, sep = '\t')

#merge all 3 dataframes
df.females <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.females, df2.females, df3.females))
df.males <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.males, df2.males, df3.males))

#save to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Tissue.dependent.AI.single.sex.analysis/"
write.table(df.females, paste0(path, "Females.tissueAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')
write.table(df.males, paste0(path, "Males.tissueAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')

