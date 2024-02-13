#load list of genes without evidence of mapping bias; only the genes in this list will be analysed further
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load output for GLM testing for parental effects on AI using data from both sexes
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Parental.Effects.on.AI.with.both.sexes/"
glm.results <- read.table(paste0(path, "Test.for.parentalEffects.in.AI.combined.sexes.tsv"), header=T, sep='\t')
#glm.results <- glm.results[glm.results$geneID %in% gene.list$geneID,]

########################################################################
## ESTIMATE % DGRP177 EXPRESSION IN MALES/FEMALES FROM EITHER CROSS ####
########################################################################

#for females, main cross
intercept <- 0
sexEffect <- 0
crossEffect <- 0
fMain.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  crossEffect <- glm.results$Cross.effect[i]
  fMain.DGRP.fractions[i] <- exp(intercept + sexEffect + crossEffect)/(1 + exp(intercept + sexEffect + crossEffect))
}

#for males, main cross
intercept <- 0
sexEffect <- 0
crossEffect <- 0
mMain.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  crossEffect <- glm.results$Cross.effect[i]
  mMain.DGRP.fractions[i] <- exp(intercept - sexEffect + crossEffect)/(1 + exp(intercept - sexEffect + crossEffect))
}

#for females, reciprocal cross
intercept <- 0
sexEffect <- 0
crossEffect <- 0
fReciprocal.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  crossEffect <- glm.results$Cross.effect[i]
  fReciprocal.DGRP.fractions[i] <- exp(intercept + sexEffect - crossEffect)/(1 + exp(intercept + sexEffect - crossEffect))
}

#for males, reciprocal cross
intercept <- 0
sexEffect <- 0
crossEffect <- 0
mReciprocal.DGRP.fractions <- c()

for (i in 1:nrow(glm.results)){
  intercept <- glm.results$Overabundance_of_DGRP177[i]
  sexEffect <- glm.results$Sex.effect[i]
  crossEffect <- glm.results$Cross.effect[i]
  mReciprocal.DGRP.fractions[i] <- exp(intercept - sexEffect - crossEffect)/(1 + exp(intercept - sexEffect - crossEffect))
}

##################################################################################################
#  make a dataframe containing fraction of DGRP-177 expression for each sex/cross combination 
df <- cbind.data.frame(geneID = glm.results$geneID, 
                       Proportion.of.DGRP177.reads.Females.Main.Cross = fMain.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Males.Main.Cross = mMain.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Females.Reciprocal.Cross = fReciprocal.DGRP.fractions,
                       Proportion.of.DGRP177.reads.Males.Reciprocal.Cross = mReciprocal.DGRP.fractions)

# also add two columns containing DGRP-177 % expression average over both crosses, for each sex
df$Avg.DGRP177.proportion.females <- (df$Proportion.of.DGRP177.reads.Females.Main.Cross + 
                                        df$Proportion.of.DGRP177.reads.Females.Reciprocal.Cross)/2
df$Avg.DGRP177.proportion.males <- (df$Proportion.of.DGRP177.reads.Males.Main.Cross + 
                                        df$Proportion.of.DGRP177.reads.Males.Reciprocal.Cross)/2


#save to file
write.table(df, paste0(path, 'Proportion.of.DGRP177.reads.tsv'), quote = F, row.names = F, col.names = T, sep = '\t')
