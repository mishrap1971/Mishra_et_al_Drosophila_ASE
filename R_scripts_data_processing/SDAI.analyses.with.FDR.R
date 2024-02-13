#estimate FDRs using p-values from GLM for testing sex effects on AI; save output to file
EstimateFDR <- function(tissue){
  glm.results <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/Sex.Effects.in.AI/",
                                   tissue, ".test.for.AI.tsv"), header = T)
  
  glm.results$FDR.intercept <- p.adjust(glm.results[,3], method = 'BH')
  glm.results$FDR.sex.term <- p.adjust(glm.results[,5], method = 'BH')
  
  write.table(glm.results, paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/FDR.Sex.Effects.in.AI/",
                                  tissue, ".FDR.test.for.AI.tsv"), sep = '\t', quote = F, row.names = F, col.names = T)
  return(glm.results)
}

EstimateFDR('Gonads')
EstimateFDR('Heads')
EstimateFDR('Whole.bodies')
EstimateFDR('Reciprocal.whole.bodies')

#####################################################################################################################
## Use code in this section to make dataframes with: geneID, log (total allelic expression) (calculated   ###########
## for male, female and combined read counts), presence/absence of AI and SDAI based on FDR cut-off,      ###########
## proportion of DGRP expression in males/females and 4 different estimates of bias in DGRP expression    ###########
#####################################################################################################################

###############################################################################################
#estimate log total allelic expression per gene
#this function outputs log (total allelic expression) for A) all male replicates, B)
#all female replicates & C) all replicates from both sexes
Log.Total.Allelic.Expression <- function(female_tissue, male_tissue){
  
  #load read count data
  setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/For.Sex.Effects.in.AI")
  f1 <- read.table(paste0(female_tissue, '_replicate_1.tsv'), header = T)
  f2 <- read.table(paste0(female_tissue, '_replicate_2.tsv'), header = T)
  f3 <- read.table(paste0(female_tissue, '_replicate_3.tsv'), header = T)
  m1 <- read.table(paste0(male_tissue, '_replicate_1.tsv'), header = T)
  m2 <- read.table(paste0(male_tissue, '_replicate_2.tsv'), header = T)
  m3 <- read.table(paste0(male_tissue, '_replicate_3.tsv'), header = T)
  
  #combine read count data from all male replicates
  combined.m <- cbind.data.frame(m1, m2$DGRP.177, m2$SP.159N, m3$DGRP.177, m3$SP.159N)
  colnames(combined.m) <- c("geneID", "DGRP.177.1", "SP.159N.1", "DGRP.177.2", "SP.159N.2", "DGRP.177.3", 
                            "SP.159N.3")
  
  #combine read count data from all female replicates
  combined.f <- cbind.data.frame(f1, f2$DGRP.177, f2$SP.159N, f3$DGRP.177, f3$SP.159N)
  colnames(combined.f) <- c("geneID", "DGRP.177.1", "SP.159N.1", "DGRP.177.2", "SP.159N.2", "DGRP.177.3", 
                            "SP.159N.3")
  
  #combine read count data from all replicates (males + females)
  combined.all <- cbind.data.frame(f1, f2$DGRP.177, f2$SP.159N, f3$DGRP.177, f3$SP.159N,
                                   m1$DGRP.177, m1$SP.159N, m2$DGRP.177, m2$SP.159N, m3$DGRP.177, m3$SP.159N)
  colnames(combined.all) <- c("geneID", "DGRP.177.1", "SP.159N.1", "DGRP.177.2", "SP.159N.2", "DGRP.177.3", 
                              "SP.159N.3", "DGRP.177.4", "SP.159N.4", "DGRP.177.5", "SP.159N.5", "DGRP.177.6", 
                              "SP.159N.6")
  
  #estimate average total read count over all replicates (3 replicates for Males/Females; 6 replicates in all)
  m.avg.read.counts <- cbind.data.frame(geneID = m1$geneID, avg.read.counts = rowSums(combined.m[,c(2:7)])/3) 
  f.avg.read.counts <- cbind.data.frame(geneID = f1$geneID, avg.read.counts = rowSums(combined.f[,c(2:7)])/3) 
  all.avg.read.counts <- cbind.data.frame(geneID = m1$geneID, avg.read.counts = rowSums(combined.all[,c(2:13)])/6) 
  
  #estimate logarithmic value, for each sex
  m.log.counts <- cbind.data.frame(geneID = m.avg.read.counts$geneID, 
                                   Log.total.allelic.expression.males = log(m.avg.read.counts$avg.read.counts))
  f.log.counts <- cbind.data.frame(geneID = f.avg.read.counts$geneID, 
                                   Log.total.allelic.expression.females = log(f.avg.read.counts$avg.read.counts))
  all.log.counts <- cbind.data.frame(geneID = all.avg.read.counts$geneID, 
                                     Log.total.allelic.expression.all = log(all.avg.read.counts$avg.read.counts))
  
  #merge dataframes of log(total allelic expression) for males and females
  df <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(m.log.counts, f.log.counts, all.log.counts))
  
  return(df)
}


##############################################################################################
#absence/presence of significant AI and significant sex-dependent AI (SDAI) from GLM outputs
presence.of.AI.and.SDAI <- function(samplename){
  #load output from GLM testing for AI
  setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/FDR.Sex.Effects.in.AI/')
  sample <- read.table(paste0(samplename, '.FDR.test.for.AI.tsv'), header = T)
  
  #generate 1s for genes with significant AI (intercept term < 0.05 OR sex term < 0.05) and 0s for genes without
  evidence.for.AI <- c()
  
  for (i in 1:nrow(sample)){
    #if (sample[,6][i] < 0.05 | sample[,7][i] < 0.05){
    if (sample[,6][i] < 0.05){
      evidence.for.AI[i] = 1
    } else {
      evidence.for.AI[i] = 0
    }
  }
  
  #generate 1s for genes with significant SDAI (intercept term < 0.05 OR sex term < 0.05) and 0s for genes without
  evidence.for.SDAI <- c()
  
  for (i in 1:nrow(sample)){
    if (sample[,7][i] < 0.05){
      evidence.for.SDAI[i] = 1
    } else {
      evidence.for.SDAI[i] = 0
    }
  }
  
  presence.of.AI.SDAI <- cbind.data.frame(geneID = sample$geneID, Presence.of.AI = evidence.for.AI, 
                                          Presence.of.sex.dependent.AI = evidence.for.SDAI)
  return(presence.of.AI.SDAI)
}

# gonads <- presence.of.AI.SDAI("Gonads")
# heads <- presence.of.AI.SDAI("Heads")
# whole.bodies <- presence.of.AI.SDAI("Whole_body")
# reciprocal <- presence.of.AI.SDAI("Reciprocal_whole_body")

########################################################################################################
#proportion of DGRP-177 reads per gene in females and males, inferred from GLM results (details in supplement)
Proportion.of.DGRP.reads <- function(samplename){
  #load output from GLM testing for AI
  setwd('C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/GLM.outputs/FDR.Sex.Effects.in.AI/')
  sample <- read.table(paste0(samplename, '.FDR.test.for.AI.tsv'), header = T)
  
  #for proportion of DGRP reads in males
  interceptTerm <- 0
  interactionTerm <- 0
  male.DGRP.fraction <- 0
  male.DGRP.fraction.list <- c()
  
  for (i in 1:nrow(sample)){
    interceptTerm <- sample$Overabundance_of_DGRP177[i]
    interactionTerm <- sample$Sex.difference[i]
    male.DGRP.fraction <- (exp(interceptTerm - interactionTerm))/(1 + exp(interceptTerm - interactionTerm))
    male.DGRP.fraction.list <- c(male.DGRP.fraction.list, male.DGRP.fraction)
  }
  
  #for proportion of DGRP reads in females
  interceptTerm <- 0
  interactionTerm <- 0
  female.DGRP.fraction <- 0
  female.DGRP.fraction.list <- c()
  
  for (i in 1:nrow(sample)){
    interceptTerm <- sample$Overabundance_of_DGRP177[i]
    interactionTerm <- sample$Sex.difference[i]
    female.DGRP.fraction <- exp(interceptTerm + interactionTerm)/(1 + exp(interceptTerm + interactionTerm))
    female.DGRP.fraction.list <- c(female.DGRP.fraction.list, female.DGRP.fraction)
  }
  
  DGRP.fractions <- cbind.data.frame(geneID = sample$geneID, 
                                     Proportion.of.DGRP177.reads.Males = male.DGRP.fraction.list, 
                                     Proportion.of.DGRP177.reads.Females = female.DGRP.fraction.list)
  
  return(DGRP.fractions)
}

# gonads <- Proportion.of.DGRP.reads("Gonads")
# heads <- Proportion.of.DGRP.reads("Heads")
# whole.bodies <- Proportion.of.DGRP.reads("Whole.bodies")
# reciprocal <- Proportion.of.DGRP.reads("Reciprocal.whole.bodies")

#########################################################################################################
#Bias in DGRP-177 expression in males, females (0.5 - DGRP-177)
#Absolute bias in DGRP-177 in males, females (|0.5 - DGRP-177|)
#Sex-averaged absolute bias ([abs. bias in males + abs. bias in female]/2)
#Sex difference in bias (abs[bias in females - bias in males])
#Input should be the output from the prev. function
#NOTE: make sure 2nd column of input is DGRP-177 proportion in males and 3rd column is proportion in females

BiasCalculations <- function(DGRP177.proportions){
  DGRP.bias.female <- c()
  DGRP.bias.male <- c()
  abs.DGRP.bias.female <- c()
  abs.DGRP.bias.male <- c()
  sex.avg.abs.bias.DGRP <- c()
  abs.sex.diff.in.bias.DGRP <- c()
  
  for (i in 1:nrow(DGRP177.proportions)){
    female.proportion <- DGRP177.proportions[i,3]
    male.proportion <- DGRP177.proportions[i,2]
    
    DGRP.bias.female[i] <- 0.5 - female.proportion
    DGRP.bias.male[i] <- 0.5 - male.proportion
    abs.DGRP.bias.female[i] <- abs(0.5 - female.proportion)
    abs.DGRP.bias.male[i] <- abs(0.5 - male.proportion)
    sex.avg.abs.bias.DGRP[i] <- (abs(0.5 - female.proportion) + abs(0.5 - male.proportion))/2
    abs.sex.diff.in.bias.DGRP[i] <- abs((0.5 - female.proportion) - (0.5 - male.proportion))
  }
  
  df.output <- cbind.data.frame(geneID = DGRP177.proportions$geneID,
                                Bias.in.DGRP.expression.Male = DGRP.bias.male,
                                Bias.in.DGRP.expression.Female = DGRP.bias.female,
                                Absolute.Bias.in.DGRP.expression.Male = abs.DGRP.bias.male,
                                Absolute.Bias.in.DGRP.expression.Female = abs.DGRP.bias.female,
                                Sex.Averaged.Absolute.Bias.in.DGRP.expression = sex.avg.abs.bias.DGRP,
                                Absolute.Sex.Difference.in.Bias.in.DGRP.expression = abs.sex.diff.in.bias.DGRP)
  
  return(df.output)
}

###########################################################################################
#compile outputs from all of the above functions to one dataframe for each tissue
#also add information about sex bias (sex bias category + estimates of log2 fold change in expression in males vs females)
#log total allele-specific expression
df1.gonads <- Log.Total.Allelic.Expression('Ovaries', 'Testes')
df1.heads <- Log.Total.Allelic.Expression('Female_heads', 'Male_heads')
df1.whole.bodies <- Log.Total.Allelic.Expression('Female_whole_bodies', 'Male_whole_bodies')
df1.reciprocal <- Log.Total.Allelic.Expression('Female_whole_bodies_reciprocal', 'Male_whole_bodies_reciprocal')

#load sex bias dataframes (specify path to dataframes before that)
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/DESeq2.outputs/"
df2.gonads <- read.table(paste0(path, "Sex.Biased.Genes.Gonads.tsv"), header = T, sep = '\t')
df2.heads <- read.table(paste0(path, "Sex.Biased.Genes.Heads.tsv"), header = T, sep = '\t')
df2.whole.bodies <- read.table(paste0(path, "Sex.Biased.Genes.Whole.bodies.tsv"), header = T, sep = '\t')
df2.reciprocal <- read.table(paste0(path, "Sex.Biased.Genes.Reciprocal.whole.bodies.tsv"), header = T, sep = '\t')

#presence or absence of AI
df3.gonads <- presence.of.AI.and.SDAI('Gonads')
df3.heads <- presence.of.AI.and.SDAI('Heads')
df3.whole.bodies <- presence.of.AI.and.SDAI('Whole.bodies')
df3.reciprocal <- presence.of.AI.and.SDAI('Reciprocal.whole.bodies')

#proportion of DGRP reads in males and females
df4.gonads <- Proportion.of.DGRP.reads('Gonads')
df4.heads <- Proportion.of.DGRP.reads('Heads')
df4.whole.bodies <- Proportion.of.DGRP.reads('Whole.bodies')
df4.reciprocal <- Proportion.of.DGRP.reads('Reciprocal.whole.bodies')

#calculations around bias in DGRP-177 expression (i.e., deviation from expected value of 0.5)
df5.gonads <- BiasCalculations(df4.gonads)
df5.heads <- BiasCalculations(df4.heads)
df5.whole.bodies <- BiasCalculations(df4.whole.bodies)
df5.reciprocal <- BiasCalculations(df4.reciprocal)

#merge all 4 dataframes
df.gonads <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.gonads, df2.gonads, df3.gonads, df4.gonads, df5.gonads))
df.heads <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.heads, df2.heads, df3.heads, df4.heads, df5.heads))
df.whole.bodies <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.whole.bodies, df2.whole.bodies, df3.whole.bodies, df4.whole.bodies, df5.whole.bodies))
df.reciprocal <- Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(df1.reciprocal, df2.reciprocal, df3.reciprocal, df4.reciprocal, df5.reciprocal))

#save to file
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI.FDR/"
write.table(df.gonads, paste0(path, "Gonads.sexAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')
write.table(df.heads, paste0(path, "Heads.sexAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')
write.table(df.whole.bodies, paste0(path, "Whole.bodies.sexAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')
write.table(df.reciprocal, paste0(path, "Reciprocal.whole.bodies.sexAI.df.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')
