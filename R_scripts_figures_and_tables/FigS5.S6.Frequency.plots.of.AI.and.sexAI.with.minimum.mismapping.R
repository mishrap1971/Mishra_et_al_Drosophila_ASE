library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(scales)
library(emmeans)

#load list of genes that had <99% reads mapped to the correct parent during competitive mapping of parental genomic data
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Mismapped.Reads.in.Parental.Genomes/"
mismapped.genes <- read.table(paste0(path, 'Below.99.percent.correct.mapping.in.both.tsv'), header = T)

#load list of genes without evidence of mapping bias; only the genes in this list will be analysed further
#remove 'mismapped genes' (above) from this list
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)
gene.list <- gene.list[!gene.list$geneID %in% mismapped.genes$geneID,]

#load dataframes derived from GLM outputs, mainly information about presence/absence of AI and sex-dependent AI
#retain genes that have no evidence of mapping bias
setwd("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/")
gonads <- read.table("Gonads.sexAI.df.tsv", header = TRUE, sep = '\t')
heads <- read.table("Heads.sexAI.df.tsv", header = TRUE, sep = '\t')
whole.bodies <- read.table("Whole.bodies.sexAI.df.tsv", header = TRUE, sep = '\t')
reciprocal <- read.table("Reciprocal.whole.bodies.sexAI.df.tsv", header = TRUE, sep = '\t')

gonads <- gonads[gonads$geneID %in% gene.list,]
heads <- heads[heads$geneID %in% gene.list,]
whole.bodies <- whole.bodies[whole.bodies$geneID %in% gene.list,]
reciprocal <- reciprocal[reciprocal$geneID %in% gene.list,]

#################################
##     frequency of AI      #####
#################################
#define functions to calculate frequency of AI and confidence intervals for a given sex bias bin
#make sure 'tissue' is a dataframe
frequency.of.AI <- function(sex.bias.bin, tissue){
  frequency <- c(length(which(tissue$Sex.Bias.Category == sex.bias.bin & tissue$Presence.of.AI == 1))/length(which(tissue$Sex.Bias.Category == sex.bias.bin)))
  return(frequency)
}

stdErr.frequency.of.AI <- function(sex.bias.bin, tissue){
  mod.AI <- glm(Presence.of.AI ~ Sex.Bias.Category + Log.total.allelic.expression.all, family = binomial, data = tissue)
  posthoc <- summary(emmeans(mod.AI, c("Sex.Bias.Category"), type = "response"))
  SE <- c(posthoc[which(posthoc$Sex.Bias.Category == sex.bias.bin),][1,3])
  return(SE)
}


# stdErr.frequency.of.AI <- function(sex.bias.bin, tissue){
#   df <- tissue[which(tissue$Sex.Bias.Category == sex.bias.bin),]
# 
#   #draw random samples of genes from a given sex bias category with replacement, estimate frequency of AI in random sample
#   #repeat 10,000 times
#   boot.df <- c()
#   frequency.AI.list <- c()
#   for (i in 1:10000){
#     boot.df <- sample(df$Presence.of.AI, nrow(df), replace = T)
#     frequency.AI <- length(which(boot.df == 1))/length(boot.df)
#     frequency.AI.list <- c(frequency.AI.list, frequency.AI)
#   }
#   #estimate 95% CIs
#   #confInts <- t.test(frequency.AI.list)$conf.int[1:2]
#   std.err <- sd(frequency.AI.list)/sqrt(length(frequency.AI.list))
#   
#   return(std.err)
# }

#apply above functions to all tissues
#estimate frequency of AI
sex.bias.bins <- c('Highly female-biased', 'Moderately female-biased', 'Unbiased', 'Moderately male-biased', 'Highly male-biased')
gonads.freq.AI <- as.numeric(sapply(sex.bias.bins, frequency.of.AI, tissue = gonads))
heads.freq.AI <- as.numeric(sapply(sex.bias.bins, frequency.of.AI, tissue = heads))
whole.bodies.freq.AI <- as.numeric(sapply(sex.bias.bins, frequency.of.AI, tissue = whole.bodies))
reciprocal.freq.AI <- as.numeric(sapply(sex.bias.bins, frequency.of.AI, tissue = reciprocal))

#estimate std. error in frequency of AI
gonads.SE.AI <- as.numeric(sapply(sex.bias.bins, stdErr.frequency.of.AI, tissue = gonads))
heads.SE.AI <- as.numeric(sapply(sex.bias.bins, stdErr.frequency.of.AI, tissue = heads))
whole.bodies.SE.AI <- as.numeric(sapply(sex.bias.bins, stdErr.frequency.of.AI, tissue = whole.bodies))
reciprocal.SE.AI <- as.numeric(sapply(sex.bias.bins, stdErr.frequency.of.AI, tissue = gonads))

#make dataframes
df.AI <- data.frame(Sample = rep(c("Gonads", "Heads", "WB (main)", "WB (reciprocal)"), each=5),
                    SexBias = rep(c("Highly female-biased", "Moderately female-biased", "Unbiased", 
                                    "Moderately male-biased", "Highly male-biased"),4), 
                    Frequency=c(gonads.freq.AI, heads.freq.AI, whole.bodies.freq.AI, reciprocal.freq.AI))

df.SE.AI <- data.frame(Sample = rep(c("Gonads", "Heads", "WB (main)", "WB (reciprocal)"), each=5),
                       SexBias = rep(c("Highly female-biased", "Moderately female-biased", "Unbiased", 
                                       "Moderately male-biased", "Highly male-biased"),4), 
                       SE =c(gonads.SE.AI, heads.SE.AI, whole.bodies.SE.AI, reciprocal.SE.AI))

#make plots
frequency1 <- df.AI$Frequency
SE1 <- df.SE.AI$SE
df.AI$SexBias <- factor(df.AI$SexBias, levels=unique(df.AI$SexBias))

#save as TIFF file with dimensions 1200x600 (width x height)
p1 <- ggplot(data=df.AI, aes(x=SexBias, y=Frequency, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge(), color = "white", lwd = 1.5, show.legend = F) + 
  theme_classic() + 
  theme(legend.position = c(.5,0.82),
        legend.box.background = element_rect(linewidth = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=25, margin = margin(0,5,0,0), color = "black"), 
        axis.title.y = element_text(size=23, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=22, face = "bold")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#0072B2", "#D55E00")) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  labs(x = "", 
       y = paste0("Frequency of genes", "\n", "with AI"), 
       tag = "A") +
  geom_errorbar(aes(ymin = frequency1 - SE1, ymax = frequency1 + SE1), width=.2, position=position_dodge(.9)) +
  ylim (0,1) 
p1


########################################################################################
#######################################################################################
########################################
#   frequency of sex-depedent AI   #####
########################################
#define functions to calculate frequency of sexAI and confidence intervals for a given sex bias bin
frequency.of.sexAI <- function(sex.bias.bin, sample){
  frequency <- c(length(which(sample$Sex.Bias.Category == sex.bias.bin & sample$Presence.of.sex.dependent.AI == 1))/length(which(sample$Sex.Bias.Category == sex.bias.bin)))
  return(frequency)
}

stderr.sexAI <- function(sex.bias.bin, sample){
  mod.sexAI <- glm(Presence.of.sex.dependent.AI ~ Sex.Bias.Category + Log.total.allelic.expression.all, family = binomial, data = sample)
  posthoc <- summary(emmeans(mod.sexAI, c("Sex.Bias.Category"), type = "response"))
  SE <- c(posthoc[which(posthoc$Sex.Bias.Category == sex.bias.bin),][1,3])
  return(SE)
}

#apply above functions to all tissues ('samples')
sex.bias.bins <- c('Highly female-biased', 'Moderately female-biased', 'Unbiased', 'Moderately male-biased', 'Highly male-biased')
gonads.freq.sexAI <- as.numeric(sapply(sex.bias.bins, frequency.of.sexAI, sample = gonads))
heads.freq.sexAI <- as.numeric(sapply(sex.bias.bins, frequency.of.sexAI, sample = heads))
whole.bodies.freq.sexAI <- as.numeric(sapply(sex.bias.bins, frequency.of.sexAI, sample = whole.bodies))
reciprocal.freq.sexAI <- as.numeric(sapply(sex.bias.bins, frequency.of.sexAI, sample = reciprocal))

gonads.SE.sexAI <- as.numeric(sapply(sex.bias.bins, stderr.sexAI, sample = gonads))
heads.SE.sexAI <- as.numeric(sapply(sex.bias.bins, stderr.sexAI, sample = heads))
whole.bodies.SE.sexAI <- as.numeric(sapply(sex.bias.bins, stderr.sexAI, sample = whole.bodies))
reciprocal.SE.sexAI <- as.numeric(sapply(sex.bias.bins, stderr.sexAI, sample = gonads))

#make dataframes
df.sexAI <- data.frame(Sample = rep(c("Gonads", "Heads", "WB (main)", "WB (reciprocal)"), each=5),
                       SexBias = rep(c("Highly female-biased", "Moderately female-biased", "Unbiased", 
                                       "Moderately male-biased", "Highly male-biased"),4), 
                       Frequency=c(gonads.freq.sexAI, heads.freq.sexAI, whole.bodies.freq.sexAI, reciprocal.freq.sexAI))

df.stderr.sexAI <- data.frame(Sample = rep(c("Gonads", "Heads", "WB (main)", "WB (reciprocal)"), each=5),
                              SexBias = rep(c("Highly female-biased", "Moderately female-biased", "Unbiased", 
                                              "Moderately male-biased", "Highly male-biased"),4), 
                              SE =c(gonads.SE.sexAI, heads.SE.sexAI, whole.bodies.SE.sexAI, reciprocal.SE.sexAI))

#make plots
frequency <- df.sexAI$Frequency
SE2 <- df.stderr.sexAI$SE
df.sexAI$SexBias <- factor(df.sexAI$SexBias, levels=unique(df.sexAI$SexBias))

#save as TIFF file with dimensions 1200x600 (width x height)
p2 <- ggplot(data=df.sexAI, aes(x=SexBias, y=Frequency, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge(), color = "white", lwd = 1.5, show.legend = T) + 
  theme_classic() + 
  theme(legend.position = "bottom",
        #legend.position = c(.5,0.82),
        legend.box.background = element_rect(linewidth = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=25, margin = margin(0,5,0,0), color = "black"), 
        axis.title.y = element_text(size=23, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=22, face = "bold")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#0072B2", "#D55E00")) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  labs(x = "", 
       y = paste0("Frequency of genes with", "\n", "sex-dependent AI"), 
       tag = "B") +
  geom_errorbar(aes(ymin = frequency - SE2, ymax = frequency + SE2), width=.2, position=position_dodge(.9)) +
  ylim (0,0.75) 
p2

#make 2-panel figure with both frequency plots; save with dimensions 1200x1000
grid.arrange(p1, p2, nrow = 2)


#######################################################################################################################
## Fig S6 histograms ##################################################################################################

gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

setwd("C:/Users/mishr/OneDrive/Desktop/ASE.Codes.Data/Data/Mismapping.of.Genomic.Reads")
low.mismapping.in.both <- read.table('Above.99.percent.correct.mapping.in.both.tsv', header = T)            #set I
low.mismapping.in.dgrp.only <- read.table('Above.99.percent.correct.mapping.in.dgrp.only.tsv', header = T)  #set II
low.mismapping.in.sp.only <- read.table('Above.99.percent.correct.mapping.in.sp.only.tsv', header = T)      #set III 
high.mismapping.in.both <- read.table('Below.99.percent.correct.mapping.in.both.tsv', header = T)           #set IV

setwd("C:/Users/mishr/OneDrive/Desktop/ASE.Codes.Data/Data/Sex.Effect.Dataframes")
whole.bodies <- read.table("Whole.bodies.main.cross.df.tsv", header = TRUE, sep = '\t')

setwd("C:/Users/mishr/OneDrive/Desktop/ASE.Codes.Data/Data/Whole.Bodies.and.SSDR/")
fraction.dgrp <- read.table("Fraction.of.DGRP.reads.both.crosses.tsv", header = T)
fraction.dgrp <- fraction.dgrp[fraction.dgrp$geneID %in% gene.list$geneID,]
avg.dgrp.fraction <- (fraction.dgrp$females.main.dgrp.fraction + fraction.dgrp$males.main.dgrp.fraction + 
                        fraction.dgrp$females.reciprocal.dgrp.fraction + fraction.dgrp$males.reciprocal.dgrp.fraction)/4
avg.dgrp.fraction <- cbind.data.frame(geneID = fraction.dgrp$geneID, average.DGRP.fraction = avg.dgrp.fraction)


#histograms of fraction of DGRP reads for genes in all 4 sets described below
#histograms include all genes tested
setI.avg.dgrp <- avg.dgrp.fraction[avg.dgrp.fraction$geneID %in% low.mismapping.in.both$geneID,]
setII.avg.dgrp <- avg.dgrp.fraction[avg.dgrp.fraction$geneID %in% low.mismapping.in.dgrp.only$geneID,]
setIII.avg.dgrp <- avg.dgrp.fraction[avg.dgrp.fraction$geneID %in% low.mismapping.in.sp.only$geneID,]
setIV.avg.dgrp <- avg.dgrp.fraction[avg.dgrp.fraction$geneID %in% high.mismapping.in.both$geneID,]


p1 <- ggplot(setI.avg.dgrp, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 80) +
  theme_minimal() +
  annotate("text", x = 0.15, y = 300, size = 7, label = paste0("N = ", nrow(setI.avg.dgrp))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "Frequency",
       tag = "A") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1)
p1

p2 <- ggplot(setII.avg.dgrp, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 83) +
  theme_minimal() +
  annotate("text", x = 0.15, y = 70, size = 7, label = paste0("N = ", nrow(setII.avg.dgrp))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "",
       tag = "B") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1)
p2

p3 <- ggplot(setIII.avg.dgrp, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 80) +
  theme_minimal() +
  annotate("text", x = 0.15, y = 55, size = 7, label = paste0("N = ", nrow(setIII.avg.dgrp))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "Fraction of DGRP-177 reads",
       y = "",
       tag = "C") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1)
p3

p4 <- ggplot(setIV.avg.dgrp, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 81) +
  theme_minimal() +
  annotate("text", x = 0.15, y = 30, size = 7, label = paste0("N = ", nrow(setI.avg.dgrp))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "",
       tag = "D") +
  geom_vline(aes(xintercept=0.5), color = "blue", linetype = "dashed", size = 1) + xlim(0.00, 1.00)
p4

#save as TIFF file with dimensions 1700 x 600
grid.arrange(p1, p2, p3, p4, nrow=2)


##histograms for genes with AI
genes.with.AI <- whole.bodies[which(whole.bodies$evidence.for.AI == 1),]$geneID
avg.dgrp.fraction.AI <- avg.dgrp.fraction[avg.dgrp.fraction$geneID %in% genes.with.AI,]

setI.avg.dgrp.AI <- avg.dgrp.fraction.AI[avg.dgrp.fraction.AI$geneID %in% low.mismapping.in.both$geneID,]
setII.avg.dgrp.AI <- avg.dgrp.fraction.AI[avg.dgrp.fraction.AI$geneID %in% low.mismapping.in.dgrp.only$geneID,]
setIII.avg.dgrp.AI <- avg.dgrp.fraction.AI[avg.dgrp.fraction.AI$geneID %in% low.mismapping.in.sp.only$geneID,]
setIV.avg.dgrp.AI <- avg.dgrp.fraction.AI[avg.dgrp.fraction.AI$geneID %in% high.mismapping.in.both$geneID,]

p1 <- ggplot(setI.avg.dgrp.AI, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 80) +
  theme_minimal() +
  annotate("text", x = 0.1, y = 135, size = 7, label = paste0("N = ", nrow(setI.avg.dgrp.AI))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "Frequency",
       tag = "A") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1) + xlim(0.00, 1.00)
p1

p2 <- ggplot(setII.avg.dgrp.AI, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 83) +
  theme_minimal() +
  annotate("text", x = 0.2, y = 35, size = 7, label = paste0("N = ", nrow(setII.avg.dgrp.AI))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "",
       tag = "B") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1) + xlim(0.00, 1.00)
p2

p3 <- ggplot(setIII.avg.dgrp.AI, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 80) +
  theme_minimal() +
  annotate("text", x = 0.2, y = 25, size = 7, label = paste0("N = ", nrow(setIII.avg.dgrp.AI))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "Fraction of DGRP-177 reads",
       y = "",
       tag = "C") +
  geom_vline(aes(xintercept=0.5), 
             color = "blue", linetype = "dashed", size = 1) + xlim(0.00, 1.00)
p3

p4 <- ggplot(setIV.avg.dgrp.AI, aes(x=average.DGRP.fraction)) + 
  geom_histogram(color="grey78", fill="grey92", bins = 81) +
  theme_minimal() +
  annotate("text", x = 0.15, y = 10, size = 7, label = paste0("N = ", nrow(setI.avg.dgrp.AI))) +
  theme(axis.text.x = element_text(size=22, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=22, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        axis.title.y = element_text(size=25, margin = margin(0,0,0,0), color = "black"), 
        plot.tag = element_text(size=25, face = "bold")) +
  labs(x = "",
       y = "",
       tag = "D") +
  geom_vline(aes(xintercept=0.5), color = "blue", linetype = "dashed", size = 1) + xlim(0.00, 1.00)
p4

#save as TIFF file with dimensions 1700 x 600
grid.arrange(p1, p2, p3, p4, nrow=2)

