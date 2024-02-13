library(ggplot2)
library(gridExtra)
library(ggpubr)

Dataframes.DGRP.Fraction.in.DNA.RNA <- function(tissue, sex){
  #load list of genes found to have no evidence of mapping bias
  gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                          header = T)
  
  #load allele-specific read count data for DNA, for given sex
  #only include genes where total read count > 30
  genomic.df <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/Genomic.read.counts/",
                          "DGRP177.SP159N.F1.", sex, ".tsv"), header = T)
  genomic.df <- genomic.df[genomic.df$Total.Count > 30,]
  
  
  #load dataframes made for analysing sex-effects on AI
  RNA.df <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/",
                              tissue, ".sexAI.df.tsv"), header = T, sep = '\t')
  
  #make a list of genes that have adequate data in both RNA and DNA; only include these genes in the final dataframe
  common.genes <- intersect(RNA.df$geneID, genomic.df$geneID)
  
  RNA.df <- RNA.df[RNA.df$geneID %in% common.genes,]
  genomic.df <- genomic.df[genomic.df$geneID %in% common.genes,]
  
  #DGRP-177 fraction in RNA
  DGRP.fraction.DNA <- genomic.df$Fraction.DGRP
  
  #DGRP-177 fraction in RNA
  DGRP.fraction.RNA <- 0
  if (sex == "male")
  {DGRP.fraction.RNA = RNA.df$Proportion.of.DGRP177.reads.Males} else {
    DGRP.fraction.RNA = RNA.df$Proportion.of.DGRP177.reads.Females
  }
  
  #presence of AI 
  RNA.df[RNA.df$Presence.of.AI == 1,]$Presence.of.AI <- "Genes with AI"
  RNA.df[RNA.df$Presence.of.AI == 0,]$Presence.of.AI <- "Genes without AI"
  
  
  
  #make a combined dataframe with geneID, DGRP-177 fraction in DNA, DGRP-177 fraction in RNA and presence/absence of AI
  df <- cbind.data.frame(geneID = RNA.df$geneID,
                         Proportion.of.DGRP177.reads.DNA = DGRP.fraction.DNA,
                         Proportion.of.DGRP177.reads.RNA = DGRP.fraction.RNA,
                         Presence.of.AI = RNA.df$Presence.of.AI)
  
  #filter above dataframe to only include genes without evidence of mapping bias
  df <- df[df$geneID %in% gene.list$geneID,]
  
  return(df)
}

f.gonads.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Gonads', 'female')
f.heads.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Heads', 'female')
f.whole.bodies.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Whole.bodies', 'female')
f.reciprocal.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Reciprocal.whole.bodies', 'female')

m.gonads.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Gonads', 'male')
m.heads.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Heads', 'male')
m.whole.bodies.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Whole.bodies', 'male')
m.reciprocal.df <- Dataframes.DGRP.Fraction.in.DNA.RNA('Reciprocal.whole.bodies', 'male')


#######################################################################################################
# make scatterplots of DGRP-177 coverage in DNA vs DGRP-177 expression for genes with and without AI ##
####################### FIGURE S2: FOR MALES #######################################################
m.gonads.Fig <- ggplot(m.gonads.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "A") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
m.gonads.Fig

m.heads.Fig <- ggplot(m.heads.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "B") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
m.heads.Fig

m.whole.bodies.Fig <- ggplot(m.whole.bodies.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "C") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
m.whole.bodies.Fig

m.reciprocal.Fig <- ggplot(m.reciprocal.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "D") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
m.reciprocal.Fig

#save as TIFF file with dimensions 900 x 570
ggarrange(m.gonads.Fig, m.heads.Fig, m.whole.bodies.Fig, m.reciprocal.Fig, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")




#######################################################################################################
# make scatterplots of DGRP-177 coverage in DNA vs DGRP-177 expression for genes with and without AI ##
####################### FIGURE S3: FOR FEMALES #######################################################
f.gonads.Fig <- ggplot(f.gonads.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "A") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
f.gonads.Fig

f.heads.Fig <- ggplot(f.heads.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "B") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
f.heads.Fig

f.whole.bodies.Fig <- ggplot(f.whole.bodies.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "C") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
f.whole.bodies.Fig

f.reciprocal.Fig <- ggplot(f.reciprocal.df, aes(Proportion.of.DGRP177.reads.DNA, Proportion.of.DGRP177.reads.RNA)) +
  geom_point(aes(color = Presence.of.AI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "D") + 
  scale_x_continuous(limits = c(0.35, 0.65)) +
  scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_abline(intercept = 0.5, slope = 0,  linewidth = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(.98,0.02),
        legend.justification = c("center", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        legend.text = element_text(color = "black", hjust = 0, size = 14),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text( face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) 
f.reciprocal.Fig

#save as TIFF file with dimensions 900 x 570
ggarrange(f.gonads.Fig, f.heads.Fig, f.whole.bodies.Fig, f.reciprocal.Fig, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


#############################################################################################
#correlations between coverage of DGRP-177 reads in DNA and DGRP-177 expression in RNA

#for genes with AI
#females
cor.test(f.gonads.df[which(f.gonads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.gonads.df[which(f.gonads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.heads.df[which(f.heads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.heads.df[which(f.heads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.whole.bodies.df[which(f.whole.bodies.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.whole.bodies.df[which(f.whole.bodies.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.reciprocal.df[which(f.reciprocal.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.reciprocal.df[which(f.reciprocal.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)

#males
cor.test(m.gonads.df[which(m.gonads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.gonads.df[which(m.gonads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.heads.df[which(m.heads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.heads.df[which(m.heads.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.whole.bodies.df[which(m.whole.bodies.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.whole.bodies.df[which(m.whole.bodies.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.reciprocal.df[which(m.reciprocal.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.reciprocal.df[which(m.reciprocal.df$Presence.of.AI == "Genes with AI"),]$Proportion.of.DGRP177.reads.RNA)


#for genes without AI
#females
cor.test(f.gonads.df[which(f.gonads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.gonads.df[which(f.gonads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.heads.df[which(f.heads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.heads.df[which(f.heads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.whole.bodies.df[which(f.whole.bodies.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.whole.bodies.df[which(f.whole.bodies.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(f.reciprocal.df[which(f.reciprocal.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         f.reciprocal.df[which(f.reciprocal.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)

#males
cor.test(m.gonads.df[which(m.gonads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.gonads.df[which(m.gonads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.heads.df[which(m.heads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.heads.df[which(m.heads.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.whole.bodies.df[which(m.whole.bodies.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.whole.bodies.df[which(m.whole.bodies.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)
cor.test(m.reciprocal.df[which(m.reciprocal.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.DNA, 
         m.reciprocal.df[which(m.reciprocal.df$Presence.of.AI == "Genes without AI"),]$Proportion.of.DGRP177.reads.RNA)


#for all genes
#females
cor.test(f.gonads.df$Proportion.of.DGRP177.reads.DNA, f.gonads.df$Proportion.of.DGRP177.reads.RNA)
cor.test(f.heads.df$Proportion.of.DGRP177.reads.DNA, f.heads.df$Proportion.of.DGRP177.reads.RNA)
cor.test(f.whole.bodies.df$Proportion.of.DGRP177.reads.DNA, f.whole.bodies.df$Proportion.of.DGRP177.reads.RNA)
cor.test(f.reciprocal.df$Proportion.of.DGRP177.reads.DNA, f.reciprocal.df$Proportion.of.DGRP177.reads.RNA)

#males
cor.test(m.gonads.df$Proportion.of.DGRP177.reads.DNA, m.gonads.df$Proportion.of.DGRP177.reads.RNA)
cor.test(m.heads.df$Proportion.of.DGRP177.reads.DNA, m.heads.df$Proportion.of.DGRP177.reads.RNA)
cor.test(m.whole.bodies.df$Proportion.of.DGRP177.reads.DNA, m.whole.bodies.df$Proportion.of.DGRP177.reads.RNA)
cor.test(m.reciprocal.df$Proportion.of.DGRP177.reads.DNA, m.reciprocal.df$Proportion.of.DGRP177.reads.RNA)

