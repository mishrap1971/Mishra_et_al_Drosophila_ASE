Dataframes.Sex.Differences.in.DGRP.fraction <- function(tissue){
  #load list of genes found to have no evidence of mapping bias
  gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                          header = T)

  #load allele-specific read counts from F1 male and female genomic sequences
  #estimate % reads that are DGRP-177 for both males and females,
  #estimate sex difference in DGRP-177 coverage (females - males)
  df.genomic <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/Genomic.read.counts/Genomic.read.counts.F1.tsv",
                           header = T)
  DGRP.fraction.female.genomic <- df.genomic$female.DGRP177/(df.genomic$female.DGRP177 + df.genomic$female.SP159N)
  DGRP.fraction.male.genomic <- df.genomic$male.DGRP177/(df.genomic$male.DGRP177 + df.genomic$male.SP159N)
  sexDiff.genomic <- DGRP.fraction.female.genomic - DGRP.fraction.male.genomic
  df1 <- cbind.data.frame(geneID = df.genomic$geneID, Sex.Difference.in.DGRP.DNA = sexDiff.genomic)
  
  #load dataframes made for analysing sex-effects on AI
  #estimate sex differences in DGRP-177 expression (females - males)
  RNA.df <- read.table(paste0("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Dataframes.for.AI.analyses/Sex.dependent.AI/",
                              tissue, ".sexAI.df.tsv"), header = T, sep = '\t')
  sexDiff.RNA <- RNA.df$Proportion.of.DGRP177.reads.Females - RNA.df$Proportion.of.DGRP177.reads.Males
  
  #label gene if it has significant AI, significant sex-dependent AI or neither
  Presence.of.AI.or.SDAI <- 0
  for (i in 1:nrow(RNA.df)){
    if (RNA.df$Presence.of.AI[i] == "1" & RNA.df$Presence.of.sex.dependent.AI[i] == "0"){
      Presence.of.AI.or.SDAI[i] = "genes with AI, no SDAI"
    } else if (RNA.df$Presence.of.AI[i] == "0" & RNA.df$Presence.of.sex.dependent.AI[i] == "1"){
      Presence.of.AI.or.SDAI[i] = "genes with SDAI, no AI"
    } else if (RNA.df$Presence.of.AI[i] == "1" & RNA.df$Presence.of.sex.dependent.AI[i] == "1"){
      Presence.of.AI.or.SDAI[i] = "genes with AI and SDAI"
    } else {Presence.of.AI.or.SDAI[i] = "no AI or SDAI"}
  }
  
  df2 <- cbind.data.frame(geneID = RNA.df$geneID, Sex.Difference.in.DGRP.RNA = sexDiff.RNA,
                          Presence.of.AI.or.SDAI = Presence.of.AI.or.SDAI)
  
  #merge df1, df2 and the list of genes without mapping bias, such that only genes common to all three vectors are retained
  df <- na.omit(Reduce(function(...) merge(..., by='geneID', all.x=TRUE), list(gene.list, df1, df2)))
  
  return(df)  
}

gonads.df <- Dataframes.Sex.Differences.in.DGRP.fraction('gonads')
heads.df <- Dataframes.Sex.Differences.in.DGRP.fraction('heads')
whole.bodies.df <- Dataframes.Sex.Differences.in.DGRP.fraction('whole.bodies')
reciprocal.df <- Dataframes.Sex.Differences.in.DGRP.fraction('reciprocal.whole.bodies')

###################################################################################################
#make scatterplots of sex difference in DGRP-177 coverage in DNA and DGRP-177 expression in RNA ###
###################################################################################################
gonads.Fig <- ggplot(gonads.df, aes(Sex.Difference.in.DGRP.DNA, Sex.Difference.in.DGRP.RNA)) +
  geom_point(aes(color = Presence.of.AI.or.SDAI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "", 
       tag = "A") + 
  scale_x_continuous(limits = c(-0.15, 0.2)) +
  scale_y_continuous(limits = c(-0.5, 0.4)) +
  geom_abline(intercept = 0, slope = 0,  size = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.text = element_text(color = "black", hjust = 0, size = 14), 
        legend.title = element_blank(),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) +
  guides(colour = guide_legend(override.aes = list(size=4)))
gonads.Fig


heads.Fig <- ggplot(heads.df, aes(Sex.Difference.in.DGRP.DNA, Sex.Difference.in.DGRP.RNA)) +
  geom_point(aes(color = Presence.of.AI.or.SDAI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "",  
       tag = "B") + 
  scale_x_continuous(limits = c(-0.15, 0.2)) +
  scale_y_continuous(limits = c(-0.5, 0.4)) +
  geom_abline(intercept = 0, slope = 0,  size = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.text = element_text(color = "black", hjust = 0, size = 14), 
        legend.title = element_blank(),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) +
  guides(colour = guide_legend(override.aes = list(size=4)))
heads.Fig


whole.bodies.Fig <- ggplot(whole.bodies.df, aes(Sex.Difference.in.DGRP.DNA, Sex.Difference.in.DGRP.RNA)) +
  geom_point(aes(color = Presence.of.AI.or.SDAI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "Sex difference in DGRP-177 genomic coverage", 
       y = "Sex difference in DGRP-177 expression", 
       tag = "C") + 
  scale_x_continuous(limits = c(-0.15, 0.2)) +
  scale_y_continuous(limits = c(-0.5, 0.4)) +
  geom_abline(intercept = 0, slope = 0,  size = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.text = element_text(color = "black", hjust = 0, size = 14), 
        legend.title = element_blank(),
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) +
  guides(colour = guide_legend(override.aes = list(size=4)))
whole.bodies.Fig 

reciprocal.Fig <- ggplot(reciprocal.df, aes(Sex.Difference.in.DGRP.DNA, Sex.Difference.in.DGRP.RNA)) +
  geom_point(aes(color = Presence.of.AI.or.SDAI), show.legend = T) + 
  scale_color_discrete() +
  labs(x = "", 
       y = "",  
       tag = "D") + 
  scale_x_continuous(limits = c(-0.15, 0.2)) +
  scale_y_continuous(limits = c(-0.5, 0.4)) +
  geom_abline(intercept = 0, slope = 0,  size = 1, linetype="solid", color = "grey") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        #legend.position = c(.98,0.02),
        legend.justification = c("right", "bottom"),
        #legend.box.just = "left",
        legend.box.background = element_rect(),
        #legend.box.margin = margin(0.25, 1, 1, 1),
        legend.text = element_text(color = "black", hjust = 0, size = 14), 
        axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 13),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 13), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) +
  guides(colour = guide_legend(override.aes = list(size=4)))
reciprocal.Fig

#save as TIFF file with dimensions 950 x 570
ggarrange(gonads.Fig, heads.Fig, whole.bodies.Fig, reciprocal.Fig, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#####################################################################################################################
#######  estimate correlations of sex difference in DGRP-177 in DNA and RNA for Table S6  ###########################
#####################################################################################################################

#for all genes
cor.test(gonads.df$Sex.Difference.in.DGRP.DNA, gonads.df$Sex.Difference.in.DGRP.RNA)
cor.test(heads.df$Sex.Difference.in.DGRP.DNA, heads.df$Sex.Difference.in.DGRP.RNA)
cor.test(whole.bodies.df$Sex.Difference.in.DGRP.DNA, whole.bodies.df$Sex.Difference.in.DGRP.RNA)
cor.test(reciprocal.df$Sex.Difference.in.DGRP.DNA, reciprocal.df$Sex.Difference.in.DGRP.RNA)

#for genes with AI, but no sex-dependent AI
cor.test(gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.DNA,
         gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.DNA,
         heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.DNA,
         whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.DNA,
         reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with AI, no SDAI'),]$Sex.Difference.in.DGRP.RNA)

#for genes with sex-dependent AI, but no AI
cor.test(gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.DNA,
         gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.DNA,
         heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.DNA,
         whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.DNA,
         reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with SDAI, no AI'),]$Sex.Difference.in.DGRP.RNA)


#for genes with AI and SDAI
cor.test(gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.DNA,
         gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.DNA,
         heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.DNA,
         whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.DNA,
         reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'genes with AI and SDAI'),]$Sex.Difference.in.DGRP.RNA)

#for genes with sex-dependent AI, and AI
cor.test(gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.DNA,
         gonads.df[which(gonads.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.DNA,
         heads.df[which(heads.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.DNA,
         whole.bodies.df[which(whole.bodies.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.RNA)
cor.test(reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.DNA,
         reciprocal.df[which(reciprocal.df$Presence.of.AI.or.SDAI == 'no AI or SDAI'),]$Sex.Difference.in.DGRP.RNA)
