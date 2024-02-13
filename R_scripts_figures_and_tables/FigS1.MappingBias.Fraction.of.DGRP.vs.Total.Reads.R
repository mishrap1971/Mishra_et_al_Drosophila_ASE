require("ggplot2")
require("gridExtra")
require("ggpubr")

setwd("C:/Users/Prashastcha Mishra/Desktop/ASE.Codes.Data/Data/Mapping.Bias/Count.Files")
read.counts.combined <- read.table('DGRP177.SP159N.F1.combined.filtered.tsv', header = T)
read.counts.male <- read.table('DGRP177.SP159N.F1.male.filtered.tsv', header = T)
read.counts.female <- read.table('DGRP177.SP159N.F1.female.filtered.tsv', header = T)

ColourBySuspicousness<-function(geneRowNumber, data){
  thisP1count<- data$DGRP.177[geneRowNumber]  ### V2 = Parent1Count
  thisTotalCount<- data$Total.Count[geneRowNumber] ### V4 = TotalCount
  if(thisP1count < thisTotalCount*0.5) x<- pbinom(thisP1count, thisTotalCount, 0.5) 
  else x<- pbinom(thisTotalCount - thisP1count, thisTotalCount, 0.5)
  y<-NA
  if(x < 0.05) y<-"red" else y<-"blue"
  return(y)
}

### Here I add a new column to my dataframe with the appropriate color code
read.counts.combined$ColorCode<-sapply(1:(dim(read.counts.combined)[1]), ColourBySuspicousness, data = read.counts.combined)
read.counts.male$ColorCode<-sapply(1:(dim(read.counts.male)[1]), ColourBySuspicousness, data = read.counts.male)
read.counts.female$ColorCode<-sapply(1:(dim(read.counts.female)[1]), ColourBySuspicousness, data = read.counts.female)

df.male <- cbind.data.frame(read.counts.male$Total.Count, read.counts.male$Fraction.DGRP, read.counts.male$ColorCode)
df.female <- cbind.data.frame(read.counts.female$Total.Count, read.counts.female$Fraction.DGRP, read.counts.female$ColorCode)
df.combined <- cbind.data.frame(read.counts.combined$Total.Count, read.counts.combined$Fraction.DGRP, read.counts.combined$ColorCode)
df.male$sex <- c(rep("Male", nrow(df.male)))
df.female$sex <- c(rep("Female", nrow(df.female)))
df.combined$sex <- c(rep("Combined", nrow(df.combined)))

colnames(df.male) <- c("TotalCount", "Fraction.DGRP.reads", "Colour", "Sex")
colnames(df.female) <- c("TotalCount", "Fraction.DGRP.reads", "Colour", "Sex")
colnames(df.combined) <- c("TotalCount", "Fraction.DGRP.reads", "Colour", "Sex")
df <- rbind.data.frame(df.male, df.female, df.combined)

#make the multi-panel plot
p <- ggplot(df, aes(x = TotalCount, y = Fraction.DGRP.reads, colour = Colour)) + 
  geom_point() + 
  scale_color_identity() +
  #theme_minimal() +
  theme(axis.text.x = element_text(size=16, margin = margin(5,0,0,0), color = "black"), 
        axis.text.y = element_text(size=16, margin = margin(0,5,0,0), color = "black"), 
        axis.title.y = element_text(size=21, margin = margin(0,0,0,0), color = "black"), 
        axis.title.x = element_text(size=21, margin = margin(0,0,0,0), color = "black"), 
        strip.text.y = element_text(size=17, face="bold"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.5) + 
  facet_grid(Sex ~ .) + 
  labs(x = paste0('\n', "Total Count"), 
       y = paste0("Fraction of DGRP-177 reads", '\n')) +
  xlim(0,4000)
p


