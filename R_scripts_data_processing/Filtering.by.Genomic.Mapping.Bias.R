#load allele-specific read counts from genomic data for F1 hybrids
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Mapping.Bias.with.F1.DNA.counts/"
genomic.data <- read.table(paste0(path, "Genomic.read.counts.F1.tsv"), header = T)

# define function to identify genes that fail to be in either of 5% tails of the binomial distribution
# under the null hypothesis that there's no bias in the genomic data being provided
# NOTE: "FALSE" genes are the ones to be filtered out eventually
GenesWithPotentialBias1 <- function(geneRowNumber, data){
  gene.DGRPCount<- data[geneRowNumber, 1]  ### DGRP177 count should be column 1 of input data; SP159N is column 2
  gene.TotalCount<- data[geneRowNumber, 1] + data[geneRowNumber, 2] 
  if(gene.DGRPCount < gene.TotalCount*0.5) x<- pbinom(gene.DGRPCount, gene.TotalCount, 0.5) 
  else x<- pbinom(gene.TotalCount - gene.DGRPCount, gene.TotalCount, 0.5)
  y<-NA
  if(x < 0.05) y<-"FALSE" else y<-"TRUE"
    return(y)
}


# define function to identify genes where the DGRP-177 read counts deviate from expected proportion of 0.5
# by more than 0.1
# NOTE: "FALSE" genes are the ones to be filtered out eventually
GenesWithPotentialBias2 <- function(geneRowNumber, data){
  gene.DGRPfraction <- data[geneRowNumber, 1]/(data[geneRowNumber, 1] + data[geneRowNumber, 2])
  if(abs(0.5 - gene.DGRPfraction) > 0.1) x <- "FALSE" else x <-"TRUE"
      return(x)
}

##############################################################################
## run functions on counts from males, females and combined data##############
##############################################################################

# add results to a new dataframe that also has geneID and the count data
df.genomic <- genomic.data

## function for Criteria 1 - to be run on male, female and combined data
df.genomic$crit1.f <- sapply(1:(nrow(df.genomic)), GenesWithPotentialBias1, data = df.genomic[,2:3])
df.genomic$crit1.m <- sapply(1:(nrow(df.genomic)), GenesWithPotentialBias1, data = df.genomic[,4:5])
df.genomic$crit1.c <- sapply(1:(nrow(df.genomic)), GenesWithPotentialBias1, data = df.genomic[,6:7])

# make column that has '1' if criteria 1 is satisfied for all 3 datasets (male, female, combined), 0 if not
df.genomic$crit1.ALL <- with(df.genomic[,8:10],  crit1.f == "TRUE" & crit1.m == "TRUE" & crit1.c == "TRUE")
df.genomic$crit1.ALL[df.genomic$crit1and2 == "TRUE"] <- 1
df.genomic$crit1.ALL[df.genomic$crit1and2 == "FALSE"] <- 0


## function for Criteria 2 - to be run on male and female data only
df.genomic$crit2.f <- sapply(1:(nrow(df.genomic)), GenesWithPotentialBias2, data = df.genomic[,2:3])
df.genomic$crit2.m <- sapply(1:(nrow(df.genomic)), GenesWithPotentialBias2, data = df.genomic[,4:5])

# make column that has '1' if criteria 2 is satisfied for male & female datasets, 0 if not
df.genomic$crit2.ALL <- with(df.genomic[,12:13],  crit2.f == "TRUE" & crit2.m == "TRUE")
df.genomic$crit2.ALL[df.genomic$crit2.ALL == "TRUE"] <- 1
df.genomic$crit2.ALL[df.genomic$crit2.ALL == "FALSE"] <- 0

#make column that has '1' if both criteria 1 and 2 are satisfied for a gene; 0 if not
df.genomic$crit1and2 <- with(df.genomic[,11:14],  crit1.ALL == "1" & crit2.ALL == "1")
df.genomic$crit1and2[df.genomic$crit1and2 == "TRUE"] <- 1
df.genomic$crit1and2[df.genomic$crit1and2 == "FALSE"] <- 0

#################################################################################################
#save df.genomic to file
write.table(df.genomic, paste0(path, "Tests.for.Bias.in.F1.genomic.data.tsv"), sep = '\t', quote = F,
            row.names = F, col.names = T)
