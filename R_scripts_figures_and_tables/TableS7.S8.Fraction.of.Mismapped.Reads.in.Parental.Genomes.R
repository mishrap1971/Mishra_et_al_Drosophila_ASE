#load list of genes found to have no evidence of mapping bias
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#load allele-specific counts obtained by competitive mapping parental genomic data to genotype-specific references
path <- "C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/HTSeq.Counts/Genomic.read.counts/"
DGRP177.counts <- read.table(paste0(path, 'DGRP177.Parent.tsv'), header = T)
SP159N.counts <- read.table(paste0(path, 'SP159N.Parent.tsv'), header = T)

###############################################################################
# percent mismapped reads in DGRP177 alignment
###############################################################################

#estimate total allele-specific reads and filter out genes where total read count < 30
DGRP177.counts$Total.Read.Count <- DGRP177.counts$DGRP.177.count + DGRP177.counts$SP.159N.count
DGRP177.counts <- DGRP177.counts[DGRP177.counts$Total.Read.Count >= 30,]

#estimate fraction of reads that are mapped to DGRP-177
DGRP177.counts$Fraction.DGRP.reads <- DGRP177.counts$DGRP.177.count/DGRP177.counts$Total.Read.Count

#percent of mismapped reads (prior to filtering out genes found to have mapping bias)
length(which(DGRP177.counts$Fraction.DGRP.reads >= 0.90))  #13974 genes; 99.8% have at least 90% reads correctly assigned
length(which(DGRP177.counts$Fraction.DGRP.reads >= 0.95))  #13910 genes; 99.3% have at least 90% reads correctly assigned
length(which(DGRP177.counts$Fraction.DGRP.reads >= 0.99))  #11136 genes; 79.5% have at least 90% reads correctly assigned

#percent of mismapped reads AFTER filtering out genes found to have mapping bias
filtered.DGRP177.counts <- DGRP177.counts[DGRP177.counts$geneID %in% gene.list$geneID,]
length(which(filtered.DGRP177.counts$Fraction.DGRP.reads >= 0.90))  #7624 genes; 99.9% have >=90% reads correctly assigned
length(which(filtered.DGRP177.counts$Fraction.DGRP.reads >= 0.95))  #7600 genes; 99.5% have >=90% reads correctly assigned
length(which(filtered.DGRP177.counts$Fraction.DGRP.reads >= 0.99))  #6117 genes; 80.1% have >=90% reads correctly assigned


###############################################################################
# percent mismapped reads in SP159N alignment
###############################################################################

#estimate total allele-specific reads and filter out genes where total read count < 30
SP159N.counts$Total.Read.Count <- SP159N.counts$DGRP.177.count + SP159N.counts$SP.159N.count
SP159N.counts <- SP159N.counts[SP159N.counts$Total.Read.Count >= 30,]

#estimate fraction of reads that are mapped to DGRP-177
SP159N.counts$Fraction.SP159N.reads <- SP159N.counts$SP.159N.count/SP159N.counts$Total.Read.Count

#percent of mismapped reads (prior to filtering out genes found to have mapping bias)
length(which(SP159N.counts$Fraction.SP159N.reads >= 0.90))  #14170 genes; 99.6% have at least 90% reads correctly assigned
length(which(SP159N.counts$Fraction.SP159N.reads >= 0.95))  #14066 genes; 99.1% have at least 90% reads correctly assigned
length(which(SP159N.counts$Fraction.SP159N.reads >= 0.99))  #11149 genes; 79.3% have at least 90% reads correctly assigned

#percent of mismapped reads AFTER filtering out genes found to have mapping bias
filtered.SP159N.counts <- SP159N.counts[SP159N.counts$geneID %in% gene.list$geneID,]
length(which(filtered.SP159N.counts$Fraction.SP159N.reads >= 0.90))  #7683 genes; 99.7% have >=90% reads correctly assigned
length(which(filtered.SP159N.counts$Fraction.SP159N.reads >= 0.95))  #7633 genes; 99.1% have >=90% reads correctly assigned
length(which(filtered.SP159N.counts$Fraction.SP159N.reads >= 0.99))  #6112 genes; 79.3% have >=90% reads correctly assigned



