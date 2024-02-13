#load list of genes without evidence of mapping bias; only the genes in this list will be analysed further
gene.list <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Genes.for.analysis.final.tsv", 
                        header = T)

#upload results from test for sex-reversals in allelic imbalance; filter out genes
SRAI.results <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Sex.Reversed.AI/Sex.Reversed.AI.test.tsv", header = T)
SRAI.results <- SRAI.results[SRAI.results$geneID %in% gene.list$geneID,]

#upload results from testing for sex*cross effects on AI in whole bodies of main and reciprocal crosses; filter out genes
#this was done to detect parental effects on AI, but now we need this to estimate false positives in sex-reversed AI
PE.results <- read.table("C:/Users/mishr/OneDrive/Desktop/ASE.New/Data/Mapping.Bias.Tests/Parental.Effects.on.AI.with.both.sexes/Test.for.parentalEffects.in.AI.combined.sexes.tsv", header=T)
PE.results <- PE.results[PE.results$geneID %in% gene.list$geneID,]

#genes with AI in males/females when testing for sex-reversed AI separately in the sexes
females.AI <- SRAI.results[which(SRAI.results$p.val.females < 0.05),] #nF = 1805
males.AI <- SRAI.results[which(SRAI.results$p.val.males < 0.05),] #nF = 1813

#genes with sex-dependent AI (SDAI) when testing for sex-dependent AI together with data from both sexes
SDAI.genes <- PE.results[PE.results$p.val_.quasibinomial._2 < 0.05,]

#genes with significant AI in females and significant SDAI
female.AI.SDAI <- females.AI[females.AI$geneID %in% SDAI.genes$geneID,] #nIntF = 824
male.AI.SDAI <- males.AI[males.AI$geneID %in% SDAI.genes$geneID,] #nIntF = 823

