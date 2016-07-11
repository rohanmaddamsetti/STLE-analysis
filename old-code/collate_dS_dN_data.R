## collate_dS_dN_data.R by Rohan Maddamsetti.
## This script makes the final csv table that is used by dS_dN_analysis.R.

## These variables should be changed when run.
#snpfile = "REL4397_dS_dN.csv"
#snpfile = "REL4398_dS_dN.csv"
#snpfile = "REL4398_dS_dN.csv"
snpfile <- "REL2545_dS_dN.csv"

possfile <- "possible_synonymous-REL606.csv"
#outfile = "REL4397_dS_dN_data.csv"
#outfile = "REL4398_dS_dN_data.csv"
#outfile = "REL4398_dS_dN_data.csv"
outfile <- "REL2545_dS_dN_data.csv"

## import data.
thetaS.data <- read.csv("Martincorena_Maddamsetti_thetaS_estimates.csv",header=T)

## this file is only useful for gene length information.
possible.synon.data <- read.csv(possfile, header=T)

gene.length.data <- subset(possible.synon.data, select=c("locus_tag","gene", "gene.length"))
snp.data <- read.csv(snpfile, header=T)

## change factors to character vectors (strings) so that they get merged
## properly.
thetaS.data$locus_tag <- sapply(thetaS.data$locus_tag,as.character)
gene.length.data$locus_tag <- sapply(gene.length.data$locus_tag,as.character)
gene.length.data$gene <- sapply(gene.length.data$gene,as.character)
snp.data$locus_tag <- sapply(snp.data$locus_tag,as.character)
snp.data$gene <- sapply(snp.data$gene,as.character)

## merge these data stepwise into a giant data frame.
murders <- merge(thetaS.data, gene.length.data, all=TRUE)

full.data <- merge(murders,snp.data,all=TRUE)
## NAs in dS and dN columns should be zeros.
full.data$dS <- sapply(full.data$dS, function(x) ifelse(is.na(x),0,x))
full.data$dN <- sapply(full.data$dN, function(x) ifelse(is.na(x),0,x))
## write these data to file.
write.csv(full.data,file=outfile, row.names=F)
