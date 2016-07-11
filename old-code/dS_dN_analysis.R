## dS_dN_analysis.R by Rohan Maddamsetti.

## This script analyses dN and dS changes in the endpoint Souza-Turner recombinants
## relative to K-12.

prepareData <- function(filename,filter=TRUE,use.maddamsetti=TRUE) {
  ## This function reads in the data from data.csv, and prepares it for
  ## analysis.

  ##import data.
  the.data <- read.csv(filename,header=T)

  if (use.maddamsetti) {
    the.data$thetaS <- the.data$Maddamsetti_thetaS
  } else {
    the.data$thetaS <- the.data$Martincorena_thetaS
  }

  if (filter) {
    ##remove entries with missing thetaS values and gene lengths.
    the.data <- the.data[!is.na(the.data$thetaS),]
    the.data <- the.data[!is.na(the.data$gene.length),]

  }
##sort the data frame by thetaS.

  the.data <- the.data[with(the.data,order(thetaS)),]
  
  return(the.data)
  
}

ks.analysis <- function (the.data) {
  ## do the following: 1) make a uniform cdf on mutation rate per base.
  ## 2) make a thetaS cdf. 3) make an empirical cdf of mutations per gene.
  ## do K-S tests for goodness of fit of the empirical cdf with the cdfs for
  ## the uniform cdf and thetaS cdf hypotheses.
  
  relevant.data <- subset(the.data, select=c("locus_tag","gene", "thetaS", "dS", "gene.length"))
  unique.relevant.data <- unique(relevant.data)

  genome.length <- sum(unique.relevant.data$gene.length)
  
  ## Calculate the empirical distribution of synonymous substitutions per gene.
  mutation.total <- sum(unique.relevant.data$dS)
  empirical.cdf <- cumsum(unique.relevant.data$dS/mutation.total)
  ## Null hypothesis: probability of a mutation per base is uniform.
  null.cdf <- cumsum(unique.relevant.data$gene.length/genome.length)
  ## Alternative hypothesis: mutation rate is proportional to thetaS.
  ## alternative 1: thetaS is the mutation rate per base pair.
  thetaS.is.per.bp.total <- sum(unique.relevant.data$thetaS*unique.relevant.data$gene.length)
  alt1.cdf <- cumsum(unique.relevant.data$thetaS*unique.relevant.data$gene.length/thetaS.is.per.bp.total)
  ## alternative 2: thetaS is the mutation rate per gene.
  thetaS.is.per.gene.total <- sum(unique.relevant.data$thetaS)
  alt2.cdf <- cumsum(unique.relevant.data$thetaS/thetaS.is.per.gene.total)
  ## Do Kolmogorov-Smirnov tests for goodness of fit.
  print(ks.test(empirical.cdf, null.cdf, simulate.p.value=TRUE))
  print(ks.test(empirical.cdf,alt1.cdf, simulate.p.value=TRUE))
  print(ks.test(empirical.cdf,alt2.cdf, simulate.p.value=TRUE))

  results.to.plot <- data.frame(locus_tag=unique.relevant.data$locus_tag, gene=unique.relevant.data$gene, thetaS=unique.relevant.data$thetaS, empirical=empirical.cdf,null=null.cdf,alt1=alt1.cdf,alt2=alt2.cdf)

  return(results.to.plot)
}

makeKSFigure <- function(the.results.to.plot) {
## This function generates the first panel that I want
## to use to illustrate how there is an association with thetaS with HGT,
## but not with synonymous substitution rates in the LTEE.
## Some post-processing with Illustrator will be needed (basically to add
## thetaS as a greek symbol on the axes).

  library(ggplot2)
  
  ## Plot thetaS across all loci.
  ## '0' in the axes is a placeholder for adding the greek letter theta in
  ## post-processing with Illustrator.
  
  ## for plotting convienence, add an index to the data frame.
  the.results.to.plot$index <- 1:length(the.results.to.plot$gene)
  
  plot <- ggplot(the.results.to.plot, aes(x=index)) +
    geom_line(aes(y=empirical), colour="red") + 
    geom_line(aes(y=null), linetype=2) + 
    geom_line(aes(y=alt1), linetype=3) +
    scale_x_continuous('Genes ranked by 0',limits=c(0,2900)) +
    scale_y_continuous('Cumulative proportion of synonymous mutations',limits=c(0,1)) +
      theme_classic() + theme(axis.title=element_text(size=18),axis.text=element_text(size=12))
  plot
  ggsave("../results/figure.pdf")
  ggsave("../results/figure.eps")
  
}

scatterplot.and.tables <- function(datafile, strain) {
data <- read.csv(datafile,header=T)

library(ggplot2)
library(reshape)

## Make scatterplots of omega values and dS.
quartz()
scat1 <- ggplot(data,aes(x=dS, y=Maddamsetti_thetaS)) + geom_point()
scat1

scat2 <- ggplot(data,aes(x=dS, y=Martincorena_thetaS)) + geom_point()
scat2

## Now make scatterplots of K-12 dS and recombinant dS (and dN).
K12.data <- read.csv("all-K-12_dS_dN.csv",header=T)
# rename the dS and dN columns.
K12.data <- rename(K12.data,c(dS="K-12 dS",dN="K-12 dN"))


comp.data <- merge(data,K12.data,all=TRUE)
## turn NAs in dS and dN to zeros.
comp.data$`K-12 dS` <- sapply(comp.data$`K-12 dS`,function(x) ifelse(is.na(x),0,x))
comp.data$`K-12 dN` <- sapply(comp.data$`K-12 dN`,function(x) ifelse(is.na(x),0,x))
comp.data$dS <- sapply(comp.data$dS,function(x) ifelse(is.na(x),0,x))
comp.data$dN <- sapply(comp.data$dN,function(x) ifelse(is.na(x),0,x))

scat3 <- ggplot(comp.data,aes(x=dS,y=`K-12 dS`)) + ylim(c(0,90)) + xlim(c(0,90)) + geom_point() 
scat3

scat4 <- ggplot(comp.data,aes(x=dN,y=`K-12 dN`)) + ylim(c(0,90)) + xlim(c(0,90)) + geom_point() 
scat4

## look at genes where K-12 dS != evolved dS, and same for dN.
comp.data$dSdiff <- comp.data$`K-12 dS` - comp.data$dS
comp.data$dNdiff <- comp.data$`K-12 dN` - comp.data$dN

dSdiff.data <- subset(comp.data,dSdiff != 0)
dSdiff.data <- dSdiff.data[order(dSdiff.data[,10]),] # sort by dS diff.

dNdiff.data <- subset(comp.data,dNdiff != 0)
dNdiff.data <- dNdiff.data[order(dNdiff.data[,11]),] # sort by dS diff.

dSfile <- paste(strain, "_dSdiff.csv",sep="") # REL4397_dSdiff.csv, for instance.
dNfile <- paste(strain, "_dNdiff.csv",sep="")

write.csv(dSdiff.data,file=dSfile, row.names=FALSE)
write.csv(dNdiff.data,file=dNfile, row.names=FALSE)


}

## analyze data for the recipient strain to give a null expectation.
scatterplot.and.tables("REL2545_dS_dN_data.csv",strain="REL2545")

scatterplot.and.tables("REL4397_dS_dN_data.csv",strain="REL4397")
scatterplot.and.tables("REL4398_dS_dN_data.csv",strain="REL4398")
## run on REL4398 as well. Note the FD equilibrium between REL4397 and REL4398,
## as well as difference in lac ability. also look at markers that differ between
## these strains, compared to K-12.


## Import the data, and get it ready for analysis.
ks.data <- prepareData(datafile,filter=TRUE,use.maddamsetti=TRUE)


## Do the K-S analysis that I did in my MBE paper.
## Note that the interpretation of these figures is subtly different.

results <- ks.analysis(ks.data)

#makeKSFigure(results)

## Do K-S test using Martincorena's thetaS estimates for robustness.
OG.data <- prepareData(datafile,filter=TRUE,use.maddamsetti=FALSE)
results2 <- ks.analysis(OG.data)
#makeKSFigure(results2)
