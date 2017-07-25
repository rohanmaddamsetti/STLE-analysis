## dissertation-analysis.R by Rohan Maddamsetti.

## This script makes figures and does statistics for the recombinant genome analysis.

## go through these imports and figure out which ones are superfluous.
library(xml2)
library(roxygen2)
library(assertthat)
library(IRanges)
library(GenomicRanges)
library(genbankr)

library(purrr)       ## consistent & safe list/vector munging
library(tidyr)       ## consistent data.frame cleaning
library(ggplot2)     ## base plots are for Coursera professors
library(scales)      ## pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   ## a helper for arranging individual ggplot objects
library(ggthemes)    ## has a clean theme for ggplot2
library(viridis)
library(DT)          ## prettier data.frame output
library(data.table)  ## faster fread()
library(dplyr)       ## consistent data.frame operations.
library(dtplyr)      ## dplyr works with data.table now.
library(ggrepel)     ## plot labeled scatterplots.
library(splines)     ## for Figure 3.
library(cowplot)     ## for Figure 6.
library(zoo)

## A reminder to what the labels are.
## 0) reference genome state.
## 1) B/K-12 marker is present that is not in reference genome (yellow)
## 2) LTEE recipient mutations (light blue)
## 3) new mutations (black)
## 4) LTEE recipient mutations that were replaced by K-12 or otherwise missing (red)
## 5) deleted markers (marker falls in a deleted region) (light purple)
## 6) REL288 specific marker
## 7) REL291 specific marker
## 8) REL296 specific marker
## 9) REL298 specific marker

#' function to rotate genome coordinates, setting oriC at the center of plots.
rotate.chr <- function(my.position,genome='REL606') {
    if (genome == 'REL606') {
        GENOME.LENGTH <- 4629812
        ORIC <- 3886105
    } else if (genome == 'K-12') {
        GENOME.LENGTH <- 4641652
        ORIC <- 3925744
    }

    midpoint <- GENOME.LENGTH/2
    L <- ORIC - midpoint
    ifelse(my.position > L,my.position-ORIC,GENOME.LENGTH-ORIC+my.position)
}

#############
## Import data.

## metadata for sequenced samples.
pop.clone.file <- file.path("../doc/Populations-and-Clones.csv")
pop.and.clone.metadata <- read.csv(pop.clone.file)

## For each dataset, add rotated genome coordinates so that oriC is at the center of the chromosome in plots.
hand.annotation <- read.csv("../results/donor_hand_annotation.csv") %>%
    mutate(rotated.position=rotate.chr(position))

##annotate auxotrophs and Hfr oriTs on Fig2 type plots.
auxotrophs <- filter(hand.annotation,annotation!="oriT"&annotation!="oriC")
hfrs <- filter(hand.annotation,annotation=="oriT")

G.score.data <- tbl_df(read.csv("../results/036806-3.csv")) %>%
    mutate(rotated.Start.position=rotate.chr(Start.position))

## clone sequencing data.
labeled.mutations <- tbl_df(read.csv("../results/labeled_mutations.csv")) %>%
    mutate(lbl=as.factor(lbl)) %>%
    mutate(rotated.position=rotate.chr(position))

## label clones as odd or even: BASED ON REL NUMBERING, NOT RM NUMBERING!
## NOTE THAT RM11734 (even) is RM3-130-1!
genome.names <- levels(labeled.mutations$genome)
is.odd <- sapply(genome.names, function(x) ifelse(strtoi(substr(x,nchar(x),nchar(x)))%%2, TRUE,FALSE))
labeled.mutations <- mutate(labeled.mutations,odd=is.odd[genome])

## evolution experiment data.
evoexp.labeled.mutations <- tbl_df(read.csv("../results/evolution-experiment/evoexp_labeled_mutations.csv")) %>%
    mutate(lbl=as.factor(lbl)) %>% mutate(rotated.position=rotate.chr(position))

## differences between K-12 and REL606.
K12.diff.data <- tbl_df(read.csv("../results/K-12-differences.csv")) %>%
    mutate(rotated.position=rotate.chr(position))

## get F-plasmid coverage data for clones and for evolution experiment.
STLE.clone.F.coverage <- tbl_df(read.csv("../results/STLE-clone-F-coverage.csv"))
STLE.evoexp.F.coverage <- tbl_df(read.csv("../results/STLE-evoexp-F-coverage.csv")) %>%
            mutate(Lineage=factor(Lineage,
                                  levels=c('Ara+1','Ara+2','Ara+3','Ara+4','Ara+5','Ara+6','Ara-1','Ara-2','Ara-3','Ara-4','Ara-5','Ara-6')))

## get labeled mutations for clones in Turner 1996 Ecology paper.
turner.clones <- tbl_df(read.csv("../results/turner_clones_labeled_mutations.csv")) %>%
    mutate(lbl=as.factor(lbl)) %>%
    mutate(rotated.position=rotate.chr(position))


######
## First, separate into odd or even clones, and
## This is used in a lot of the following code (Fig S2, Fig 2, Fig 4 or something like that.)
odd.genomes <- filter(labeled.mutations, odd==TRUE)
even.genomes <- filter(labeled.mutations, odd==FALSE)

## This x-axis label is re-used a lot.
oriC.xlab <- expression(paste("Distance from ",italic("oriC")))

##################################################################################
## Make Figures and Tables.

## Figure S1. Density of differences between K-12 and REL606.

FigS1 <- ggplot(K12.diff.data, aes(x=rotated.position)) + geom_histogram(bins=556) + theme_tufte() + xlab(oriC.xlab) + ylab("Differences between K-12 and REL606")
ggsave("/Users/Rohandinho/Desktop/FigS1.pdf",FigS1,height=3,width=4)


#############################
## Table 2: Probable beneficial mutations found in STLE recipients.

## Calculation that 32 genes that are in Figure 1 and in Table 2.
## Take top 57 G-scoring genes with 2 or more dN in non-mutators (Tenaillon et al. 2016).
top.G.score.genes <- filter(G.score.data,Observed.nonsynonymous.mutation>1)
top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name)
recipient.top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4')
genes.to.label <- filter(top.G.score.genes,Gene.name %in% recipient.top.hits$gene.annotation)
## gene.to.label contains 32 genes.

Table2.df <- mutate(genes.to.label,Gene = Gene.name) %>%
    select(Gene,Start.position,Coding.length,G.score)

write.csv(file="~/Desktop/Table2.csv",Table2.df)

#####################################################################################
## Make Figure 1 and Figure 8 (LCA for STLE continuation experiment analysis).

Fig1.function <- function(mut.df, genes.to.label, analysis.type='Fig1') {
    ## REL606 markers are negative space on the plot.
    ## also, reverse levels to get Ara+1 genomes on top of plot.
    ## ylim from 0 to 130. plot 12 genomes/lineages on 10,20,30,40,50,60,70,80,90,110,120. so map genome/lineage to a center position.

    ## for max ylim.
    y.pad <- 10
    y.max <- 13 * y.pad ## 130

    no.B.genomes <- filter(mut.df,lbl!='0') %>%
        mutate(lineage=factor(lineage,
                              levels=c('Ara+1','Ara+2','Ara+3','Ara+4','Ara+5','Ara+6','Ara-1','Ara-2','Ara-3','Ara-4','Ara-5','Ara-6'))) %>%
        arrange(lineage,position) %>%
        mutate(y.center=y.max-match(lineage,levels(lineage))*y.pad,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3) %>%
        ## turn lbls for donor-specific markers into K-12 lbls.
        mutate(lbl = replace(lbl,lbl %in% c('6','7','8','9'),'1')) %>%
        ## then turn into a factor.
        mutate(lbl=factor(lbl))

    ## label the top 32 beneficial LTEE genes with mutations in the recipients with symbols.
    top.hit.labels <- filter(no.B.genomes, gene.annotation %in% genes.to.label$Gene.name) %>%
        filter(lbl %in% c(2,3,4)) %>%
        mutate(y.center=y.max-match(lineage,unique(lineage))*y.pad,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3)

    ## reorder top.hit.labels$mut.annotation factor levels for nice plotting.
    top.hit.labels$mut.annotation <- factor(top.hit.labels$mut.annotation, levels = c("dN", "dS", "indel", "IS-insertion", "base-substitution", "non-coding", "non-point"))

    Fig1.xlab <- expression(paste("Distance from ",italic("oriC")))
    ## all donor mutations should be labeled as yellow.
    panel <- ggplot(no.B.genomes, aes(x=x1,xend=x2,y=y1,yend=y2,colour=lbl)) +
        geom_segment(size=0.05) +
        theme_tufte() +
        xlab(Fig1.xlab) +
        ylim(5,y.max-5) +
        geom_point(data=top.hit.labels,aes(x=rotated.position,
                                          y=y.center,
                                          color=lbl,
                                          shape=mut.annotation)) +
        ## set the symbols this way:
        ## dN:open circle, dS:open square, indel:open triangle, IS-insertion: an 'x'.
        scale_shape_manual(values = c(1,0,2,4),name="symbol") +
        theme(legend.position = "top") +
        geom_label_repel(data=filter(top.hit.labels,lbl %in% c(3,4)),
                         aes(x=rotated.position,
                             y=y.center,
                             fill=lbl,
                             label=gene.annotation),
                         fontface='italic',
                         color='white',
                         box.padding = unit(0.35, "lines"),
                         point.padding = unit(0.5, "lines"),
                         segment.color = 'grey50',
                         size=2.5) +
        scale_fill_manual(values=c('black','#d7191c'),
                          labels=c("new","replaced")) +
        guides(fill=FALSE,shape=FALSE,color=FALSE) +
        scale_colour_manual(values=c('#ffffbf', '#abd9e9', 'black','#d7191c','#f1b6da'),
                            name='color',
                            labels=c("K-12","LTEE","new","replaced","deleted")) +
        scale_y_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100,110,120),labels=rev(levels(no.B.genomes$lineage)))

    if (analysis.type == 'Fig1') {
        panel <- panel + ylab("Odd-numbered Clone")
    } else if (analysis.type == 'FigS2') {
        panel <- panel + ylab("Even-numbered Clone")
    } else if (analysis.type == 'evoexpLCA') {
        panel <- panel + ylab("Inferred LCA of Population")
    }
    return(panel)
}


## make Figure 1 and S2.
Fig1 <- Fig1.function(odd.genomes, genes.to.label, analysis.type='Fig1')
ggsave("/Users/Rohandinho/Desktop/Fig1.pdf", Fig1,width=8,height=10)

FigS2 <- Fig1.function(even.genomes, genes.to.label, analysis.type='FigS2')
ggsave("/Users/Rohandinho/Desktop/FigS2.pdf", FigS2,width=8,height=10)


##############################################################################
## Figures 2, 6, and S3.
donor.mutations <- filter(labeled.mutations,lbl == 6 | lbl==7 | lbl == 8 | lbl == 9) %>%
    mutate(genome=factor(paste(lineage,genome,sep=': '))) %>%
    ## reorder factor for plotting.
    mutate(genome=factor(genome,levels=levels(genome)[c(12:21,1:11)])) %>%
    select(-frequency) %>%
    arrange(genome,position) %>%
    distinct(genome, position, .keep_all = TRUE)

score.introgression <- function(independent.genomes) {
    ##initialize introgression.score column.
    independent.genomes <- mutate(independent.genomes, introgression.score = 0)

    ##This is a little helper to sum up labels (which is a factor)
    label.levels <- levels(independent.genomes$lbl)
    K12.labels <- c(1,4,6,7,8,9)
    labels.to.int <- sapply(label.levels, function(x) ifelse(x %in% K12.labels,1,0))

    ##take all positions labeled with 1, and sum them up to get the introgression score.
    ## omit Ara-3 (mostly K-12).
    introgression.scores <- filter(independent.genomes,lineage != "Ara-3") %>%
        group_by(rotated.position) %>%
        summarize(introgression.score = sum(labels.to.int[lbl])) %>%
        left_join(distinct(select(independent.genomes,rotated.position,gene.annotation))) %>%
        arrange(desc(introgression.score))

return(introgression.scores)
}

get.introgressed.genes <- function(introgression.scores) {
    select(introgression.scores,-rotated.position) %>% distinct()
}

makeFig2A <- function(scores, auxotrophs, hfrs) {

    oriC.xlab <- expression(paste("Distance from ",italic("oriC")))

    Fig2A <- ggplot(scores, aes(x=rotated.position,y=introgression.score)) +
        geom_line(size=0.05) +
        ## fit a natural cubic spline with 100 degrees of freedom.
        geom_smooth(method='lm', formula = y ~ ns(x,100)) +
        theme_tufte() +
        xlab(oriC.xlab) +
        ylab("Introgression Score") +
        ylim(-1,11.5) +
        ## add auxotroph lines
        geom_vline(data=auxotrophs,
                   aes(xintercept=rotated.position,
                       color=Donor.strain),
                   size=0.5,
                   linetype="dashed") +
        geom_label_repel(data=auxotrophs,
                         aes(x=rotated.position,
                             y=10,
                             label=annotation,
                             fill=Donor.strain,
                             fontface="italic"),
                         inherit.aes=FALSE,
                         box.padding = unit(0.35, "lines"),
                         point.padding = unit(0.25, "lines"),
                         arrow = arrow(length = unit(0.01, 'npc')),
                         segment.color = 'grey50',
                         color='white',
                         size=2.5) +
        ## add annotation of Hfr and oriC on chromosome.
        geom_label_repel(data=hfrs,
                         aes(x = rotated.position,
                             y = -0.65,
                             label = Hfr.orientation,
                             fill = Donor.strain),
                         inherit.aes = FALSE,
                         size = 2.5,
                         point.padding = unit(0, "lines"),
                         color = 'white') +
        guides(fill=FALSE,color=FALSE)
}

odd.introgression.scores <- score.introgression(odd.genomes)

## Figure 2A.
## Visual comparisons of parallel recombination across lineages
## with the G-score of mutations over the genome, and with the occurrence of
## new mutations over the genome.
## This is only made with odd genomes at the moment.

Fig2A <- makeFig2A(odd.introgression.scores,auxotrophs,hfrs)
ggsave("/Users/Rohandinho/Desktop/Fig3A.pdf",Fig2A,height=2.5,width=6)
## NOTE: The delta does not print properly but can fix in Illustrator.

## which genes get introgressed the most? is this a sign of positive selection?
odd.introgressed.genes <- get.introgressed.genes(odd.introgression.scores)

############# Fig. 2B: location of donor specific mutations on the chromosome.

donor.lbl.map <- function (vec) sapply(vec, function(lbl) {
    if (lbl == 6) {
        return('REL288')
    } else if (lbl == 7) {
        return('REL291')
    } else if (lbl == 8) {
        return('REL296')
    } else if (lbl == 9) {
        return('REL298')
    } else {
        return('NA')
    }
} )

only.donor.mutations <- donor.mutations %>%
    mutate(y.center=5,y1=y.center-3,y2=y.center+3) %>%
    droplevels() %>% filter(odd == TRUE) %>% mutate(Donor.strain=donor.lbl.map(lbl))

Fig2B <- ggplot(only.donor.mutations,aes(x=rotated.position,xend=rotated.position,y=y1,yend=y2,colour=Donor.strain)) +
    geom_segment(size=0.05) +
    theme_tufte() +
    xlab(oriC.xlab) +
    theme(axis.title.y = element_text(color = "white"),
          axis.text.y = element_text(color='white'),
          axis.ticks.y=element_blank()) +
    ## add auxotroph lines.
    geom_vline(data=auxotrophs,aes(xintercept=rotated.position,color=Donor.strain),size=0.5,linetype="dashed") +
    guides(fill=FALSE,color=FALSE) +
    ## annotate Hfr oriT on chromosome.
    geom_label_repel(data=hfrs,aes(x=rotated.position,y=-2,label=Hfr.orientation,
                                   fill=Donor.strain),inherit.aes=FALSE,size=2.5,
                    point.padding = unit(0,"lines"),color = 'white')

ggsave("/Users/Rohandinho/Desktop/Fig2B.pdf",Fig2B,height=1,width=6)

########################################
## Figure S3: 'Figure 1' for the Turner clones.

FigS3.function <- function(mut.df, genes.to.label) {
    ## REL606 markers are negative space on the plot.
    ## also, reverse levels to get Ara+1 genomes on top of plot.
    ## ylim from 0 to 30. plot 2 genomes/lineages on 10,20. so map genome/lineage to a center position.

    ## for max ylim.
    y.pad <- 10
    y.max <- 3 * y.pad ## 30

    no.B.genomes <- filter(mut.df,lbl!='0') %>%
                mutate(genome=factor(genome,
                              levels=c('REL4397','REL4398'))) %>%
        arrange(genome,position) %>%
        mutate(y.center=y.max-match(genome,levels(genome))*y.pad,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3) %>%
        ## turn lbls for donor-specific markers into K-12 lbls.
        mutate(lbl = replace(lbl,lbl %in% c('6','7','8','9'),'1')) %>%
        ## then turn into a factor.
        mutate(lbl=factor(lbl))

    ## label the top 32 beneficial LTEE genes with mutations in the recipients with symbols.
    top.hit.labels <- filter(no.B.genomes, gene.annotation %in% genes.to.label$Gene.name) %>%
        filter(lbl %in% c(2,3,4)) %>%
        mutate(y.center=y.max-match(genome,unique(genome))*y.pad,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3)

    ## reorder top.hit.labels$mut.annotation factor levels for nice plotting.
    top.hit.labels$mut.annotation <- factor(top.hit.labels$mut.annotation, levels = c("dN", "dS", "indel", "IS-insertion", "base-substitution", "non-coding", "non-point"))

    Fig.xlab <- expression(paste("Distance from ",italic("oriC")))

    ## all donor mutations should be labeled as yellow.
    panel <- ggplot(no.B.genomes, aes(x=x1,xend=x2,y=y1,yend=y2,colour=lbl)) +
        geom_segment(size=0.05) +
        theme_tufte() +
        xlab(Fig.xlab) +
        ylim(5,y.max-5) +
        geom_point(data=top.hit.labels,aes(x=rotated.position,
                                          y=y.center,
                                          color=lbl,
                                          shape=mut.annotation)) +
        ## set the symbols this way:
        ## dN:open circle, dS:open square, indel:open triangle, IS-insertion: an 'x'.
        scale_shape_manual(values = c(1,0,2,4),name="symbol") +
        theme(legend.position = "top") +
        geom_label_repel(data=filter(top.hit.labels,lbl %in% c(3,4)),
                         aes(x=rotated.position,
                             y=y.center,
                             fill=lbl,
                             label=gene.annotation),
                         fontface='italic',
                         color='white',
                         box.padding = unit(0.35, "lines"),
                         point.padding = unit(0.5, "lines"),
                         segment.color = 'grey50',
                         size=2.5) +
        scale_fill_manual(values=c('#d7191c'),
                          labels=c("replaced")) +
        guides(fill=FALSE,color=FALSE,shape=FALSE) +
        scale_colour_manual(values=c('#ffffbf', 'black','#d7191c','#f1b6da'),
                            name='color',
                            labels=c("K-12","new","replaced","deleted")) +
        scale_y_continuous(breaks=c(10,20),labels=rev(levels(no.B.genomes$genome))) +
        ylab("Ara-3 Clone")
    return(panel)
}

## make Figure S3. Turner clones.
FigS3 <- FigS3.function(turner.clones, genes.to.label)
ggsave("/Users/Rohandinho/Desktop/FigS3.pdf", FigS3,width=8,height=3)

REL4397.data <- turner.clones %>% select(-frequency) %>% filter(genome=='REL4397') %>%
    select(-genome,-lineage,-reference)
REL4398.data <- turner.clones %>% select(-frequency) %>% filter(genome=='REL4398') %>%
    select(-genome,-lineage,-reference)

## find mutations specific to each clone.
REL4397.diffs <- setdiff(REL4397.data,REL4398.data)
REL4398.diffs <- setdiff(REL4398.data,REL4397.data)

REL4397.diffs2 <- REL4397.diffs %>% filter(mut.annotation!='dS',mut.annotation!='non-coding')
REL4398.diffs2 <- REL4398.diffs %>% filter(mut.annotation!='dS',mut.annotation!='non-coding')
########################################
## Fig. S4: Does sequence divergence/conservation predict location of markers? NO.
## REMEMBER: divergence and introgression is correlated-- since when there's no divergence,
## all introgression events are undetectable! So, more divergence = better resolution of introgression.
## this is extra clear when regressing sum(introgression) against divergence in the windows.

REL606.LENGTH <- 4629812
## use 556 bins. Each is 8327 bp long (integer divisors)
bin.length <- REL606.LENGTH/556

## BUG: yegZ dS mutation has NA introgression score; there are two mutations at the same position in yegZ in annotated_K-12.gd,
## probably due to conflicting donor-specific mutations. This is very minor (only seems to affect one mutation)
## so fix this later, if at all.

introgression.vs.divergence.data <- left_join(K12.diff.data,odd.introgression.scores) %>%
    ## bin mutations across the genome.
    ## NOTE: equal size bins, NOT equal numbers of mutations!
    mutate(my_bin=ceiling(position/bin.length)) %>% group_by(my_bin)

two.classed.introgression.data <- introgression.vs.divergence.data %>%
    summarize(introgression=ifelse(mean(introgression.score)>0,1,0),divergence=n())

## no difference in divergence between regions with and without introgression. p = 0.8756.
kruskal.test(divergence ~ introgression,data=two.classed.introgression.data)

median.introgression.data <- introgression.vs.divergence.data %>%
    summarize(median.introgression=median(introgression.score),divergence=n()) %>% filter(median.introgression > 0)

introgression.divergence.model <- lm(median.introgression~divergence,data=median.introgression.data)
confint(introgression.divergence.model)

FigS4 <- ggplot(median.introgression.data,aes(x=divergence, y=median.introgression)) +
    geom_jitter() + ylab("Median Introgression within Bin") + xlab("K-12 Differences within Bin") +
    theme_tufte()

ggsave("/Users/Rohandinho/Desktop/FigS4.pdf",FigS4,width=5,height=4)


####################################################################################################
## Replaced mutations. Table 3 and Figure 3 (Fig. 3 is made separately).

#### Make in-depth alignments of all replaced mutations, with special attention to Ara+1 and Ara-4.
#### Print out a csv of genes for align_replaced_mutations.py to align,
#### and annotate as 0) reversion to pre-LTEE state, 1) K-12 state, 3) new allele.

#### Only look at dN mutations in odd REL clones.
replaced.gene.list <- filter(labeled.mutations,lbl==4,mut.annotation=='dN',odd==TRUE) %>%
    #select(gene.annotation,lineage,genome) %>%
    group_by(lineage,genome) %>% distinct(gene.annotation)
write.csv(replaced.gene.list,"../results/align_these.csv",row.names=FALSE,quote=FALSE)

## Run pythonw align_replaced_mutations.py to get results for Table 3 in the paper,
## as well as the numbers in this section of the manuscript.
## 31 of the 61 replaced dN are in non-mutators, and these are shown in Table 3.

################################
## Table 4: Putative gene conversion events.
## Started by looking for parallelism in new mutations: super strong parallelism (multiple new mutations in the same
## gene is probably gene conversion or something.

## omit mutator lineages (Ara+6,Ara+3,Ara-2) +6 is mutT, +3 is mutS, -2 is mutL mutator.
new.odd.mutations <- filter(odd.genomes,lineage != "Ara+6" & lineage != "Ara+3" & lineage != "Ara-2") %>% filter(lbl=='3')

## First: filter out cases when 3 or more new mutations occur in the same gene in the same odd-numbered genome
## (most likely not a new mutation).

gene.convs <- group_by(new.odd.mutations,lineage) %>% group_by(gene.annotation) %>% filter(n()>=3)
rest.new.odd.mutations <- group_by(new.odd.mutations,lineage) %>% group_by(gene.annotation) %>% filter(n()<3)

gene.convs.summary <- group_by(gene.convs,gene.annotation) %>% summarize(uniq.lineage=length(unique(lineage)),mutcount=length(gene.annotation) ,total.pos=length(unique(position)))

## There seems to be a relationship between number of new mutations...and number of deleted LTEE markers?
## remember that this is only looking at non-mutators.
rest.new.muts.summary <- group_by(rest.new.odd.mutations,gene.annotation) %>% summarize(lineages=length(unique(lineage)),mutcount=length(gene.annotation),total.pos=length(unique(position)))

rest.new.muts.lineage.summary <- group_by(rest.new.odd.mutations,lineage) %>% summarize(mutcount=n())

deleted.odd.mutations <- filter(odd.genomes,lineage != "Ara+6" & lineage != "Ara+3" & lineage != "Ara-2") %>% filter(lbl=='4')

deleted.muts.lineage.summary <- group_by(deleted.odd.mutations,lineage) %>% summarize(deletedcount=n())

## can't conclude that new mutations due to mutagenesis--strong signal of selection!!
rest.new.odd.G.score <- filter(G.score.data, Gene.name %in% rest.new.odd.mutations$gene.annotation)
rest.new.odd.G.score.summary <- summarize(rest.new.odd.G.score,mean.G.score=mean(G.score))
## mean G-score is 13!
gene.convs.G.score <- filter(G.score.data, Gene.name %in% gene.convs$gene.annotation)
gene.convs.G.score.summary <- summarize(gene.convs.G.score,mean.G.score=mean(G.score))
## while mean G-score is 0.42 for the 'gene conversion events'.

t.test(rest.new.odd.G.score$G.score,gene.convs.G.score$G.score)

## how many new mutations occur in top G.score genes?
top.G.score.genes <- filter(G.score.data,G.score>8.9)
top.new.odd.mutations <- filter(new.odd.mutations,gene.annotation %in% top.G.score.genes$Gene.name)
## All in Ara+1 or Ara-4! Does this mean that because these had mutations removed, they get bigger beneficial mutations?
## also note that Ara-4 lost a pykF allele and picked up a different pykF allele!
## lets rank rest.new.odd.mutations by G.score.
## there's a better join that this so that I don't have to set all the columns to NULL but whatever.
G.score.temp <- mutate(G.score.data,gene.annotation=Gene.name,Gene.name=NULL,Observed.nonsynonymous.mutation=NULL,
                       Expected.nonsynonymous.mutation=NULL,Synonymous.mutation=NULL, Intergenic.point.mutation=NULL,
                       Point.mutation.in.pseudogene.or.noncoding.gene=NULL,IS.insertion=NULL,Short.indel=NULL,
                       Large.deletion=NULL, Long.duplication=NULL,
                       Total.mutations=NULL, Excluding.nonsynonymous..and.synonymous=NULL)
G.score.rest.new.odd.mutations <- left_join(rest.new.odd.mutations,G.score.temp) %>% mutate(odd=NULL,Coding.length=NULL,Start.position=NULL,introgression.score=NULL,mutation=NULL,mut.type=NULL,lbl=NULL,reference=NULL,Gene.order=NULL) %>% na.omit() %>% arrange(desc(G.score))

## get average G.score of genes mutated in Ara+1, Ara-3, Ara-4.
special.G.score <- filter(G.score.rest.new.odd.mutations, lineage == 'Ara+1'|lineage=='Ara-3'|lineage=='Ara-4') %>% summarize(mean.G.score=mean(G.score))

#######################################################
## Map recombination breakpoints. Currently breakpoints occur ON
## mutations, NOT between mutations.
## Main goal of this analysis is to make Figure 4,
## a plot of chunk length distributions.

## This function maps a vector of labels to a vector of transitions between
## chunks from K-12 and chunks from LTEE recipient.
## '1-2' is a start of a K-12 chunk,
## '2-1' is the start of a LTEE chunk,
## '0' marks positions in between breakpoints.
## NOTE: in reality breakpoints lie between the i-1 and ith markers,
## whereas this code places breakpoints at the ith marker,
## so these lengths are an approximation at best.

labels.to.chunks <- function(labelz) {
    chunks <- rep('0',length(labelz))
    ## in.K12.chunk might not be FALSE at the first marker in the genome.
    ## the final if statement in the function checks this assumption.
    in.K12.chunk <- FALSE
    last.transition.index <- 1
    for (i in 1:length(labelz)) {
        # handle transitions from K-12 to LTEE and vice-versa.
        if (!in.K12.chunk & (labelz[i] == '1'|labelz[i] == '4')) {
            in.K12.chunk <- TRUE
            chunks[i] <- '1-2'
            last.transition.index <- i
        } else if (in.K12.chunk & (labelz[i] == '0'|labelz[i] == '2')) {
            in.K12.chunk <- FALSE
            chunks[i] <- '2-1'
            last.transition.index <- i
        }
    }
    ## if the last transition in chunks is '1-2', then the first index
    ## should be '0', not '1-2' (in.K12.chunk was TRUE!).
    if (chunks[1] == '1-2' & chunks[last.transition.index] == '1-2')
        chunks[1] <- '0'

    return(chunks)
}

## label each mutation as TRUE if in.K12.chunk and FALSE if not in.K12.chunk.
label.segments <- function(labelz) {
    segment.label <- rep(FALSE,length(labelz))
    ## in.K12.chunk might not be FALSE at the first marker in the genome.
    ## the final if statement in the function checks this assumption.
    in.K12.chunk <- FALSE
    last.transition.index <- 1
    for (i in 1:length(labelz)) {
        # handle transitions from K-12 to LTEE and vice-versa.
        if (!in.K12.chunk & (labelz[i] == '1'|labelz[i] == '4')) {
            in.K12.chunk <- TRUE
            last.transition.index <- i
        } else if (in.K12.chunk & (labelz[i] == '0'|labelz[i] == '2')) {
            in.K12.chunk <- FALSE
            last.transition.index <- i
        }
        segment.label[i] <- in.K12.chunk
    }
    ## if the last transition in chunks is '1-2', then the first index
    ## should be '0', not '1-2' (in.K12.chunk was TRUE!).
    if (segment.label[1] == FALSE & segment.label[last.transition.index] == TRUE)
        segment.label[1] <- TRUE

    return(segment.label)
}

## index segments in order to group sites in the same segment.
## segment.labels is the vector of TRUE or FALSE returned by label.segments.
index.segments <- function(segment.labels) {
    segment.indexes <- rep(0,length(segment.labels))
    cur.index = 1
    segment.indexes[cur.index] = 1
    for (i in 2:length(segment.labels)) {
        if (segment.labels[i] != segment.labels[i-1])
            cur.index <- cur.index + 1
        segment.indexes[i] <- cur.index
    }
    ## check if last segment is actually part of first segment.
    ## if so, replace the last segment label (cur.index) with 1.
    if (segment.labels[1] == segment.labels[length(segment.labels)])
        segment.indexes <- sapply(segment.indexes, function(x) ifelse(x==cur.index,1,x))

    return(segment.indexes)
}

## This function labels transitions, then labels and indexes segments,
## and then calculates changes to segment length due to indels.

## It is important that the input gets grouped by variable 'group', so that
## different genomes are processed separately.
calc.indel.change <- function(genomes.df) {
    df <- genomes.df %>%
        group_by(lineage) %>%
        mutate(chunk.transitions=labels.to.chunks(lbl)) %>%
        mutate(K12.chunk=label.segments(lbl)) %>%
        mutate(chunk.index=index.segments(K12.chunk)) %>%
        group_by(chunk.index)

    ## add up insertions and subtract deletions in each chunk.
    insertions <- filter(df,mut.type=='INS')
    deletions <- filter(df,mut.type=='DEL')

    del.sizes <- summarize(deletions, deleted = -sum(as.numeric(mutation)))
    ins.sizes <- summarize(insertions, inserted = sum(as.numeric(mutation)))

    indels <- full_join(del.sizes,ins.sizes) %>%
        mutate(deleted=ifelse(is.na(deleted),0,deleted)) %>%
        mutate(inserted=ifelse(is.na(inserted),0,inserted)) %>%
        mutate(delta.length=inserted-deleted)

    df <- left_join(df,indels) %>%
        mutate(deleted=ifelse(is.na(deleted),0,deleted)) %>%
        mutate(inserted=ifelse(is.na(inserted),0,inserted)) %>%
        mutate(delta.length=ifelse(is.na(delta.length),0,delta.length)) %>%
        ungroup() ## remove grouping to avoid problems downstream.

    return(df)

}

## works on both even and odd.genomes.
calc.chunks <- function(odd.genomes) {
    REL606.GENOME.LENGTH <- 4629812
    genome.chunks <- calc.indel.change(odd.genomes) %>%
        filter(chunk.transitions != '0') %>%
        ## position - lag(position) gives the distance between the previous transition
        ## marker. So, there are N-1 chunk lengths for N markers (first entry is NA).
        group_by(lineage) %>% ## make sure groups are by genome.
        mutate(chunk.length=position-lag(position)) %>%
        ## since genomes are circular, we add up the chunks at the beginning and
        ## end of each genome and assign to the first entry of chunk.length.
        ##(start of first chunk + genome.length-start of last chunk.
        mutate(chunk.length=replace(chunk.length,is.na(chunk.length),
                                    position[1]+REL606.GENOME.LENGTH-position[n()])) %>%
        ## finally, correct for indels.
        mutate(corrected.chunk.length=chunk.length+delta.length)

    return(genome.chunks)
}

## The group_by call is important, so that different genomes are processed
## separately.
odd.genome.chunks <- calc.chunks(odd.genomes)
even.genome.chunks <- calc.chunks(even.genomes)

## If a chunk.transition is '1-2', the corresponding chunk length is LTEE (1).
## If a chunk.transition is '2-1', the corresponding chunk length is K-12 (2).
## This is because of the position-lag(position) code.

odd.K12.chunks <- filter(odd.genome.chunks,chunk.transitions=='2-1') %>% mutate(segment.type='K-12')
odd.LTEE.chunks <- filter(odd.genome.chunks,chunk.transitions=='1-2') %>% mutate(segment.type='REL606')

### Do evens.
even.K12.chunks <- filter(even.genome.chunks,chunk.transitions=='2-1') %>% mutate(segment.type='K-12')
even.LTEE.chunks <- filter(even.genome.chunks,chunk.transitions=='1-2') %>% mutate(segment.type='REL606')

## join even and odd chunks.
all.K12.chunks <- full_join(even.K12.chunks,odd.K12.chunks)
K12.chunk.summary <- group_by(all.K12.chunks,genome) %>% summarize(donor.DNA=sum(chunk.length))
all.LTEE.chunks <- full_join(even.LTEE.chunks,odd.LTEE.chunks)
LTEE.chunk.summary <- group_by(all.LTEE.chunks,genome) %>% summarize(LTEE.DNA=sum(chunk.length))
all.chunk.summary <- full_join(K12.chunk.summary,LTEE.chunk.summary) %>% mutate(percent.donor=donor.DNA/(donor.DNA+LTEE.DNA))

## find the number of greens in K-12 chunks out of total number of greens and reds in each
## genome.
genome.segments <- labeled.mutations %>% group_by(genome) %>%
    mutate(chunk.transitions=labels.to.chunks(lbl))

total.LTEE.markers <- labeled.mutations %>% group_by(genome) %>% filter(lbl=='2'|lbl=='4') %>%
    filter(gene.annotation %in% top.G.score.genes$Gene.name) %>%
    filter(mut.annotation != 'dS') %>% summarize(LTEE.marker.count=n())

green.markers <- labeled.mutations %>% group_by(genome) %>% filter(lbl=='4') %>%
    filter(gene.annotation %in% top.G.score.genes$Gene.name) %>%
    filter(mut.annotation != 'dS') %>% summarize(green.LTEE.marker.count=n())

## find the number of greens in K-12 chunks out of total number of greens and reds in each
## genome.
K12.deleted.markers <- labeled.mutations %>% group_by(genome) %>%
    mutate(in.K12=label.segments(lbl)) %>%
    filter(lbl=='4' & in.K12==TRUE) %>%
    filter(gene.annotation %in% top.G.score.genes$Gene.name) %>%
    filter(mut.annotation != 'dS') %>% summarize(K12.deleted.LTEE.marker.count=n())

deleted.marker.summary <- full_join(total.LTEE.markers,K12.deleted.markers) %>% mutate(percent.green=K12.deleted.LTEE.marker.count/LTEE.marker.count)

donor.and.deleted.markers <- full_join(all.chunk.summary,deleted.marker.summary)
## write table to file for Rich to look at.
write.csv(donor.and.deleted.markers,"/Users/Rohandinho/Desktop/percent-greens-in-donor.csv")

##############################
## Make Figure 4 (chunk length distributions)

all.odd.chunks <- full_join(odd.K12.chunks,odd.LTEE.chunks) %>%
    ## change segment.type to Recipient or Donor
    mutate(segment.type = ifelse(segment.type == 'K-12','K-12 Donor','Recipient')) %>%
    ## reorder factor for plotting.
    ungroup() %>%
    mutate(lineage=factor(lineage,levels=levels(lineage)[c(7:12,1:6)]))

Fig4 <- ggplot(all.odd.chunks, aes(x=log10(chunk.length))) + geom_histogram(bins=35) +
    facet_grid(lineage ~ segment.type, scales="free_y") +
    theme_classic() +
    xlab(expression("log"[10]*"(Segment Length)")) +
    ylab("Count") +
    theme(text=element_text(family="serif")) +
    theme(strip.background=element_blank()) +
    theme(panel.grid.minor.x=element_line(color='grey90',linetype="dashed"))

ggsave("/Users/Rohandinho/Desktop/Fig4.pdf",Fig4,width=4,height=7)

## STATISTICAL TEST:
## are the distributions of K-12 chunks (or LTEE chunks)
## identical across replicate recombinant lines? Answer: similar, but not draws
## from the same distribution.

## at least one lineage has a different distribution of K-12 chunk lengths
kruskal.test(chunk.length ~ lineage, data=odd.K12.chunks)

## distribution of LTEE chunk lengths are somewhat similar across lineages.
kruskal.test(chunk.length ~ lineage, data=odd.LTEE.chunks)

## omit mutators and weird ones.
skip.me <- c("Ara+2","Ara-3","Ara-2", "Ara+3", "Ara+6")
kruskal.test(chunk.length ~ lineage, data=filter(odd.K12.chunks,!lineage %in% skip.me))
kruskal.test(chunk.length ~ lineage, data=filter(odd.LTEE.chunks,!lineage %in% skip.me))

###################################
## Table 5. New synonymous mutations in evolved clones.
## Tabulate dS and number of recombination chunks.
##calculate ratio of dS by recombination/dS by mutation
##empirically for each clone.

all.K12.chunk.count <- all.K12.chunks %>% group_by(genome) %>%
    summarize(recombination.events=n())

recomb.mut.ratio.table <- filter(labeled.mutations,mut.annotation=='dS') %>%
    group_by(genome) %>% mutate(recomb.dS=ifelse(lbl %in% c(1,4,6,7,8,9),1,0)) %>%
    mutate(new.dS=ifelse(lbl==3,1,0)) %>%
    summarise(tot.new.dS=sum(new.dS),tot.recomb.dS=sum(recomb.dS)) %>%
    left_join(all.K12.chunk.count) %>%
    mutate(dS.r.over.m.ratio=tot.recomb.dS/tot.new.dS) %>%
    mutate(event.r.over.m.ratio=recombination.events/tot.new.dS)

write.csv(recomb.mut.ratio.table,"../results/Table5_recomb_mut_ratio.csv",row.names=FALSE,quote=FALSE)

####################################################################################################
## Figure 5. Ara-3 clones have the F-plasmid.

Fig5.data <- STLE.clone.F.coverage %>% filter(Strain.Type != 'Recipient')

Fig5A.data <- Fig5.data %>% filter (Lineage == 'Ara-3')
Fig5B.data <- Fig5.data %>% filter (Lineage == 'Donor')
Fig5C.data <- Fig5.data %>% filter (Lineage != 'Ara-3' & Lineage != 'Donor')

## plot the Ara-3 clone F coverage in shades of orange.
Fig5A <- ggplot(data=Fig5A.data,aes(x=Position,y=Coverage,color=Clone,group=Clone)) +
    geom_line() +
    theme_tufte() +
    guides(color=FALSE) +
    scale_color_manual(values=c('#feedde','#fdbe85','#fd8d3c','#d94701'))

## keep default colors so that donor colors match with other figures.
Fig5B <- ggplot(data=Fig5B.data,aes(x=Position,y=Coverage,color=Clone,group=Clone)) +
    geom_line() +
    theme_tufte() +
    guides(color=FALSE)

## plot the rest of the recombinant clones in shades of grey.
Fig5C <- ggplot(data=Fig6C.data,aes(x=Position,y=Coverage,color=Lineage,group=Clone)) +
    geom_line() +
    theme_tufte() +
    guides(color=FALSE) +
    scale_colour_grey()

## arrange panes with cowplot.
Fig5 <- plot_grid(Fig5A, Fig5B, Fig5C, labels = c("A", "B", "C"), ncol = 1)
save_plot("~/Desktop/Fig5.pdf", Fig5, ncol = 1, nrow = 3, base_aspect_ratio = 3,base_height=2)

## clean up memory.
rm(Fig5.data,Fig5A.data,Fig5B.data,Fig5C.data)

############## Fig. 6: Number of donor specific mutations in each clone.
## NOTE: donor.mutations is defined around Fig. 2B code.
donor.mutations.summary <- donor.mutations %>% group_by(genome,lbl) %>% summarise(count=n())

Fig6 <- ggplot(donor.mutations, aes(x=genome,fill=lbl)) + geom_bar() + theme_tufte() + ylab("Number of donor-specific markers") + xlab("Clone") + scale_fill_discrete(name='Donor',labels=c('REL288','REL291','REL296','REL298')) + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(text=element_text(family="serif")) + guides(fill=FALSE)
ggsave("/Users/Rohandinho/Desktop/Fig6.pdf",Fig6,width=6,height=4)

###############################################################################################
## Analyze evolution experiment results.

## A reminder to what the labels are.
## 0) reference genome state.
## 1) B/K-12 marker is present that is not in reference genome (yellow)
## 2) LTEE recipient mutations (light blue)
## 3) new mutations (black)
## 4) LTEE recipient mutations that were replaced by K-12 or otherwise missing (red)
## 5) deleted markers (marker falls in a deleted region). (light purple)
## 6) REL288 specific marker
## 7) REL291 specific marker
## 8) REL296 specific marker
## 9) REL298 specific marker

evoexp.data <- evoexp.labeled.mutations %>%
    mutate(generation=ifelse(grepl('149',genome),1000,1200)) %>%
    ## if a mutation has frequency 'NA', it is either 0,4,5.
    ## therefore, it is a fixed B marker or a fixed replaced or deleted mutation.
    mutate(frequency=ifelse(is.na(frequency),1,frequency)) %>%
    arrange(lineage,position,generation) %>% group_by(lineage,position) %>%
    mutate(initial.freq = frequency) %>% mutate(final.freq = lead(frequency)) %>%
    mutate(delta.freq=lead(frequency)-frequency) %>%
    na.omit()  %>% select(-generation,-genome) %>% ungroup()

################
## Figure 7.
## Do K-12 mutation tend to decay over time? Not in this picture.
## However, selection is already acting on deleterious K-12 mutations.
## See Good and Desai (2014) for discussion of deleterious hitchhikers.
K12.evoexp.data <- filter(evoexp.data,lbl==1)

## cluster mutations based on initial frequency, final frequency, and genomic position.
make.Fig7 <- function(K12.evoexp.data) {

    Fig7.cluster.plot <- ggplot(K12.evoexp.data,
                                aes(x=initial.freq,y=final.freq,color=rotated.position)) +
        geom_point(size=0.3) +
        theme_classic() +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
        theme(text=element_text(family='serif')) +
        theme(strip.background=element_blank()) +
        geom_abline(slope=1,intercept=0,size=0.5,linetype="dashed") +
        scale_color_viridis(option="magma") +
        facet_wrap(~lineage , ncol=4) +
        xlim(0,1) +
        ylim(0,1) +
        ylab('Frequency of K-12 alleles at generation 1200') +
        xlab('Frequency of K-12 alleles at generation 1000') +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        guides(color=FALSE)

    return(Fig7.cluster.plot)
}

Fig7 <- make.Fig7(K12.evoexp.data)
ggsave("/Users/Rohandinho/Desktop/Fig7.pdf",Fig7,width=7.5,height=7.5)

## make the same figure, but for 'new' mutations (lbl==3).
new.evoexp.data <- filter(evoexp.data,lbl==3)
new.plot <- make.Fig7(new.evoexp.data)
ggsave("/Users/Rohandinho/Desktop/new.cluster.pdf",new.plot,width=7.5,height=7.5)
## This plot shows an extremely strong correlation between K-12 mutations and 'new'
## mutations. The best interpretation is that 'new' mutations are introduced by
## recombination, either by gene conversion, or by K-12 reads mapping to
## non-homologous E. coli B regions.

################

## How many replaced mutations fixed in the STLE and continuation?
replaced.evoexp.mutations <- evoexp.data %>% filter(lbl==4)
initial.erasures <- replaced.evoexp.mutations %>% filter(initial.freq == 1)
final.erasures <- replaced.evoexp.mutations %>% filter(final.freq == 1)
fixed.erasures <- replaced.evoexp.mutations %>% filter(final.freq == 1 & initial.freq == 1)
non.fixed.erasures <- replaced.evoexp.mutations %>% filter(final.freq != 1 & initial.freq == 1)

## Is there evidence of new mutations compensating for replaced beneficial LTEE mutations?
top.fixed.erasures <- filter(fixed.erasures,gene.annotation %in% genes.to.label$Gene.name)
new.replacements <- filter(evoexp.data,gene.annotation %in% top.fixed.erasures$gene.annotation)
new.replacements2 <- filter(new.replacements,final.freq==1,mut.annotation=='dN')
################################################################TODO!!!! FINISH THIS ANALYSIS!!!

## filter out mutations also in the recombinant clones.

##    Do K-12/new mutations that reached very high frequency (near fixation)
##    or very low frequency (near extinction) during the STLE continuation reject
##    the neutral expectation that P(fixation) = initial.freq? No.

neutrality.test.data <- filter(evoexp.data,final.freq==1 |final.freq==0) %>%
    filter(lbl==1 | lbl==3)

##    Compare K-12 to 'new' mutations in STLE continuation.
##    Consider a mutation fixed, if it is at frequency 1
##    at both 1000 gen and 1200 gen timepoints.
evoexp.fixations <- evoexp.data %>% filter(final.freq==1 & initial.freq==1) %>%
    filter(lbl==1 | lbl==3)

make.parallel.table <- function(grouped_df,K12=TRUE, dS=TRUE) {
    my.lbl <- ifelse(K12,1,3)
    my.mut.type <- ifelse(dS,'dS','dN')
    parallel.table <- filter(grouped_df,lbl==my.lbl) %>%
        filter(mut.annotation==my.mut.type) %>%
        summarise(parallelism=n_distinct(lineage)) %>% filter(parallelism > 1) %>% arrange(desc(parallelism))
    return(parallel.table)
}

new.parallel.evoexp.dN.fixations.by.allele <- evoexp.fixations %>%
    group_by(position,gene.annotation) %>% make.parallel.table(K12=FALSE, dS=FALSE)
new.parallel.evoexp.dS.fixations.by.allele <- evoexp.fixations %>%
    group_by(position,gene.annotation) %>% make.parallel.table(K12=FALSE, dS=TRUE)
new.parallel.evoexp.dN.fixations.by.gene <- evoexp.fixations %>%
    group_by(gene.annotation) %>% make.parallel.table(K12=FALSE, dS=FALSE)
new.parallel.evoexp.dS.fixations.by.gene <- evoexp.fixations %>%
    group_by(gene.annotation) %>% make.parallel.table(K12=FALSE, dS=TRUE)

K12.parallel.evoexp.dN.fixations.by.allele <- evoexp.fixations %>%
    group_by(position,gene.annotation) %>% make.parallel.table(K12=TRUE, dS=FALSE)
K12.parallel.evoexp.dS.fixations.by.allele <- evoexp.fixations %>%
    group_by(position,gene.annotation) %>% make.parallel.table(K12=TRUE, dS=TRUE)
K12.parallel.evoexp.dN.fixations.by.gene <- evoexp.fixations %>%
    group_by(gene.annotation) %>% make.parallel.table(K12=TRUE, dS=FALSE)
K12.parallel.evoexp.dS.fixations.by.gene <- evoexp.fixations %>%
    group_by(gene.annotation) %>% make.parallel.table(K12=TRUE, dS=TRUE)

## New mutations look like false positives due to gene conversion (or bad mapping).
## That's because dN and dS cluster on position, and dS 'new mutations' occur at the
## allele level in multiple lineages.
new.parallel.evoexp.dS.fixations.by.allele
new.parallel.evoexp.dN.fixations.by.allele
new.parallel.evoexp.dS.fixations.by.gene
new.parallel.evoexp.dN.fixations.by.gene

K12.parallel.evoexp.dS.fixations.by.allele
K12.parallel.evoexp.dN.fixations.by.allele
## 9408 K12 dS alleles to 2967 dN alleles.

K12.parallel.evoexp.dS.fixations.by.gene
K12.parallel.evoexp.dN.fixations.by.gene
## 914 genes with K12 dS alleles to 1031 genes with K12 dN alleles.

## Compare to the total number of K12 dS and dN differences from REL606.
all.K12.dS <- K12.diff.data %>% filter(mut.annotation=='dS')
all.K12.dN <- K12.diff.data %>% filter(mut.annotation=='dN')
## 23827 K-12 dS, 7969 K-12 dN.
binom.test(2967,9408,p=7969/23827)
## This test seems to indicate purifying selection against dN.
## the difference is significant, but this test is a hack. Might pattern go away if done
## more precisely?

## Compare to  number of non-fixed differences in evoexp data.
obs.evoexp.dS <- evoexp.data %>% filter(mut.annotation=='dS',lbl=='1')
obs.evoexp.dN <- evoexp.data %>% filter(mut.annotation=='dN',lbl=='1')
## 64037 dS, 19730 dN.
## what about those that didn't fix?
non.fixed.evoexp.dS <- obs.evoexp.dS %>% filter(final.freq != 1 | initial.freq != 1)
non.fixed.evoexp.dN <- obs.evoexp.dN %>% filter(final.freq != 1 | initial.freq != 1)
## 21551 dS, 6726 dN.
binom.test(2967,9408,p=19730/64037)
binom.test(2967,9408,p=6726/21551)
## no difference. Therefore, I don't conclude that there is a difference in the numbers of
## dS and dN that went to fixation compare to those observed in the STLE continuation,
## and therefore there is no evidence of purifying selection on dN. Probably a result of the
## dynamics being driven by hidden driver mutations. AND dN have already been pre-screened by
## prior selection on the K-12 background.

################
## Figures 8 and 9. Respectively these are
## versions of Figs. 1 and 3, using mutations that fixed in the STLE to estimate
## the LCA of each STLE population.

LCA.evoexp.data <- filter(evoexp.data,initial.freq==1 & final.freq == 1)

## take intersection of even and odd clones.
clone.intersection <- intersect(select(odd.genomes,-genome,-odd,-frequency),select(even.genomes,-genome,-odd,-frequency))
## to prevent different factor levels to get in the way of comparisons.
clone.intersection[] <- lapply(clone.intersection, as.character)

## now, calculate Jaccard index with LCA.evoexp.data to see how good
    ## the LCA.evoexp.data inference is.

evoexp.precomp <- select(LCA.evoexp.data,-frequency,-initial.freq,-final.freq,-delta.freq)
## to prevent different factor levels to get in the way of comparisons.
evoexp.precomp[] <- lapply(evoexp.precomp,as.character)

evoexp.clone.intersection <- intersect(clone.intersection,evoexp.precomp)
evoexp.clone.union <- union(clone.intersection,evoexp.precomp)

my.jaccard.index <- nrow(evoexp.clone.intersection)/nrow(evoexp.clone.union)
## 0.82. OK but not great.

## convert data frame to datatypes expected by Fig1.function.
inferred.LCA <- mutate(evoexp.clone.intersection,
                       lineage=factor(lineage),
                       mut.type=factor(mut.type),
                       reference=factor(reference),
                       position=as.integer(position),
                       mutation=factor(mutation),
                       mut.annotation=factor(mut.annotation),
                       gene.annotation=factor(gene.annotation),
                       lbl=factor(lbl),
                       rotated.position=as.double(rotated.position))

## For the LCA picture, only plot mutations in A & B (found in both clones and initial.freq=1 & final.freq=1)

Fig8 <- Fig1.function(inferred.LCA,genes.to.label,analysis.type='evoexpLCA')
ggsave("/Users/Rohandinho/Desktop/Fig8.pdf", Fig8,width=8,height=10)

scored.LCA <- score.introgression(inferred.LCA)
LCA.introgressed.genes <- get.introgressed.genes(scored.LCA)

Fig9 <- makeFig3(scored.LCA,auxotrophs,hfrs)
ggsave("/Users/Rohandinho/Desktop/Fig9.pdf",Fig9,height=2.5,width=6)
## NOTE: fix delta symbol in Illustrator.

########################################
## Figure S5: F-plasmid coverage in STLE continuation experiment
## at generation 1000 and 1200.

FigS5 <- ggplot(data=STLE.evoexp.F.coverage,aes(x=Position,y=log10(Coverage + 1),color=Lineage)) +
    geom_line() +
    theme_tufte() +
    ylab(expression("log"[10]*"(Coverage + 1)")) +
    scale_color_manual(values=c("red",rep('black',6),'red',rep('black',3))) +
    guides(color=FALSE) +
    facet_grid(Lineage ~ Generation,labeller=labeller(Generation = c('1000'='Generation 1000','1200'='Generation 1200')))

ggsave("/Users/Rohandinho/Desktop/FigS5.pdf",FigS5)
