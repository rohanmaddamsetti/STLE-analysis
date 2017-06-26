## dissertation-analysis.R by Rohan Maddamsetti.

## This script makes figures and does statistics for the recombinant genome analysis.

## go through these imports and figure out which ones are superfluous.
library(xml2)
library(roxygen2)
library(assertthat)
library(IRanges)
library(GenomicRanges)
library(genbankr)

library(purrr)       # consistent & safe list/vector munging
library(tidyr)       # consistent data.frame cleaning
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(DT)          # prettier data.frame output
library(data.table)  # faster fread()
library(dplyr)       # consistent data.frame operations.
library(dtplyr)      # dplyr works with data.table now.
library(ggrepel)     # plot labeled scatterplots.
library(zoo)


#' function to rotate genome coordinates, setting oriC at the center of plots.
rotate.chr <- function(my.position) {
    REL606.LENGTH <- 4629812
    ORIC <- 3886105
    midpoint <- REL606.LENGTH/2
    L <- ORIC - midpoint
    ifelse(my.position > L,my.position-ORIC,REL606.LENGTH-ORIC+my.position)
}

#############
## Import data.

## metadata for sequenced samples.
pop.clone.file <- file.path("../doc/Populations-and-Clones.csv")
pop.and.clone.metadata <- read.csv(pop.clone.file)

## For each dataset, add rotated genome coordinates so that oriC is at the center of the chromosome in plots.
hand.annotation <- tbl_df(read.csv("../results/donor_hand_annotation.csv")) %>%
    mutate(rotated.position=rotate.chr(position))

##annotate auxotrophs and Hfr oriTs on Fig2 type plots.
auxotrophs <- filter(hand.annotation,annotation!="Hfr"&annotation!="oriC")
hfrs <- filter(hand.annotation,annotation=="Hfr")

G.score.data <- tbl_df(read.csv("../results/036806-3.csv")) %>%
    mutate(rotated.Start.position=rotate.chr(Start.position))

## clone sequencing data.
labeled.mutations <- tbl_df(read.csv("../results/labeled_mutations.csv")) %>% mutate(lbl=as.factor(lbl)) %>%
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

##################################################################################
## Make Figures and Tables.

######
## Separate into odd or even clones, and
## This is used in a lot of the following code (Fig S2, Fig 2, Fig 4 or something like that.)
odd.genomes <- filter(labeled.mutations, odd==TRUE)
even.genomes <- filter(labeled.mutations, odd==FALSE)

#initialize introgression.score column.
odd.genomes <- mutate(odd.genomes, introgression.score = 0)

##This is a little helper to sum up labels (which is a factor)
label.levels <- levels(labeled.mutations$lbl)
labels.to.int <- sapply(label.levels, function(x) strtoi(x))

##take all positions labeled with 1, and sum them up to get the introgression score.
## omit Ara-3 (all K-12).
scored.odd.genomes <- filter(odd.genomes,lineage != "Ara-3") %>%
    filter(lbl=='1' | lbl=='0') %>%
    group_by(position) %>%
    summarize(introgression.score = sum(labels.to.int[lbl])) %>%
    mutate(rotated.position=rotate.chr(position))


############
## Figure S1. Density of differences between K-12 and REL606.

oriC.xlab <- expression(paste("Distance from ",italic("oriC")))
FigS1 <- ggplot(K12.diff.data, aes(x=rotated.position)) + geom_histogram(bins=400) + theme_classic() + xlab(oriC.xlab) + ylab("Differences between K-12 and REL606")
ggsave("/Users/Rohandinho/Desktop/FigS1.pdf",FigS1)

############################################################################################
## make Figure 1.
## break Figure 1 into 4 panels, three lineages per page (both clones on a page).

##IDEA: add a border to points where a new mutation occurs in a gene where an old mutation was deleted.

makeFig1panel <- function(labeled.mutations,G.score.data,l1,l2,l3) {
    ## REL606 markers are negative space on the plot.
    ## also, reverse levels to get Ara+1 genomes on top of plot.
    ## ylim from 0 to 70. plot 6 genomes on 10,20,30,40,50. so map genome to a center position.

    no.B.genomes <- filter(labeled.mutations,lbl!='0') %>%
        mutate(genome=factor(paste(lineage,genome,sep=': '))) %>%
        mutate(genome=factor(genome,levels=rev(levels(genome)))) %>%
        filter(lineage==l1 | lineage==l2 | lineage==l3) %>%
        mutate(y.center=70-match(genome,unique(genome))*10,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3) %>%
        mutate(lbl=as.factor(lbl))

 ## label the top 54 G-scoring genes with symbols.
    ## below this threshold, very short genes with 1 dN have a high G-score.
    top.G.score.genes <- filter(G.score.data,G.score>8.9)

    top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name)

    ## label these mutations as points in the plot.
    LTEE.top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4' | lbl=='3') %>%
        mutate(y.center=70-match(genome,unique(genome))*10,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3)

    ## draw labels on the plot for high G-score mutations.
    genes.to.label <- filter(top.G.score.genes,Gene.name %in% top.hits$gene.annotation)

    panel.labels <- select(LTEE.top.hits,lineage,genome,rotated.position,gene.annotation,y.center) %>%
        group_by(lineage) %>% mutate(ypos=mean(y.center))

    Fig1.xlab <- expression(paste("Distance from ",italic("oriC")))
    ## all donor mutations should be labeled as yellow.
    panel <- ggplot(no.B.genomes, aes(x=x1,xend=x2,y=y1,yend=y2,colour=lbl)) +
        geom_segment(size=0.05) +
        scale_colour_manual(values=c('yellow', 'red', 'black','green','blue','tan1','tan2','tan3','tan4')) +
        theme_classic() +
        xlab(Fig1.xlab) +
        ylab("Clone") + ylim(5,70) + guides(colour=FALSE) +
        scale_y_continuous(breaks=c(10,20,30,40,50,60),labels=rev(unique(no.B.genomes$genome)))+
        geom_point(data=LTEE.top.hits,aes(x=rotated.position,y=y.center,color=lbl,shape=mut.annotation)) + guides(shape=FALSE) +
        #geom_segment(data=genes.to.label, aes(x=Start.position,xend=Start.position,y=1,yend=65),size=0.1,inherit.aes=FALSE) +
        geom_segment(data=filter(genes.to.label, Gene.name %in% c('ybaL', 'fabF', 'topA', 'pykF', 'mreC/B', 'malT', 'spoT', 'hslU', 'iclR', 'nadR')), aes(x=rotated.Start.position,xend=rotated.Start.position,y=1,yend=65),size=0.1,inherit.aes=FALSE) +
        ## label genes that are in the text at the top.
        geom_text(data=filter(genes.to.label, Gene.name %in% c('ybaL', 'fabF', 'topA', 'pykF', 'mreC/B', 'malT', 'spoT', 'hslU', 'iclR', 'nadR')),aes(label=Gene.name,x=rotated.Start.position, y=67),size=3,angle=45,inherit.aes=FALSE)

    return(panel)
}

panel1 <- makeFig1panel(labeled.mutations,G.score.data,"Ara+1","Ara+2","Ara+3")
ggsave("/Users/Rohandinho/Desktop/Fig1A.pdf", panel1,width=9.5,height=8)

panel2 <- makeFig1panel(labeled.mutations,G.score.data,"Ara+4","Ara+5","Ara+6")
ggsave("/Users/Rohandinho/Desktop/Fig1B.pdf", panel2,width=9.5,height=8)

panel3 <- makeFig1panel(labeled.mutations,G.score.data,"Ara-1","Ara-2","Ara-3")
ggsave("/Users/Rohandinho/Desktop/Fig1C.pdf", panel3,width=9.5,height=8)

panel4 <- makeFig1panel(labeled.mutations,G.score.data,"Ara-4","Ara-5","Ara-6")
ggsave("/Users/Rohandinho/Desktop/Fig1D.pdf", panel4,width=9.5,height=8)

##############################################################################
############## Fig. 2: Number of donor specific mutations in each clone.
donor.mutations <- filter(labeled.mutations,lbl == 6 | lbl==7 | lbl == 8 | lbl == 9) %>%
    mutate(genome=factor(paste(lineage,genome,sep=': '))) %>%
    ## reorder factor for plotting.
    mutate(genome=factor(genome,levels=levels(genome)[c(12:21,1:11)])) %>%
    select(-odd,-frequency) %>%
    arrange(genome,position) %>%
    distinct(genome, position, .keep_all = TRUE)


donor.mutations.summary <- donor.mutations %>% group_by(genome,lbl) %>% summarise(count=n())

Fig2 <- ggplot(donor.mutations, aes(x=genome,fill=lbl)) + geom_bar() + theme_classic() + ylab("Number of donor-specific markers") + xlab("Clone") + scale_fill_discrete(name='Donor',labels=c('REL288','REL291','REL296','REL298')) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/Users/Rohandinho/Desktop/Fig2.pdf",Fig2,width=11,height=8)

############# Fig. S2: location of donor specific mutations on the chromosome.
## TODO: Make FigS2 look much better.
only.donor.mutations <- donor.mutations %>%
    mutate(y.center=35-match(lbl,sort(unique(lbl)))*10,y1=y.center-3,y2=y.center+3)


FigS2 <- ggplot(only.donor.mutations,aes(x=rotated.position,xend=rotated.position,y=y1,yend=y2,colour=lbl)) +
    geom_segment(size=0.05) + theme_classic() +
    scale_colour_discrete(name='Donor',labels=c('REL288','REL291','REL296','REL298')) +
    xlab(oriC.xlab) +
    ## add auxotroph lines
    geom_vline(data=auxotrophs,aes(xintercept=rotated.position),size=0.05,linetype="dashed") +
    geom_text_repel(data=auxotrophs,aes(x=rotated.position,y=8.2,label=annotation),inherit.aes=FALSE, size=3) +
    ## add annotation of Hfr and oriC on chromosome.
   geom_text_repel(data=hfrs,aes(x=rotated.position,y=-0.3,label=Hfr.orientation),inherit.aes=FALSE,size=3,nudge_y=-0.3)

ggsave("/Users/Rohandinho/Desktop/FigS2.pdf",FigS2,width=11,height=8)

#############################################################
## Figure 3. Visual comparisons of parallel recombination across lineages
## with the G-score of mutations over the genome, and with the occurrence of
## new mutations over the genome.
## This is only made with odd genomes at the moment.

## which genes get introgressed the most? is this a sign of positive selection?
gene.test <- filter(odd.genomes,lineage != "Ara-3") %>%
    filter(lbl=='1' | lbl=='0') %>%
    group_by(gene.annotation) %>%
    summarize(introgression.score = sum(labels.to.int[lbl])) %>% arrange(desc(introgression.score))

## omit mutator lineages (Ara+6,Ara+3,Ara-2) +6 is mutT, +3 is mutS, -2 is mutL mutator.
new.odd.mutations <- filter(odd.genomes,lineage != "Ara+6" & lineage != "Ara+3" & lineage != "Ara-2") %>% filter(lbl=='3')

Fig3 <- ggplot(scored.odd.genomes, aes(x=rotated.position,y=introgression.score)) + geom_line(size=0.05) + theme_classic() +
    xlab(oriC.xlab) + ylab("Introgression Score") +
    ## add auxotroph lines
    geom_vline(data=auxotrophs,aes(xintercept=rotated.position),size=0.05,linetype="dashed") +
    geom_text_repel(data=auxotrophs,aes(x=rotated.position,y=8.2,label=annotation),inherit.aes=FALSE, size=3) +
    ## add annotation of Hfr and oriC on chromosome.
   geom_text_repel(data=hfrs,aes(x=rotated.position,y=-0.3,label=Hfr.orientation),inherit.aes=FALSE,size=3,nudge_y=-0.3)

ggsave("/Users/Rohandinho/Desktop/Fig3.pdf",Fig3)


#############################
## Probable beneficial mutations.

## Calculation that 31 genes that are in Figure 1 and in Table 2.
## the top 54 G-scoring genes: below this threshold, very short genes with 1 dN have a high G-score.
top.G.score.genes <- filter(G.score.data,G.score>8.9)
top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name)
LTEE.top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4')
genes.to.label <- filter(top.G.score.genes,Gene.name %in% LTEE.top.hits$gene.annotation)
## gene.to.label contains 31 genes.


###########
## Fig. S3: Does sequence divergence/conservation predict location of markers? NO.
## Perhaps do this differently after showing to Rich.
## make a scatterplot of Fig. S1 and Fig. 3.

## REMEMBER: divergence and introgression is correlated-- since when there's no divergence,
## all introgression events are undetectable! So, more divergence = better resolution of introgression.
## this is extra clear when regressing sum(introgression) against divergence in the windows.

REL606.LENGTH <- 4629812
## use 556 bins. Each is 8327 bp long (integer divisors)
bin.length <- REL606.LENGTH/556

introgression.vs.divergence.data <- right_join(K12.diff.data,scored.odd.genomes) %>%
    ## bin mutations across the genome.
    ## NOTE: equal size bins, NOT equal numbers of mutations!
    mutate(my_bin=ceiling(position/bin.length)) %>% group_by(my_bin)

two.classed.introgression.data <- introgression.vs.divergence.data %>%
    summarize(introgression=ifelse(mean(introgression.score)>0,1,0),divergence=n())

## no difference in divergence between regions with and without introgression. p = 0.8691.
kruskal.test(divergence ~ introgression,data=two.classed.introgression.data)

median.introgression.data <- introgression.vs.divergence.data %>%
    summarize(introgression=median(introgression.score),divergence=n()) %>% filter(introgression > 0)

introgression.divergence.model <- lm(introgression~divergence,data=median.introgression.data)
confint(introgression.divergence.model)

##ggsave("/Users/Rohandinho/Desktop/FigS3.pdf",FigS3,width=11,height=8)


####################################################################################################
## Erased mutations. Table 3 and Figure 3 (Fig. 3 is made separately).

#### Make in-depth alignments of all erased mutations, with special attention to Ara+1 and Ara-4.
#### Print out a csv of genes for align_erased_mutations.py to align,
#### and annotate as 0) reversion to pre-LTEE state, 1) K-12 state, 3) new allele.

#### Only look at dN mutations in odd REL clones.
erased.gene.list <- filter(labeled.mutations,lbl==4,mut.annotation=='dN',odd==TRUE) %>%
    #select(gene.annotation,lineage,genome) %>%
    group_by(lineage,genome) %>% distinct(gene.annotation)
write.csv(erased.gene.list,"../results/align_these.csv",row.names=FALSE,quote=FALSE)

## Run pythonw align_erased_mutations.py to get results for Table 3 in the paper,
## as well as the numbers in this section of the manuscript.
## 31 of the 61 erased dN are in non-mutators, and these are shown in Table 3.

################################
## Table 4: Putative gene conversion events.
## Started by looking for parallelism in new mutations: super strong parallelism (multiple new mutations in the same
## gene is probably gene conversion or something.

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

#######################################################################################################################
## K-12 specific genes analysis.

## positions have already been filtered by a coverage threshold (past 5% probability
## under REL606 1X coverage distribution).

#' Find K-12 specific segments. return boundaries for each segment.
#' @export
find.K12.segments <- function(filtered.coverage) {

    ## calculate intervals.
    boundaries <- filtered.coverage %>% mutate(left.diff=K12_position - lag(K12_position)) %>%
        mutate(right.diff=lead(K12_position) - K12_position) %>%
        ## corner case: check for the NA values at the endpoints and set them as boundaries.
        mutate(is.right.boundary=is.na(right.diff)|ifelse(right.diff>1,TRUE,FALSE)) %>%
        mutate(is.left.boundary=is.na(left.diff)|ifelse(left.diff>1,TRUE,FALSE)) %>%
        filter(is.left.boundary==TRUE | is.right.boundary==TRUE)

    left.boundaries <- filter(boundaries,is.left.boundary==TRUE) %>% arrange(K12_position)
    right.boundaries <- filter(boundaries,is.right.boundary==TRUE) %>% arrange(K12_position)
    assert_that(nrow(left.boundaries) == nrow(right.boundaries))

    K12.segments <- data.frame(left.boundary=left.boundaries$K12_position,
                               right.boundary=right.boundaries$K12_position) %>%
        ## filter out intervals less than the length of a read (250 bp for these data).
        mutate(len=right.boundary-left.boundary) %>% filter(len>250) %>%
        mutate(segment.index=row_number())
    return(K12.segments)
}

## this helper function is a wrapper around functions from genbankr package.
parse.ref.genes <- function(ref.gbk) {
    ref.genes <- ref.gbk %>% parseGenBank %>% make_gbrecord %>% genes
    return(ref.genes)
}

## input: ref.gbk: file.path of the reference genome,
##    : data.frame returned by find.K12.segments.
annotate.segments <- function(segments,ref.genes) {

    ## create the IRanges object.
    amp.ranges <- IRanges(segments$left.boundary,segments$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with ref genes.
    g.amp.ranges <- GRanges("K-12",ranges=amp.ranges)
    ## and add the data.frame of segments as metadata.
    mcols(g.amp.ranges) <- segments

    ## find overlaps between ref genes and segments.
    hits <- findOverlaps(ref.genes,g.amp.ranges,ignore.strand=FALSE)
    ## take the hits, the ref annotation, and the segments,
    ## and produce a table of genes found in each segment.

    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))

    query.df <- data.frame(query.index=1:length(ref.genes),
                           gene=ref.genes$gene,locus_tag=ref.genes$locus_tag,
                           start=start(ranges(ref.genes)),end=end(ranges(ref.genes)))

    subject.df <- bind_cols(data.frame(subject.index=1:length(g.amp.ranges)),data.frame(mcols(g.amp.ranges)))

    segment.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))

    return(segment.genes.df)
}

projdir <- "/Users/Rohandinho/Desktop/Projects/STLE-analysis"
K12.specific.dir <- file.path(projdir,"results/K12-specific-genes")
breseq.out.dir <- file.path(projdir,"breseq-assemblies/REL606-polymorphism")
REL606.gb <- file.path(projdir,"references/REL606.7.gbk")
K12.gb <- file.path(projdir,"references/K-12.1.gbk")

## parse K12.genes (this is sloooow.)
K12.genes <- parse.ref.genes(K12.gb)

## pop.and.clone.metadata contains info on which samples are mixed pops. etc.
K12.LENGTH <- 4641652

K12.cov.file <- file.path(K12.specific.dir,"K12-coverage.csv")

K12.specific.coverage <- fread(K12.cov.file)
pop.and.clone.metadata <- mutate(pop.and.clone.metadata,genome=Name)
K12.specific.data <- tbl_df(left_join(data.frame(K12.specific.coverage),pop.and.clone.metadata))

K12.clone.data <- filter(K12.specific.data,is.Clone==1,Generation==1000) %>%
    select(-Name) %>% mutate(Name=factor(paste(Population,REL.Name,sep=': '))) %>%
    group_by(Name)

## label clones as even or odd based on REL numbering for parallelism analysis.

clone.genome.names <- levels(K12.clone.data$REL.Name)
is.odd <- sapply(clone.genome.names, function(x) ifelse(strtoi(substr(x,nchar(x),nchar(x)))%%2, TRUE,FALSE))

odd.K12.clone.data <- K12.clone.data %>% mutate(odd=is.odd[REL.Name]) %>%
    filter(odd==TRUE)

## group_by mixed.pop.data and then filter for positions
## that passed the threshold in both T0 and T1 samples.
K12.mixed.pop.data <- filter(K12.specific.data,is.Clone==0,Generation %in% c(1000,1200)) %>%
    group_by(Population,K12_position) %>% summarize(row.num=n()) %>% filter(row.num==2)

## Split-Apply-Combine to get K-12 specific segments.

fixed.K12.segments <- K12.mixed.pop.data %>% do(find.K12.segments(.))
clone.K12.segments <- K12.clone.data %>% do(find.K12.segments(.))

odd.clone.K12.segments <- odd.K12.clone.data %>% do(find.K12.segments(.))

## Now, annotate those segments.
annotated.fixed.K12.segments <- annotate.segments(fixed.K12.segments,K12.genes)
annotated.clone.K12.segments <- annotate.segments(clone.K12.segments,K12.genes)
odd.annotated.clone.K12.segments <- annotate.segments(odd.clone.K12.segments,K12.genes)

write.csv(x=annotated.clone.K12.segments,file=file.path(K12.specific.dir,"clone.K12.segs.csv"))
write.csv(x=annotated.fixed.K12.segments,file=file.path(K12.specific.dir,"fixed.K12.segs.csv"))

## Make tables and figures.

## This shows that IS5 transposons found in K-12 but not in REL606,
## as shown by a BLAST search, have introgressed into evolved genomes
## (can't tell if fixed since occur at multiple copies)
fixed.K12.seg.parallelism <- annotated.fixed.K12.segments %>% group_by(gene,locus_tag) %>% summarise(count=n()) %>% arrange(desc(count),locus_tag)


## TODO: DOUBLE CHECK THESE NUMBERS BY HAND!
clone.K12.seg.parallelism <- odd.annotated.clone.K12.segments %>% group_by(gene,locus_tag) %>% summarise(count=n()) %>% arrange(desc(count),locus_tag)

## TODO: generalize the code for making Fig. 1 and Fig. 2 to avoid code duplication,
## and then use here for the K-12 specific introgression, fixed mutations
## in the STLE continuation experiment, and for the clones.


####TODO:  ####          ######                     ##########      ##  ###  #      #######        ##########       ######    ####    ####
##              MAKE             K-12 SPECIFIC GENES           FIGURES          AND         TABLES!!!!
## Figure S3 A,B,C,D: K-12 specific genes based on K12-reference numbering.
## Figure S4: Parallel introgression of K-12 specific genes.


#######################################################
## Map recombination breakpoints. Currently breakpoints occur ON
## mutations, NOT between mutations.
## Make Figure 5 (chunk length distributions) and do two statistical tests:
## Distance between new mutations and breakpoints (are breakpoints mutagenic?)
## recombination hot and cold spots

##global variable used in Fig. 5 code.
REL606.GENOME.LENGTH <- 4629812

## This function maps a vector of labels to a vector of transitions between
## chunks from K-12 and chunks from LTEE recipient.
## '1-2' is a start of a K-12 chunk,
## '2-1' is the start of a LTEE chunk,
## '0' marks positions in between breakpoints.
## NOTE: in reality breakpoints lie between the i-1 and ith markers,
## whereas this code places breakpoints at the ith marker,
## so these lengths are an approximation at best.

## It would be GREAT if I could write this without a for loop,
## this is a 15 second speed bottleneck.


############################  IMPORTANT BUG TO FIX!!!!
############################  THIS CODE DOES NOT HANDLE DELETIONS OR INSERTIONS.


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

#########################################

## The group_by call is important, so that different genomes are processed
## separately.
odd.genome.chunks <- group_by(odd.genomes, lineage) %>%
    mutate(chunk.transitions=labels.to.chunks(lbl)) %>%
    filter(chunk.transitions =='1-2'|chunk.transitions=='2-1') %>%
## position - lag(position) gives the distance between the previous transition
## marker. So, there are N-1 chunk lengths for N markers (first entry is NA).
    mutate(chunk.length=position-lag(position)) %>%
## since genomes are circular, we add up the chunks at the beginning and
## end of each genome and assign to the first entry of chunk.length.
##(start of first chunk + genome.length-start of last chunk.
    mutate(chunk.length=replace(chunk.length,is.na(chunk.length),
                                position[1]+REL606.GENOME.LENGTH-position[n()]))

## If a chunk.transition is '1-2', the corresponding chunk length is LTEE (1).
## If a chunk.transition is '2-1', the corresponding chunk length is K-12 (2).
## This is because of the position-lag(position) code.

odd.K12.chunks <- filter(odd.genome.chunks,chunk.transitions=='2-1')
odd.LTEE.chunks <- filter(odd.genome.chunks,chunk.transitions=='1-2')

### Do evens.

even.genome.chunks <- group_by(even.genomes, lineage) %>%
    mutate(chunk.transitions=labels.to.chunks(lbl)) %>%
    filter(chunk.transitions =='1-2'|chunk.transitions=='2-1') %>%
## position - lag(position) gives the distance between the previous transition
## marker. So, there are N-1 chunk lengths for N markers (first entry is NA).
    mutate(chunk.length=position-lag(position)) %>%
## since genomes are circular, we add up the chunks at the beginning and
## end of each genome and assign to the first entry of chunk.length.
##(start of first chunk + genome.length-start of last chunk.
    mutate(chunk.length=replace(chunk.length,is.na(chunk.length),
                                position[1]+REL606.GENOME.LENGTH-position[n()]))

## If a chunk.transition is '1-2', the corresponding chunk length is LTEE (1).
## If a chunk.transition is '2-1', the corresponding chunk length is K-12 (2).
## This is because of the position-lag(position) code.

even.K12.chunks <- filter(even.genome.chunks,chunk.transitions=='2-1')
even.LTEE.chunks <- filter(even.genome.chunks,chunk.transitions=='1-2')

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

#########################################################################################################
## I put Ara-3 'recipient' chunks into donor chunk plot to make figure 5.
## This is what plot.K12.chunks is.

old.Fig4 <- ggplot(odd.K12.chunks, aes(x=log10(chunk.length))) + geom_histogram() + facet_grid(lineage ~ .) + theme_classic() + xlab("log10(donor segment length)") + ylab("Count")

plot.K12.chunks <- select(odd.K12.chunks,-chunk.transitions) %>% filter(lineage!='Ara-3') %>% full_join(filter(select(odd.LTEE.chunks,-chunk.transitions),lineage=='Ara-3'))

Fig5 <- ggplot(plot.K12.chunks, aes(x=log10(chunk.length))) + geom_histogram() + facet_grid(lineage ~ .) + theme_classic() + xlab("log10(donor segment length)") + ylab("Count")

ggsave("/Users/Rohandinho/Desktop/Fig5.pdf",Fig5)

#Fig4B <- ggplot(LTEE.chunks,aes(x=log10(chunk.length),fill=lineage)) + geom_histogram() + facet_grid(lineage ~ .) + theme_classic() + xlab("log10(recipient chunk length)") + ylab("Count")

## NOTE! This test was done by inputting odd.genomes2 in the group_by call
## for odd.genome.chunks to remove mutators and obviously different genomes.
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

## STATISTICAL TEST:
## Is there evidence that recombination breakpoints are mutagenic?
## TODO: find out what the average minimum distance is between new mutations
## labeled '3' and breakpoints.
## Then, drop mutations at random, and calculate average minimum distance to breakpoints.
## repeat 10000 times to calculate a p-value.
odd.breaks.and.muts <- filter(odd.genome.chunks, chunk.transitions=='2-1' | chunk.transitions=='1-2' | lbl == '3')
odd.breaks.and.muts <- select(odd.genome.chunks,lineage,genome,mut.type,position,mut.annotation,chunk.transitions,lbl)


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


write.csv(recomb.mut.ratio.table,"../results/recomb_mut_ratio.csv",row.names=FALSE,quote=FALSE)

####################################################################################################
## Analyze evolution experiment results.

evoexp.data <- evoexp.labeled.mutations %>%
    mutate(generation=ifelse(grepl('149',genome),1000,1200)) %>%
    mutate(frequency=ifelse(is.na(frequency),0,frequency))

## cluster mutations based on initial frequency, final frequency, and genomic position.
evoexp.data2 <- evoexp.data %>%
    arrange(lineage,position,generation) %>% group_by(lineage,position) %>%
    mutate(initial.freq = frequency) %>% mutate(final.freq = lead(frequency)) %>%
    mutate(delta.freq=lead(frequency)-frequency) %>%
    na.omit()  %>% select(-generation) %>% ungroup()

## look at erased mutations (lbl == 4).
evoexp.erased.mutations <- evoexp.data2 %>% filter(lbl==4,final.freq>0)
## Almost all of these mutations are in the hypermutators Ara+6,Ara+3.
## Interestingly enough, the initial.freq for all these erased mutations is 0.

################
## Figure 7.
## Do K-12 mutation tend to decay over time? Not in this picture.
## However, selection is already acting on deleterious K-12 mutations.
## See Good and Desai (2014) for discussion of deleterious hitchhikers.
K12.evoexp.data2 <- filter(evoexp.data2,lbl==1)
Fig6.cluster.plot <- ggplot(K12.evoexp.data2, aes(x=initial.freq,y=final.freq,color=position)) +
    geom_point(size=0.3) + theme_classic() + geom_abline(slope=1,intercept=0) +
    scale_color_viridis(option="plasma") +
    facet_wrap( ~ lineage , ncol=4)
ggsave("/Users/Rohandinho/Desktop/Fig6.pdf",Fig6.cluster.plot)

dN.K12.evoexp.data2 <- K12.evoexp.data2 %>% filter(mut.annotation=='dN')
dN.cluster.plot <- ggplot(dN.K12.evoexp.data2, aes(x=initial.freq,y=final.freq,color=position)) +
    geom_point(size=0.3) + theme_classic() + geom_abline(slope=1,intercept=0) +
    scale_color_viridis(option="plasma") +
    facet_wrap( ~ lineage , ncol=4)
ggsave("/Users/Rohandinho/Desktop/dNcluster.pdf",dN.cluster.plot)

dS.K12.evoexp.data2 <- K12.evoexp.data2 %>% filter(mut.annotation=='dS')
dS.cluster.plot <- ggplot(dS.K12.evoexp.data2, aes(x=initial.freq,y=final.freq,color=position)) +
    geom_point(size=0.3) + theme_classic() + geom_abline(slope=1,intercept=0) +
    scale_color_viridis(option="plasma") +
    facet_wrap( ~ lineage , ncol=4)
ggsave("/Users/Rohandinho/Desktop/dScluster.pdf",dS.cluster.plot)

## How many erased mutations fixed in the STLE and continuation?
erased.evoexp.mutations <- evoexp.data2 %>% filter(lbl==4)
initial.erasures <- erased.evoexp.mutations %>% filter(initial.freq == 1)
final.erasures <- erased.evoexp.mutations %>% filter(final.freq == 1)
fixed.erasures <- erased.evoexp.mutations %>% filter(final.freq == 1 & initial.freq == 1)

##    Do  K-12/new mutations that reached very high frequency (near fixation)
##    or very low frequency (near extinction) during the STLE continuation reject
##    the neutral expectation that P(fixation) = initial.freq? No.

neutrality.test.data <- filter(evoexp.data2,final.freq==1 |final.freq==0) %>%
    filter(lbl==1 | lbl==3)

neutrality.plot <- ggplot(neutrality.test.data,aes(x=initial.freq)) + geom_density() +
    theme_classic() + facet_wrap(final.freq ~ lbl)

## Compare K-12 to new mutations in STLE continuation.
##    Consider a mutation fixed, if it is at frequency 1
##    at both 1000 gen and 1200 gen timepoints.
evoexp.fixations <- evoexp.data2 %>% filter(final.freq==1 & initial.freq==1) %>%
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
obs.evoexp.dS <- evoexp.data2 %>% filter(mut.annotation=='dS',lbl=='1')
obs.evoexp.dN <- evoexp.data2 %>% filter(mut.annotation=='dN',lbl=='1')
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

########
## Make Figures 7 and 8. Respectively these are
## versions of Figs. 1 and 2, using mutations that fixed in the STLE to estimate
## the LCA of each STLE population.

make_LCA_Fig1panel <- function(labeled.mutations,G.score.data,l1,l2,l3,l4,l5,l6) {
    ## REL606 markers are negative space on the plot.
    ## also, reverse levels to get Ara+1 genomes on top of plot.
    ## ylim from 0 to 70. plot 6 genomes on 10,20,30,40,50. so map genome to a center position.

    no.B.genomes <- filter(labeled.mutations,lbl!='0') %>%
        mutate(lineage=factor(lineage,levels=rev(levels(lineage)))) %>%
        filter(lineage==l1 | lineage==l2 | lineage==l3 | lineage == l4 | lineage == l5 | lineage == l6) %>%
        mutate(y.center=70-match(lineage,unique(lineage))*10,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3) %>%
        mutate(lbl=as.factor(lbl))

 ## label the top 54 G-scoring genes with symbols.
    ## below this threshold, very short genes with 1 dN have a high G-score.
    top.G.score.genes <- filter(G.score.data,G.score>8.9)

    top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name)

    ## label these mutations as points in the plot.
    LTEE.top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4' | lbl=='3') %>%
        mutate(y.center=70-match(genome,unique(genome))*10,
               x1=rotated.position,x2=rotated.position,y1=y.center-3,y2=y.center+3)

    ## draw labels on the plot for high G-score mutations.
    genes.to.label <- filter(top.G.score.genes,Gene.name %in% top.hits$gene.annotation)

    panel.labels <- select(LTEE.top.hits,lineage,genome,rotated.position,gene.annotation,y.center) %>%
        group_by(lineage) %>% mutate(ypos=mean(y.center))

    Fig1.xlab <- expression(paste("Distance from ",italic("oriC")))
    ## all donor mutations should be labeled as yellow.
    panel <- ggplot(no.B.genomes, aes(x=x1,xend=x2,y=y1,yend=y2,colour=lbl)) +
        geom_segment(size=0.05) +
        ## NOTE: NO GREEN MUTATIONS BECAUSE NO ERASED MUTATIONS. SHOULD MAP LBL TO COLOR IN A BETTER WAY.
        scale_colour_manual(values=c('yellow', 'red', 'black','blue','tan1','tan2','tan3','tan4')) +
        theme_classic() +
        xlab(Fig1.xlab) +
        ylab("Inferred LCA of Population") + ylim(5,70) + guides(colour=FALSE) +
        scale_y_continuous(breaks=c(10,20,30,40,50,60),labels=rev(unique(no.B.genomes$lineage)))+
        geom_point(data=LTEE.top.hits,aes(x=rotated.position,y=y.center,color=lbl,shape=mut.annotation)) + guides(shape=FALSE) +
        #geom_segment(data=genes.to.label, aes(x=Start.position,xend=Start.position,y=1,yend=65),size=0.1,inherit.aes=FALSE) +
        geom_segment(data=filter(genes.to.label, Gene.name %in% c('ybaL', 'fabF', 'topA', 'pykF', 'mreC/B', 'malT', 'spoT', 'hslU', 'iclR', 'nadR')), aes(x=rotated.Start.position,xend=rotated.Start.position,y=1,yend=65),size=0.1,inherit.aes=FALSE) +
        ## label genes that are in the text at the top.
        geom_text(data=filter(genes.to.label, Gene.name %in% c('ybaL', 'fabF', 'topA', 'pykF', 'mreC/B', 'malT', 'spoT', 'hslU', 'iclR', 'nadR')),aes(label=Gene.name,x=rotated.Start.position, y=67),size=3,angle=45,inherit.aes=FALSE)

    return(panel)
}

LCA.evoexp.data <- filter(evoexp.data2,initial.freq==1 & final.freq == 1)

LCApanel1 <- make_LCA_Fig1panel(LCA.evoexp.data,G.score.data,"Ara+1","Ara+2","Ara+3", "Ara+4", "Ara+5", "Ara+6")
ggsave("/Users/Rohandinho/Desktop/LCA_Fig1A.pdf", LCApanel1,width=9.5,height=8)

LCApanel2 <- make_LCA_Fig1panel(LCA.evoexp.data,G.score.data,"Ara-1","Ara-2","Ara-3", "Ara-4", "Ara-5", "Ara-6")
ggsave("/Users/Rohandinho/Desktop/LCA_Fig1B.pdf", LCApanel2,width=9.5,height=8)


makeLCA_Fig2 <- function(LCA.evoexp.data,hand.annotation) {
    ## no need to worry about even or odd clones.

    ##annotate auxotrophs and Hfr oriTs on Fig2 type plots.
    auxotrophs <- filter(hand.annotation,annotation!="Hfr"&annotation!="oriC")
    hfrs <- filter(hand.annotation,annotation=="Hfr")

    ##initialize introgression.score column.
    LCA.evoexp.data <- mutate(LCA.evoexp.data, introgression.score = 0)

    ##This is a little helper to sum up labels (which is a factor)
    label.levels <- levels(LCA.evoexp.data$lbl)
    labels.to.int <- sapply(label.levels, function(x) strtoi(x))

    ##take all positions labeled with 1, and sum them up to get the introgression score.
    ## omit Ara-3 (all K-12).
    scored.LCA.evoexp.data <- filter(LCA.evoexp.data,lineage != "Ara-3") %>%
        filter(lbl=='1' | lbl=='0') %>%
        group_by(position) %>%
        summarize(introgression.score = sum(labels.to.int[lbl])) %>%
        mutate(rotated.position=rotate.chr(position))

    ## which genes get introgressed the most? is this a sign of positive selection?
    gene.test <- filter(LCA.evoexp.data,lineage != "Ara-3") %>%
        filter(lbl=='1' | lbl=='0') %>%
        group_by(gene.annotation) %>%
        summarize(introgression.score = sum(labels.to.int[lbl])) %>% arrange(desc(introgression.score))

    Fig2.xlab <- expression(paste("Distance from ",italic("oriC")))
    LCA_Fig2 <- ggplot(scored.LCA.evoexp.data, aes(x=rotated.position,y=introgression.score)) + geom_line(size=0.05) + theme_classic() +
        xlab(Fig2.xlab) + ylab("Introgression Score") +
        ## add auxotroph lines
        geom_vline(data=auxotrophs,aes(xintercept=rotated.position),size=0.05,linetype="dashed") +
        geom_text_repel(data=auxotrophs,aes(x=rotated.position,y=8.2,label=annotation),inherit.aes=FALSE, size=3) +
        ## add annotation of Hfr and oriC on chromosome.
        geom_text_repel(data=hfrs,aes(x=rotated.position,y=-0.3,label=Hfr.orientation),inherit.aes=FALSE,size=3,nudge_y=-0.3)

    return(LCA_Fig2)

}


LCA_Fig2 <- makeLCA_Fig2(LCA.evoexp.data,hand.annotation)
ggsave("/Users/Rohandinho/Desktop/LCA_Fig2.pdf",LCA_Fig2)

## A reminder to what the labels are.
## 0) reference genome state.
## 1) B/K-12 marker is present that is not in reference genome (yellow)
## 2) LTEE recipient mutations (red)
## 3) new mutations (black)
## 4) LTEE recipient mutations that were erased by K-12 or otherwise missing (green)
## 5) deleted markers (marker falls in a deleted region).
## 6) REL288 specific marker
## 7) REL291 specific marker
## 8) REL296 specific marker
## 9) REL298 specific marker

## first, plot each lineage on a different panel. On each panel, plot both initial and final timepoint.
## second, plot (final - initial) allele frequency for each lineage.
for (l in levels(evoexp.data$lineage)) {
    ## for now filter out 'new' mutations.
    l.data <- filter(evoexp.data,lineage == l,lbl %in% c(0,1,2,6,7,8,9))
    l.panel <- ggplot(l.data,aes(x=position,y=frequency)) + geom_line() + facet_grid(generation ~ .) + theme_classic() + ggtitle(l)
    #ggsave(l.panel,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/",l,"_evolexp.pdf"))

    ## get final - initial allele frequency, and drop generation column as meaningless (b/c took the difference)
    l.delta.data <- l.data %>% arrange(position,generation) %>% group_by(position) %>%
        mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation)
    delta.panel <- ggplot(l.delta.data,aes(x=position,y=delta.freq)) + geom_line() + theme_classic() + ggtitle(l)
    ggsave(delta.panel,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/",l,"_delta.pdf"))
}

###  1) sum data across LTEE lines to look for selection (only look at dN)
parallel.delta.freq <- evoexp.data %>% filter(mut.annotation=='dN') %>%
    arrange(lineage,position,generation) %>% group_by(lineage,position) %>%
    mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation) %>%
    group_by(position) %>% summarise(sum.delta.freq=sum(delta.freq))

sum.delta.plot <- ggplot(parallel.delta.freq,aes(x=position,y=sum.delta.freq)) + geom_line() + theme_classic() + ggtitle('WHAT?')
ggsave(sum.delta.plot,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/","sumplot.pdf"))

## now omit Ara+1 and Ara-3 (only look at dN)
parallel.delta.freq2 <- evoexp.data %>% filter(mut.annotation=='dN') %>%
    filter(lineage !='Ara+1') %>%
    filter(lineage!='Ara-3') %>%
    arrange(lineage,position,generation) %>% group_by(lineage,position) %>%
    mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation) %>%
    group_by(position) %>% summarise(sum.delta.freq=sum(delta.freq))

sum.delta.plot2 <- ggplot(parallel.delta.freq2,aes(x=position,y=sum.delta.freq)) + geom_line() + theme_classic() + ggtitle('Change in K-12 allele frequency over time')
ggsave(sum.delta.plot2,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/","sumplot2.pdf"))

## now, look at mutations where the magnitude of sum delta freq is >=2.
interesting.delta.pos <- filter(parallel.delta.freq,abs(sum.delta.freq)>=2)
interesting.delta.muts <- filter(evoexp.data,position %in% interesting.delta.pos$position)
unique(interesting.delta.muts$gene.annotation)

### 2) look at delta freq in mutations found in the clones: what is the fate of these
###    lineages? these mutations should be more or less correlated.
parallel.delta.freq3 <- evoexp.data %>% filter(mut.annotation=='dN') %>%
    arrange(lineage,position,generation) %>% group_by(lineage,position) %>%
    mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation)

## use facet_grid to make a small multiple.
pop.multiple <- ggplot(data=parallel.delta.freq3,aes(x=position,y=delta.freq)) +
    geom_line() + theme_classic() + facet_grid(.  ~ lineage)
ggsave(pop.multiple,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/","delta_small_multiple.pdf"))

## plot each panel separately.
for (g in levels(labeled.mutations$genome)) {
    g.data <- filter(labeled.mutations,genome==g)
    ## This next line doesn't strictly test for equality but probably good enough.
    clone.lineage <- filter(evoexp.data,lineage %in% g.data$lineage,position %in% g.data$position)
    clone.delta.freq <- clone.lineage %>%
        arrange(position,generation) %>% group_by(position) %>%
        mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation)

    clone.panel <- ggplot(clone.delta.freq,aes(x=position,y=delta.freq)) + geom_line() + theme_classic() + ggtitle(g)
    ggsave(clone.panel,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/clone_plots/",g,"_delta.pdf"))

    ## also make a plot of mutations that are in the pop. but not in this clone.
    clone.complement <- filter(evoexp.data,lineage %in% g.data$lineage,
                               !position %in% g.data$position)

    complement.delta.freq <- clone.complement %>%
        arrange(position,generation) %>% group_by(position) %>%
        mutate(delta.freq=lead(frequency)-frequency) %>% na.omit() %>% select(-generation)

    complement.panel <- ggplot(complement.delta.freq,aes(x=position,y=delta.freq)) + geom_line() + theme_classic() + ggtitle(g)
    ggsave(complement.panel,file=paste("/Users/Rohandinho/Desktop/evolexp_plots/clone_complement_plots/",g,"_delta.pdf"))

}
