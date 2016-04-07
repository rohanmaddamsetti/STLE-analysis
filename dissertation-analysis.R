## dissertation-analysis.R by Rohan Maddamsetti.

## This script makes figures and does statistics for the recombinant genome analysis.

library(ggplot2)
library(ggrepel)
library(dplyr)

labeled.mutations <- tbl_df(read.csv("../results/labeled_mutations.csv"))
labeled.mutations <- mutate(labeled.mutations,lbl=as.factor(lbl))
G.score.data <- tbl_df(read.csv("../results/036806-3.csv"))
hand.annotation <- tbl_df(read.csv("../results/donor_hand_annotation.csv"))

###############################################
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
               x1=position,x2=position,y1=y.center-3,y2=y.center+3)

 ## label the top 54 G-scoring genes with symbols.
    ## below this threshold, very short genes with 1 dN have a high G-score.
    top.G.score.genes <- filter(G.score.data,G.score>8.9)

    top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name)

    ## draw and label vertical lines on the plot for high G-score mutations.
    genes.to.label <- filter(top.G.score.genes,Gene.name %in% top.hits$gene.annotation) %>%
        mutate(Gene.index=row_number())

    ## label these mutations as points in the plot.
    LTEE.top.hits <- filter(no.B.genomes, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4' | lbl=='3') %>%
        mutate(y.center=70-match(genome,unique(genome))*10,
               x1=position,x2=position,y1=y.center-3,y2=y.center+3)

    ## this is for plotting the labels of genes mentioned in the text.
    #genes.to.label <- mutate(genes.to.label,label.name=Gene.name) %>%
    #    mutate(label.name=plyr::revalue(label.name,c('mreC'='mreC/B')))

    panel <- ggplot(no.B.genomes, aes(x=x1,xend=x2,y=y1,yend=y2,colour=lbl)) + geom_segment(size=0.05) +
        scale_colour_manual(values=c('yellow', 'red', 'black','green')) + theme_classic() + xlab("Position") + ylab("Clone") + ylim(5,70) + guides(colour=FALSE) + scale_y_continuous(breaks=c(10,20,30,40,50,60),labels=rev(unique(no.B.genomes$genome))) +
        geom_point(data=LTEE.top.hits,aes(x=position,y=y.center,color=lbl,shape=mut.annotation)) + geom_segment(data=genes.to.label, aes(x=Start.position,xend=Start.position,y=1,yend=65),size=0.1,inherit.aes=FALSE) + guides(shape=FALSE) +
        ## number genes near at the bottom.
        geom_text(data=filter(genes.to.label,Gene.index %in% c(1,6,11,16,21,26,31)),aes(label=Gene.index,x=Start.position,y=0),size=2,inherit.aes=FALSE) #+
        ## label genes that are in the text at the top.
        #geom_text(data=filter(genes.to.label, label.name %in% c('ybaL', 'fabF', 'topA', 'pykF', 'mreC/B', 'malT', 'spoT', 'hslU', 'iclR', 'nadR', 'arcA')),aes(label=Gene.name,x=Start.position, y=67),size=3,angle=45,inherit.aes=FALSE)

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


## Calculation that 31 genes that are in Figure 1 and in Table 2.
## the top 54 G-scoring genes: below this threshold, very short genes with 1 dN have a high G-score.
top.G.score.genes <- filter(G.score.data,G.score>8.9)
top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name)
LTEE.top.hits <- filter(labeled.mutations, gene.annotation %in% top.G.score.genes$Gene.name) %>%
        filter(lbl=='2' | lbl=='4')
genes.to.label <- filter(top.G.score.genes,Gene.name %in% LTEE.top.hits$gene.annotation)
## gene.to.label contains 31 genes.



#############################################################
## Figures 2. Visual comparisons of parallel recombination across lineages
## with the G-score of mutations over the genome, and with the occurrence of
## new mutations over the genome.
## This is only made with odd genomes at the moment.

##annotate auxotrophs and Hfr oriTs on the plot.
auxotrophs <- filter(hand.annotation,annotation!="Hfr"&annotation!="oriC")
hfrs.and.origin <- filter(hand.annotation,annotation=="Hfr"|annotation=="oriC")

## label genomes as odd or even.
genome.names <- levels(labeled.mutations$genome)
is.odd <- sapply(genome.names, function(x) ifelse(strtoi(substr(x,nchar(x),nchar(x)))%%2, TRUE,FALSE))
labeled.mutations <- mutate(labeled.mutations,odd=is.odd[genome])
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
    summarize(introgression.score = sum(labels.to.int[lbl]))

## omit mutator lineages (Ara+6,Ara+3,Ara-2) +6 is mutT, +3 is mutS, -2 is mutL mutator.
new.odd.mutations <- filter(odd.genomes,lineage != "Ara+6" & lineage != "Ara+3" & lineage != "Ara-2") %>%
    filter(lbl=='3')

################################
## Look for parallelism in new mutations: super strong parallelism (multiple new mutations in the same
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


#######################################
Fig2 <- ggplot(scored.odd.genomes, aes(x=position,y=introgression.score)) + geom_line() + theme_classic() + xlab("Position") + ylab("Introgression Score") +
    ## add auxotroph lines
    geom_vline(data=auxotrophs,aes(xintercept=position),size=0.1,linetype="dashed") +
    geom_text_repel(data=auxotrophs,aes(x=position,y=8.2,label=annotation),inherit.aes=FALSE, size=2) +
    ## add annotation of Hfr and oriC on chromosome.
   geom_text_repel(data=hfrs.and.origin,aes(x=position,y=-0.3,label=Hfr.orientation),inherit.aes=FALSE,size=2,nudge_y=-0.3)

ggsave("/Users/Rohandinho/Desktop/Fig2.pdf",Fig2)


#######################################################
## Map recombination breakpoints. Currently breakpoints occur ON
## mutations, NOT between mutations.
## Make Figure 3 (chunk length distributions) and do two statistical tests:
## Distance between new mutations and breakpoints (are breakpoints mutagenic?)
## recombination hot and cold spots

##global variable used in Fig. 3 code.
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

K12.chunks <- filter(odd.genome.chunks,chunk.transitions=='2-1')
LTEE.chunks <- filter(odd.genome.chunks,chunk.transitions=='1-2')

Fig3A <- ggplot(K12.chunks, aes(x=log10(chunk.length),fill=lineage)) + geom_histogram() + facet_grid(lineage ~ .) + theme_classic() + xlab("log10(donor chunk length)") + ylab("Count")

ggsave("/Users/Rohandinho/Desktop/Fig3A.pdf",Fig3A)

Fig3B <- ggplot(LTEE.chunks,aes(x=log10(chunk.length),fill=lineage)) + geom_histogram() + facet_grid(lineage ~ .) + theme_classic() + xlab("log10(recipient chunk length)") + ylab("Count")

ggsave("/Users/Rohandinho/Desktop/Fig3B.pdf",Fig3B)

## NOTE! This test was done by inputting odd.genomes2 in the group_by call
## for odd.genome.chunks to remove mutators and obviously different genomes.
## STATISTICAL TEST:
## are the distributions of K-12 chunks (or LTEE chunks)
## identical across replicate recombinant lines? Answer: similar, but not draws
## from the same distribution.

## at least one lineage has a different distribution of K-12 chunk lengths
kruskal.test(chunk.length ~ lineage, data=K12.chunks)

## distribution of LTEE chunk lengths are somewhat similar across lineages.
kruskal.test(chunk.length ~ lineage, data=LTEE.chunks)

## omit mutators and weird ones.
skip.me <- c("Ara+2","Ara-3","Ara-2", "Ara+3", "Ara+6")
kruskal.test(chunk.length ~ lineage, data=filter(K12.chunks,!lineage %in% skip.me))
kruskal.test(chunk.length ~ lineage, data=filter(LTEE.chunks,!lineage %in% skip.me))

## STATISTICAL TEST:
## Is there evidence that recombination breakpoints are mutagenic?
## TODO: find out what the average minimum distance is between new mutations
## labeled '3' and breakpoints.
## Then, drop mutations at random, and calculate average minimum distance to breakpoints.
## repeat 10000 times to calculate a p-value.
odd.breaks.and.muts <- filter(odd.genome.chunks, chunk.transitions=='2-1' | chunk.transitions=='1-2' | lbl == '3')
odd.breaks.and.muts <- select(odd.genome.chunks,lineage,genome,mut.type,position,mut.annotation,chunk.transitions,lbl)


############################
## Additional figures.

## Figure 4. Density of differences between K-12 and REL606.

K12.diff.data <- tbl_df(read.csv("../results/K-12-differences.csv"))

Fig4 <- ggplot(K12.diff.data, aes(x=position)) + geom_histogram(bins=400) + theme_classic() + xlab("Position") + ylab("Count")

## Figure S1. Density of LTEE mutations in recipient clones.

LTEE.marker.data <- tbl_df(read.csv("../results/LTEE-recipient-markers.csv"))
## omit mutator lineages (Ara+6,Ara+3,Ara-2)
LTEE.marker.data <- filter(LTEE.marker.data,lineage != "Ara+6" & lineage != "Ara+3" & lineage != "Ara-2")
FigS1 <- ggplot(LTEE.marker.data, aes(x=position),colour='purple') + geom_histogram(bins=400) + theme_classic() + xlab("Position") + ylab("Count")

## Figure S2. plot G-scores of mutations over the genome.
FigS2 <- ggplot(G.score.data, aes(x=Start.position,y=G.score)) + geom_point(size=0.5) + theme_classic() + geom_text(data=filter(G.score.data,G.score>45),aes(label=Gene.name),size=2,hjust=-0.2) + xlab("Position") + ylab("Count")


ggsave("/Users/Rohandinho/Desktop/Fig4.pdf",Fig4)
ggsave("/Users/Rohandinho/Desktop/FigS1.pdf",FigS1)
ggsave("/Users/Rohandinho/Desktop/FigS2.pdf",FigS2)
