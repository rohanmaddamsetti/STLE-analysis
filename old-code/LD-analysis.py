#!/usr/bin/env python3

## LD-analysis.py by Rohan Maddamsetti.

## By default, this can only be used with K-12 mutations mapped onto REL606,
## since it reads annotated_K-12.gd in hardwritten code.

## Usage: python3 LD-analysis.py ../annotated-diffs/K-12_minus_REL4397.gd

import fileinput
from Bio import SeqIO

diffloc = "/Users/Rohandinho/Desktop/Projects/recombinant-assembly/annotated-diffs/"
refgenome = "/Users/Rohandinho/Desktop/REL606.6.gbk"

def main():

    absence_dict = {}
    allK12file = open(diffloc + "annotated_K-12.gd","r")
    for line in allK12file:
        if line.startswith('#'):
            continue
        data = line.split()
        pos = int(data[4])
        absence_dict[pos] = 0
    for line in fileinput.input():
        if line.startswith('#'):
            continue
        data = line.split()
        pos = int(data[4])
        absence_dict[pos] = 1

    ## now, find windows. the starts and ends are K-12 markers that DID make it in.
    ## TODO: rewrite to use midpoint method.

    windows = []
    start = 0
    end = 0
    for x in sorted(absence_dict.keys()):
        if absence_dict[x] == 1:
            end = x
        else:
            if end > start: # spit out the current window.
                end = x
                # empty lists are for locus and gene data to be added later.
                #windows.append([start,end,[],[]])
                windows.append([start,end])
                start = x
            else: ## continue
                start = x
                end = start

    ## print output.
    print("Start, End")
    for x in windows:
        print(','.join([str(field) for field in x]))
    quit() # don't deal with genes and loci for the time being.

    ## Use Biopython to get loci and genes within each window.
    gene_starts = {}
    gene_ends = {}
    for record in SeqIO.parse(open(refgenome), "genbank"):
        for f in record.features:
            if f.type != 'CDS':
                continue
            locus_tag = f.qualifiers['locus_tag'][0]
            if 'gene' in f.qualifiers:
                gene = f.qualifiers['gene'][0]
            else:
                gene = locus_tag
            g_start, g_end = int(f.location.start), int(f.location.end)
            gene_starts[g_start] = (locus_tag,gene)            
            gene_ends[g_end] = (locus_tag,gene)
    ## Annotate windows with locus_tag and gene data.
    for win in windows:
        for j in sorted(gene_ends.keys()):
            if j < win[0] or j > win[1]:
                continue
            else:
                if gene_ends[j][0] not in win[2]: # then add locus and gene.
                    win[2].append(gene_ends[j][0])
                    win[3].append(gene_ends[j][1])
        for i in sorted(gene_starts.keys()):
            if i < win[0] or i > win[1]:
                continue
            else:
                if gene_starts[i][0] not in win[2]: # then add locus and gene.
                    win[2].append(gene_starts[i][0])
                    win[3].append(gene_starts[i][1])
    ## print output.
    print("Start, End, Loci, Genes")
    for x in windows:
        print(','.join([str(field) for field in x]))

main()

