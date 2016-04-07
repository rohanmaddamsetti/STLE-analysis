#!/usr/bin/env python

## label_mutations.py by Rohan Maddamsetti.

## This script generates a csv file from annotated genome diff files.
## It labels mutations in recombinant genomes from the end of the Souza-Turner
## experiment as being
## 0) missing markers between REL606 and K-12 (white/no color)
## 1) present K-12 markers (yellow)
## 2) LTEE recipient mutations (red)
## 3) new mutations (black)
## 4) LTEE recipient mutations that were erased by K-12 or otherwise missing (green)
## 5) deleted markers (marker falls in a deleted region).

## Usage1: python label_mutations.py 1 > ../results/labeled_mutations.csv

## Usage2: python label_mutations 2 > ../results/K-12-differences.csv

## Usage3: python label_mutations 3 > ../results/LTEE-recipient-markers.csv

## An R script called dissertation_analysis.R makes figures and does stats
## using these three csv files (output of this script)

## TODO: work with K-12 reference rather than REL606 reference.

## TODO: add options that do different stages of the analysis.
## TODO: automate the whole pipeline with a makefile.

## NOTA BENE: this is a rough analysis, and the input gd files from breseq resequencing aren't standardized yet. In particular, K-12 assemblies used breseq 0.26, and assemblies may or may not have used an additional Tn10 reference sequence.

## This usage is a wrapper for gdtools and parallyze to look for parallel
## positive selection on new mutations in the Souza-Turner experiment.
## Usage4: python label_mutations 4

import argparse
from os.path import join, basename, exists
from os import listdir, makedirs, chdir, getcwd
import csv
import pprint
import sys
from Bio import SeqIO
from subprocess import call, Popen

## one line function to get name of genome.
get_genome_name = lambda gdfile: basename(gdfile).split('.')[0].split('_')[1]

def get_annotation_from_field(d,field_key):
    '''helper function to get the value associated with a given field
    in list containing the fields in one line of an annotated genome diff file.'''
    for f in d:
        if f.startswith(field_key):
            return f.split('=')[1]

def parse_annotated_gd(gdfile):

    gd_dict = {}
    gd_handle = open(gdfile, 'r')
    for line in gd_handle:
        if "gene_name=araA" in line and "aa_position=92" in line and "snp_type=nonsynonymous" in line:
            continue
        elif "gene_name=recD" in line and "aa_position=10" in line and "snp_type=nonsynonymous" in line:
            continue
        elif line.startswith('#'):
            continue
        line = line.strip()
        data = line.split()
        mut_type = data[0]
        ref_genome = data[3]
        ref_position = data[4]
        mutation = data[5]
        ## we want to get the annotation of which genes (or intergenic region)
        ## mutated, and we want to classify mutations as either dN, dS, non-coding,
        ## or non-point mutation.
        mut_annotation = ''
        gene_annotation = ''
        if mut_type == "SNP":
            gene_annotation = get_annotation_from_field(data, "gene_name")
            if "snp_type=synonymous" in line:
                mut_annotation = "dS"
            ## treat nonsense mutations as nonsynonymous mutations.
            elif "snp_type=nonsynonymous" or "snp_type=nonsense" in line:
                mut_annotation = "dN"
            elif "snp_type=intergenic" in line:
                mut_annotation = "non-coding"
            elif "snp_type=pseudogene" in line:
                mut_annotation = "non-coding"
        else:
            gene_annotation = get_annotation_from_field(data, "gene_name")
            mut_annotation = "non-point"
        ## I discovered that some of the K-12 donor genomes have mutations at
        ## the same coordinate, probably due to hypermutability (such as slippage).
        ## to prevent overwriting mutations with the same position, each position
        ## maps to a list of mutations at that position.
        if int(ref_position) not in gd_dict:
            gd_dict[int(ref_position)] = []
        gd_dict[int(ref_position)].append(','.join([mut_type,ref_genome,ref_position,mutation, mut_annotation, gene_annotation]))
    return gd_dict

def label_mutations(donor_dict, recipient_dict, recombinant_dict, recombinant_name,lineage):
    relevant_coords = sorted(donor_dict.keys()|recipient_dict.keys()|recombinant_dict.keys())
    line_buffer = [] ## this is to handle deletions that span the end and beginning.
    in_deletion = False
    ## in_deletion is false by assumption-- code at end handles the case that
    ## in_deletion is actually true at the beginning of the loop, in which case
    ## the labels for the first mutations should be set to the deleted state '5'.

    for i in relevant_coords:
        label = ''
        if i in recombinant_dict:
            if i in donor_dict and recombinant_dict[i][0] in donor_dict[i]:
                label = '1' ## K-12 mutation in recombinant
            elif i in recipient_dict and recombinant_dict[i][0] in recipient_dict[i]:
                label = '2' ## LTEE marker in recombinant
            else:
                label = '3' ## new mutation in recombinant
        elif i in donor_dict:
            label = '0' ## this K-12 marker did not make it, site in REL606 state.
        elif i in recipient_dict:
            label = '4' ## this LTEE marker did not make it, erased by K-12?

        if label == '0':
            ## NOTE: an odd feature is that some K-12 markers have the same position
            ## since K-12 annotation comes from 4 donor strains, so different markers
            ## can occur at the same position.
            for x in donor_dict[i]:
                line_buffer.append([lineage,recombinant_name,x,label])
        elif label == '1' or label == '2' or label == '3':
            line_buffer.append([lineage,recombinant_name,recombinant_dict[i][0],label])
        elif label == '4': ##LTEE markers that are missing from the recombinant.
            line_buffer.append([lineage,recombinant_name,recipient_dict[i][0],label])

    for l in line_buffer:
        print(','.join(l))

def make_REL606_ref_labeled_muts_csv(data_dir):
    ''' data_dir contains folders called Ara+1,...,Ara+6,...,Ara-1,...,Ara-6.
    each of those folders contains an annotated diff of the recipient strains
    and annotated diffs of the recombinant strains. data_dir also contains
    an annotated diff of all differences between K-12 and REL606.
    '''

    ## print header. lbl is used instead of label to avoid confusion with keyword in ggplot2.
    print("lineage,genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation,lbl")

    donor_gd = join(data_dir,'annotated_K-12.gd')
    donor_dict = parse_annotated_gd(donor_gd)

    lineages = ['Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6',
                'Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6']

    ## to convert my IDs to REL IDs.
    rel_name = {'RM3-130-1':'REL11734','RM3-130-2':'REL11735',
                'RM3-130-3':'REL11736','RM3-130-4':'REL11737',
                'RM3-130-5':'REL11738','RM3-130-6':'REL11739',
                'RM3-130-7':'REL11740','RM3-130-8':'REL11741',
                'RM3-130-9':'REL11742','RM3-130-10':'REL11743',
                'RM3-130-11':'REL11744','RM3-130-12':'REL11745',
                'RM3-130-13':'REL11746','RM3-130-14':'REL11747',
                'RM3-130-15':'REL11748','RM3-130-16':'REL11749',
                'RM3-130-17':'REL11750','RM3-130-18':'REL11751',
                'RM3-130-19':'REL11752','RM3-130-20':'REL11753',
                'RM3-130-21':'REL11754','RM3-130-22':'REL11755',
                'RM3-130-23':'REL11756','RM3-130-24':'REL11757'}

    for l in lineages:
        lineage_dir = join(data_dir,l)
        lineage_diffs = listdir(lineage_dir)
        for x in lineage_diffs:
            if 'REL25' in x :
                recipient_gd = join(lineage_dir,x)
            elif int(x[-4]) % 2:
                odd_recombinant_gd = join(lineage_dir,x)
            else:
                even_recombinant_gd = join(lineage_dir,x)
        recipient_dict = parse_annotated_gd(recipient_gd)
        odd_recombinant_dict = parse_annotated_gd(odd_recombinant_gd)
        odd_name = rel_name[get_genome_name(odd_recombinant_gd)]
        even_recombinant_dict = parse_annotated_gd(even_recombinant_gd)
        even_name = rel_name[get_genome_name(even_recombinant_gd)]
        label_mutations(donor_dict, recipient_dict,odd_recombinant_dict, odd_name,l)
        label_mutations(donor_dict, recipient_dict,even_recombinant_dict, even_name,l)

def make_K12_diff_csv(data_dir):
    donor_gd = join(data_dir,'annotated_K-12.gd')
    donor_name = get_genome_name(donor_gd)
    donor_dict = parse_annotated_gd(donor_gd)
    ## print header.
    print("genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation")
    for i in sorted(donor_dict.keys()):
        for x in donor_dict[i]:
            print(','.join([donor_name,x]))

def make_LTEE_marker_csv(data_dir):

    lineages = ['Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6',
                'Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6']
    ## print header.
    print("lineage,genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation")
    for l in lineages:
        lineage_dir = join(data_dir,l)
        lineage_diffs = listdir(lineage_dir)
        for x in lineage_diffs:
            if 'REL25' in x :
                recipient_gd = join(lineage_dir,x)
                recipient_dict = parse_annotated_gd(recipient_gd)
                recipient_name = get_genome_name(recipient_gd)
                for i in sorted(recipient_dict.keys()):
                    for x in recipient_dict[i]:
                        print(','.join([l,recipient_name,x]))

def run_parallyze_on_new_mutations(proj_dir,gd_dir):
    '''this function is a wrapper for gdtools and parallyze--it uses
    gdtools SUBTRACT to filter only new mutations in the recombinant
    genome diffs, writes the filtered diffs to results/new-mutation-diffs/,
    writes an input file for parallyze, and then runs parallyze.
    The point of this analysis is to ask whether there is evidence of parallel
    evolution among new mutations in the recombinant clones.'''

    data_dir = join(proj_dir,gd_dir)
    output_path = join(proj_dir, 'results/new_mutation_diffs/REL606-ref')
    if not exists(output_path):
        makedirs(output_path)

    ####### NOTE: even recombinants are commented out, because odds OR evens
    ###### are independent genomes. this wouldn't matter if the parallyze option
    ###### that allows related genomes within a lineage to be included worked.
    ###### I could see what Jeff and Olivier do in the 264 genome project
    ###### to get this option working.
    donor_gd = join(proj_dir,gd_dir,'annotated_K-12.gd')
    lineages = ['Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6',
                'Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6']
    for l in lineages:
        lineage_dir = join(data_dir,l)
        lineage_diffs = listdir(lineage_dir)
        for x in lineage_diffs:
            if 'REL25' in x :
                recipient_gd = join(lineage_dir,x)
            elif int(x[-4]) % 2:
                odd_recombinant_gd = join(lineage_dir,x)
            else:
                even_recombinant_gd = join(lineage_dir,x)
        odd_name = get_genome_name(odd_recombinant_gd)
        ##even_name = get_genome_name(even_recombinant_gd)
        odd_new_mut_output = join(output_path,odd_name+'_new_muts.gd')
        ##even_new_mut_output = join(output_path,even_name+'_new_muts.gd')
        if not exists(odd_new_mut_output):
            call(["gdtools", "SUBTRACT", odd_recombinant_gd, donor_gd, recipient_gd, "-o", odd_new_mut_output])
        ##if not exists(even_new_mut_output):
        ##    call(["gdtools", "SUBTRACT", even_recombinant_gd, donor_gd, recipient_gd, "-o", even_new_mut_output])

    ## do a sanity check to make sure that the number of new mutations in the gd files
    ## matches the number of new mutations labeled in /results/labeled_mutations.csv.
    ##new_mutation_sanity_check(proj_dir,lineages, output_path)

    ## now, generate an configuration file for parallyze.
    conf_path = join(output_path,"parallyze.conf")
    with open(conf_path,'w') as pf:
        print("# parallyze.conf for new mutation analysis",file=pf)
        print("# of odd recombinant genomes in Souza-Turner experiment",file=pf)
        print("REF_GENOME ",join(proj_dir, "references/REL606.6.gbk"),file=pf)
        print("GENOMEDIFF_FILES",file=pf)
        for f in listdir(output_path):
            if f.endswith(".gd"):
                print(join(output_path,f),file=pf)
        print("NONSYNONYMOUS True",file=pf)
        print("SYNONYMOUS False",file=pf)
        print("NONCODING False",file=pf)
        print("PSEUDOGENE False",file=pf)
        print("INTERGENIC False",file=pf)
        print("GENE_PRODUCT False",file=pf)
        print("REPLICATES 1000",file=pf)
        print("GENES_TO_DISPLAY 10",file=pf)
    pf.close()
    parallyze_path = "/Users/Rohandinho/Desktop/Projects/parallyze/parallyze.py"
    call([parallyze_path, "--config", conf_path, "--procedure", "mutationTally"])

def new_mutation_sanity_check(proj_dir,lineages,output_path,count_only_odd=True):
    for f in listdir(output_path):
        call(["wc", join(output_path,f)])
    labeled_muts = open(join(proj_dir,"results/labeled_mutations.csv"))
    new_mut_count = {x:0 for x in lineages}
    for i,l in enumerate(labeled_muts):
        if i == 0: ## skip the header
            continue
        data = l.strip().split(',')
        my_lineage = data[0]
        my_clone = data[1]
        my_label = data[-1]
        if my_label == '3': ## new mutation
            ## in this mod expression, even clones evaluate to false
            if count_only_odd and not int(my_clone[-1]) % 2:
                continue
            #print(data)
            new_mut_count[my_lineage] += 1
    print(new_mut_count)

def main():
    proj_dir = '/Users/Rohandinho/Desktop/Projects/recombinant-assembly/'
    REL606_gd_dir = 'annotated-diffs/REL606-ref-runs'

    parser = argparse.ArgumentParser(description='options run different sections of analysis for genomes from Souza-Turner experiment')
    parser.add_argument('analysis_stage', type=int, help='number representing which function to do', action='store')
    args = parser.parse_args()
    if args.analysis_stage == 1:
        make_REL606_ref_labeled_muts_csv(join(proj_dir,REL606_gd_dir))
    elif args.analysis_stage == 2:
        make_K12_diff_csv(join(proj_dir,REL606_gd_dir))
    elif args.analysis_stage == 3:
        make_LTEE_marker_csv(join(proj_dir,REL606_gd_dir))
    elif args.analysis_stage == 4:
        run_parallyze_on_new_mutations(proj_dir,REL606_gd_dir)

main()
