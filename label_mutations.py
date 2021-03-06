#!/usr/bin/env python

'''
label_mutations.py by Rohan Maddamsetti.

This script generates a csv file from annotated genome diff files.
It labels mutations in recombinant genomes from the end of the Souza-Turner
experiment as being
0) reference genome state.
1) B/K-12 marker is present that is not in reference genome (yellow)
2) LTEE recipient mutations (red)
3) new mutations (black)
4) LTEE recipient mutations that were erased by K-12 or otherwise missing (green)
5) deleted markers (marker falls in a deleted region).
6) REL288 specific marker
7) REL291 specific marker
8) REL296 specific marker
9) REL298 specific marker

Usage1: python label_mutations.py 1 > ../results/labeled_mutations.csv

Usage2: python label_mutations.py 2 > ../results/K-12-differences.csv

Usage3: python label_mutations.py 3 > ../results/LTEE-recipient-markers.csv

Usage4: python label_mutations.py 4 > ../results/evoexp_labeled_mutations.csv

Usage5: python label_mutations.py 5 > ../results/turner_clones_labeled_mutations.csv

An R script called dissertation_analysis.R makes figures and does stats
using these csv files (output of this script for 1,2,3,4,5)

'''

import argparse
from os.path import join, basename, exists
from os import listdir, makedirs, chdir, getcwd
import csv
import pprint
import sys
from Bio import SeqIO
from subprocess import call, Popen

## global variable.
REL606_GENOME_LENGTH = 4629812

## one line function to get name of genome.
get_genome_name = lambda gdfile: basename(gdfile).split('.')[0].split('_')[1]

def get_annotation_from_field(d,field_key):
    '''helper function to get the value associated with a given field
    in list containing the fields in one line of an annotated genome diff file.'''
    for f in d:
        if f.startswith(field_key):
            return f.split('=')[1]

def parse_annotated_gd(gdfile):
    gd_dict = {} # dict of integer:string pairs.
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

        ## skip mutations that don't correspond to REL606 reference (such as mutations
        ## in the Tn10 transposon or the F-plasmid.)
        if ref_genome != 'REL606':
            continue
        ## skip 'UN' mutations.
        if mut_type == 'UN':
            continue

        ## frequency column is printed even in consensus mode.
        frequency = get_annotation_from_field(data,"frequency")
        ## we want to get the annotation of which genes (or intergenic region)
        ## mutated, and we want to classify mutations as either dN, dS, non-coding,
        ## or non-point mutation.
        mut_annotation = ''
        gene_annotation = get_annotation_from_field(data, "gene_name")
        if mut_type == "SNP":
            snp_type = get_annotation_from_field(data,"snp_type")
            if snp_type == 'synonymous':
                mut_annotation = "dS"
            ## treat nonsense mutations as nonsynonymous mutations.
            elif snp_type == 'nonsynonymous' or snp_type== 'nonsense':
                mut_annotation = "dN"
            elif snp_type == 'intergenic':
                mut_annotation = "non-coding"
            elif snp_type == 'pseudogene':
                mut_annotation = "non-coding"
            elif snp_type == 'noncoding':
                mut_annotation = "non-coding"
        elif mut_type == "DEL":
            mut_annotation = "indel"
        elif mut_type == "INS":
            mut_annotation = "indel"
        elif mut_type == 'SUB':
            mut_annotation = 'base-substitution'
        elif mut_type == 'MOB':
            mut_annotation = 'IS-insertion'
        else:
            mut_annotation = "non-point"

        all_annotation = [mut_type,ref_genome,ref_position,mutation, mut_annotation, gene_annotation,frequency]
        ## Turn any missing annotation into NA (i.e. if frequency not donor diffs).
        all_annotation = [x if x is not None else 'NA' for x in all_annotation]
        gd_dict[int(ref_position)] = ','.join(all_annotation)
    return gd_dict

def muts_equal(mut1,mut2):
    ''' hack method, I could rewrite mutations as a class and
    overload the == operator to compare mutations in the future.
    '''
    ## cut off the last field in the string (which should be mutation frequency)
    ## and compare the two strings.
    i = mut1.rfind(',')
    j = mut2.rfind(',')
    return mut1[:i] == mut2[:j]

def label_mutations(donor_dictz, recipient_dict, recombinant_dict, recombinant_name,lineage):

    REL288_dict = donor_dictz['REL288']
    REL291_dict = donor_dictz['REL291']
    REL296_dict = donor_dictz['REL296']
    REL298_dict = donor_dictz['REL298']
    all_donor_dict = donor_dictz['donor_intersection']
    donor_coords = {k for k in set.union(*[set(donor_dictz[x].keys()) for x in donor_dictz])}

    ## I discovered that some of the K-12 donor genomes have mutations at
    ## the same coordinate, probably due to hypermutability (such as slippage).
    ## to prevent overwriting mutations with the same position, each position
    ## maps to a list of mutations at that position.

    relevant_coords = sorted(donor_coords|recipient_dict.keys()|recombinant_dict.keys())
    line_buffer = [] ## this slows down the code, but is useful in debugging.
    last_deletion_end = -1
    for i in relevant_coords:
        label = ''
        if i < last_deletion_end: ## this marker is within a deleted region in the recombinant.
            label = '5'
        elif i in recombinant_dict:
            if 'DEL' in recombinant_dict[i]:
                del_data = recombinant_dict[i].split(',')
                last_deletion_end = int(del_data[3]) + i
            if i in all_donor_dict and muts_equal(recombinant_dict[i], all_donor_dict[i]):
                label = '1' ## K-12 mutation in recombinant
            elif i in REL288_dict and muts_equal(recombinant_dict[i], REL288_dict[i]):
                label = '6' ## REL288 specific marker
            elif i in REL291_dict and muts_equal(recombinant_dict[i], REL291_dict[i]):
                label = '7' ## REL291 specific marker
            elif i in REL296_dict and muts_equal(recombinant_dict[i], REL296_dict[i]):
                label = '8' ## REL296 specific marker
            elif i in REL298_dict and muts_equal(recombinant_dict[i], REL298_dict[i]):
                label = '9' ## REL298 specific marker
            elif i in recipient_dict and muts_equal(recombinant_dict[i], recipient_dict[i]):
                label = '2' ## LTEE recipient marker in recombinant
            else:
                label = '3' ## new mutation in recombinant
        elif i in all_donor_dict or i in REL288_dict or i in REL291_dict or i in REL296_dict or i in REL298_dict:
            label = '0' ## this K-12 marker did not make it, site in REL606 state.
        elif i in recipient_dict:
            label = '4' ## this LTEE marker did not make it, erased by K-12?

        if label == '0':
            ## NOTE: an odd feature is that some K-12 markers have the same position
            ## since K-12 annotation comes from 4 donor strains, so different markers
            ## can occur at the same position. Greedily choose the donor marker at this position.
            if i in all_donor_dict:
                line_buffer.append([lineage,recombinant_name,all_donor_dict[i],label])
            elif i in REL288_dict:
                line_buffer.append([lineage,recombinant_name,REL288_dict[i],label])
            elif i in REL291_dict:
                line_buffer.append([lineage,recombinant_name,REL291_dict[i],label])
            elif i in REL296_dict:
                line_buffer.append([lineage,recombinant_name,REL296_dict[i],label])
            elif i in REL298_dict:
                line_buffer.append([lineage,recombinant_name,REL298_dict[i],label])
        elif label == '1' or label == '2' or label == '3' or label == '6' or label == '7' or label == '8' or label == '9':
            line_buffer.append([lineage,recombinant_name,recombinant_dict[i],label])
        elif label == '4': ##LTEE markers that are missing from the recombinant.
            line_buffer.append([lineage,recombinant_name,recipient_dict[i],label])
        elif label == '5': ## find the dict that the deleted marker is in, and add to the file.
            all_dicts = [recombinant_dict,recipient_dict,all_donor_dict,REL288_dict,REL291_dict,REL296_dict,REL298_dict]
            for some_dict in all_dicts:
                if i in some_dict:
                    line_buffer.append([lineage,recombinant_name,some_dict[i],label])
                    break

    if last_deletion_end > REL606_GENOME_LENGTH:
        print('ERROR: DELETION SPANS BEGINNING AND END OF GENOME FILE!')
        quit()

    for l in line_buffer:
        print(','.join(l))

def make_ref_labeled_muts_csv(data_dir,analysis='recombinants'):
    ''' data_dir contains folders called Ara+1,...,Ara+6,...,Ara-1,...,Ara-6.
    each of those folders contains an annotated diff of the recipient strains
    and annotated diffs of the recombinant strains. We want to track markers specific
    to each of the 4 Hfr K-12 donor strains.
    '''

    if analysis not in ('recombinants','evoexp','turner_clones'):
        print('ERROR: UNKNOWN ANALYSIS:', analysis)
        quit()

    ## print header. lbl is used instead of label to avoid confusion with this keyword in ggplot2.
    print("lineage,genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation,frequency,lbl")

    donor_marker_dir = "../donor-specific-markers/"

    REL288_gd = join(data_dir,donor_marker_dir,'annotated_REL288_specific_markers.gd')
    REL291_gd = join(data_dir,donor_marker_dir,'annotated_REL291_specific_markers.gd')
    REL296_gd = join(data_dir,donor_marker_dir,'annotated_REL296_specific_markers.gd')
    REL298_gd = join(data_dir,donor_marker_dir,'annotated_REL298_specific_markers.gd')
    donor_intersection_gd = join(data_dir,donor_marker_dir,'annotated_donor_intersection.gd')

    REL288_dict = parse_annotated_gd(REL288_gd)
    REL291_dict = parse_annotated_gd(REL291_gd)
    REL296_dict = parse_annotated_gd(REL296_gd)
    REL298_dict = parse_annotated_gd(REL298_gd)
    donor_dict = parse_annotated_gd(donor_intersection_gd)

    donor_dictz = {'REL288':REL288_dict, 'REL291':REL291_dict,'REL296':REL296_dict,'REL298':REL298_dict,'donor_intersection':donor_dict}

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

    ## note that for the evolution experiment the endpoints don't have an REL designation,
    ## and instead of even or odd clones, there's the T=1000 and T=1200 pop. samples.
    if analysis == 'recombinants':
        for l in lineages:
            lineage_dir = join(data_dir,l)
            lineage_diffs = [x for x in listdir(lineage_dir) if x.endswith('.gd')]
            for x in lineage_diffs:
                if 'REL25' in x:
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

            label_mutations(donor_dictz, recipient_dict,odd_recombinant_dict, odd_name,l)
            label_mutations(donor_dictz, recipient_dict,even_recombinant_dict, even_name,l)

    elif analysis == 'evoexp':
        for l in lineages:
            lineage_dir = join(data_dir,l)
            lineage_diffs = [x for x in listdir(lineage_dir) if x.endswith('.gd')]
            ## VERY CHEAP HACK FOR NOW TO GET RECIPIENT FROM REL606-ref-runs:
            my_recipient = [x for x in listdir(join(data_dir,"../REL606-ref-runs",l)) if x.endswith('.gd') and 'REL25' in x][0]
            recipient_gd = join(data_dir,"../REL606-ref-runs",l,my_recipient)
            lineage_diffs = [x for x in listdir(lineage_dir) if x.endswith('.gd')]

            for x in lineage_diffs:
                if 'RM3-149' in x:
                    initial_gd = join(lineage_dir,x)
                elif 'RM3-153' in x:
                    final_gd = join(lineage_dir,x)
                else:
                    print('ERROR IN EVOEXP INPUT')
                    print(x)
                    quit()
            recipient_dict = parse_annotated_gd(recipient_gd)
            initial_dict = parse_annotated_gd(initial_gd)
            initial_name = get_genome_name(initial_gd)
            final_dict = parse_annotated_gd(final_gd)
            final_name = get_genome_name(final_gd)

            label_mutations(donor_dictz, recipient_dict, initial_dict, initial_name,l)
            label_mutations(donor_dictz, recipient_dict, final_dict, final_name,l)

    elif analysis == 'turner_clones':
        REL4397_gd = join(data_dir,'turner_clones','annotated_REL4397.gd')
        REL4398_gd = join(data_dir,'turner_clones','annotated_REL4398.gd')
        recipient_gd = join(data_dir,'Ara-3/annotated_REL2545.gd')
        REL4397_dict = parse_annotated_gd(REL4397_gd)
        REL4398_dict = parse_annotated_gd(REL4398_gd)
        recipient_dict = parse_annotated_gd(recipient_gd)

        REL4397_name = get_genome_name(REL4397_gd)
        REL4398_name = get_genome_name(REL4398_gd)

        label_mutations(donor_dictz, recipient_dict, REL4397_dict, REL4397_name,'Ara-3')
        label_mutations(donor_dictz, recipient_dict, REL4398_dict, REL4398_name,'Ara-3')


def make_K12_diff_csv(data_dir):
    donor_gd = join(data_dir,'annotated_donor_union.gd')
    donor_name = get_genome_name(donor_gd)
    donor_dict = parse_annotated_gd(donor_gd)
    ## print header.
    print("genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation,frequency")
    for i in sorted(donor_dict.keys()):
        x = donor_name + ',' + donor_dict[i]
        print(x)

def make_LTEE_marker_csv(data_dir):
    lineages = ['Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6',
                'Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6']
    ## print header.
    print("lineage,genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation,frequency")
    for l in lineages:
        lineage_dir = join(data_dir,l)
        lineage_diffs = [x for x in listdir(lineage_dir) if x.endswith('.gd')]
        for x in lineage_diffs:
            if 'REL25' in x :
                recipient_gd = join(lineage_dir,x)
                recipient_dict = parse_annotated_gd(recipient_gd)
                recipient_name = get_genome_name(recipient_gd)
                for i in sorted(recipient_dict.keys()):
                    print(recipient_name+','+recipient_dict[i])

def make_donor_csv(data_dir):
    print("genome,mut.type,reference,position,mutation,mut.annotation,gene.annotation,frequency")
    donor_diffs = [x for x in listdir(data_dir) if x.endswith('.gd')]
    for x in donor_diffs:
        donor_gd = join(data_dir,x)
        donor_dict = parse_annotated_gd(donor_gd)
        donor_name = get_genome_name(donor_gd)
        for i in sorted(donor_dict.keys()):
            print(donor_name+','+donor_dict[i])

def main():
    REL606_gd_dir = '../annotated-diffs/REL606-ref-runs'
    evo_exp_gd_dir = '../annotated-diffs/evoexp'
    donor_gd_dir = '../annotated-diffs/donor-specific-markers'

    parser = argparse.ArgumentParser(description='options run different sections of analysis for genomes from STLE' + \
    'that have been run through breseq using REL606 as reference genome')

    parser.add_argument('analysis_stage', type=int, help='number representing which function to do', action='store')
    args = parser.parse_args()
    if args.analysis_stage == 1:
        make_ref_labeled_muts_csv(REL606_gd_dir,analysis='recombinants')
    elif args.analysis_stage == 2:
        make_K12_diff_csv(donor_gd_dir)
    elif args.analysis_stage == 3:
            make_LTEE_marker_csv(REL606_gd_dir)
    elif args.analysis_stage == 4:
        make_ref_labeled_muts_csv(evo_exp_gd_dir,analysis='evoexp')
    elif args.analysis_stage == 5:
        make_ref_labeled_muts_csv(REL606_gd_dir,analysis='turner_clones')


main()
