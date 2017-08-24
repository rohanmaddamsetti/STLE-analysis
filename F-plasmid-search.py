#!/usr/bin/env python

'''
F-plasmid-search.py by Rohan Maddamsetti.

Run this on orchestra! source activate evcouplings_env first.

This script writes out csv files of F plasmid coverage for clones and
for evolution experiment samples. Graphs of these data are made in
dissertation-analysis.R.

'''

from os import makedirs, listdir, chdir, getcwd
from os.path import exists, join, basename, isdir
import pandas as pd

def make_STLE_clone_coverage_df():
    '''
    This function works on the STLE clones. A separate function works on the
    evolution experiment data, due to the different directory stucture.
    '''

    stle_dir = "../breseq-assemblies/REL606-ref-runs/"

    lineage_dict = {'REL288':'Donor', 'REL291':'Donor','REL296':'Donor','REL298':'Donor',
                    'REL2537':'Ara+1', 'REL2538':'Ara+2',
                    'REL2539':'Ara+3', 'REL2540':'Ara+4',
                    'REL2541':'Ara+5', 'REL2542':'Ara+6',
                    'REL2543':'Ara-1', 'REL2544':'Ara-2',
                    'REL2545':'Ara-3', 'REL2546':'Ara-4',
                    'REL2547':'Ara-5', 'REL2548':'Ara-6',
                    'RM3-130-1':'Ara+1', 'RM3-130-2':'Ara+1',
                    'RM3-130-3':'Ara+2', 'RM3-130-4':'Ara+2',
                    'RM3-130-5':'Ara+3', 'RM3-130-6':'Ara+3',
                    'RM3-130-7':'Ara+4', 'RM3-130-8':'Ara+4',
                    'RM3-130-9':'Ara+5', 'RM3-130-10':'Ara+5',
                    'RM3-130-11':'Ara+6', 'RM3-130-12':'Ara+6',
                    'RM3-130-13':'Ara-1', 'RM3-130-14':'Ara-1',
                    'RM3-130-15':'Ara-2', 'RM3-130-16':'Ara-2',
                    'RM3-130-17':'Ara-3', 'RM3-130-18':'Ara-3',
                    'REL4397':'Ara-3', 'REL4398':'Ara-3',
                    'RM3-130-19':'Ara-4', 'RM3-130-20':'Ara-4',
                    'RM3-130-21':'Ara-5', 'RM3-130-22':'Ara-5',
                    'RM3-130-23':'Ara-6', 'RM3-130-24':'Ara-6' }

    clone_col = []
    strain_type_col = []
    lineage_col = []
    unique_cov_col = []
    position_col = []

    for l1 in [x for x in listdir(stle_dir) if not x.startswith(".")]:
        if l1 == 'donors':
            strain_type = "Donor"
        elif l1 == 'recipients':
            strain_type = "Recipient"
        elif l1 == 'recombinants':
            strain_type = "Recombinant"
        for l2 in listdir(join(stle_dir,l1)):
            if l2 not in lineage_dict:
                continue
            else:
                lineage = lineage_dict[l2]
                coverage_f = join(stle_dir,l1,l2,"08_mutation_identification/NC_002483.coverage.tab")
                coverage_fh = open(coverage_f)
                for i,line in enumerate(coverage_fh):
                    if i == 0:
                        continue
                    line_data = line.split()
                    unique_coverage = int(line_data[0]) + int(line_data[1])
                    position = int(line_data[-1])
                    ## Now add row data into each column list.
                    strain_type_col.append(strain_type)
                    lineage_col.append(lineage)
                    unique_cov_col.append(unique_coverage)
                    position_col.append(position)
                    clone_col.append(l2)

    df = pd.DataFrame( { 'Clone':clone_col,
        'Strain Type' : strain_type_col,
        'Lineage' : lineage_col,
        'Coverage' :  unique_cov_col,
        'Position' : position_col })
    return df

def make_evolexp_coverage_df():

    stle_dir = "../breseq-assemblies/evoexp-polymorphism/"

    lineage_dict = {
        'RM3-149-1':'Ara+1',
        'RM3-149-2':'Ara+2',
        'RM3-149-3':'Ara+3',
        'RM3-149-4':'Ara+4',
        'RM3-149-5':'Ara+5',
        'RM3-149-6':'Ara+6',
        'RM3-149-7':'Ara-1',
        'RM3-149-8':'Ara-2',
        'RM3-149-9':'Ara-3',
        'RM3-149-10':'Ara-4',
        'RM3-149-11':'Ara-5',
        'RM3-149-12':'Ara-6',
        'RM3-153-1':'Ara+1',
        'RM3-153-2':'Ara+2',
         'RM3-153-3':'Ara+3',
        'RM3-153-4':'Ara+4',
        'RM3-153-5':'Ara+5',
        'RM3-153-6':'Ara+6',
        'RM3-153-7':'Ara-1',
        'RM3-153-8':'Ara-2',
        'RM3-153-9':'Ara-3',
        'RM3-153-10':'Ara-4',
        'RM3-153-11':'Ara-5',
        'RM3-153-12':'Ara-6'
    }

    generation_dict = {
        'RM3-149-1':1000,
        'RM3-149-2':1000,
        'RM3-149-3':1000,
        'RM3-149-4':1000,
        'RM3-149-5':1000,
        'RM3-149-6':1000,
        'RM3-149-7':1000,
        'RM3-149-8':1000,
        'RM3-149-9':1000,
        'RM3-149-10':1000,
        'RM3-149-11':1000,
        'RM3-149-12':1000,
        'RM3-153-1':1200,
        'RM3-153-2':1200,
        'RM3-153-3':1200,
        'RM3-153-4':1200,
        'RM3-153-5':1200,
        'RM3-153-6':1200,
        'RM3-153-7':1200,
        'RM3-153-8':1200,
        'RM3-153-9':1200,
        'RM3-153-10':1200,
        'RM3-153-11':1200,
        'RM3-153-12':1200,
    }

    lineage_col = []
    generation_col = []
    unique_cov_col = []
    position_col = []

    for l in [x for x in listdir(stle_dir) if not x.startswith(".")]:
        my_strain = l.split('_')[0] ## cheap string hack
        my_lineage = lineage_dict[my_strain]
        my_generation = generation_dict[my_strain]
        coverage_f = join(stle_dir,l,"08_mutation_identification/NC_002483.coverage.tab")
        coverage_fh = open(coverage_f)
        for i, line in enumerate(coverage_fh):
            if i == 0:
                continue
            line_data = line.split()
            unique_coverage = int(line_data[0]) + int(line_data[1])
            position = int(line_data[-1])
            ## Now add row data into each column list.
            lineage_col.append(my_lineage)
            generation_col.append(my_generation)
            unique_cov_col.append(unique_coverage)
            position_col.append(position)

    evolexp_df = pd.DataFrame( { 'Generation' : generation_col,
                         'Lineage' : lineage_col,
                         'Coverage' :  unique_cov_col,
                         'Position' : position_col })
    return evolexp_df

def main():

    ## assert that current working directory is STLE-analysis/src.
    assert getcwd().endswith('STLE-analysis/src')
    projdir = ".."
    outdir = join(projdir,"results/F-plasmid-search")

    df = make_STLE_clone_coverage_df()
    df2 = make_evolexp_coverage_df()

    ## write clone coverage to csv.
    df.to_csv(join(projdir,"results/STLE-clone-F-coverage.csv"))
    ## write evoexp coverage to csv.
    df2.to_csv(join(projdir,"results/STLE-evoexp-F-coverage.csv"))

main()

