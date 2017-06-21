#!/usr/bin/env python

'''
F-plasmid-search.py by Rohan Maddamsetti.

This script makes 3 figures, when examining breseq assemblies against
the F-plasmid reference genome.
1) makes violin plots of the coverage distribution for clones from each population.
2) plot ribbon plots of the coverage over the F-plasmid reference genome.
3) compares F-plasmid coverage to chromosome coverage in the donors versus
   in the 4 Ara-3 endpoint clones. By eye I see a 2-fold excess of the plasmid.

This script repeats 1) and 2) on the endpoint samples of my STLE continuation
evolution experiment, in order to show that the F-plasmid is apparently present
in the Ara+1 population as well as the Ara-3 population.

'''

from os import makedirs, listdir, chdir, getcwd
from os.path import exists, join, basename, isdir
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ks_2samp, mannwhitneyu
##import statsmodels.api as sm
##import statsmodels.formula.api as smf

def make_STLE_clone_coverage_df():
    '''
    This function works on the STLE clones. A separate function works on the
    evolution experiment data, due to the different directory stucture.
    '''

    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"
    stle_dir = join(projdir,"breseq-assemblies/F-plasmid-ref-runs/")
    lineage_dict = {'REL288':'REL288', 'REL291':'REL291','REL296':'REL296','REL298':'REL298',
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

    df = pd.DataFrame( { 'Strain Type' : strain_type_col,
                         'Lineage' : lineage_col,
                         'Unique Coverage' :  unique_cov_col,
                         'Position' : position_col })

    return df

def make_evolexp_coverage_df():
    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"
    coverage_dir = join(projdir,"breseq-assemblies/evolution-experiment/F-plasmid-coverage/")
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
        'RM3-149-1':0,
        'RM3-149-2':0,
        'RM3-149-3':0,
        'RM3-149-4':0,
        'RM3-149-5':0,
        'RM3-149-6':0,
        'RM3-149-7':0,
        'RM3-149-8':0,
        'RM3-149-9':0,
        'RM3-149-10':0,
        'RM3-149-11':0,
        'RM3-149-12':0,
        'RM3-153-1':1000,
        'RM3-153-2':1000,
        'RM3-153-3':1000,
        'RM3-153-4':1000,
        'RM3-153-5':1000,
        'RM3-153-6':1000,
        'RM3-153-7':1000,
        'RM3-153-8':1000,
        'RM3-153-9':1000,
        'RM3-153-10':1000,
        'RM3-153-11':1000,
        'RM3-153-12':1000,
    }

    lineage_col = []
    generation_col = []
    unique_cov_col = []
    position_col = []

    for l in [x for x in listdir(coverage_dir) if x.endswith(".tab")]:
        my_strain = l.split('_')[0] ## cheap string hack
        my_lineage = lineage_dict[my_strain]
        my_generation = generation_dict[my_strain]
        coverage_f = join(coverage_dir,l)
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
                         'Unique Coverage' :  unique_cov_col,
                         'Position' : position_col })
    return evolexp_df

def make_violin_plot(df,evolexp=False):
    '''
    make violin plot with seaborn.
    '''
    sns.set(style="ticks")
    plt.figure()

    if evolexp:
        ax = sns.violinplot(x="Lineage", y="Unique Coverage",hue='Generation',split=True, data=df)
        plt.savefig("/Users/Rohandinho/Desktop/evolexp_test.pdf", format='PDF')
    else:
        ax = sns.violinplot(x="Lineage", y="Unique Coverage", data=df)
        plt.savefig("/Users/Rohandinho/Desktop/test.pdf", format='PDF')

def make_coverage_small_multiple(df,evolexp=False):
    '''
    make small multiples of coverage over F-plasmid reference.
    '''
    sns.set(style="ticks")
    plt.figure()

    if evolexp:
        g = sns.FacetGrid(df, row="Generation", col="Lineage")
        g.map_dataframe(plt.plot, "Position", "Unique Coverage")
        ## TODO: fix the x-axis. try to use scientific notation on x-axis.
        #g.set(xaxis.get_major_formatter().set_powerlimits((0, 1))
        g.set(xticklabels=[])
        plt.savefig("/Users/Rohandinho/Desktop/test1.pdf", format='PDF')

    else:
        ## sum coverage by position (summing over clones in the same lineage)
        df2 = df.groupby(["Lineage","Position"], as_index=False).sum()
        g = sns.FacetGrid(df2, col="Lineage", col_wrap=4, size=2)
        g.map_dataframe(plt.plot, "Position", "Unique Coverage")
        ## TODO: fix the x-axis. try to use scientific notation on x-axis.
        #g.set(xaxis.get_major_formatter().set_powerlimits((0, 1))
        g.set(xticklabels=[])
        plt.savefig("/Users/Rohandinho/Desktop/test2.pdf", format='PDF')

def evol_exp_coverage_statistical_tests():
    '''
    Two functions due to different directory structures.

    Make a pandas dataframe containing coverage information
    among the Ara-3 STLE continuation populations,
    and for the Ara+1 Day 30 STLE continuation population.

    To avoid auto-correlation in the coverage sampling, use coverage every 500 bp
    (since short reads at most gives 350 bp of information).
'''

    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"

    lineage_col = []
    unique_cov_col = []
    position_col = []
    refseq_col = [] ## values are either "CHR" or "F"
    sample_col = []

    lineage_dict = {"RM3-149-9":'Ara-3',"RM3-153-1":'Ara+1',"RM3-153-9":'Ara-3'}
    samples = lineage_dict.keys()
    chr_coverage_dir = join(projdir,"breseq-assemblies/evolution-experiment/coverage")
    F_coverage_dir = join(projdir,"breseq-assemblies/evolution-experiment/F-plasmid-coverage")

    for xf in [y for y in listdir(chr_coverage_dir) if y.endswith(".tab")]:
        for sample in samples:
            ## make sure we have the right sample.
            if sample in xf:
                this_cov_fh = open(join(chr_coverage_dir,xf))
                for i,line in enumerate(this_cov_fh):
                    ## skip header, and sample once every 500 positions.
                    if i == 0 or i % 500 != 1:
                        continue
                    line_data = line.split()
                    unique_coverage = int(line_data[0]) + int(line_data[1])
                    position = int(line_data[-1])
                    lineage = lineage_dict[sample]
                    ## Now add row data into each column list.
                    sample_col.append(sample)
                    refseq_col.append("CHR")
                    lineage_col.append(lineage)
                    unique_cov_col.append(unique_coverage)
                    position_col.append(position)

    for xf in [y for y in listdir(F_coverage_dir) if y.endswith(".tab")]:
        for sample in samples:
            ## make sure we have the right sample.
            if sample in xf:
                this_cov_fh = open(join(F_coverage_dir,xf))
                for i,line in enumerate(this_cov_fh):
                    ## skip header, and sample once every 500 positions.
                    if i == 0 or i % 500 != 1:
                        continue
                    line_data = line.split()
                    unique_coverage = int(line_data[0]) + int(line_data[1])
                    position = int(line_data[-1])
                    lineage = lineage_dict[sample]
                    ## Now add row data into each column list.
                    sample_col.append(sample)
                    refseq_col.append("F")
                    lineage_col.append(lineage)
                    unique_cov_col.append(unique_coverage)
                    position_col.append(position)

    df = pd.DataFrame({
        'Sample' : sample_col,
        'Reference' : refseq_col,
        'Lineage' : lineage_col,
        'Unique Coverage' :  unique_cov_col,
        'Position' : position_col
    })

    ## Make a small multiple comparing coverage (ref in different colors)
    sns.set(style="ticks")
    plt.figure()
    g = sns.FacetGrid(df, row="Sample", col="Reference",sharex=False,sharey=False)
    g.map_dataframe(plt.plot, "Position", "Unique Coverage")
    axes = g.axes
    ## change x axis for the F-plasmid since it's so much smaller than E. coli
    ##chromosome, and make sure y-axes are on the same scale
    axes[0,1].set_xlim(0,100000)
    axes[1,1].set_xlim(0,100000)

    axes[1,1].set_ylim(0,1800)

    axes[1,1].set_xlim(0,100000)
    axes[2,1].set_xlim(0,100000)
    plt.savefig("/Users/Rohandinho/Desktop/evolexp-F-plasmid-coverage.pdf", format='PDF')

def coverage_statistical_tests():
    '''
    In the populations with evidence of the F-plasmid,
    does coverage indicate elevated copy-number, i.e. existence
    as a plasmid rather than being chromosomally integrated?

    Make a pandas dataframe containing coverage information
    among the Ara-3 clones.

    As a control, do these tests on the donor clones as well.

    Then, do either a Mann-Whitney U test or a 2-sample Kolmogorov-Smirnov test
    to see whether coverage on the F-plasmid reference is higher than the coverage
    on the reference genome.

    On the Ara-3 clones, use the K-12 reference genome.

    To avoid auto-correlation in the coverage sampling, use coverage every 500 bp
    (since short reads at most gives 350 bp of information).
'''

    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"

    lineage_col = []
    unique_cov_col = []
    position_col = []
    refseq_col = [] ## values are either "CHR" or "F"
    sample_col = []

    lineage_dict = {"REL288":'Donor',"REL291":'Donor',"REL296":'Donor',
                        "REL298":'Donor',"RM3-130-17":"Ara-3",'RM3-130-18':'Ara-3',
                        'REL4397':'Ara-3','REL4398':'Ara-3'}
    samples = list(lineage_dict.keys())
    assembly_dir = join(projdir,"breseq-assemblies")
    chr_coverage_base = join(assembly_dir,"K-12-ref-runs")
    for yf in listdir(chr_coverage_base):
        if not isdir(join(chr_coverage_base,yf)):
            continue
        for xf in listdir(join(chr_coverage_base,yf)):
            for sample in samples:
                ## make sure we have the right sample.
                if sample in xf:
                    this_cov_file = join(chr_coverage_base,yf,xf,"08_mutation_identification/NC_000913.coverage.tab")
                    this_cov_fh = open(this_cov_file)
                    for i,line in enumerate(this_cov_fh):
                        ## skip header, and sample once every 500 positions.
                        if i == 0 or i % 500 != 1:
                            continue
                        line_data = line.split()
                        unique_coverage = int(line_data[0]) + int(line_data[1])
                        position = int(line_data[-1])
                        lineage = lineage_dict[sample]
                        ## Now add row data into each column list.
                        sample_col.append(sample)
                        refseq_col.append("CHR")
                        lineage_col.append(lineage)
                        unique_cov_col.append(unique_coverage)
                        position_col.append(position)

    F_coverage_base = join(assembly_dir,"F-plasmid-ref-runs")
    for yf in listdir(F_coverage_base):
        if not isdir(join(F_coverage_base,yf)):
            continue
        for xf in listdir(join(F_coverage_base,yf)):
            for sample in samples:
                ## make sure we have the right sample.
                if sample in xf:
                    this_cov_file = join(F_coverage_base,yf,xf,"08_mutation_identification/NC_002483.coverage.tab")
                    this_cov_fh = open(this_cov_file)
                    for i,line in enumerate(this_cov_fh):
                        ## skip header, and sample once every 500 positions.
                        if i == 0 or i % 500 != 1:
                            continue
                        line_data = line.split()
                        unique_coverage = int(line_data[0]) + int(line_data[1])
                        position = int(line_data[-1])
                        lineage = lineage_dict[sample]
                        ## Now add row data into each column list.
                        sample_col.append(sample)
                        refseq_col.append("F")
                        lineage_col.append(lineage)
                        unique_cov_col.append(unique_coverage)
                        position_col.append(position)

    df = pd.DataFrame({
        'Sample' : sample_col,
        'Reference' : refseq_col,
        'Lineage' : lineage_col,
        'Unique Coverage' :  unique_cov_col,
        'Position' : position_col
    })

    ## Now, take subsets of the data frame by sample_col, and do KS-tests by refseq_col.

    for sample in samples:
        testdf = df.loc[df['Sample'] == sample]
        chr_cov_vec = testdf[testdf['Reference']=="CHR"]['Unique Coverage']
        F_cov_vec = testdf[testdf['Reference']=="F"]['Unique Coverage']
        this_KS_test = ks_2samp(chr_cov_vec,F_cov_vec)
        print(sample)
        print(this_KS_test)
        this_U_test = mannwhitneyu(chr_cov_vec,F_cov_vec)
        print(this_U_test)

    ## Make a small multiple comparing coverage (ref in different colors)
    sns.set(style="ticks")
    plt.figure()
    g = sns.FacetGrid(df, row="Sample", col="Reference",sharex=False,sharey=True)
    g.map_dataframe(plt.plot, "Position", "Unique Coverage")
    axes = g.axes
    ## change x axis for the F-plasmid since it's so much smaller than E. coli
    ##chromosome, and make sure y-axes are on the same scale
    axes[0,1].set_xlim(0,100000)
    axes[1,1].set_xlim(0,100000)
    axes[1,1].set_xlim(0,100000)
    axes[2,1].set_xlim(0,100000)
    plt.savefig("/Users/Rohandinho/Desktop/coverage_tests.pdf", format='PDF')

    ## try to make a violin plot as well.
    plt.figure()
    v = sns.violinplot(x="Sample", y="Unique Coverage",hue="Reference", data=df)
    plt.savefig("/Users/Rohandinho/Desktop/violin.pdf", format='PDF')

    return df

def main():
    coverage_statistical_tests()
    ##evol_exp_coverage_statistical_tests()
    ##df = make_STLE_clone_coverage_df()
    ##make_violin_plot(df)
    ##make_coverage_small_multiple(df)

    ##df2 = make_evolexp_coverage_df()
    #make_violin_plot(df2,evolexp=True)
    ##make_coverage_small_multiple(df2,evolexp=True)



    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"
    stle_dir = join(projdir,"breseq-assemblies/F-plasmid-ref-runs/")
    outdir = join(projdir,"results/F-plasmid-search")

main()

