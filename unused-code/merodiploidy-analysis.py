#!/usr/bin/env python

'''
merodiploidy-analysis.py by Rohan Maddamsetti

 decision conditions for merodiploidy:
  1) intermediate frequency in breseq -p output for clones, 0.2 < x < 0.8
  2) better fits 2X coverage model than 1X coverage model (2X model: Poisson(2*mu))
  3) not found in any donor genomes (w/ REL606 reference in breseq -p output).

'''

import os
from lxml import html
from scipy.stats import nbinom, poisson
from math import isclose

def REL606_coverage_from_summary(assemblydir):
    '''
    parse the summary.html breseq output file, and return the mean and dispersion
    of the negative binomial fit to the read coverage distribution, returned as a
    dict with keys {mean, dispersion, variance}.
    NOTE: this code has only been tested on the summary file
    output by breseq 0.30.0. It might fail on earlier or later versions.
    '''
    coverage_dict = {}
    for my_genome in [x for x in os.listdir(assemblydir) if not x.startswith('.')]:
        summaryf = os.path.join(assemblydir,my_genome,"output/summary.html")
        tree = html.parse(summaryf)
        ## print text in the table 'Reference Sequence Information.
        query = '//table[./tr/th[contains(text(),"fit dispersion")]]'
        table = tree.xpath(query)[0]
        table_data = table.xpath('./tr/td')
        avg = table_data[4].text_content()
        disp = table_data[5].text_content()
        var = str(float(avg)*float(disp))
        coverage_dict[my_genome] = {'mean':avg,'dispersion':disp,'variance':var}

    return coverage_dict

def calc_2X_coverage_threshold(cov_dict):
    '''
    calculate coverage threshold for each key in cov_dict, based on a likelihood ratio 
    between empirical Nbinom(mu,disp) 1X coverage distribution, and a theoretical 
    Poisson(2*mu) 2X coverage distribution.
    see end of 'alternative parameterization' section of Negative binomial page
    and scipy negative binomial documentation for details of calculation.

    choose coverage threshold s.t. log likelihood ratio > 10.

    '''

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
                'RM3-130-23':'REL11756','RM3-130-24':'REL11757',
                'REL4397':'REL4397', 'REL4398':'REL4398',
                'REL288':'REL288','REL291':'REL291','REL296':'REL296','REL298':'REL298'}

    
    threshold_dict = {}
    for g in cov_dict:
        mean = float(cov_dict[g]['mean'])
        var = float(cov_dict[g]['variance'])
        q = (var-mean)/var
        n = mean**2/(var-mean)
        p = 1 - q
        
        ## assert that I did the math correctly.
        assert(isclose(nbinom.mean(n,p), mean))
        assert(isclose(nbinom.var(n,p), var))

        ## find the integer threshold that includes ~95% of REL606 distribution,
        ## excluding 5% on the left hand side.
        for x in range(int(mean),int(2*mean)):
            p0 = nbinom.pmf(x,n,p)
            p1 = poisson.pmf(x,2*mean)
            lratio = p1/p0
            if lratio > 10:
                my_threshold = x
                my_threshold_p0 = p0
                my_threshold_p1 = p1
                my_lratio = lratio
                break    
        threshold_dict[rel_name[g]] = {'threshold':str(my_threshold),
                             'threshold_p0':str(my_threshold_p0),
                             'threshold_p1':str(my_threshold_p1),
                             'lratio':str(lratio)}
    return threshold_dict

def make_highcov_dict(assemblydir, threshold_dict):
    ''' dictionary of genome => set of positions with high coverage.'''
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
                'RM3-130-23':'REL11756','RM3-130-24':'REL11757',
                'REL4397':'REL4397', 'REL4398':'REL4398',
                'REL288':'REL288','REL291':'REL291','REL296':'REL296','REL298':'REL298'}

    
    hicov_dict = {rel_name[x]:set() for x in os.listdir(assemblydir) if not x.startswith('.')}
    for my_genome in [x for x in os.listdir(assemblydir) if not x.startswith('.')]:
        breseq_outdir = os.path.join(assemblydir,my_genome)
        covfile = os.path.join(breseq_outdir,"08_mutation_identification","REL606.coverage.tab")
        covfh = open(covfile)
        for i,l in enumerate(covfh):
            if i == 0:
                continue
            dat = l.split()
            pos = dat[-1]
            pos_coverage = int(dat[0]) + int(dat[1])
            if int(threshold_dict[rel_name[my_genome]]['threshold']) < pos_coverage:
                hicov_dict[rel_name[my_genome]].add(pos)
    return hicov_dict

def get_annotation_from_field(d,field_key):
    '''helper function to get the value associated with a given field
    in list containing the fields in one line of an annotated genome diff file.'''
    for f in d:
        if f.startswith(field_key):
            return f.split('=')[1]

def find_donor_false_positives(highcov):
    ''' 
    get positions with intermediate frequencies in the donor annotated gds, that
    are also fit the 2X coverage model better than the 1X coverage model. 
    '''
    false_positive_pos = set()
    for my_donor in ['REL288','REL291','REL296', 'REL298']:
        my_gd = '../annotated-diffs/REL606-poly-runs/Donor/annotated-poly_'+my_donor+'.gd'
        for line in open(my_gd):
            if line.startswith('#'):
                continue
            fields = line.split()
            pos = fields[4]
            freq = get_annotation_from_field(fields,'frequency')
            if freq is None:
                continue
            freq = float(freq)
            #if pos in highcov[my_donor] and 0.2 < freq and freq < 0.8:
            if 0.2 < freq and freq < 0.8:
                false_positive_pos.add(pos)
    return false_positive_pos

def main():
    ## run source activate evcouplings_env if needed
    ## and module load bowtie2 and samtools if on orchestra.

    ## assume script is called from STLE-analysis/src
    assert os.getcwd().endswith('/src')
    projdir = os.path.dirname(os.getcwd())

    assemblydir = os.path.join(projdir, 'breseq-assemblies/REL606-polymorphism')

    ## get REL606 ref coverage statistics, and
    ## calculate threshold coverage for 10X likelihood ratio.
    cov_dict = REL606_coverage_from_summary(assemblydir)
    threshold_dict = calc_2X_coverage_threshold(cov_dict)
    #assert cov_dict.keys() == threshold_dict.keys()
    print(threshold_dict)

    ## get positions in each genome that are more likely to be 2X copy number than 1X.
    highcov = make_highcov_dict(assemblydir, threshold_dict)

    donor_false_positives = find_donor_false_positives(highcov)
    
    ## filter for intermediate frequency mutations, and filter out false positives
    ## in the donor genomes.
    poly_labeled_muts = os.path.join(projdir,"results/poly_labeled_mutations.csv")
    for i, line in enumerate(open(poly_labeled_muts)):
        line = line.strip()
        if i == 0:
            print(line)
            continue
        linedata = line.split(',')
        freq = linedata[-2]
        genome = linedata[1]
        pos = linedata[4]
        if freq == 'NA':
            continue
        elif 0.2 < float(freq) and float(freq) < 0.8 and pos in highcov[genome] and pos not in donor_false_positives:
            print(line)


main()
