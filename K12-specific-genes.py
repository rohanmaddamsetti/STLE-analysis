#!/usr/bin/env python

'''
K12-specific-genes.py by Rohan Maddamsetti.

Examine introgression of K-12 genes that
are not found in the REL606 reference genome.

1) map unaligned reads (those that don't map to REL606) to K-12
   using bowtie2.
2) use samtools depth or samtools mpileup to print coverage per
   position.
3) use K-12 reference and biopython to print out a text file with
   the following header:
   evolved clone, pop, region_start, region_end, gene_list,
   coverage. (OR do this in R as in my copy number analysis).
'''

from os import listdir, getcwd, makedirs
from os.path import join, exists, dirname
import subprocess
from lxml import html
from scipy.stats import nbinom
from math import isclose

def run_bowtie2_build(ref_genome,outdir,indexname='K12_index'):
    bowtie2_buildargs = ['bowtie2-build', ref_genome, indexname]
    return subprocess.run(bowtie2_buildargs,cwd=outdir)

def get_unaligned_readfiles(breseq_outdir):
    data_dir = join(breseq_outdir,"data")
    flist = [join(data_dir,x) for x in listdir(data_dir) if x.endswith('.unmatched.fastq')]
    return flist

def make_sorted_bam_aln(resultdir,genome_name,readfile_list,indexname='K12_index'):
    '''
    Map unaligned reads (those that don't map to REL606) to K-12
    using bowtie2. Then use samtools view to compress to bam, and samtools sort
    to sort the bam file.

    Bowtie2 settings:
    ignore paired-end information in reads, and do end-to-end alignment of reads.

    '''
    samfile = genome_name + '.sam'
    bowtie2_args = ['bowtie2', '--time', '-x', indexname, '-U', ','.join(readfile_list),
                    '-S', samfile]
    subprocess.run(bowtie2_args,cwd=resultdir)

    bamfile = genome_name + '.bam'
    bamfile_h = open(join(resultdir,bamfile),'w')
    samtools_view_args = ['samtools', 'view', '-bS', samfile]
    subprocess.run(samtools_view_args,cwd=resultdir,stdout=bamfile_h)
    bamfile_h.close()

    sorted_bamfile = genome_name + '.sorted.bam'
    samtools_sort_args = ['samtools', 'sort', bamfile, '-o', sorted_bamfile]
    subprocess.run(samtools_sort_args,cwd=resultdir)

    ## now clean up by removing .sam and .bam files.
    rm_sam_args = ['rm', samfile]
    subprocess.run(rm_sam_args,cwd=resultdir)
    rm_bam_args = ['rm', bamfile]
    subprocess.run(rm_bam_args,cwd=resultdir)

def REL606_coverage_from_summary(assemblydir,resultdir):
    '''
    parse the summary.html breseq output file, and return the mean and dispersion
    of the negative binomial fit to the read coverage distribution, returned as a
    dict with keys {mean, dispersion, variance}.
    NOTE: this code has only been tested on the summary file
    output by breseq 0.30.0. It might fail on earlier or later versions.
    '''
    coverage_dict = {}
    for my_genome in [x for x in listdir(assemblydir) if not x.startswith('.')]:
        ## skip K-12 donor genomes.
        if my_genome in ['REL288', 'REL291', 'REL296', 'REL298']:
            continue
        summaryf = join(assemblydir,my_genome,"output/summary.html")
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

def calc_coverage_threshold(cov_dict):
    '''
    calculate minimum coverage threshold for each key in cov_dict.
    see end of 'alternative parameterization' section of Negative binomial page
    and scipy negative binomial documentation for details of calculation.
    '''
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
        my_threshold = nbinom.ppf(0.05,n,p)
        my_threshold_p = nbinom.cdf(my_threshold,n,p)
        threshold_dict[g] = {'threshold':str(my_threshold),
                             'threshold_p':str(my_threshold_p)}
    return threshold_dict

def write_coverage(outfh, my_genome, my_sorted_bamfile, threshold):
    '''
    see snippet that Benni emailed me.
    '''
    samtools_depth_args = ['samtools', 'depth', my_sorted_bamfile]
    depth_out = subprocess.Popen(samtools_depth_args,stdout=subprocess.PIPE,universal_newlines=True)
    while True:
        line = depth_out.stdout.readline()
        if line == '' and depth_out.poll() is not None:
            break
        if line:
            line = line.strip()
            ## fields for samtools depth: reference header, position, coverage
            fields = line.split('\t')
            position = fields[1]
            coverage = fields[2]
            if coverage > threshold:
                csv_line = ','.join([my_genome, position, coverage])
                print(csv_line,file=outfh)
    rc = depth_out.poll()
    return rc

def generate_K12_coverage_file(assemblydir,resultdir,threshold_dict):
    ## open a fh to write coverage.
    coverage_outfile = join(resultdir,"K12-coverage.csv")
    coverage_out_fh = open(coverage_outfile, 'w')
    coverage_header = "Name, K12_position, coverage"
    print(coverage_header, file=coverage_out_fh)

    for my_genome in [x for x in listdir(assemblydir) if not x.startswith('.')]:
        ## skip K-12 donor genomes.
        if my_genome in ['REL288', 'REL291', 'REL296', 'REL298']:
            continue
        fullf = join(assemblydir,my_genome)
        my_unaligned_readfiles = get_unaligned_readfiles(fullf)
        ## Run bowtie2 and samtools view and sort if sorted.bam files don't exist.
        my_sorted_bamfile = join(resultdir,my_genome+'.sorted.bam')
        if not exists(my_sorted_bamfile):
            ## Ignore paired-end information in reads,
            ## and do end-to-end alignments of reads.
            make_sorted_bam_aln(resultdir,my_genome, my_unaligned_readfiles)

        ## 2) use samtools depth to get coverage per
        ## position, and write out to the open csv fh to analyze with R.
        print(my_genome)
        threshold = threshold_dict[my_genome]['threshold']
        write_coverage(coverage_out_fh,my_genome, my_sorted_bamfile, threshold)

def main():

    ## run source activate evcouplings_env if needed
    ## and module load bowtie2 and samtools if on orchestra.

    ## assume script is called from STLE-analysis/src
    assert getcwd().endswith('/src')
    projdir = dirname(getcwd())

    resultdir = join(projdir, 'results/K12-specific-genes')
    if not exists(resultdir):
        makedirs(resultdir)

    assemblydir = join(projdir, 'breseq-assemblies/REL606-polymorphism')

    K12fasta = join(projdir,'references/K-12.fasta')

    ## get REL606 ref coverage, calculate a threshold, and write to file.
    cov_dict = REL606_coverage_from_summary(assemblydir,resultdir)
    threshold_dict = calc_coverage_threshold(cov_dict)
    assert cov_dict.keys() == threshold_dict.keys()

    outfile = join(resultdir,"REL606-ref-coverage-summary.csv")
    out_fh = open(outfile, 'w')
    header = "Name, mean, dispersion, variance, threshold.coverage, threshold.probability"
    print(header,file=out_fh)
    for g in cov_dict:
        print(','.join([g, cov_dict[g]['mean'],
                               cov_dict[g]['dispersion'],
                               cov_dict[g]['variance'],
                               threshold_dict[g]['threshold'],
                               threshold_dict[g]['threshold_p']]),
              file=out_fh)


    ## make a bowtie2 index for K-12 genome if not found.
    if not len([x for x in listdir(resultdir) if x.startswith('K12_index')]):
        run_bowtie2_build(K12fasta,resultdir)

    ## filter coverage based on the threshold.
    generate_K12_coverage_file(assemblydir,resultdir,threshold_dict)

main()
