#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/breseq/0.31.0.

import subprocess
import os

## run breseq 0.30.2 on all samples.

def run_breseq(myout,reads1,reads2,polymorphism=False):
    ref_seq = '../references/REL606.7.gbk'
    Tn10_ref = '../references/Tn10.gbk'
    F_plasmid_seq = '../references/F-plasmid.1.gbk'

    prefix_args = ['bsub', '-q', 'medium', '-W', '72:0', '-R','rusage[mem=30000]' ,'breseq']
    suffix_args = ['-r', ref_seq, '-r', Tn10_ref, '-r', F_plasmid_seq, '-o', myout, reads1, reads2]

    if polymorphism:
        breseq_args = prefix_args + ['-p'] + suffix_args
    else:
        breseq_args = prefix_args + suffix_args    
    output_done = os.path.join(myout,'output/output.done')
    if os.path.exists(output_done):
        print("output exists: skipping",myout)
    else:
        print(' '.join(breseq_args))
        subprocess.run(breseq_args)

data_fh = open("../doc/Populations-and-Clones.csv")
for i,l in enumerate(data_fh):
    l = l.strip()
    #print(l)
    ## skip header
    if i == 0:
        continue
    ldata = l.split(',')
    name, population, is_clone, generation, rel_name = ldata

    ## Now set up each breseq run.
    if population == 'Donor':
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','donors',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        read_dir = "../read-data/combined_reads/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_out,reads1,reads2)
    elif is_clone == '1' and generation == '0': ## recipient clone.
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','recipients',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        read_dir = "../read-data/combined_reads/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_out,reads1,reads2)
    elif is_clone == '1' and generation == '1000': ## recombinant clone.
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','recombinants',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        if name in ['REL4397','REL4398']: ## corner case: handle turner clones.
            read_dir = "../read-data/combined_reads/"
        else:
            read_dir = "../read-data/20160120_DNASeq_PE/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_out,reads1,reads2)
    elif is_clone == '0': ## STLE continuation mixed pop.
        ref_out = os.path.join('../breseq-assemblies/','evoexp-polymorphism',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        read_dir = "../read-data/20160902_DNASeq_PE/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_out,reads1,reads2,polymorphism=True)
