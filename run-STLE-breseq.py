#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/breseq/0.30.2

import subprocess
import os

## run breseq 0.30.2 on all samples.

def run_breseq(reference_seq,myout,reads1,reads2,polymorphism=False):
    if polymorphism:
        breseq_args = ['bsub', '-q', 'medium', '-W', '36:0', 'breseq', '-p', '-r', reference_seq, '-o', myout, reads1, reads2]
    else:
        breseq_args = ['bsub', '-q', 'short', '-W', '12:0', 'breseq', '-r', reference_seq, '-o', myout, reads1, reads2]
    output_done = os.path.join(myout,'output/output.done')
    if not os.path.exists(output_done):
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

    ref_seq = '../references/REL606.7.gbk'
    F_plasmid_seq = '../references/F-plasmid.1.gbk'

    ## Now set up each breseq run.
    if population == 'Donor':
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','donors',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        F_plasmid_out = os.path.join('../breseq-assemblies/','F-plasmid-ref-runs','donors',name)
        if not os.path.exists(F_plasmid_out):
            os.makedirs(F_plasmid_out)
        read_dir = "../read-data/combined_reads/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_seq,ref_out,reads1,reads2)
        run_breseq(F_plasmid_seq,F_plasmid_out,reads1,reads2)
    elif is_clone == '1' and generation == '0': ## recipient clone.
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','recipients',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        read_dir = "../read-data/combined_reads/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_seq,ref_out,reads1,reads2)

    elif is_clone == '1' and generation == '1000': ## recombinant clone.
        ref_out = os.path.join('../breseq-assemblies/','REL606-ref-runs','recombinants',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        F_plasmid_out = os.path.join('../breseq-assemblies/','F-plasmid-ref-runs','recombinants',name)
        if not os.path.exists(F_plasmid_out):
            os.makedirs(F_plasmid_out)
        ## handle corner case for REL4397 and REL4398 clones (data in combined_reads).
        if name in ('REL4397','REL4398'):
            read_dir = "../read-data/combined_reads/"
        else:
            read_dir = "../read-data/20160120_DNASeq_PE/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_seq,ref_out,reads1,reads2)
        run_breseq(F_plasmid_seq,F_plasmid_out,reads1,reads2)
    elif is_clone == '0': ## STLE continuation mixed pop.
        ref_out = os.path.join('../breseq-assemblies/','evoexp-polymorphism',name)
        if not os.path.exists(ref_out):
            os.makedirs(ref_out)
        F_plasmid_out = os.path.join('../breseq-assemblies/','F-plasmid-ref-runs','evoexp',name)
        if not os.path.exists(F_plasmid_out):
            os.makedirs(F_plasmid_out)
        read_dir = "../read-data/20160902_DNASeq_PE/"
        reads1, reads2 = [os.path.join(read_dir,x) for x in os.listdir(read_dir) if name+'_' in x]
        run_breseq(ref_seq,ref_out,reads1,reads2,polymorphism=True)
        run_breseq(F_plasmid_seq,F_plasmid_out,reads1,reads2)
