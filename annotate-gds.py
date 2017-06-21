#!/usr/bin/env python

'''
annotate-gds.py by Rohan Maddamsetti.

## run this script before running label_mutations.py.

If running on orchestra, make sure to module load seq/breseq/X and source activate evcouplings_env.

This script annotates output.gd files and puts in a nice directory structure for
label_mutations.py to work on.

Basically a wrapper around gdtools ANNOTATE that makes this particular directory structure:

REL606-poly-runs (outdir)
  -donors
  -Ara-1
  -Ara-2
  ...
  -Ara+5
  -Ara+6
   -annotated-poly_REL2542.gd
   -annotated-poly_RM3-130-11.gd
   -annotated-poly_RM3-130-12.gd

'''

import os
import subprocess

def main():

    ## make sure runnning from src dir.
    assert os.getcwd().endswith('/src')

    ## hardcode input arguments for now.
    outdir = '../annotated-diffs/REL606-poly-runs'
    assemblydir  = '../breseq-assemblies/REL606-polymorphism'

    ## get strain info.
    strain_csv = '../doc/Populations-and-Clones.csv'
    strain_fh = open(strain_csv)

    for i,l in enumerate(strain_fh):
        l = l.strip()
        if i == 0:
            print(l)
            continue
        ldata = l.split(',')
        name, population, is_clone, generation, rel_name = ldata
        if is_clone == '0': ## skip mixed pop. samples.
            continue
        infile = os.path.join(assemblydir,name,'output', 'output.gd')
        outpath = os.path.join(outdir,population)
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outfile = os.path.join(outpath,'annotated-poly_'+name+'.gd')
        ##print(infile)
        ##print(outfile)
        ##annotate_args = ['gdtools', 'ANNOTATE', '-f', 'GD', '-r', '../references/REL606.7.gbk', '-o', outfile, infile]
        ##subprocess.run(annotate_args)
        orch_annotate_args = ['bsub', '-q', 'short','-W','6:0', 'gdtools', 'ANNOTATE', '-f', 'GD', '-r', '../references/REL606.7.gbk', '-o', outfile, infile]
        subprocess.run(orch_annotate_args)

main()
