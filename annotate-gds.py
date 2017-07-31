#!/usr/bin/env python

'''
annotate-gds.py by Rohan Maddamsetti.

Run this script before running label_mutations.py.

If running on orchestra, make sure to module load seq/breseq/0.30.2 and source activate evcouplings_env.

This script annotates output.gd files and puts in a nice directory structure for
label_mutations.py to work on.

Basically a wrapper around gdtools ANNOTATE that makes this particular directory structure:

REL606-ref-runs (outdir)
  -Ara-1
  -Ara-2
  ...
  -Ara+5
  -Ara+6
   -annotated-poly_REL2542.gd
   -annotated-poly_RM3-130-11.gd
   -annotated-poly_RM3-130-12.gd

ALSO: handle case for turner clones (REL4397 and REL4398),
and output donors in donor-specific-markers.

A script called annotate-donors-specific-markers.py will do gdtools INTERSECT etc.
upon the gd output of this script.

'''

import os
import subprocess


def annotate_gds(indir,outdir,name,population,is_donor=False):
    infile = os.path.join(indir,name,'output', 'output.gd')
    outpath = os.path.join(outdir,population)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = os.path.join(outpath,'annotated_'+name+'.gd')
    if not os.path.exists(infile):
        print("ERROR: output.gd file not found.")
        quit()
    if os.path.exists(outfile):
        print("skipping: ",outfile,"exists")
        return
    if is_donor: ## give Tn10 reference as well.
        orch_annotate_args = ['bsub', '-q', 'short','-W','6:0', 'gdtools', 'ANNOTATE', '-f', 'GD', '-r', '../references/REL606.7.gbk', '-r', '../references/Tn10.gbk', '-o', outfile, infile]
    else:
        orch_annotate_args = ['bsub', '-q', 'short','-W','6:0', 'gdtools', 'ANNOTATE', '-f', 'GD', '-r', '../references/REL606.7.gbk', '-o', outfile, infile]
    subprocess.run(orch_annotate_args)
    print(' '.join(orch_annotate_args))

def main():

    ## make sure runnning from src dir.
    assert os.getcwd().endswith('/src')

    indir  = '../breseq-assemblies/REL606-ref-runs'
    outdir = '../annotated-diffs/REL606-ref-runs'
    evoexp_indir = '../breseq-assemblies/evoexp-polymorphism'
    evoexp_outdir = '../annotated-diffs/evoexp'
    donor_outdir = '../annotated-diffs/donor-specific-markers'

    ## get strain info.
    strain_csv = '../doc/Populations-and-Clones.csv'
    strain_fh = open(strain_csv)

    for i,l in enumerate(strain_fh):
        l = l.strip()
        if i == 0:
            ##print(l)
            continue
        ldata = l.split(',')
        name, population, is_clone, generation, rel_name = ldata
        ##print(ldata)
        if is_clone == '0': ## mixed pop. samples for evolution experiment.
            annotate_gds(evoexp_indir,evoexp_outdir,name,population)
        else:
            if population == 'Donor': ## write into donor-specific-markers.
                my_indir = os.path.join(indir,'donors')
                annotate_gds(my_indir,donor_outdir,name,population,is_donor=True)
            else:
                ## handle corner case for REL4397 and REL4398 clones (write to 'turner_clones').
                if name in ('REL4397','REL4398'):
                    my_indir = os.path.join(indir,'recombinants')
                    annotate_gds(my_indir,outdir,name,"turner_clones")
                else:
                    if name.startswith('REL'):
                        my_indir = os.path.join(indir,'recipients')
                        annotate_gds(my_indir,outdir,name,population)
                    elif name.startswith('RM'):
                        my_indir = os.path.join(indir,'recombinants')
                        annotate_gds(my_indir,outdir,name,population)

main()
