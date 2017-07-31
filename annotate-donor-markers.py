#!/usr/bin/env python

'''
annotate-donor-markers.py by Rohan Maddamsetti.

Takes the donor gd output of annotate-gds.py, which are in 
../../STLE-analysis/annotated-diffs/donor-specific-markers/Donor.

Calculates markers specific to each donor strain using the breseq utility gdtools.

Calculation of donor-specific markers using gdtools UNION and SUBTRACT.

REL288: A
REL291: B
REL296: C
REL298: D

annotated_donor_union.gd is union(A,B,C,D).
annotated_donorBCD.gd is union(B,C,D).
annotated_donorACD.gd is union(A,C,D).
annotated_donorABD.gd is union(A,B,D).
annotated_donorABC.gd is union(A,B,C).

annotated_REL288_specific_markers.gd is A' = A - union(B,C,D)
annotated_REL291_specific_markers.gd is B' = B - union(A,C,D)
annotated_REL296_specific_markers.gd is C' = C - union(A,B,D)
annotated_REL298_specific_markers.gd is D' = D - union(A,B,C)

annotated_donor_intersection.gd = union(A,B,C,D) - A' - B' - C' - D'.


'''

import os
import subprocess

def gd_union(gdlist,outfile):
    if os.path.exists(outfile):
        print("skipping: ",outfile,"exists")
        return
    union_args = ['gdtools','UNION'] + gdlist + ['-o',outfile]
    subprocess.run(union_args)
    
def main():
    diffdir = '../annotated-diffs/donor-specific-markers'
    indir = os.path.join(diffdir,'Donor')

    donors = ['REL288','REL291','REL296','REL298']
    all_donor_gds = [os.path.join(indir,f) for f in os.listdir(indir)]
    
    ## check that the input looks right in indir.
    for f in [os.path.basename(x) for x in all_donor_gds]:
        prefix,rest = f.split('_')
        my_donor,suffix = rest.split('.')
        assert my_donor in donors
        donors.remove(my_donor)
    assert len(donors) == 0

    donorsABCD = os.path.join(diffdir,'annotated_donor_union.gd')
    gd_union(all_donor_gds,donorsABCD)
    
    
            

main()
