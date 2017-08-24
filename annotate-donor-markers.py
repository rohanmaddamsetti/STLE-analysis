#!/usr/bin/env python

'''
annotate-donor-markers.py by Rohan Maddamsetti.

If running on orchestra, make sure to module load seq/breseq/0.31.0 and source activate evcouplings_env.

Takes the donor gd output of annotate-gds.py, which are in 
../../STLE-analysis/annotated-diffs/donor-specific-markers/Donor.

Calculates markers specific to each donor strain using the breseq utility gdtools.

Calculation of donor-specific markers using gdtools UNION and SUBTRACT.

REL288: A
REL291: B
REL296: C
REL298: D

annotated_donor_union.gd is union(A,B,C,D).
annotated_donorABC.gd is union(A,B,C).
annotated_donorABD.gd is union(A,B,D).
annotated_donorACD.gd is union(A,C,D).
annotated_donorBCD.gd is union(B,C,D).

annotated_REL288_specific_markers.gd is A' = A - union(B,C,D)
annotated_REL291_specific_markers.gd is B' = B - union(A,C,D)
annotated_REL296_specific_markers.gd is C' = C - union(A,B,D)
annotated_REL298_specific_markers.gd is D' = D - union(A,B,C)

annotated_donor_intersection.gd = union(A,B,C,D) - A' - B' - C' - D'.


'''

import os
import subprocess
import itertools

def run_gdtools(gdlist,outfile,call='UNION'):
    if call not in ['UNION','SUBTRACT']:
        print('ERROR: CAN\'T RUN gdtools',call)
        quit()
    if os.path.exists(outfile):
        print("skipping: ",outfile,"exists")
        return
    orch_args = ['bsub', '-K', '-q', 'short', '-W', '6:0','gdtools',call,'-o',outfile] + list(gdlist)
    subprocess.run(orch_args)
    print(' '.join(orch_args))
    
def main():
    diffdir = '../annotated-diffs/donor-specific-markers'
    indir = os.path.join(diffdir,'Donor')

    donors = ['REL288','REL291','REL296','REL298']
    ## sort all_donor_gds to get output files to match up right.
    all_donor_gds = sorted([os.path.join(indir,f) for f in os.listdir(indir)])
    
    ## check that the input looks right in indir.
    for f in [os.path.basename(x) for x in all_donor_gds]:
        prefix,rest = f.split('_')
        my_donor,suffix = rest.split('.')
        assert my_donor in donors
        donors.remove(my_donor)
    assert len(donors) == 0

    donorsABCD = os.path.join(diffdir,'annotated_donor_union.gd')
    run_gdtools(all_donor_gds,donorsABCD,call='UNION')

    donor_outfiles = [os.path.join(diffdir,'annotated_donorABC.gd'),
                      os.path.join(diffdir,'annotated_donorABD.gd'),
                      os.path.join(diffdir,'annotated_donorACD.gd'),
                      os.path.join(diffdir,'annotated_donorBCD.gd')]
    
    donor_infiles = itertools.combinations(all_donor_gds,3)
    for my_infiles, my_outfile in zip(donor_infiles,donor_outfiles):
        run_gdtools(my_infiles, my_outfile,call='UNION')

    '''
    Now make the donor specific marker files.

    annotated_REL288_specific_markers.gd is A' = A - union(B,C,D)
    annotated_REL291_specific_markers.gd is B' = B - union(A,C,D)
    annotated_REL296_specific_markers.gd is C' = C - union(A,B,D)
    annotated_REL298_specific_markers.gd is D' = D - union(A,B,C)
    '''
    donor_specific_gdfiles = []
    for donor_gd, subtract_me in zip(all_donor_gds,reversed(donor_outfiles)):
        donor_specific_out = os.path.basename(donor_gd).split('.')[0] + '_specific_markers.gd'
        full_donor_specific_out = os.path.join(diffdir,donor_specific_out)
        run_gdtools([donor_gd,subtract_me],full_donor_specific_out,call='SUBTRACT')
        donor_specific_gdfiles.append(full_donor_specific_out)
    
    '''annotated_donor_intersection.gd = union(A,B,C,D) - A' - B' - C' - D'.'''
    run_gdtools([donorsABCD]+donor_specific_gdfiles,os.path.join(diffdir,'annotated_donor_intersection.gd'),call='SUBTRACT')
        
main()
