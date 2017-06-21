#!/usr/bin/python

'''

BLAST-F-plasmid-form.py by Rohan Maddamsetti.

PROBLEM: REL296 doesn't have the 3' F-plasmid end in the denovo assembly!!!
LOOK INTO THIS.

Determine whether F-plasmid is integrated in chromosome or is in a
plasmid form by BLASTing oriT F-plasmid against de novo assemblies.

ALSO: Write out the true REL606 coordinate and orientation
for the location of the F-plasmid in the donors.

Procedure:
- make BLAST dbs for each de novo assembly, the K-12 chromosome, and F-plasmid.
- BLAST de novo assemblies with F-plasmid as the query.
- use query length to get the first BLAST hit that matches 3' end of query.
- get the flanking region in the de novo scaffold, call this 'Flank'.
- Blast Flank against K-12 and F-plasmid dbs.
- Assign F-plasmid as integrated in chromosome or not based on score.
- (possible TODO: replace assignment by a formal hypothesis test, 
   where H1: match to plasmid
   and H2: match to K-12 chromosome).

RESULT: all flanking regions map to repeat sequences that match both the 5'
        end of the F-plasmid, as well as several homologous regions of the K-12
        chromosome. So, these data could equally well come from free or integrated
        plasmid.

'''

from subprocess import run
from os import listdir, environ, makedirs
from os.path import join, exists, dirname, split, basename
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord
from pprint import pprint

##put the path from bash into the path in this script.
mypath = environ['PATH']
sys.path.append(mypath)

def cp_scaffolds(dna_file,blastdir):
    ## copy the dna_file into the results directory for looking up flanking sequences.
    args = ['cp', dna_file, blastdir]
    run(args,cwd=blastdir)

def makeBLASTdb(dna_file,blastdir):
    ## cheap hack to get database name right.
    ## if case: dna_file is a de novo assembly scaffolds.fasta.
    if basename(dna_file).startswith('scaffold'):
        my_genome = split(dirname(dna_file))[1]
    else: ## dna_file is either F-plasmid.fasta or K-12.fasta.
        my_genome = basename(dna_file).split('.')[0]
    args = ['makeblastdb', '-in', dna_file, '-dbtype', 'nucl',
                        '-input_type', 'fasta', '-out', join('.',my_genome), '-title', my_genome]
    run(args,cwd=blastdir)

def runBLAST(results_dir,queryfile):
    my_genome = basename(results_dir)
    blast_cmd = ["blastn", "-task", "megablast", "-db", my_genome,
                 "-outfmt", "5", "-max_target_seqs", "1",
                 "-query", queryfile, "-out",
                 "./results.xml"]
    run(blast_cmd, cwd=results_dir)


def get_flank(alignment,hsp, subject_file):
    subject_name = alignment.title.split(' ')[1]
    subject_start = hsp.sbjct_start
    subject_end = hsp.sbjct_end
    orientation = 1
    ##print(subject_name)
    ##print(subject_start)
    ##print(subject_end)
    ##print(orientation)
    if subject_end < subject_start:
        orientation = -1

    with open(subject_file, "rU") as subj_handle:
        for record in SeqIO.parse(subj_handle, "fasta"):
            if record.id != subject_name:
                continue
            else:
                if orientation == 1:
                    return record.seq[subject_end:]
                elif orientation == -1:
                    return record.seq[:subject_end]
    
def write_flanks(rbase,flanksfile):
    '''
    Parse the results from BLASTing the F-plasmid against the de novo assemblies.
    get the query length, get the first BLAST hit that matches the 3'-end of the
    query, and write the flanking region to file. 
    
    '''
    flank_record_list = []
    ## iterate over BLASTs against de novo assemblies.
    denovo_dirs = [x for x in listdir(rbase) if x.startswith('REL') or x.startswith('RM')] 
    for mygenome in denovo_dirs:
        myfulldir = join(rbase, mygenome)
        ##print(myfulldir)
        result_f = join(myfulldir,"results.xml")
        result_h = open(result_f)
        blast_record = NCBIXML.read(result_h)
        query_length = int(blast_record.query_letters)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect > 0.0000000001:
                    ## skip bad hits.
                    continue
                if hsp.query_end != query_length:
                    ## skip hits that don't match 3' end of F-plasmid query.
                    continue
                subject_seq = join(myfulldir,"scaffolds.fasta")
                ##print(mygenome)
                my_flank_seq = get_flank(alignment, hsp, subject_seq)
                flank_record_list.append(SeqRecord(seq=my_flank_seq, id=mygenome+'_flank'))
    with open(flanksfile,'w') as flanks_outhandle:               
        SeqIO.write(flank_record_list,flanks_outhandle, format="fasta")

def main():
    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis"
    denovodir = join(projdir, "denovo-assemblies")
    Fplasmid_file = join(projdir, "references/F-plasmid.fasta")
    K12ref_file = join(projdir, "references/K-12.fasta")
    results_basedir = join(projdir,"results/BLAST-F-plasmid-form")

    ## set up BLAST databases and results directories for K-12 and F-plasmid refs.
    F_resultsdir = join(results_basedir, "F-plasmid")
    if not exists(F_resultsdir):
        makedirs(F_resultsdir)
    K12_resultsdir = join(results_basedir, "K-12")
    if not exists(K12_resultsdir):
        makedirs(K12_resultsdir)        
    ##makeBLASTdb(Fplasmid_file,F_resultsdir)
    ##makeBLASTdb(K12ref_file,K12_resultsdir)

    denovos = [f for f in listdir(denovodir) if f.startswith('REL') or f.startswith('RM')]
    for f in denovos:
        full_f = join(denovodir,f)
        this_scaffolds = join(full_f,'scaffolds.fasta')
        BLASTresultdir = join(results_basedir,f)
        if not exists(BLASTresultdir):
            makedirs(BLASTresultdir)
        ## copy scaffolds file into BLASTresultdir for looking up flanking region.
        cp_scaffolds(this_scaffolds, BLASTresultdir)
        ## set up BLAST databases and results directories for denovo assemblies.
        ##makeBLASTdb(this_scaffolds,BLASTresultdir)
        ## BLAST F-plasmid against de novo assemblies.
        ##runBLAST(BLASTresultdir,Fplasmid_file)

    flanksfile = join(results_basedir, "flanks.fasta")
    write_flanks(results_basedir, flanksfile)

    ## BLAST flanks against K-12 and F-plasmid.
    runBLAST(F_resultsdir, flanksfile)
    runBLAST(K12_resultsdir, flanksfile)

main()
