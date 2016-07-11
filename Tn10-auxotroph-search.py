
#!/usr/bin/env python

## Tn10-auxotroph-search.py by Rohan Maddamsetti
## blast all de novo assemblies for Tn10.
## check to see if Tn10 maps to known auxotroph mutations in donors.

from os.path import join, basename, exists, isdir
from os import listdir, makedirs, chdir, getcwd
import subprocess

from Bio.Blast.Applications import NcbiblastnCommandline as BLAST
from Bio.Blast import NCBIXML

def make_blast_dbs_from_de_novo_assemblies():
    proj_dir = "/Users/Rohandinho/Desktop/Projects/recombinant-assembly/"

    result_dir = join(proj_dir,'results','Tn10-auxotroph-search')

    assembly_dir = join(proj_dir,'denovo-assemblies')
    denovo_assemblies = [p for p in listdir(assembly_dir)]
    ## the next line is a stupid filter hack but it works.
    denovo_assemblies = [p for p in denovo_assemblies if p.startswith('R')]

    working_dir = getcwd()

    for cur_dir in denovo_assemblies:
        cur_path = join(assembly_dir,cur_dir)
        blastdb_infile = join(cur_path,'scaffolds.fasta')

        cur_out_dir = join(result_dir,cur_dir,'blast_dbs')
        if not exists(cur_out_dir):
            makedirs(cur_out_dir)

        outfile = join(cur_out_dir,'scaffolds.fasta')
        args = ['makeblastdb', '-in', blastdb_infile, '-dbtype', 'nucl', '-out', outfile]
        #print(args)
        subprocess.call(args)

def blast_Tn10_against_de_novo_assemblies():
    proj_dir = "/Users/Rohandinho/Desktop/Projects/recombinant-assembly/"
    ## use Tn10 reference sequences as a query against denovo assembly blastdbs.
    Tn10_query = join(proj_dir,'references','Tn10.fasta')
    result_dir = join(proj_dir,'results','Tn10-auxotroph-search')
    genomes = [x for x in listdir(result_dir)]
    genomes = [x for x in genomes if isdir(join(result_dir,x)) and x.startswith('R')] ## quik and dirty hack.
    for cur_dir in genomes:
        this_db = join(result_dir,cur_dir,'blast_dbs','scaffolds.fasta')
        cur_out_dir = join(result_dir,cur_dir,'blast_results')
        if not exists(cur_out_dir):
            makedirs(cur_out_dir)
        this_blast_out = join(cur_out_dir,'Tn10_hits.xml')
        this_blast = BLAST(cmd='blastn',
                           query=Tn10_query,
                           db=this_db,
                           out=this_blast_out,
                           outfmt=5)

        #print(this_blast)
        this_blast()

def blast_auxotrophy_genes_against_donors():
    proj_dir = "/Users/Rohandinho/Desktop/Projects/recombinant-assembly/"
    result_dir = join(proj_dir,'results','Tn10-auxotroph-search')
    argA_query = join(result_dir,'queries','argA.fasta')
    argE_query = join(result_dir,'queries','argE.fasta')
    ilvD_query = join(result_dir,'queries','ilvD.fasta')
    leuABCD_query = join(result_dir,'queries','leuABCD.fasta')

    donors = ['REL288', 'REL291', 'REL296', 'REL298']
    for cur_dir in donors:
        this_db = join(result_dir,cur_dir,'blast_dbs','scaffolds.fasta')
        cur_out_dir = join(result_dir,cur_dir,'blast_results')
        if not exists(cur_out_dir):
            makedirs(cur_out_dir)
        if cur_dir == 'REL288':
            this_query = ilvD_query
        elif cur_dir == 'REL291':
            this_query = argE_query
        elif cur_dir == 'REL296':
            this_query = leuABCD_query
        elif cur_dir == 'REL298':
            this_query = argA_query
        this_blast_out = join(cur_out_dir,'auxotroph_hits.xml')
        this_blast = BLAST(cmd='blastn',
                           query=this_query,
                           db=this_db,
                           out=this_blast_out,
                           outfmt=5)

        #print(this_blast)
        this_blast()

def main():
    make_blast_dbs_from_de_novo_assemblies()
    blast_Tn10_against_de_novo_assemblies()
    blast_auxotrophy_genes_against_donors()

main()
