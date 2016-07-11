#!/usr/bin/env python

## align_special_cases.py by Rohan Maddamsetti.

from os.path import join, basename, exists
from os import listdir, makedirs, chdir, getcwd
import sys
from pprint import pprint

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline as BLAST
from Bio.Blast import NCBIXML
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline as MAFFT
from io import StringIO
from Bio import AlignIO
from copy import deepcopy

def get_gene_and_201bp_upstream(genefeature,genomeseq):
    mystart = genefeature.location.start
    myend = genefeature.location.end
    mystrand = genefeature.location.strand
    if mystrand == 1:
        newfeature = SeqFeature(FeatureLocation(mystart-201,myend),strand=mystrand)
    elif mystrand == -1:
        newfeature = SeqFeature(FeatureLocation(mystart,myend+201),strand=mystrand)
    return newfeature.extract(genomeseq)

def get_K12_seqs(proj_dir, a_minus4_genes,a_plus1_genes,include_upstream):
    ## this code is a bit more complicated to check for name differences between REL606 and K12 genes,
    ## that is I check gene synonyms as well.
    K12_seqs = {}
    K12_genome = next(SeqIO.parse(open(join(proj_dir,"references/K-12.1.gbk")),"genbank"))
    for feat in K12_genome.features:
        if feat.type == 'CDS' and 'gene' in feat.qualifiers:
            gene_name = feat.qualifiers['gene'][0]
            gene_synonyms = [x.strip() for x in feat.qualifiers['gene_synonym'][0].split(';')]
            all_gene_names = set([gene_name] + gene_synonyms)
            if all_gene_names & set(a_minus4_genes):
                this_gene1 = set(all_gene_names & set(a_minus4_genes)).pop()
                if include_upstream:
                    K12_seqs['Ara-4_'+this_gene1] = get_gene_and_201bp_upstream(feat,K12_genome.seq)
                else:
                    K12_seqs['Ara-4_'+this_gene1] = feat.extract(K12_genome.seq)
            if all_gene_names & set(a_plus1_genes):
                this_gene2 = set(all_gene_names & set(a_plus1_genes)).pop()
                if include_upstream:
                    K12_seqs['Ara+1_'+this_gene2] = get_gene_and_201bp_upstream(feat,K12_genome.seq)
                else:
                    K12_seqs['Ara+1_'+this_gene2] = feat.extract(K12_genome.seq)
    return K12_seqs

def get_REL606_seqs(proj_dir, a_minus4_genes,a_plus1_genes,include_upstream):
    REL606_seqs = {}
    REL606_genome = next(SeqIO.parse(open(join(proj_dir,"references/REL606.6.gbk")),"genbank"))
    for feat in REL606_genome.features:
        if feat.type == 'CDS' and 'gene' in feat.qualifiers:
            this_gene = feat.qualifiers['gene'][0]
            if this_gene in a_minus4_genes:
                if include_upstream:
                    REL606_seqs['Ara-4_'+this_gene] = get_gene_and_201bp_upstream(feat,REL606_genome.seq)
                else:
                    REL606_seqs['Ara-4_'+this_gene] = feat.extract(REL606_genome.seq)
            if this_gene in a_plus1_genes:
                if include_upstream:
                    REL606_seqs['Ara+1_'+this_gene] = get_gene_and_201bp_upstream(feat,REL606_genome.seq)
                else:
                    REL606_seqs['Ara+1_'+this_gene] = feat.extract(REL606_genome.seq)
    return REL606_seqs

def print_alns(msa_dict,print_type):
    ''' function depends on value of print_type.
    1) print full alignment.
    2) print alignment, filtering consensus columns.
    3) print alignment, masking consensus columns as dots. '''
    ## define an sort order for sequences in the alignment.
    aln_order = {'REL606':0, 'recipient':1,'recombinant':2,'K12':3}
    ## print out a nice human readable alignment (mask consensus columns).
    for k,v in sorted(msa_dict.items()):
        v.sort(key=lambda record: aln_order[record.id])
        print(k)

        ## filter identical sites from the alignment.
        keep_cols = []
        for i in range(v.get_alignment_length()):
            column_nucs = [sr.seq[i] for sr in v]
            if len(set(column_nucs)) > 1:
                keep_cols.append(i)

        v2 = deepcopy(v)
        for x in range(len(v)):
            v2[x].seq = Seq(''.join([v[x].seq[i] for i in keep_cols]))

        v3 = deepcopy(v)
        for x in range(len(v)):
            v3[x].seq = Seq(''.join( [v[x].seq[i] if i in keep_cols else '.' for i in range(v.get_alignment_length())]))

        print()
        print('full seq length:',v.get_alignment_length())
        if print_type == 1:
            print(v.format('phylip'))
        elif print_type == 2:
            print(v2.format('phylip'))
        elif print_type == 3:
            print(v3.format('clustal'))
        print('*****************************************')


def make_alignments(upstream_bool=False,make_protein_aln=False):

    ## Never translate upstream regions!
    assert not (make_protein_aln and upstream_bool)

    proj_dir = "/Users/Rohandinho/Desktop/Projects/recombinant-assembly/"
    align_f = open(join(proj_dir,"results/align_these.csv"))

    genes_to_align = []
    for i,line in enumerate(align_f):
        if i == 0:
            continue
        if '/' in line:
            continue
        if 'ECB' in line:
            continue
        if '[' in line:
            continue
        if 'ins' in line:
            continue
        line = line.strip()
        data = line.split(',')
        data.pop(0)
        ## cheap way to remove inside quotes:
        data = [x[1:-1] for x in data]
        genes_to_align.append(data)

    a_minus4_genes = [x[0] for x in genes_to_align if x[1] == 'Ara-4']
    a_plus1_genes = [x[0] for x in genes_to_align if x[1] == 'Ara+1']

    ## initialize the data structure that holds the sequences to align.
    ## an example of a key is 'Ara-1_nadR'.
    alignment_dict = {x[1]+'_'+x[0]: {'K12':"", 'REL606':"", 'recipient':"", 'recombinant':""} for x in genes_to_align}

    ## get the sequences in the K-12 genome.
    K12_seqs = get_K12_seqs(proj_dir, a_minus4_genes,a_plus1_genes,include_upstream=upstream_bool)

    ## put K-12 sequences into alignment_dict.
    for x in K12_seqs:
        alignment_dict[x]['K12'] = str(K12_seqs[x])

    ## get the sequences in the B genome.
    REL606_seqs = get_REL606_seqs(proj_dir, a_minus4_genes,a_plus1_genes,include_upstream=upstream_bool)

    ## put REL606 sequences into alignment_dict.
    for x in REL606_seqs:
        alignment_dict[x]['REL606'] = str(REL606_seqs[x])

    fasta_base = join(proj_dir,"results/special-alignments/")

    ## write BLAST queries to file.
    REL606_plus1_query_file = join(fasta_base,'Ara+1/REL606_queries.fasta')
    REL606_plus1_queries = [SeqRecord(v,id=k,description='') for k,v in REL606_seqs.items() if k.startswith('Ara+1')]
    SeqIO.write(REL606_plus1_queries, REL606_plus1_query_file,'fasta')

    REL606_minus4_query_file = join(fasta_base,'Ara-4/REL606_queries.fasta')
    REL606_minus4_queries = [SeqRecord(v,id=k,description='') for k,v in REL606_seqs.items() if k.startswith('Ara-4')]
    SeqIO.write(REL606_minus4_queries, REL606_minus4_query_file,'fasta')

    K12_plus1_query_file = join(fasta_base,'Ara+1/K12_queries.fasta')
    K12_plus1_queries = [SeqRecord(v,id=k,description='') for k,v in K12_seqs.items() if k.startswith('Ara+1')]
    SeqIO.write(K12_plus1_queries, K12_plus1_query_file,'fasta')

    K12_minus4_query_file = join(fasta_base,'Ara-4/REL606_queries.fasta')
    K12_minus4_queries = [SeqRecord(v,id=k,description='') for k,v in K12_seqs.items() if k.startswith('Ara-4')]
    SeqIO.write(K12_minus4_queries, K12_minus4_query_file,'fasta')

    ## I made BLAST databases containing the recipient and the recombinant FASTA files
    ## for Ara+1 and Ara-4 (these are full genomes).
    for ref in ['REL606','K12']:
        for lineage in ['Ara+1','Ara-4']:
            for db in ['recipient','recombinant']:
                if ref == 'REL606':
                    if lineage == 'Ara+1':
                        this_query_file = REL606_plus1_query_file
                    elif lineage == 'Ara-4':
                        this_query_file = REL606_minus4_query_file
                elif ref == 'K12':
                    if lineage == 'Ara+1':
                        this_query_file = K12_plus1_query_file
                    elif lineage == 'Ara-4':
                        this_query_file = K12_minus4_query_file
                this_db = join(fasta_base, lineage+'/blast_dbs/'+db+'.fasta')
                this_blast_out = join(fasta_base,
                                      lineage+'/blast_results/'+ref+'_'+db+'_results.xml')
                this_blast = BLAST(cmd='blastn',
                                   query=this_query_file,
                                   db=this_db,
                                   out=this_blast_out,
                                   outfmt=5,
                                   max_target_seqs=1,)
                ## run BLASTN.
                this_blast()
                ## Now parse the output to finish populating the fields in alignment_dict.
                result_handle = open(this_blast_out)
                blast_records = NCBIXML.parse(result_handle)
                for rec in blast_records:
                    ## we only want the first hit for each query.
                    y = rec.alignments[0].hsps[0].sbjct
                    ## strip gap characters from the sequence before aligning with MAFFT.
                    alignment_dict[rec.query][db] = y.replace("-","")

    ## make MSAs with mafft using the sequences stored in alignment_dict.
    ## write the alignments to file for mafft to use, run mafft, and store the MSAs.
    msa_dict = {}
    for k,v in alignment_dict.items():
        this_alignment_out = join(fasta_base, 'mafft-input/'+k+'.fasta')
        if make_protein_aln:
            these_align_seqs = [SeqRecord(Seq(dna).translate(),id=nm,description='') for nm,dna in v.items()]
        else:
            these_align_seqs = [SeqRecord(Seq(dna),id=nm,description='') for nm,dna in v.items()]
        SeqIO.write(these_align_seqs, this_alignment_out, 'fasta')
        ## run mafft.
        run_mafft = MAFFT(input=this_alignment_out)
        stdout, stderr = run_mafft()
        ## parse and save mafft output in a data structure.
        msa_dict[k] = AlignIO.read(StringIO(stdout), "fasta")

    ## I checked strange alignments (all in Ara+1) by hand, and I found that nadR, cycA,
    ## and yfcTU have IS150 insertions in the recipient clones. gltL has an IS150 insertion
    ##in REL606 that is preserved in the recipient, but not in K12 or the recombinant clone.

    return msa_dict

def main():

    ## Usage: python align_special_cases.py > special-alignments.txt
    prot_dict = make_alignments(upstream_bool=False,make_protein_aln=True)
    print_alns(prot_dict,3)

    ## as a sanity check, compare the upstream run to the run that simply extracts the
    ## ORF in a different text file.

    ## Usage: python align_special_cases.py > special-alignments2.txt
    #msa_dict = make_alignments(upstream_bool=True)
    #msa_dict = make_alignments(upstream_bool=False)
    #print_alns(msa_dict,3)

main()
