#!/usr/bin/env python

## align_replaced_mutations.py by Rohan Maddamsetti.

'''
1) Automatically make blast dbs.
2) annotate alignments as 0) REL606 state, 1) K-12 state, 2) new state.

For each lineage (using odd recombinant genomes) count the number of replaced
mutations in each class, and make a figure.
'''

from os.path import join, basename, exists
from os import listdir, makedirs, chdir, getcwd
import sys
from pprint import pprint
from copy import deepcopy
import subprocess
import argparse


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

def get_gene_and_201bp_upstream(genefeature,genomeseq):
    mystart = genefeature.location.start
    myend = genefeature.location.end
    mystrand = genefeature.location.strand
    if mystrand == 1:
        newfeature = SeqFeature(FeatureLocation(mystart-201,myend),strand=mystrand)
    elif mystrand == -1:
        newfeature = SeqFeature(FeatureLocation(mystart,myend+201),strand=mystrand)
    return newfeature.extract(genomeseq)

def get_K12_seqs(proj_dir,p, p_genes,include_upstream):
    '''
    more complicated to check for name differences between REL606 and K12 genes;
    that is I check gene synonyms as well.

    p is the name of the population (e.g. 'Ara+3').
    '''

    K12_seqs = {}
    K12_genome = next(SeqIO.parse(open(join(proj_dir,"references/K-12.1.gbk")),"genbank"))
    for feat in K12_genome.features:
        if feat.type == 'CDS' and 'gene' in feat.qualifiers:
            gene_name = feat.qualifiers['gene'][0]
            gene_synonyms = [x.strip() for x in feat.qualifiers['gene_synonym'][0].split(';')]
            all_gene_names = set([gene_name] + gene_synonyms)
            if all_gene_names & set(p_genes):
                this_gene = set(all_gene_names & set(p_genes)).pop()
                p_gene_key = p + '_' + this_gene
                if include_upstream:
                    K12_seqs[p_gene_key] = get_gene_and_201bp_upstream(feat,K12_genome.seq)
                else:
                    K12_seqs[p_gene_key] = feat.extract(K12_genome.seq)
    return K12_seqs

def get_REL606_seqs(proj_dir, p, p_genes,include_upstream):
    '''
        p is the name of the population (e.g. 'Ara+3').
    '''
    REL606_seqs = {}
    REL606_genome = next(SeqIO.parse(open(join(proj_dir,"references/REL606.7.gbk")),"genbank"))
    for feat in REL606_genome.features:
        if feat.type == 'CDS' and 'gene' in feat.qualifiers:
            this_gene = feat.qualifiers['gene'][0]
            if this_gene in p_genes:
                p_gene_key = p + '_' + this_gene
                if include_upstream:
                    REL606_seqs[p_gene_key] = get_gene_and_201bp_upstream(feat,REL606_genome.seq)
                else:
                    REL606_seqs[p_gene_key] = feat.extract(REL606_genome.seq)
    return REL606_seqs

def setup_blastdbs():
    '''
    this only blasts the odd REL-numbered recombinants.
    some lines in the Ara-3 recombinant gd files cause problems, comment out those lines for now.
    '''

    projdir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"
    basedir = join(projdir,"results/replaced-alignments/")

    pops = ['Ara+1','Ara+2','Ara+3','Ara+4','Ara+5','Ara+6',
            'Ara-1','Ara-2','Ara-3','Ara-4','Ara-5','Ara-6']

    recipient_dict = {'Ara+1':'REL2537', 'Ara+2':'REL2538', 'Ara+3':'REL2539',
                      'Ara+4':'REL2540', 'Ara+5':'REL2541', 'Ara+6':'REL2542',
                      'Ara-1':'REL2543', 'Ara-2':'REL2544', 'Ara-3':'REL2545',
                      'Ara-4':'REL2546', 'Ara-5':'REL2547', 'Ara-6':'REL2548'}

    ## to convert my IDs to REL ID for the recombinant clones.
    ## REMEMBER: we want ODD REL clones (which have EVEN RM numbers),
    ## so choose the EVEN RM clones.
    odd_recombinant_dict = {p:'RM3-130-'+str(i) for p,i in zip(pops,range(2,24+1,2))}

    for p in pops:
        cur_dir = join(basedir,p)
        if not exists(cur_dir):
            makedirs(cur_dir)
        if not exists(join(cur_dir,'blast_dbs')):
            makedirs(join(cur_dir,'blast_dbs'))
        if not exists(join(cur_dir,'blast_results')):
            makedirs(join(cur_dir,'blast_results'))
        chdir(join(cur_dir,'blast_dbs'))
        ## first, run gdtools ANNOTATE to make recipient.fasta and recombinant.fasta.
        refgenome = join(projdir,'references/REL606.7.gbk')
        refF = join(projdir,'references/F-plasmid.1.gbk')
        refTn10 = join(projdir,'references/Tn10.gbk')
        cur_gdpath = join(projdir,'annotated-diffs/REL606-ref-runs/',p)
        cur_recipient = join(cur_gdpath,'annotated_' + recipient_dict[p] +'.gd')
        cur_recombinant = join(cur_gdpath,'annotated_' + odd_recombinant_dict[p] + '.gd')
        partial_args = ['gdtools','APPLY','-r',refgenome,'-r',refF, '-r',refTn10,'-f','FASTA','-o']
        if not exists("recipient.fasta"):
            subprocess.run(partial_args + ["recipient.fasta", cur_recipient])
        if not exists("recombinant.fasta"):
            subprocess.run(partial_args + ["recombinant.fasta", cur_recombinant])
        ## second, if blastdbs haven't been made, make blast dbs.
        if not exists("recipient.fasta.nhr"):
            recipient_args = ['makeblastdb','-in','recipient.fasta',
                              '-dbtype','nucl' ]
            subprocess.run(recipient_args)
        if not exists("recombinant.fasta.nhr"):
            recombinant_args = ['makeblastdb','-in','recombinant.fasta',
                              '-dbtype','nucl' ]
            subprocess.run(recombinant_args)


def make_alignments(upstream_bool=False,make_protein_aln=False):
    '''
    Make both protein and DNA alignments in order to compare both.
'''

    ## Never translate upstream regions!
    assert not (make_protein_aln and upstream_bool)

    proj_dir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/"
    align_f = open(join(proj_dir,"results/align_these.csv"))
    fasta_base = join(proj_dir,"results/replaced-alignments/")

    if make_protein_aln:
        mafft_indir = join(fasta_base,'protein-mafft-input')
        mafft_outdir = join(fasta_base,'protein-mafft-output')
    else:
        mafft_indir = join(fasta_base,'dna-mafft-input')
        mafft_outdir = join(fasta_base,'dna-mafft-output')

    pops = ['Ara+1','Ara+2','Ara+3','Ara+4','Ara+5','Ara+6',
            'Ara-1','Ara-2','Ara-3','Ara-4','Ara-5','Ara-6']

    genes_to_align = []
    for i,line in enumerate(align_f):
        if i == 0:
            continue
        if 'ECB' in line:
            continue
        line = line.strip()
        data = line.split(',')
        genes_to_align.append(data)


    ## initialize the data structure that holds the sequences to align.
    ## an example of a key is 'Ara-1_nadR'.
    alignment_dict = {x[1]+'_'+x[0]: {'K12':"", 'REL606':"", 'recipient':"", 'recombinant':""} for x in genes_to_align}

    for p in pops:
        p_genes = [x[0] for x in genes_to_align if x[1] == p]
        if not p_genes: # skip if no alignments to examine.
            continue
        ## get the sequences in the K-12 genome.
        K12_seqs = get_K12_seqs(proj_dir, p, p_genes,include_upstream=upstream_bool)

        ## put K-12 sequences into alignment_dict.
        for x in K12_seqs:
            alignment_dict[x]['K12'] = str(K12_seqs[x])

        ## get the sequences in the B genome.
        REL606_seqs = get_REL606_seqs(proj_dir, p, p_genes,include_upstream=upstream_bool)

        ## put REL606 sequences into alignment_dict.
        for x in REL606_seqs:
            alignment_dict[x]['REL606'] = str(REL606_seqs[x])

        ## write BLAST queries to file.
        REL606_query_file = join(fasta_base,p,'REL606_queries.fasta')
        REL606_queries = [SeqRecord(v,id=k,description='') for k,v in REL606_seqs.items() if k.startswith(p)]
        SeqIO.write(REL606_queries, REL606_query_file,'fasta')

        K12_query_file = join(fasta_base,p,'K12_queries.fasta')
        K12_queries = [SeqRecord(v,id=k,description='') for k,v in K12_seqs.items() if k.startswith(p)]
        SeqIO.write(K12_queries, K12_query_file,'fasta')

        for ref in ['REL606','K12']:
            for db in ['recipient','recombinant']:
                if ref == 'REL606':
                    this_query_file = REL606_query_file
                elif ref == 'K12':
                    this_query_file = K12_query_file
                this_db = join(fasta_base, p,'blast_dbs',db+'.fasta')
                this_blast_out = join(fasta_base,p,
                              'blast_results',ref+'_'+db+'_results.xml')
                this_blast = BLAST(cmd='blastn',
                           query=this_query_file,
                           db=this_db,
                           out=this_blast_out,
                           outfmt=5,
                           max_target_seqs=1)
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
    if not exists(mafft_indir):
        makedirs(mafft_indir)
    if not exists(mafft_outdir):
        makedirs(mafft_outdir)
    msa_dict = {}
    for k,v in alignment_dict.items():
        this_alignment_in = join(mafft_indir, k+'.fasta')
        if make_protein_aln:
            these_align_seqs = [SeqRecord(Seq(dna).translate(),id=nm,description='') for nm,dna in v.items()]
        else:
            these_align_seqs = [SeqRecord(Seq(dna),id=nm,description='') for nm,dna in v.items()]
        SeqIO.write(these_align_seqs, this_alignment_in, 'fasta')
        ## run mafft.
        run_mafft = MAFFT(input=this_alignment_in)
        stdout, stderr = run_mafft()
        this_aln_out = join(mafft_outdir, k+'.fasta')
        ## parse and save mafft output in a data structure.
        msa_dict[k] = AlignIO.read(StringIO(stdout), "fasta")
        with open(this_aln_out, 'w') as handle:
            SeqIO.write(msa_dict[k], handle, "fasta")
## I checked strange alignments (all in Ara+1) by hand, and I found that nadR, cycA,
    ## and yfcTU have IS150 insertions in the recipient clones. gltL has an IS150 insertion
    ##in REL606 that is preserved in the recipient, but not in K12 or the recombinant clone.
    return msa_dict

def make_msa_dict(proj_dir,protein_aln=False):
    '''
    input: a directory containing aligned sequences.
    output: a dictionary containing the MSA.
'''
    msa_dict = {}
    if protein_aln:
        mafft_out = join(proj_dir,"results/replaced-alignments/protein-mafft-output")
    else:
        mafft_out = join(proj_dir,"results/replaced-alignments/dna-mafft-output")
    for f in [x for x in listdir(mafft_out) if x.endswith('fasta')]:
        k = basename(f).split('.')[0]
        full_f = join(mafft_out,f)
        msa_dict[k] = AlignIO.read(full_f,'fasta')
    return msa_dict

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

def count_aln_type(msa_dict):
    '''
    annotate alignments as 0) REL606 state, 1) recipient state, 2) K-12 state 3) New state.
    (in protein world, not DNA world).

    if the sequence in REL606 and K-12 is identical, mark as REL606 state (ancestral state).

'''
    aln_type_count = {'REL606':0,'Recipient':0,'K-12':0,'New':0,'REL606/K-12':0}

    ## make a dict of dicts to turn into a pandas dataframe.

   ## define an sort order for sequences in the alignment.
    aln_order = {'REL606':0, 'recipient':1,'recombinant':2,'K12':3, 'REL606/K-12':4}

    ## lists for DataFrame columns.
    lineage_col = []
    gene_col = []
    allele_type_col = []

    for k,v in sorted(msa_dict.items()):
        lin,gen = k.split('_')
        lineage_col.append(lin)
        gene_col.append(gen)

        v.sort(key=lambda record: aln_order[record.id])
        if v[2].seq == v[0].seq == v[3].seq:
            aln_type_count['REL606/K-12'] += 1
            print(' '.join([k,'REL606/K-12']))
            #print(' '.join([k,'4']))
            ## note: mark as REL606 allele in this case.
            allele_type_col.append('REL606 allele')
        elif v[2].seq == v[0].seq:
            aln_type_count['REL606'] += 1
            print(' '.join([k,'REL606']))
            #print(' '.join([k,'0']))
            allele_type_col.append('REL606 allele')
        elif v[2].seq == v[1].seq:
            aln_type_count['Recipient'] += 1
            print(' '.join([k,'Recipient']))
            #print(' '.join([k,'2']))
            allele_type_col.append('Recipient allele')
        elif v[2].seq == v[3].seq:
            aln_type_count['K-12'] += 1
            print(' '.join([k,'K-12']))
            #print(' '.join([k,'1']))
            allele_type_col.append('K-12 allele')
        elif v[2].seq != v[0].seq and v[2].seq != v[1].seq and v[2].seq != v[3].seq:
            aln_type_count['New'] += 1
            print(' '.join([k,'New']))
            #print(' '.join([k,'2']))
            allele_type_col.append('New allele')
        else:
            print('ERROR: unknown allele type')
            print(k)
            quit()

    print(aln_type_count)

def main():

    ## Usage: python align_replaced_mutations.py 3
    parser = argparse.ArgumentParser(description='print alignments for replaced mutations.')
    parser.add_argument('aln_print_mode', type=int,help='number representing how to print alignment')
    args = parser.parse_args()
    if args.aln_print_mode not in (1,2,3):
        print("error:print mode must be 1, 2, or 3")
        quit()

    proj_dir = "../"
    fasta_base = join(proj_dir,"results/replaced-alignments/")

    setup_blastdbs()

    prot_mafft_outdir = join(fasta_base,'protein-mafft-output')
    if exists(prot_mafft_outdir) and listdir(prot_mafft_outdir):
        prot_dict = make_msa_dict(proj_dir,protein_aln=True)
    else:
        prot_dict = make_alignments(upstream_bool=False,make_protein_aln=True)

    print_alns(prot_dict,args.aln_print_mode)
    count_aln_type(prot_dict)

    dna_mafft_outdir = join(fasta_base,'dna-mafft-output')
    if exists(dna_mafft_outdir) and listdir(dna_mafft_outdir):
        dna_dict = make_msa_dict(proj_dir, protein_aln=False)
    else:
        dna_dict = make_alignments(upstream_bool=True,make_protein_aln=False)

    ## now for DNA world as well. Include bp upstream.
    print_alns(dna_dict,args.aln_print_mode)
    count_aln_type(dna_dict)


main()
