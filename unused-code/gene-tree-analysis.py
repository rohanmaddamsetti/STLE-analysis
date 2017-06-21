#!/usr/bin/env python

'''
gene-tree-analysis.py by Rohan Maddamsetti

Usage: pythonw gene-tree-analysis.py (need python framework for treeCl to work)

## TODO: might be faster to run this as an ipython notebook.

-test my synonymous genetic variation model using
  experimental data.
  Treat Clones and LCA as data, and K-12 as hidden
  variable. How many gene trees have elongated branches, but same
  topology? basic premise: if recombination + selection preserves the clonal frame of a
                  gene, that should be sufficient to explain a lot of
                  synonymous variation. How many genes have preserved
   clonal frame? Do these genes have elongated branches but same
   topology?

How has bacterial recombination affected gene trees?

the gene tree-lengths for non-synonymous (pN)
 and synonymous sites (pS)
 were estimated using CODEML from the PAML package (version 4.2)
'''

from Bio import SeqIO
from subprocess import run
import os
import treeCl
from treeCl.clustering import spectral, methods
import pickle
from pandas import DataFrame
from numpy import savetxt

def get_reference_CDS(rel606path):
    ''' get sequence, start and end of every CDS in REL606'''
    gene_dict = {}
    ref_genome = SeqIO.read(rel606path, "genbank")
    for feature in ref_genome.features:
        if feature.type != "CDS":
            continue
        try:
            locus_tag = feature.qualifiers["locus_tag"][0]
        except:
            continue
        start = feature.location.start.position
        end = feature.location.end.position
        strand = feature.strand
        sequence = ref_genome.seq[start:end]
        gene_dict[locus_tag] = {'seq':sequence,'start':start,'end':end,'strand':strand}
    return gene_dict

def filter_for_SNPs(gddir,filterdir):
    clone_files = [(x,os.path.join(gddir,x,"output","output.gd")) for x in os.listdir(gddir) if 'RM3' in x]
    for f_tup in clone_files:
        outf = os.path.join(filterdir,f_tup[0]+'.gd')
        outfh = open(outf,'w')
        for line in open(f_tup[1]):
            line = line.strip()
            if line.startswith('#') or line.startswith('SNP'):
                print(line,file=outfh)
            else:
                print('#'+line,file=outfh)

def run_gdtools_apply(filterdir,ref_file, fastadir):
    filtered_gds = [os.path.join(filterdir,x) for x in os.listdir(filterdir) if 'RM3' in x]
    for f in filtered_gds:
        outf = os.path.join(fastadir,os.path.basename(f))[:-3] + ".fasta"
        gdtools_args = ['gdtools', 'APPLY','-r', ref_file, '-o', outf, f]
        print(' '.join(gdtools_args))
        ## This step is slow. Only run if output doesn't already exist.
        if not os.path.exists(outf):
            run(gdtools_args)

def make_gene_aln_inputs(fastadir,gdict,aln_indir):

    fastadict = {}
    for f in [x for x in os.listdir(fastadir) if x.endswith('.fasta')]:
        fullf = os.path.join(fastadir,f)
        ##print(fullf)
        my_name = os.path.basename(fullf)[:-6]
        x = SeqIO.read(fullf,"fasta").seq
        fastadict[my_name] = x

    for g in sorted(gdict):
        start = gdict[g]['start']
        end = gdict[g]['end']
        strand = gdict[g]['strand']
        my_outf = os.path.join(aln_indir,g+".fasta")
        ## skip if the output file already exists.
        if os.path.exists(my_outf):
            continue
        cur_gene_outf = open(my_outf,'w')
        for f in fastadict:
            print('>'+f,file=cur_gene_outf)
            if strand == 1:
                my_gene = fastadict[f][start:end]
            elif strand == -1:
                my_gene = fastadict[f][start:end].reverse_complement()
            print(my_gene,file=cur_gene_outf)

def score_clustering_results(c, partition, resultdir, nclusters):
    my_cache_dir = os.path.join(resultdir,"mdsclust") + str(nclusters)
    raxml = treeCl.tasks.RaxmlTaskInterface()
    sc = treeCl.Scorer(c, cache_dir=my_cache_dir, task_interface=raxml)
    sc.write_partition(partition)
    results = sc.analyse_cache_dir(executable='raxml', threads=8)
    full_results = sc.get_partition_results(partition)
    return full_results

def main():

    ## if running on orchestra, module load breseq and source activate evcouplings_env.

    ## assume script is called from STLE-analysis/src
    assert os.getcwd().endswith('/src')
    projdir = os.path.dirname(os.getcwd())

    ref_path = os.path.join(projdir,"references/REL606.7.gbk")
    ## get start and end coordinates of REL606 genes.
    ##gdict = get_reference_CDS(ref_path)

    resultdir = os.path.join(projdir,'results/gene-tree-analysis')
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)

    gddir = os.path.join(projdir,"breseq-assemblies/REL606-ref-runs/recombinants")
    ## CHEAP HACK: only allow SNP mutations, so that gdtools APPLY won't screw up the
    ## frame and thus the genome coordinates. might fix later if results are promising.

    filterdir = os.path.join(resultdir,"filtered-gds")
    if not os.path.exists(filterdir):
        os.makedirs(filterdir)

    ##filter_for_SNPs(gddir,filterdir)
    fastadir = os.path.join(resultdir,'evolved-fasta')
    if not os.path.exists(fastadir):
        os.makedirs(fastadir)
    ##run_gdtools_apply(filterdir,ref_path,fastadir)
    ## write out gene sequences to file for alignment.

    aln_indir = os.path.join(resultdir,"alninput")
    if not os.path.exists(aln_indir):
        os.makedirs(aln_indir)
    ##make_gene_aln_inputs(fastadir,gdict,aln_indir)

    '''
    Use TreeCl module to make and cluster trees. The following follows the
    TreeCl documentation on github.
    '''

    ## trees are expensive to calculate. load pickled or cached results if possible.
    treespicklefile = os.path.join(resultdir,"tree_collection_pickle.p")
    cachedir = os.path.join(projdir,'results/gene-tree-analysis/gene-trees')
    if os.path.exists(treespicklefile):
        c = pickle.load(open(treespicklefile,'rb'))
    elif os.path.exists(cachedir):
        c = treeCl.Collection(input_dir=aln_indir, param_dir=cachedir,file_format='fasta')
        pickle.dump(c,open(treespicklefile,'wb'))
    else:
        c = treeCl.Collection(input_dir=aln_indir, file_format='fasta')
        ## use RAxML to calculate trees, default parameters.
        c.calc_trees(executable='raxml')
        ## write results to disk.
        c.write_parameters(cachedir)

    ## use the geodesic distance: uses both branch lengths and topology.
    ## this is super expensive, so pickle.
    treedist_pickle_file = os.path.join(resultdir,"treedist_pickle.p")
    if os.path.exists(treedist_pickle_file):
        distmatrix = pickle.load(open(treedist_pickle_file, 'rb'))
    else:
        ## with pure python code, it is better to use processpools to parallelise for speed.
        processes = treeCl.parutils.ProcesspoolJobHandler(8)
        ## jobs are done in batches to reduce overhead.
        distmatrix = c.get_inter_tree_distances('geo',
                                        jobhandler=processes,
                                        batchsize=1000)
        pickle.dump(dm, open(treedist_pickle_file, 'wb'))


    ## rather than using a likelihood ratio criterion
    ## to find the number of clusters, just choose 20 clusters for now.
    scorer_pickle_file = os.path.join(resultdir, "scorer_pickle.p")
    partition_pickle_file = os.path.join(resultdir, "partition_pickle.p")
    if not os.path.exists(scorer_pickle_file) or not os.path.exists(partition_pickle_file):
        mdsclust = treeCl.MultidimensionalScaling(distmatrix)
        partition =  mdsclust.cluster(20)
        ##Score the result likelihood
        my_cache_dir = os.path.join(resultdir,"mdsclust") + str(20)
        raxml = treeCl.tasks.RaxmlTaskInterface()
        sc = treeCl.Scorer(c, cache_dir=my_cache_dir, task_interface=raxml)
        sc.write_partition(partition)
        results = sc.analyse_cache_dir(executable='raxml', threads=8)
        full_results = sc.get_partition_results(partition)
        pickle.dump(sc, open(scorer_pickle_file, 'wb'))
        pickle.dump(partition, open(partition_pickle_file, 'wb'))
    else:
        sc = pickle.load(open(scorer_pickle_file,'rb'))
        partition = pickle.load(open(partition_pickle_file,'rb'))
        ## Get concatenated sequence alignments for each group
        concats = [c.concatenate(grp) for grp in partition.get_membership()]
        alignments = [conc.alignment for conc in concats]

        ## make trees from the concatenated alignments.
        concat_aln_dir = os.path.join(resultdir,'concatenated-alignments')
        if not os.path.exists(concat_aln_dir):
            os.makedirs(concat_aln_dir)
        for i,aln in enumerate(alignments):
            clustname = 'cluster' + str(i+1) + '.fasta'
            aln.write_alignment(os.path.join(concat_aln_dir,clustname), "fasta")

        concat_c = treeCl.Collection(input_dir=concat_aln_dir, file_format='fasta')
        concat_c.calc_trees(executable='raxml')
        ## write results to disk.
        concat_cachedir = os.path.join(resultdir,'concat-gene-trees')
        concat_c.write_parameters(concat_cachedir)

        ## Get a list of the loci in each group
        loci = sc.get_partition_members(partition)

        ## Get trees for each group
        trees = sc.get_partition_trees(partition)

        ## Get full model parameters for each group
        full_results = sc.get_partition_results(partition)

main()
