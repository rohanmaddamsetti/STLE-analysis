#!/usr/bin/env python3

import fileinput

## Usage: python3 dS_dN_data.py ../annotated-diffs/REL4397.gd

def main():

    locus_list = [] # for a sorted order.
    dS_dict = {}
    dN_dict = {}
    gene_dict = {}

    for line in fileinput.input():
        snp_type = None
        if "snp_type=synonymous" in line:
            snp_type = 'synonymous'
        elif "snp_type=nonsynonymous" in line:
            snp_type = 'nonsynonymous'
        else:
            continue
        data = line.split()
        gene = data[13][10:]
        locus_tag = data[-3][10:]
        if locus_tag not in locus_list:
            locus_list.append(locus_tag)        
        if locus_tag not in gene_dict:
            gene_dict[locus_tag] = gene
        if snp_type == 'synonymous':
            if locus_tag in dS_dict:
                dS_dict[locus_tag] += 1
            else:
                dS_dict[locus_tag] = 1
                dN_dict[locus_tag] = 0
        elif snp_type == 'nonsynonymous':
            if locus_tag in dN_dict:
                dN_dict[locus_tag] += 1
            else:
                dN_dict[locus_tag] = 1
                dS_dict[locus_tag] = 0
        else:
            raise AssertionError()
    print("locus_tag,gene,dS,dN")
    for locus in locus_list:
        print(locus + ',' + gene_dict[locus] + ',' + str(dS_dict[locus]) + ',' + str(dN_dict[locus]))
        
main()
