#!/usr/bin/env python

'''
plot-trees.py by Rohan Maddamsetti

hack to plot trees made from the concatenated tree clusters.
'''


import json
from pprint import pprint
import os
from ete3 import Tree


def main():

    tree_list = []
    treedir = "/Users/Rohandinho/Desktop/Projects/STLE-analysis/results/gene-tree-analysis/concat-gene-trees"
    for f in os.listdir(treedir):
        if not f.endswith('.json'):
            continue
        with open(os.path.join(treedir,f)) as data_file:
            data = json.load(data_file)
            jsontree = data['ml_tree']
            t = Tree(jsontree)
            print(jsontree)
            print()


main()
