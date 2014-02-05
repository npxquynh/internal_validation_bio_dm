import os
import csv
import pdb

import helper

def read_lgn(filepath):
    edges_in_lgn = []
    with open(filepath) as lgn_file:
        lgn_reader = csv.reader(lgn_file)
        lgn_reader.next() # skip the header row

        edges_in_lgn = []
        for row in lgn_reader:
            gene_1, gene_2 = helper.swap(row[0], row[1])
            edges_in_lgn.append([gene_1, gene_2])

    genes_in_lgn = helper.genes_from_edges(edges_in_lgn)

    return genes_in_lgn, edges_in_lgn

