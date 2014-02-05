import math
import numpy
import operator
from collections import defaultdict

# TODO: remove pdb
import pdb

import helper

class InternalValidationRls:
    def __init__(self, genes_in_tile, blocks, genes_in_lgn, edges_in_lgn):
        self.genes_in_tile = genes_in_tile
        self.blocks = blocks
        self.genes_in_lgn = genes_in_lgn
        self.edges_in_lgn = edges_in_lgn
        self.list_intra = dict()
        self.list_extra = dict()

        self.__generate_zero_matrices()

    def __find_genes_connected_with_LGN(self, block):
        """
        Find both (i) extra genes and (ii) LGN genes connected with LGN gene

        Output: [set(extra_genes_connected_with_LGN),
        set(LGN_gene_connected_with_LGN)]
        """
        # connected_nodes[i] = set() = nodes connected with LGN[i]
        connected_nodes = dict()
        for gene in self.genes_in_lgn:
            nodes = set()
            for edge in block:
                try:
                    if edge[0] == gene:
                        nodes.add(edge[1])
                    if edge[1] == gene:
                        nodes.add(edge[0])
                except IndexError:
                    print "Index error"
            extra_nodes = set(nodes) - set(self.genes_in_lgn)
            LGN_nodes = set(nodes).intersection(set(self.genes_in_lgn))

            connected_nodes[gene] = [extra_nodes, LGN_nodes]

        return connected_nodes

    def __find_extra_genes_connected_with_LGN(self, connected_nodes):
        """
        connected_nodes: = output of __find_genes_connected_with_LGN
        """
        extra_genes = set()
        for key in connected_nodes:
            extra_genes = extra_genes.union(connected_nodes[key][0])
        return extra_genes

    def __count_appearance_of_extra_genes(self, extra_gene):
        count = 0
        for set_of_genes in self.genes_in_tile:
            if extra_gene in set_of_genes:
                count += 1
        return count

    def __refine_list_extra(self):
        extra_genes = list(self.list_extra.keys())
        for key in self.list_extra:
            [_, extra_gene] = helper.decode_edge_key(key)
            count = self.__count_appearance_of_extra_genes(extra_gene)
            self.list_extra[key][1] = count

    def create_list_intra_extra(self):
        print 'create_list_intra_extra'

        for block in self.blocks:
            connected_nodes = self.__find_genes_connected_with_LGN(block)

            # create list_intra
            # TODO: remove duplicated edge between LGN
            # since A, B belongs to LGN
            # then connected_nodes[A] = {[], [B]}
            # connected_nodes[B] = {[], [A]}
            # But they're actually the same connection
            for key in connected_nodes:
                for LGN_connected_gene in connected_nodes[key][1]:
                    ordered_key = helper.swap(key, LGN_connected_gene)
                    temp_key = helper.encode_edge_key(ordered_key[0], ordered_key[1])
                    if temp_key in self.list_intra:
                        self.list_intra[temp_key] += 1
                    else:
                        self.list_intra[temp_key] = 1

            # create list_extra
            for key in connected_nodes:
                for extra_connected_gene in connected_nodes[key][0]:
                    # TODO: without doing swap, can the code run faster
                    """ I think by placing the LGN gene first, and the
                    extra_connected_gene later, we can always ensure the
                    unique of connection
                    """
                    temp_key = helper.encode_edge_key(key, extra_connected_gene)
                    if temp_key in self.list_extra:
                        self.list_extra[temp_key][0] += 1
                    else:
                        self.list_extra[temp_key] = [1, 0]

        self.__refine_list_extra()

    def __calculate_the_frequency(self):
        n = len(self.blocks)
        for key in self.list_intra:
            self.list_intra[key] /= float(n)
        for key in self.list_extra:
            self.list_extra[key][0] /= float(self.list_extra[key][1])

    def __generate_zero_matrices(self):
        f = [x / 100.0 for x in range(100)]
        ncols = 100
        self.TP = [ 0 for i in range(ncols) ]
        self.FN = [ 0 for i in range(ncols) ]
        self.FP = [ 0 for i in range(ncols) ]

    def __connection_is_in_lgn(self, gene_1, gene_2):
        """
        connection in 'edges_in_lgn' is already sorted
        """
        for connection in self.edges_in_lgn:
            if gene_1 == connection[0] and gene_2 == connection[1]:
                return True
        return False

    def calculate_TPFNFP(self):
        ncols = 100
        # TP & FN
        for key in self.list_intra:
            gene_1, gene_2 = helper.decode_edge_key(key)
            if self.__connection_is_in_lgn(gene_1, gene_2):
                equal_point = int(math.floor(self.list_intra[key] * 100))
                print equal_point
                for i in range(equal_point, ncols):
                    self.TP[i] += 1
                for i in range(0, equal_point):
                    self.FN[i] += 1
            else:
                for i in range(0, ncols):
                    self.FP[i] += 1

        for key in self.list_extra:
            equal_point = int(math.floor(self.list_extra[key][0] * 100))
            for i in range(equal_point, ncols):
                    self.FP[i] += 1

    def statistical_result(self):
        self.PPV = numpy.true_divide(self.TP, numpy.add(self.TP, self.FP))
        self.SE = numpy.true_divide(self.TP,
            numpy.add(self.TP, self.FN))

    def calculate_cutoff_frequency(self):
        ones = [ 1 for i in range(100) ]
        dist = numpy.sqrt(numpy.square(numpy.subtract(self.PPV, ones)) +
            numpy.square(numpy.subtract(self.SE, ones)))

        self.cutoff_frequency = min(dist)

    def expansion_list(self):
        self.create_list_intra_extra()
        self.__calculate_the_frequency()
        self.calculate_TPFNFP()
        self.statistical_result()
        self.calculate_cutoff_frequency()

        expansion_list = []
        for key in self.list_extra:
            freq = self.list_extra[key][0]
            if freq >= self.cutoff_frequency:
                gene_1, gene_2 = helper.decode_edge_key(key)
                expansion_list.append([gene_1, gene_2, freq])

        return expansion_list



