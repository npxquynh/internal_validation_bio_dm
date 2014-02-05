import math
import numpy
import operator

# TODO: remove pdb
import pdb

import helper

class InternalValidation:
    def __init__(self, genes_in_tile, blocks, genes_in_lgn, edges_in_lgn):
        self.genes_in_tile = genes_in_tile
        self.blocks = blocks
        self.genes_in_lgn = genes_in_lgn
        self.edges_in_lgn = edges_in_lgn
        self.list_intra = dict()
        self.list_extra = dict()
        self.__generate_zero_matrices()

    def __find_extra_genes_connected_with_LGN(self, block):
        # connected_nodes[i] = set() = nodes connected with LGN[i]
        connected_nodes = []
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
            # exclude those connection between LGN gene & another LGN gene
            nodes = set(nodes) - set(self.genes_in_lgn)
            connected_nodes.append(nodes)
        return connected_nodes

    # Input:
    # LGN = [1, 2, 3]
    #   [(1, 10), (1, 2, 3), (2)]
    #   0
    # Output: (1, 2, 3)
    def __find_connected_nodes_in_SubLGN_1(self, connected_nodes, index):
        # connected_nodes = list[set()]
        result = set()
        for (idx, nodes) in enumerate(connected_nodes):
            if idx != index:
                result = result.union(nodes)

        return result

    def __refine_list_intra(self):
        number_of_blocks = float(len(self.blocks))
        # TODO: handle division by zero
        for key in self.list_intra:
            self.list_intra[key] /= number_of_blocks

    # Input:
    # LGN = [1, 2, 3]
    #   [(1, 10), (1, 2, 3), (2)]
    #   0
    # Output: [1, 2, 3, 2]
    def __find_connected_nodes_in_SubLGN_2(self, connected_nodes, index):
        # connected_nodes = list[set()]
        result = list()
        for (idx, nodes) in enumerate(connected_nodes):
            if idx != index:
                for gene in nodes:
                    result.append(gene)

        return result

    def __connection_in_extra_genes_step_0(self, block):
        """
        calculate c[0]: whether or not this 'extra gene' appear in a specific block
        """
        control = dict()
        max_length = len(self.genes_in_lgn) + 2

        genes_in_block = helper.genes_from_edges(block)
        extra_genes = set(genes_in_block) - set(self.genes_in_lgn)

        """ since we do not care how many time extra genes appear in block
        we ONLY pay attention to whether it appears or not
        """
        for gene in extra_genes:
            control[gene] = [0] * max_length
            control[gene][0] = 1

            # number of times it's connected with the LGN - c[1]
            # number of times it's connected with the subLGN - c[2]
            # will be calculated later

        return control

    # TODO: connected_nodes = list[set()] so it might remove duplicated node that
    # are connected to the same gene in LGN.
    # But another thinking suggests to me that this case will never happen,
    # since if it happens, we will have duplicated connection
    # Like At1 <-> At2, and At1 <-> At2
    # So I wonder whether my thinking is correct?
    def __connection_in_extra_genes_step_1(self, control, connected_nodes):
        """
        calculate c[1]: whether 'extra-gene' is directly connected to the LGN

        connected_nodes = list[set()] ==> remember to 'flatten' them first
        """
        # use set() since we do not pay attention to how many times
        # extra genes connect with LGN
        nodes = set(helper.list_of_sets_to_list(connected_nodes))
        for gene in nodes:
            if gene in control:
                control[gene][1] += 1
            else:
                # In the same block, this thing should never happen!!!
                print "Problem!!! internal validation 1"

        return control

    def __connection_in_extra_genes_step_2(self, control, connected_nodes):
        """
        calculate c[2] ... c[n+1]: number of times it's connected to SLGN(0), ...,
        whether 'extra-gene' is connected to SLGN(i), ... ,
        whether 'extra-gene' is connected to SLGN(n-1))

        connected_nodes_in_subLGN = list[]
        """
        number_of_LGN_genes = len(self.genes_in_lgn)
        list_of_connected_nodes = helper.list_of_sets_to_list(connected_nodes)
        calculated_node = set()
        set_of_connected_nodes = set(list_of_connected_nodes)
        for (index, set_of_genes) in enumerate(connected_nodes):
            for gene in set_of_genes:
                if gene in control:
                    for k in range(0, number_of_LGN_genes):
                        control[gene][k + 2] = 1
                    if gene not in calculated_node:
                        control[gene][index + 2] = 0
                else:
                    # In the same block, this thing should never happen!!!
                    print "Problem!!! internal validation 2"
        return control

    def __merge_list_extra(self, control):
        for key in control:
            if key in self.list_extra:
                self.list_extra[key] = map(operator.add, self.list_extra[key], control[key])
            else:
                self.list_extra[key] = control[key]

    def __refine_list_extra(self):
        max_length = len(self.genes_in_lgn) + 2

        keys_to_be_deleted = []

        # TODO: refactor constant number of LGN genes
        for key in self.list_extra:
            if self.list_extra[key][1] == 0:
                keys_to_be_deleted.append(key)
            else:
                for k in range(1, max_length):
                    self.list_extra[key][k] /= float(self.list_extra[key][0])

        for key in keys_to_be_deleted:
            self.list_extra.__delitem__(key)

    def create_list_intra_extra(self):
        print 'create_list_intra_extra'

        for gene in self.genes_in_lgn:
            self.list_intra[gene] = 0

        for block in self.blocks:
            connected_nodes = self.__find_extra_genes_connected_with_LGN(block)

            """ calculate list_intra """
            for (index, gene) in enumerate(self.genes_in_lgn):
                connected_nodes_in_subLGN = self.__find_connected_nodes_in_SubLGN_1(connected_nodes, index)
                if connected_nodes_in_subLGN.__contains__(gene):
                    self.list_intra[gene] += 1

            """ calculate list_extra, we just simply do it for every genes in
            the block, and then we can remove those genes in LGN
            """
            control = self.__connection_in_extra_genes_step_0(block)
            control = self.__connection_in_extra_genes_step_1(control, connected_nodes)
            control = self.__connection_in_extra_genes_step_2(control, connected_nodes)

            self.__merge_list_extra(control)

        self.__refine_list_intra()
        self.__refine_list_extra()

    def __generate_zero_matrices(self):
        f = [x / 100.0 for x in range(100)]
        ncols = 100
        nrows = len(self.genes_in_lgn)
        self.TP = [ [0] * ncols for i in range(nrows) ]
        self.FN = [ [1] * ncols for i in range(nrows) ]
        self.FP = [ [0] * ncols for i in range(nrows) ]

    def largest_frequency_list_extra(self, index_of_subLGN):
        max_freq = 0
        for key in self.list_extra:
            if self.list_extra[key][index_of_subLGN + 2] > max_freq:
                # print self.list_extra[key][index_of_subLGN + 2]
                max_freq = self.list_extra[key][index_of_subLGN + 2]

        return max_freq

    def calculate_TPFNFP(self):
        ncols = 100
        nrows = len(self.genes_in_lgn)

        # True Positive & False Negative
        for (i, key) in enumerate(self.list_intra):
            # finding an index in which a remaining part
            # of a row is filled with 1
            # Ex:
            # 0 0 0 1
            # 0 1 1 1
            # 0 0 0 1
            equal_point = int(math.floor(self.list_intra[key] * 100))
            for j in range(equal_point, ncols):
                self.TP[i][j] = 1
                self.FN[i][j] = 0

        # False Positive
        f = [x / 100.0 for x in range(100)]
        for key in self.list_extra:
            for i in range(0, nrows):
                for j in range(ncols):
                    try:
                        if self.list_extra[key][i+2] >= f[j]:
                            self.FP[i][j] += 1
                    except KeyError:
                        pdb.set_trace()

    def statistical_result(self):
        TP = numpy.array(self.TP, dtype = numpy.float)
        FP = numpy.array(self.FP, dtype = numpy.float)
        FN = numpy.array(self.FN, dtype = numpy.float)

        self.PPV = TP / (TP + FP)
        self.SE = TP / (TP + FN)

    def calculate_cutoff_frequency(self):
        Y = numpy.mean(self.PPV, axis=0)
        X = numpy.mean(self.SE, axis=0)
        ones = [ 1 for i in range(100) ]

        dist = numpy.sqrt(numpy.square(numpy.subtract(X, ones)) +
            numpy.square(numpy.subtract(Y, ones)))

        self.cutoff_frequency = min(dist)
        print "Cutoff Freq 1: %f" % self.cutoff_frequency

    def expansion_list(self):
        self.create_list_intra_extra()
        self.calculate_TPFNFP()
        self.statistical_result()
        self.calculate_cutoff_frequency()

        expansion_list = []
        for key in self.list_extra:
            freq = self.list_extra[key][1]
            if freq >= self.cutoff_frequency:
                expansion_list.append([key, freq])

        return expansion_list



