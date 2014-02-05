import os
from itertools import chain

def generate_filepath(dirname, filename):
    return os.path.join(dirname, filename)

def genes_from_edges(edges):
    return set(chain.from_iterable(edges))

def list_of_sets_to_list(l):
    result = []
    for s in l:
        for element in s:
            result.append(element)
    return result

def count_element_in_list_of_sets(list_of_sets):
    elements_list = list_of_sets_to_list(list_of_sets)
    elements_set = set(elements_list)
    result = dict()
    for element in elements_set:
        result[e] = elements_list.count(element)
    return result

def add_two_dim_array(array_1, array_2):
    nrows = len(array_1)
    ncols = len(array_1[0])

    for row in range(nrows):
        for col in range(ncols):
            array_1[row][col] += array_2[row][col]
    return array_1

def divide_two_dim_array(array_1, array_2):
    nrows = len(array_1)
    ncols = len(array_1[0])

    for row in range(nrows):
        for col in range(ncols):
            array_1[row][col] /= float(array_2[row][col])
    return array_1

def mean_of_columns(matrix):
    nrows = len(matrix)
    ncols = len(matrix[0])

def swap(a, b):
    if a > b:
        return b,a
    return a,b

def encode_edge_key(a, b):
    """
    [gene_1, gene_2] ==> 'gene_1*gene_2'
    """
    return '%s*%s' % (a, b)

def decode_edge_key(code):
    return code.split('*')

