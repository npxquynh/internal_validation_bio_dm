import os
import validate_datastructure as vds
import pdb

def merge_different_pc_run(code, folder):
    """ List all 'expanded network files' from different run of PC-algorithm
    for the same network & same observation

    code: all PC runs for the same network & same observation shares the same
    prefix ~ xml_filename
    """
    list_of_expanded_network_file = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.startswith(code):
                list_of_expanded_network_file.append(os.path.join(root, file))

    return list_of_expanded_network_file

def parse_first_line(line):
    """
    line = '#abc gene_1, gene_2, gene_3\n'
    """
    genes_in_expanded_network = set()

    line = line.strip()
    # remove the #abc at the beginning of line
    line = line[(line.find(' ') + 1):]
    genes_in_expanded_network = set(line.split(', ')[1:])

    return genes_in_expanded_network

def read_expanded_network(list_of_expanded_network_file):
    # TODO: should I use set() here?
    genes_in_tiles = []
    blocks = [] # ~ edges_in_tiles

    for filepath in list_of_expanded_network_file:
        with open(filepath, 'r') as f:

            count = -1

            for line in f:
                if line.startswith('#'):
                    count = count + 1
                    genes_in_tiles.append(parse_first_line(line))
                    blocks.append([])
                else:
                    blocks[count].append(line.strip().split(','))

    return genes_in_tiles, blocks


