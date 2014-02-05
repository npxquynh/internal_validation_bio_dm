def genes_in_tile(set_of_genes):
    for gene in set_of_genes:
        # no trailing whitespace and
        if gene != gene.strip():
            return False
        if gene.find('#') != -1:
            return False

    return True

def list_intra(list_intra):
    for key in list_intra:
        if list_intra[key] > 1:
            print "Error with frequency in list intra"