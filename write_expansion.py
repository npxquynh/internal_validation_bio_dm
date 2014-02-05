import os
import csv

def write_expansion_list(filepath, file_content):
    """
    file_content = [
        [gene_name, frequency_1],
        [gene_name_2, frequency_2]
    ]
    """
    with open(filepath, 'w') as output_file:
        w = csv.writer(output_file, delimiter=',')
        for line in file_content:
            w.writerow(line)
