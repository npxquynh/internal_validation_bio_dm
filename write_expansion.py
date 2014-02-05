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

def write_statistical_result(filepath, result):
    """
    result = [ [], [], ... , []]
    """
    with open(filepath, 'w') as output_file:
        for item in result:
            output_file.write('%s\n' % ' '.join(str(i) for i in item))