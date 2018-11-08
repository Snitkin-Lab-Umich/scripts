__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict


parser = argparse.ArgumentParser(description='Parsing All position with label file and investigating positions to determine the reason why it was filtered out from the final list')
#All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-positions_file_dir', action='store', dest="positions_file_dir", help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-positions_filenames', action='store', dest="positions_filenames", help='Names of filter2 only SNP vcf files with name per line.')
parser.add_argument('-unique_positions', action='store', dest="unique_positions", help='Names of unique_positions_file')

args = parser.parse_args()

All_position_file = args.positions_filenames
position_label = OrderedDict()


filename = args.unique_positions
unique_position = []
with open(filename) as fp:
    for line in fp:
        line = line.strip()
        unique_position.append(line)



with open(All_position_file, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #position = row[0]
            #print row
            if row[0] in unique_position:
                print row
            #position_label[row[0]] = row[1:]

##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def generate_heatmap_position():
    f1=open("All_label_raw_only_closely_related_positions", 'w+')
    for key in position_label:
        if key in unique_position:
            strr = key + "\n" + value + "\n"
            f1.write(strr)

#generate_heatmap_position()

