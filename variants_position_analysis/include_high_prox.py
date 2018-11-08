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


with open(All_position_file, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            for row in csv_reader:
                position_label[row[0]] = row[1:]
# unique_position = []
# with open(args.unique_positions) as fp:
#         for line in fp:
#             line = line.strip()
#             unique_position.append(line)
##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def include_high_prox_positions():

    f1=open("High_FQ_proximate_positions", 'w+')
    high_prox_array = []
    #print "here"
    high_prox = ['7','7']
    lll = ['0', '2', '3', '4', '5', '6']
    for value in position_label:


        prox_but_filtered_out = 0
        #print position_label[value]
        if set(high_prox) & set(position_label[value]):
            # if set(lll) not in position_label[value]:
            if set(lll) & set(position_label[value]):
                print "bakwas:%s" % value
            else:
                high_prox_array.append(value)
    for i in high_prox_array:
        strr = i + "\n"
	f1.write(strr)


include_high_prox_positions()
