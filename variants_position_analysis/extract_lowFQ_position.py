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
        #reader = csv.reader(current_fp, delimiter="\t")
        # file_name = file + "_positions"
        # f1=open(file_name, 'w+')
	print "reading position file"
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #position = row[0]
            #print row
            position_label[row[0]] = row[1:]

##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def generate_heatmap_position():
    print "generate heatmap matrix"
    f1=open("Only_lowFQ_positions_for_three_samples", 'w+')
    #any(item > 0 or item < 1 for item in position_label[value])
    for value in position_label:
	any(item > 0 or item < 1 for item in position_label[value])
        #ref_var = ['1', '1']
	#any(item > 0 or item < 1 for item in position_label[value])
        #for i in position_label[value]:
	    #print i

generate_heatmap_position()
