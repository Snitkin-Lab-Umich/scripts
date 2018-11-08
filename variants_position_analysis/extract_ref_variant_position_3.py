from __future__ import division
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
        print "reading position file"
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #position = row[0]
            #print row
            position_label[row[0]] = row[1:]

##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def generate_FQ_25perc_heatmap_position():
    print "generate heatmap matrix"
    f1=open("All_FQ_less_than_25_perc", 'w+')
    f2=open("All_FQ_less_than_25_perc_2", 'w+')
    f3=open("All_FQ_more_than_25_perc_2", 'w+')
    f4=open("All_FQ_less_than_10_perc_2", 'w+')
    for value in position_label:
        lll = ['0', '2', '3', '4', '5', '6', '7']
        other_filter_positions = ['0', '2', '3', '4', '6', '7']
        unmapped = ['0', '0']
        ref_var = ['1', '1']
        fq = ['5', '5']
        if set(ref_var) & set(position_label[value]):
            if set(lll) & set(position_label[value]):
                if set(unmapped) & set(position_label[value]):
                    sfsf = "bakwaas"
                elif set(fq) & set(position_label[value]):
                    occurence = position_label[value].count('5')
                    perc = float(occurence / len(position_label[value])) * 100
                    if perc > 25.00:
                        print perc
                        STRR1 = value + "\t" + str(position_label[value]) + "\n"
                        f3.write(STRR1)
                    #print position_label[value]
                    #STRR2 = value + "\t" + str(position_label[value]) + "\n"
                    #f1.write(STRR2)
		# STRR2 = value + "\t" + str(position_label[value]) + "\n"
		# f3.write(STRR2)

generate_FQ_25perc_heatmap_position()




