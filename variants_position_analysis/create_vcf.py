__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict


parser = argparse.ArgumentParser(description='Create VCF from Position list: Remove all those positions from vcf file that are not in the list. This list is obtained after removing those positions that contain atlease one filtered out variant because of not passing the final filter2')
#All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-filter2_file_dir', action='store', dest="filter2_file_dir", help='Filter2 file directory')
parser.add_argument('-filter2_filenames', action='store', dest="filter2_filenames", help='Filter2 filenames')
parser.add_argument('-position_list', action='store', dest="position_list", help='List of positions: This positions are from heatmap i.e: proximate, unmapped, filtered-out')
args = parser.parse_args()

filter2_only_snp_vcf_filenames = args.filter2_filenames

vcf_filenames = []
with open(filter2_only_snp_vcf_filenames) as fp:
    for line in fp:
        line = line.strip()
        line = args.filter2_file_dir + line
        #tree_filenames = line
        vcf_filenames.append(line)


#This positions are from heatmap i.e: proximate, unmapped, filtered-out
#This positions are from heatmap i.e: proximate, unmapped, filtered-out
#This positions are from heatmap i.e: proximate, unmapped, filtered-out
position_array = []
with open(args.position_list) as fp:
    for line in fp:
        line = line.strip()
        position_array.append(line)


filtered_out_vcf_files = []

for i in vcf_filenames:
    #print_array = i
    print_array =[]
    with open(i) as file_open:
         for line in file_open:
            line = line.strip()
            if line.startswith("#"):
                print_array.append(line)
            else:
                #line.split('	')
                split_array = re.split(r'\t+', line)
                if split_array[1] in position_array:
                    print_array.append(line)
    file_name = i + "_no_unmapped_no_proximate_no_filtered_out_position.vcf"
    filtered_out_vcf_files.append(file_name)
    f1 = open(file_name, 'w+')
    for ios in print_array:
        print_string = str(ios) + "\n"
        f1.write(print_string)


