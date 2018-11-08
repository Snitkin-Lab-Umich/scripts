__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser(description='Parsing All position with label file and investigating positions to determine the reason why it was filtered out from the final list')
#All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-label_file_names', action='store', dest="label_file_names", help='label_file_names')
parser.add_argument('-positions', action='store', dest="positions", help='position file name')
parser.add_argument('-dir', action='store', dest="dir", help='directory')
args = parser.parse_args()

All_position_file = args.positions

#full filenames array
label_filenames = args.label_file_names
label_file_names = []
with open(label_filenames) as fp:
    for line in fp:
        line = line.strip()
        line = args.dir + line
        label_file_names.append(line)
    fp.close()
#full filenames array

#full filenames array
unique_position_array = []
with open(All_position_file) as fpp:
    for line in fpp:
        line = line.strip()
        unique_position_array.append(line)
    fpp.close()
#full filenames array





num_cores = multiprocessing.cpu_count()
for i in label_file_names:
#def reason(i):
    label_files = i
    out_file_name = str(label_files) + "_updated"
    f1=open(out_file_name, 'w+')
    array_name = os.path.basename(label_files)
    current_unmapped_file = label_files.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "unmapped.bed_positions")
    current_proximate_file = label_files.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf_positions_array")
    current_variant_position_file = label_files.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf_positions")
    variant_position_array = "variant_" + str(array_name)
    variant_position_array = []
    with open(current_variant_position_file, 'rU') as fpq:
        for line in fpq:
            line = line.strip()
            variant_position_array.append(line)
    fpq.close()
    unmapped_array = "unmapped_" + str(array_name)
    unmapped_array = []
    with open(current_unmapped_file, 'rU') as fp1:
        for line in fp1:
            line = line.strip()
            unmapped_array.append(line)
    fp1.close()
    proximate_array = "proximate_" + str(array_name)
    proximate_array = []
    with open(current_proximate_file, 'rU') as fp2:
        for liness in fp2:
            liness = liness.strip()
            proximate_array.append(liness)
    fp2.close()
    label_array = "label_array_" + str(array_name)
    label_array = []
    with open(label_files, 'rU') as fp3:
        for l in fp3:
            l = l.strip()
            label_array.append(l)
    fp3.close()
    for i in range(len(unique_position_array)):
	print label_files
	print label_array[i]
        if label_array[i] == "reference_unmapped_position":
            if unique_position_array[i] not in unmapped_array:
                print_string = "reference_allele"
            else:
                print_string = label_array[i]
        elif "_proximate_SNP" in label_array[i]:
            if unique_position_array[i] not in proximate_array:
                print_string = label_array[i].replace("_proximate_SNP", "")
            else:
                print_string = label_array[i]
        elif label_array[i] == "VARIANT":
            if unique_position_array[i] not in variant_position_array:
                if unique_position_array[i] in unmapped_array:
                    print_string = "reference_unmapped_position"
                else:
                    print_string = "reference_allele"
            else:
                print_string = label_array[i]
        else:
            print_string = label_array[i]
        #print_final = unique_position_array[i] + "\t" + print_string + "\n"
        print_final = print_string + "\n"
        f1.write(print_final)



#results = Parallel(n_jobs=num_cores)(delayed(reason)(i) for i in label_file_names)
