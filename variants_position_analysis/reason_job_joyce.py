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
from pyfasta import Fasta

parser = argparse.ArgumentParser(description='Creating Label files individual jobs')
# All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-filter2_only_snp_vcf_file', action='store', dest="filter2_only_snp_vcf_file",
                    help='Names of filter2 only SNP vcf file')
args = parser.parse_args()
dir = args.filter2_only_snp_vcf_dir
unique_positions_file = dir + "unique_positions_file"

############################################ Generate unique positions array###########################################
position_array_sort = []
f = open(unique_positions_file, 'r+')
for line in f:
    line = line.strip()
    position_array_sort.append(line)
f.close()
############################################ End: Generate unique positions array######################################

file = args.filter2_only_snp_vcf_file
print "Processing %s" % file
out_file_name = file + "_positions_label"


#initialize array
array_name = os.path.basename(out_file_name)

reference = "/home/apirani/bin/reference/Ecoli_CD306/Ecoli_CD306.fasta"

################# Generate label files and check why the positions were filtered out from the final vcf file #########
f1=open(out_file_name, 'w+')
for j in position_array_sort:
    # Changing Grep command 25 jan
    #cmd = "grep -wP \'\s+" + j + "\s+\' " + file
    cmd = "grep -v \'^#\' %s | awk -F\'\t\' \'{print $2}\' | grep -w \'%s\'" % (file, j)
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    if not out:
        f = Fasta(reference)
        if len(f.keys()) == 1:
            ref_id = str(f.keys())

        fasta_string = ""
        extract_base = "tr -d \'\\n\' < %s | cut -b%s" % (reference, j)
        #print extract_base
        # proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
        # (out, err) = proc.communicate()
        # out = out.strip()
        # fasta_string = fasta_string + out
        # if not out:
        #     print "Error extracting reference allele"
        #out = str(f.sequence({'chr': str(f.keys()[0]), 'start': int(lines), 'stop': int(lines)}))
        fasta_string = fasta_string + str(f.sequence({'chr': str(f.keys()[0]), 'start': int(j), 'stop': int(j)}))

        pattern = re.compile(r'\s+')
        fasta_string = re.sub(pattern, '', fasta_string)

        st = fasta_string + fasta_string + "\n"
        f1.write(st)
    else:
        cmd2 = "grep -P \'\s+" + j + "\s+\' " + args.filter2_only_snp_vcf_file
        #cmd2 =  "grep -v \'^#\' %s | awk -F\'\t\' \'{print $2}\' | grep -w \'%s\'" % (final_file, j)
        proc = subprocess.Popen([cmd2], stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc.communicate()
        line_string_array = out2.split('\t')
        print line_string_array
        ref_allele = line_string_array[3]
        alt_allel = line_string_array[4]
        st = ref_allele + alt_allel + "\n"
        f1.write(st)
f1.close()
############# End: Generate label files and check why the positions were filtered out from the final vcf file #########
