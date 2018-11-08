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
parser = argparse.ArgumentParser(description='Creating Label files individual jobs')
# All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-filter2_only_snp_vcf_file', action='store', dest="filter2_only_snp_vcf_file",
                    help='Names of filter2 only SNP vcf file')
args = parser.parse_args()


dir = args.filter2_only_snp_vcf_dir

proximate_positions_file = dir + "proximate_positions_file"
unmapped_positions_file = dir + "unmapped_positions_file"
unique_positions_file = dir + "unique_positions_file"

position_array_sort = []
f = open(unique_positions_file, 'r+')
for line in f:
    line = line.strip()
    position_array_sort.append(line)
f.close()


unmapped_position_array_unique = []
f = open(unmapped_positions_file, 'r+')
for line in f:
    line = line.strip()
    unmapped_position_array_unique.append(line)
f.close()

proximate_position_array_unique = []
f = open(proximate_positions_file, 'r+')
for line in f:
    line = line.strip()
    proximate_position_array_unique.append(line)
f.close()




file = args.filter2_only_snp_vcf_file
print "Processing %s" % file
# with open(file, 'rU') as fp2:
out_file_name = file + "_positions_label"
f1=open(out_file_name, 'w+')
for j in position_array_sort:
    cmd = "grep \'" + j + "\' " + file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    if not out:
        # filter2_final.vcf_no_proximate_snp.vcf
        # aln_mpileup_raw.vcf_5bp_indel_removed.vcf
        mpileup_file = file.replace("filter2_final.vcf_no_proximate_snp.vcf", "aln_mpileup_raw.vcf_5bp_indel_removed.vcf")
        final_file = mpileup_file
        cmd2 = "grep \'" + j + "\' " + final_file
        proc = subprocess.Popen([cmd2], stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc.communicate()
        if not out2:
            if j in unmapped_position_array_unique:
                st = "reference_unmapped_position\n"
                f1.write(st)
            else:
                st = "reference_allele\n"
                f1.write(st)
        else:
            #print out2
            line_string_array = out2.split('\t')
            # if line_string_array[1] in proximate_position_array_unique:
            #     pst = "_proximate_SNP"
            if line_string_array[1] in proximate_position_array_unique:
                pst = "_proximate_SNP"
            else:
                pst = ""
            format_string_array = line_string_array[7].split(';')
            for i in format_string_array:
                if i.startswith('FQ='):
                    m = re.match("FQ=(\W*\d+)", i)
                    FQ_value = m.group(1)
                    if FQ_value.startswith('-'):
                        st = "HighFQ"
                        #st = "high_qualty_filtered_out\n"
                        qual_string = line_string_array[5]
                        if float(qual_string) < 100.00:
                            st = st + "_QUAL"
                        for i in format_string_array:
                            if i.startswith('DP='):
                                m = re.match("DP=(\W*\d+)", i)
                                DP_value = m.group(1)
                                if int(DP_value) <= 15:
                                    st = st + "_DP"

                    else:
                        #st = "Uncertain_low_quality\n"
                        st = "LowFQ"
                        qual_string = line_string_array[5]
                        if float(qual_string) < 100.00:
                            st = st + "_QUAL"
                        for i in format_string_array:
                            if i.startswith('DP='):
                                m = re.match("DP=(\W*\d+)", i)
                                DP_value = m.group(1)
                                if int(DP_value) < 15:
                                    st = st + "_DP"
            #st = st + "\n"
            st = st + pst + "\n"
            f1.write(st)
    else:
        st = "VARIANT\n"
        f1.write(st)
f1.close()

