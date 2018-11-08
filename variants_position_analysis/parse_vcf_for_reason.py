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

parser = argparse.ArgumentParser(
    description='Parsing filtered VCF files and investigating Variants to determine the reason why it was filtered out from the final list')
# All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-filter2_only_snp_vcf_filenames', action='store', dest="filter2_only_snp_vcf_filenames",
                    help='Names of filter2 only SNP vcf files with name per line.')
parser.add_argument('-full', action='store', dest="full",
                    help='Names of filter2 only SNP vcf files with name per line.')
parser.add_argument('-proximate_position_file', action='store', dest="proximate_position_file",
                    help='Names of proximate position files with name per line.')
parser.add_argument('-unmapped_position_file', action='store', dest="unmapped_position_file",
                    help='Names of unmapped position files with name per line.')
args = parser.parse_args()

filter2_only_snp_vcf_dir = args.filter2_only_snp_vcf_dir

############################################prepare file array###################################################

#full filenames array
full_filenames = args.full
full_names_array = []
with open(full_filenames) as fp:
    for line in fp:
        line = line.strip()
        line = args.filter2_only_snp_vcf_dir + line
        full_names_array.append(line)
    fp.close()
#full filenames array


#proximate file array
proximate_position_file = args.proximate_position_file
proximate_position_file_array = []
with open(proximate_position_file) as fp:
    for line in fp:
        line = line.strip()
        line = args.filter2_only_snp_vcf_dir + line
        proximate_position_file_array.append(line)
    fp.close()
#proximate file array


#unmapped position file array
unmapped_position_file = args.unmapped_position_file
unmapped_position_file_array = []
with open(unmapped_position_file) as fp:
    for line in fp:
        line = line.strip()
        line = args.filter2_only_snp_vcf_dir + line
        unmapped_position_file_array.append(line)
    fp.close()
#unmapped position file array

#filter2 vcf filenames array
filter2_only_snp_vcf_filenames = args.filter2_only_snp_vcf_filenames
vcf_filenames = []
with open(filter2_only_snp_vcf_filenames) as fp:
    for line in fp:
        line = line.strip()
        line = args.filter2_only_snp_vcf_dir + line
        vcf_filenames.append(line)
    fp.close()
#filter2 vcf filenames array

#print str(full_names_array) + "\n"
#print str(proximate_position_file_array) + "\n"
#print str(unmapped_position_file_array) + "\n"
#print str(vcf_filenames) + "\n"
############################################End: prepare file array###################################################

############################################Create _positions files###################################################
#From this filtered vcf files get all the positions and write to files with entension _position
#print "creating _positions file step\n"
filter2_only_snp_position_files_array = args.filter2_only_snp_vcf_dir + "position_filenames"
for file in vcf_filenames:
    file_base_name = os.path.basename(file)
    with open(filter2_only_snp_position_files_array, 'a') as fc:
        #print "Processing file: %s\n" % file
        with open(file, 'rU') as csv_file:
            file_name = str(file) + "_positions"
            addpositionfilenametoarray = file_name + "\n"
            #filter2_only_snp_position_files_array.append(addpositionfilenametoarray)
            fc.write(addpositionfilenametoarray)
            f1 = open(file_name, 'w+')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                position = row[0]
                #print row
                if position.startswith('gi'):
                    p_string = row[1] + "\n"
                    f1.write(p_string)
            f1.close()
        csv_file.close()
#print "End of creating _positions file step\n"
#print "The position files are:\n" + str(filter2_only_snp_position_files_array) + "\n"
############################################Create _positions files###################################################


####################Create position array containing unique positiones from positions file############################
#Get unique positions from all the _position extension files and save to an array position_array
position_filenames = []
position_file = filter2_only_snp_position_files_array
with open(position_file) as fp:
    for line in fp:
        line = line.strip()
        #print line
        position_filenames.append(line)
    fp.close()

position_array = []
for filess in position_filenames:
    #print filess
    f = open(filess, 'r+')
    for line in f:
        line = line.strip()
        position_array.append(line)
    f.close()

#print position_array
position_array_unique = set(position_array)
position_array_sort = sorted(position_array_unique)
#print "\nThe number of unique positions:\n" + str(len(position_array_sort)) + "\n"
for i in position_array_sort:
    print i
#del position_filenames[:]
#del position_array_unique[:]
####################END: Create position array containing unique positiones from positions file#######################



####################Creating Proximate Position Array and unmapped position array#####################################
#Create a proximate_position_array containing 10 bp snps positions that were filtered out
proximate_position_array = []
for file in proximate_position_file_array:
    #print "Processing proximate position array for file: %s" % file
    with open(file, 'rU') as fp1:
        #reader = csv.reader(current_fp, delimiter="\t")
        for line in fp1:
            line = line.strip()
            proximate_position_array.append(line)
    fp1.close()
    #print "END: Processing proximate position array for file: %s" % file
proximate_position_array_unique = set(proximate_position_array)
#print "\nThe number of proximate positions:\n" + str(len(proximate_position_array_unique)) + "\n"
#del proximate_position_array[:]



#Create a unmapped position array containing unmapped positions
unmapped_position_array = []
for file in unmapped_position_file_array:
    #print "Processing unmapped position array for file: %s" % file
    with open(file, 'rU') as fp1:
        #reader = csv.reader(current_fp, delimiter="\t")
        for line in fp1:
            line = line.strip()
            unmapped_position_array.append(line)
    fp1.close()
    #print "END: Processing unmapped position array for file: %s" % file

unmapped_position_array_unique = set(unmapped_position_array)
#print "\nThe number of unmapped positions:\n" + str(len(unmapped_position_array_unique)) + "\n"
#del unmapped_position_array[:]
########################################################################################################################

#For each of these position:
#check if it has been filtered out from the raw mpileup vcf files
#if not its a TRUE variant which satisfies all the filter parameters
#else:
#check if its a reference allele
#if not check the reason why was it filtered out: HIGH FQ but low DP/low QUAL or LOW FQ/low DP/lowQUAL
position_label_file_array = []
num_cores = multiprocessing.cpu_count()
#print "The number of cores: %s\n" % num_cores
#def reason(i):
for i in full_names_array:
    file = i
    print "Processing %s" % file
    # with open(file, 'rU') as fp2:
    out_file_name = file + "_positions_label"
    position_label_file_array.append(out_file_name)
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
                print out2
                line_string_array = out2.split('\t')
                # if line_string_array[1] in proximate_position_array_unique:
                #     pst = "_proximate_SNP"
                if line_string_array[1] in proximate_position_array_unique:
                    pst = "_proximate_SNP"
                else:
                    pst = ""
                format_string_array = line_string_array[7].split(';')
                # DP = format_string_array[0]
                # FQ = format_string_array[5]
                # MQ = format_string_array[8]
                # for i in format_string_array:
                #     if i.startswith('FQ='):
                #         m = re.match("FQ=(\W*\d+)", i)
                #         FQ_value = m.group(1)
                #         if FQ_value.startswith('-'):
                #             st = "high_qualty_filtered_out\n"
                #
                #             f1.write(st)
                #         else:
                #             st = "Uncertain_low_quality\n"
                #             f1.write(st)
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

#results = Parallel(n_jobs=4)(delayed(reason)(i) for i in full_names_array)
#results = Parallel(n_jobs=num_cores)(delayed(paircomp)(i) for i in tree_filenames)

# All_file_strring = ' '.join(position_label_file_array)
# paste_command = "paste " + All_file_strring + " > All_position_label_final_file"
# os.system(paste_command)
#Get the positions from position_array in order to make a matrix for heatmap
#pending
#I am adding this position column manually
# position_column = ""
# for i in position_array_sort:
#     position_column = position_column + i + "\n"
# position_column_file = args.filter2_only_snp_vcf_dir + "position_column"
# f1=open(position_column_file, 'w+')
#f1.write(position_column)

