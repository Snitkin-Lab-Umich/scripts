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
    #if filess != "/nfs/esnitkin/Ali/Project_MRSA_analysis/2016_03_03_MRSA_position_analysis/MRSA_USA_300/KP_mrsa_56_filter2_final.vcf_no_proximate_snp.vcf_positions":

    f = open(filess, 'r+')
    for line in f:
        line = line.strip()
        position_array.append(line)
    f.close()
#print position_array
position_array_unique = set(position_array)
position_array_sort = sorted(position_array_unique)
#print "\nThe number of unique positions:\n" + str(len(position_array_sort)) + "\n"
unique_position_file = "unique_positions_file"
f=open(unique_position_file, 'w+')
for i in position_array_sort:
    f.write(i + "\n")
f.close()
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
proximate_position_file = "proximate_positions_file"
f1=open(proximate_position_file, 'w+')
for i in proximate_position_array_unique:
    f1.write(i + "\n")
f1.close()
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
unmapped_position_file = "unmapped_positions_file"
f2=open(unmapped_position_file, 'w+')
for i in unmapped_position_array_unique:
    f2.write(i + "\n")
f2.close()

def create_job():
    for i in vcf_filenames:
        job_name = os.path.basename(i)
        job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m abe\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=72:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n/home/apirani/anaconda/bin/python %s/reason_job.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, filter2_only_snp_vcf_dir, filter2_only_snp_vcf_dir, i)
        job_file_name = "%s.pbs" % i
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)




create_job()
