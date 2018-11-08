__author__ = 'alipirani'
# This script assigns phylogroup based on the combination of presence and absence of certain primer target sequences such as: arpA  chuA  tspe4  yjaA C_specific  E_specific
# The primer sequences for these targets can be found in the following paper: http://onlinelibrary.wiley.com/doi/10.1111/1758-2229.12019/abstract
# Various combinations of presence/absence of these genes used for phylogroup assignment can also be found in these paper.
# The script take the name of genome fasta file(one genome per line) and the directory where these fasta file are located.
# The script also requires path to pre-computed Blast db of each primer target sequences.
# These genes are assumed to be named in the following fashion(arpA.fasta  chuA.fasta  C_specific.fasta  E_specific.fasta  tspe4.fasta  yjaA.fasta)
# Make sure the Blast db is also named similarly to fasta sequence file names


import argparse
import re
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict


parser = argparse.ArgumentParser(description='E. coli Phylogroup Determination')
parser.add_argument('-db', action='store', dest="db_directory",
                    help='Path to Database directory')
parser.add_argument('-genomes', action='store', dest="genome_name",
                    help='File containing genome fasta file name(One name per line). Assuming genome fasta file contains nucleotide sequences.')
parser.add_argument('-dir', action='store', dest="directory",
                    help='Directory of genome fasta files')
args = parser.parse_args()




# Read Genome fasta filename file
genome_names_array = []
with open(args.genome_name) as fp:
    for line in fp:
        line = line.strip()
        line = args.directory + "/" + line
        genome_names_array.append(line)
    fp.close()




# Check if the database files exists
primer_target_seq = ['arpA.fasta', 'chuA.fasta', 'yjaA.fasta', 'tspe4.fasta','C_specific.fasta', 'E_specific.fasta']
print "\nAssuming the primer target sequence fasta files are:\narpA.fasta\nchuA.fasta\nC_specific.fasta\nE_specific.fasta\ntspe4.fasta\nyjaA.fasta\n"
for i in primer_target_seq:
    full_path = args.db_directory + "/" + i
    full_path_suffix = args.db_directory + "/" + i + ".nhr"
    if not os.path.isfile(full_path):
        file_basename = os.path.basename(full_path)
        print "The primer target sequence file " + file_basename + " does not exists.\n Please create Blast db for this fasta sequence and place the db files along with fasta sequence in -db directory\n"
        exit()
    if not os.path.isfile(full_path_suffix):
        file_basename = os.path.basename(full_path_suffix)
        print "The blast database files does not exists.\n Please create Blast db for your primer target sequence fasta files and place the it in -db directory\n"
        exit()
    if os.path.isfile(full_path) and os.path.isfile(full_path_suffix):
        file_basename = os.path.basename(full_path)
        print "The primer target fasta file " + file_basename + " and associated Blast db files exists.\n"

print "Proceeding...\n"

#Blast
temp_combination_file = "combination_file_tmp"
f1=open(temp_combination_file, 'w+')
header = "genome"
for gene in primer_target_seq:
    header = header + "\t" + gene
f1.write(header + '\n')
for i in genome_names_array:
    print_string = str(i)
    if not os.path.isfile(i):
        file_basename = os.path.basename(i)
        print "The genome fasta file " + file_basename + " does not exists.\n"
    else:
        for gene in primer_target_seq:
            output = i + "_" + gene + "_tmp.xml"
            db = args.db_directory + "/" + gene
            blastn_cline = NcbiblastnCommandline(query=i, db=db, outfmt=6)
            #print(str(blastn_cline))
            #os.system(str(blastn_cline))
            stdout, stderr = blastn_cline()
            if stdout == "":
                print_string = print_string + "\t" + "-"
            else:
                print_string = print_string + "\t" + "+"
    f1.write(print_string + "\n")



def assign_phylogroup():
    temp_combination_file = "combination_file_tmp"
    with open(temp_combination_file, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            print row[1]

assign_phylogroup()

























