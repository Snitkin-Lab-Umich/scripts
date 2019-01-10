__author__ = 'alipirani'
import os
import argparse
import re
import subprocess
import statistics
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser(description='Generate Plasmid Assembly Report')
parser.add_argument('-filename', action='store', dest="filename", help='filename of assembly in fasta format')
parser.add_argument('-AMR_gene', action='store', dest="AMR_gene", help='AMR gene blast database')
args = parser.parse_args()
sample_name = args.filename

# Run Bioawk on the Assembly fasta file to get contig names and its corresponsing length
command = "~/bin/bioawk-master/bioawk -c fastx '{ print $name, length($seq) }' < %s" % sample_name
output = subprocess.check_output(command, shell=True)

# Save 'Main' plasmid names(component_*) into uniq_plasmid_names array. Refer Spades Plasmid Assemblt output page for info.
uniq_plasmid_names = []
# Parse Bioawk results output and get information such as unique plasmid names and number of plasmids found by Spades assembler
for line in output.splitlines():
    linesplit = line.split('\t')
    plasmid_name = linesplit[0]
    plasmid_name_split = plasmid_name.split('_')
    plasmid_name_uniq = str(plasmid_name_split[6]) + "_" + str(plasmid_name_split[7])
    #plasmid_name_uniq = str(plasmid_name_split[2]) + "_" + str(plasmid_name_split[3])
    if plasmid_name_uniq not in uniq_plasmid_names:
        uniq_plasmid_names.append(plasmid_name_uniq)
plasmid_number = len(uniq_plasmid_names)


# Get number of total contigs in spades assembly output
command_2 = "grep \'>\' %s | wc -l" % sample_name
no_of_contigs = subprocess.check_output(command_2, shell=True)

# Initiate dictionaries for various information
plasmid_length = OrderedDict()  # Save length of all the contigs corresponding to respective unique plasmid
median_coverage = OrderedDict() # Calculate and Save median of all the contigs corresponding to respective unique plasmid
longest_contig = defaultdict(list) # extract the length of longest contig for respective unique plasmid

# Again parse the Bioawk output for extracting information about individual plasmids
for uniq_plasmid in uniq_plasmid_names:
    uniq_plasmid_length = 0
    uniq_plasmid_coverage = 0
    median_array = []
    for line in output.splitlines():
        linesplit = line.split('\t')
        plasmid_name = linesplit[0]
        plasmid_name_split = plasmid_name.split('_')
        plasmid_name_uniq = str(plasmid_name_split[6]) + "_" + str(plasmid_name_split[7])
        if uniq_plasmid == plasmid_name_uniq:
            uniq_plasmid_length = uniq_plasmid_length + int(linesplit[1])
            uniq_plasmid_coverage = uniq_plasmid_coverage + float(plasmid_name_split[5])
            median_array.append(uniq_plasmid_coverage)
            longest_contig[uniq_plasmid].append(linesplit[1])

    median_coverage[uniq_plasmid] = statistics.median(median_array)

    # Calculate Mean Coverage
    command_3 = "grep \'%s\' %s | wc -l" % (uniq_plasmid, sample_name)
    no_of_uniq_plasmid_contigs = subprocess.check_output(command_3, shell=True)
    mean_coverage = float(uniq_plasmid_coverage) / int(no_of_uniq_plasmid_contigs)

    # Save length of each unique plasmid in dictionary
    plasmid_length[uniq_plasmid] = uniq_plasmid_length


# Calculate The entire length in bp in an assembly
plasm_comp_size = sum(plasmid_length.values())


# Blast KPC2 gene against each assembly and find out which plasmid contains the KPC2 gene
command_4 = "blastn -db %s -query %s -outfmt 6" % (args.AMR_gene, sample_name)
blast_output = subprocess.check_output(command_4, shell=True)
kpc_positive = OrderedDict()
for line in blast_output.splitlines():
    linesplit = line.split('\t')
    query_id = linesplit[0]
    align_length = linesplit[3]
    if align_length > 800:
        query_id_split = query_id.split('_')
        query_id_split_uniq = str(query_id_split[6]) + "_" + str(query_id_split[7])
        #sample_name = "%s"
        kpc_positive[query_id_split_uniq] = "KPC+"

# Prepare Final printable output string
final_string = ""
no_of_plasmids = len(uniq_plasmid_names)
total_no_of_contigs = no_of_contigs.strip()

plasmid_name_length_string = ""
for i in plasmid_length:
    plasmid_name_length_string = plasmid_name_length_string + i + ":" + str(plasmid_length[i]) + ","

plasmid_median_string = ""
for i in median_coverage:
    plasmid_median_string = plasmid_median_string + i + ":" + str(median_coverage[i]) + ","

longest_contig_string = ""
for i in longest_contig:
    longest_contig_string = longest_contig_string + i + ":" + str(longest_contig[i][0]) + ","

length_median_longest_string = "\t"
for i in plasmid_length:
    if i in kpc_positive.keys():
        kpc_string = "KPC+,%s" % sample_name
    else:
        kpc_string = ""
    length_median_longest_string = length_median_longest_string + ",,,,," + i + "," + str(plasmid_length[i]) + "," + str(median_coverage[i]) + "," + str(longest_contig[i][0]) + "," + kpc_string + "\n"

header = "sample\tTotal # of contigs\t# of plasmids\ttotal_contigs_length(bp)\tplasmid_length(bp)\tplasmid_median_cov\tmedian_coverage\tlongest_contigs_in_plasmid"
#print header
#final_string = os.path.basename(sample_name) + "\t" + str(total_no_of_contigs) + "\t" + str(no_of_plasmids) + "\t" + plasmid_name_length_string + "\t" + plasmid_median_string + "\t" + str(plasm_comp_size) + "\t" + longest_contig_string
print sample_name + "\t" + str(total_no_of_contigs) + "\t" + str(no_of_plasmids) + "\t" + str(plasm_comp_size) + "\t\t\t" + length_median_longest_string






















