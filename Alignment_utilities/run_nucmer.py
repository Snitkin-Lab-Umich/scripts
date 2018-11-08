from __future__ import division
import os
import argparse
import re
import subprocess
import statistics
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser(description='Run Nucmer on Sample 2008 All-vs-All')
parser.add_argument('-filename', action='store', dest="filename", help='filename of pseudomolecule assembly fasta file')
parser.add_argument('-out', action='store', dest="out", help='out_directory')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fasta files')
args = parser.parse_args()

filenames_array = []

with open(args.filename) as fp:
    for line in fp:
        line = line.strip()
        line = args.dir + "/" + line
        filenames_array.append(line)


def run_nucmer():
    for file in filenames_array:
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
        os.system(mkdir_command)
        for file_again in filenames_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=%s %s %s && show-coords -r -c -l -T %s.delta > %s.coords && mv %s.delta %s.coords %s/%s" % (prefix, file, file_again, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
            os.system(command)

run_nucmer()
def parse_coord():
    for folder in filenames_array:
        folder_base = os.path.basename(folder)
        list_command = "ls %s/%s/*.coords" % (args.out, folder_base.replace('.fasta', ''))
        list_of_coords_files = subprocess.check_output(list_command, shell=True)
        comparison_file = args.out + (os.path.basename(folder)) + ".score"
        print comparison_file
        f_out = open(comparison_file, 'w+')
        for lines in list_of_coords_files.splitlines():
            aligned_bases = 0
            f = open(lines,'r')
            r_line = f.readlines()[4:]
            bed_out_file = lines + ".bed"
            f1 = open(bed_out_file, 'w+')
            cmd1 = "tail -n1 %s | awk -F\'\t\' \'{print $8}\'" % lines
            r_length = subprocess.check_output(cmd1, shell=True)
            cmd2 = "tail -n1 %s | awk -F\'\t\' \'{print $9}\'" % lines
            q_length = subprocess.check_output(cmd2, shell=True)
            for i in r_line:
                i = i.strip()
                i = i.split('\t')
                name = os.path.basename(lines)
                write_string = name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                #r_length = i[7]
                #q_length = i[8]
                f1.write(write_string)
            f1.close()
            bed_cmd = "bedtools merge -i %s" % bed_out_file
            merged_bed_positions = subprocess.check_output(bed_cmd, shell=True)
            for rows in merged_bed_positions.splitlines():
                rows_split = rows.split('\t')
                sum = int(rows_split[2]) - int(rows_split[1]) + 1
                aligned_bases = aligned_bases + sum
            uniq_ref = int(r_length.strip()) - aligned_bases
            uniq_query = int(q_length.strip()) - aligned_bases
            deno = aligned_bases + uniq_ref + uniq_query
            score = float(aligned_bases / deno)
            #score = aligned_bases/(aligned_bases + uniq_ref + uniq_query)
            print score
            f_out.write(str(score) + '\n')
#parse_coord()

# for i in *.fa*; do cat $i | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | cat <( echo ">$i") - > All_pseudomolecule/$i
# for i in */*.delta; do show-coords -r -c -l -T $i > $i.coords; done

# Create pseudomolecule from contigs and save it in a directory
# for i in *.fa*; do cat $i | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | cat <( echo ">$i") - > All_pseudomolecule/$i
# Run this script from the directory where pseudomolecule is saved
# It will create score files in the output directory.
# example: ~/anaconda/bin/python run_nucmer.py -filename filenames -out alignment_delta/ -dir ./All_pseudomolecule

