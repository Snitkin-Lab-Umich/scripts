from __future__ import division
import os
import argparse
import re
import subprocess
import readline
from joblib import Parallel, delayed
import multiprocessing
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser(description='Run Nucmer on Sample 2008 All-vs-All')
parser.add_argument('-folder', action='store', dest="folder", help='filename of pseudomolecule assembly fasta file')
parser.add_argument('-species', action='store', dest="species", help='species')
parser.add_argument('-patient', action='store', dest="patient", help='patient')
args = parser.parse_args()

filenames_array = []
species_array = []


with open(args.folder) as fp:
    for line in fp:
        line = line.strip()
        #line = args.dir + "/" + line
        filenames_array.append(line)

with open(args.species) as fp:
    for line in fp:
        line = line.strip()
        #line = args.dir + "/" + line
        species_array.append(line)



patient_array = []

with open(args.patient) as fp:
    for line in fp:
        line = line.strip()
        #line = args.dir + "/" + line
        patient_array.append(line)

uniq_array = []
reverse_string_array = []

# for file in filenames_array:
#     for file_again in filenames_array:
#         if file != file_again:
#             generate_string = str(file) + "_" + str(file_again) + ".coords.annotations.txt"
#             reverse_string = str(file_again) + "_" + str(file) + ".coords.annotations.txt"
#             if reverse_string not in uniq_array:
#                 uniq_array.append(generate_string)
                #reverse_string_array.append(reverse_string)
#print uniq_array
#print len(uniq_array)
# f = open("list_uniq_files.txt", 'w+')
# for file in uniq_array:
#     f.write(file + '\n')





species_uniq_array = []
for file in species_array:
    for file_again in species_array:
        if file != file_again:
            generate_string = str(file) + "_" + str(file_again)
            reverse_string = str(file_again) + "_" + str(file)
            if reverse_string not in species_uniq_array:
                species_uniq_array.append(generate_string)
                #reverse_string_array.append(reverse_string)

for i in species_uniq_array:
    print i
#print species_uniq_array



copy_commands = []
for i in patient_array:
    for species_pair in species_uniq_array:
        species_1 = species_pair.split('_')[0]
        species_2 = species_pair.split('_')[1]
        # generate_copy = "cp *-" + i + "-*" + species_1 + "*-" + species_2
        # generate_copy_2 = "cp *-" + i + "-*" + species_2 + "*-" + species_1

        generate_copy = "cp *-" + i + "-*" + species_1 + "*-" + i + "-*" + species_2 + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"
        generate_copy_2 = "cp *-" + i + "-*" + species_2 + "*-" + i + "-*" + species_1 + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"

        copy_commands.append(generate_copy)
        copy_commands.append(generate_copy_2)
        if "coli" in species_1:
            # generate_copy_3 = "cp *-" + i + "-*" + + "E_coli" + "*-" + species_2
            # generate_copy_4 = "cp *-" + i + "-*" + species_2 + "*-" + "species_2"

            generate_copy_3 = "cp *-" + i + "-*" + "E_coli" + "*-" + i + "-*" + species_2 + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"
            generate_copy_4 = "cp *-" + i + "-*" + species_2 + "*-" + i + "-*" + "species_2" + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"


            copy_commands.append(generate_copy_3)
            copy_commands.append(generate_copy_4)
        elif "coli" in species_2:
            # generate_copy_5 = "cp *-" + i + "-*" + species_1 + "*-" + "Ecoli"
            # generate_copy_6 = "cp *-" + i + "-*" + "Ecoli" + "*-" + species_1

            generate_copy_5 = "cp *-" + i + "-*" + species_1 + "*-" + i + "-*" + "Ecoli" + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"
            generate_copy_6 = "cp *-" + i + "-*" + "Ecoli" + "*-" + i + "-*" + species_1 + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"

            copy_commands.append(generate_copy_5)
            copy_commands.append(generate_copy_6)

        if "bauman" in species_1:
            generate_copy_7 = "cp *-" + i + "-*" + "baumannii" + "*-" + i + "-*" + species_2 + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"
            generate_copy_8 = "cp *-" + i + "-*" + species_2 + "*-" + i + "-*" + "baumannii" + "*" + " /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_09_22_HGT_Mummer_analysis/All_Aim2_assembly_pseudomolecule/analysis/facility/*/" + i + "/" + species_pair + "/"
            copy_commands.append(generate_copy_7)
            copy_commands.append(generate_copy_8)



for i in  copy_commands:
    print i
