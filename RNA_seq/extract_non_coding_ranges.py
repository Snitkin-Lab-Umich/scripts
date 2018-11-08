# Use this script to extract all the non-coding genomic ranges not present in the given gff file

import sys
import argparse
import os
import csv
from collections import OrderedDict
from collections import defaultdict

parser = argparse.ArgumentParser(description='Extract non-coding genomic ranges')
parser.add_argument('-gff', action='store', dest="gff", help='GFF file')
parser.add_argument('-bed', action='store', dest="bed", help='BED file')
parser.add_argument('-count', action='store', dest="count", help='count file')
args = parser.parse_args()



def extract_intergenic_regions_from_gff():
    gff_start = []
    gff_end = []
    genome_dict = defaultdict(list)
    with open(args.gff, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #print row[1][0]
            #gff_start.append(line_split[3])
            #gff_end.append(line_split[4])
            for row in csv_reader:

                gff_start.append(row[3])
                gff_end.append(row[4])

    basename = os.path.basename(args.gff)
    genome_name = os.path.splitext(basename)[0]

    first_range = genome_name + "\tFIG\tnc\t" + "1\t%s" % str(int(gff_start[0])-1) + "\t.\t+\t0\tName=Non-coding Region 0"
    print first_range
    for i in range(1, len(gff_start)):
        range_start = i - 1
        range_end = i
        if int(gff_end[range_start]) < int(gff_start[i]):
            #print gff_end[range_start] + "\t" + gff_start[i]
            print genome_name + "\tFIG\tnc\t" + str(int(gff_end[range_start])+1) + "\t" + str(int(gff_start[i])-1) + "\t.\t+\t0\tName=Non-coding Region " + str(i)
        else:
            continue

# this method will work on only one chromosome/genome/single large combined supercontig
#extract_intergenic_regions_from_gff()

def bed_to_gff():
    with open(args.bed, 'rU') as csv_file:
        new_nc_gff_file = str(args.bed).replace('.bed', '.gff')
        new_nc_gff = open(new_nc_gff_file, 'w+')
	new_nc_gff.write('##gff-version 3\n')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
	    if int(row[1]) == 0:
                    row[1] = "1"
	    print_string = row[0] + "\tFIG\tnc\t" + str(row[1]) + "\t" + row[2] + "\t.\t.\t.\tName=Intergenic region 0\n"
            new_nc_gff.write(print_string)
	    count = 1
            for row in csv_reader:
		print row
		if int(row[1]) == 0:
		    row[1] = "1"
                print_string = row[0] + "\tFIG\tnc\t" + str(row[1]) + "\t" + row[2] + "\t.\t.\t.\tName=Intergenic region %s\n" % str(count)
		#print print_string
                new_nc_gff.write(print_string)
                count = count + 1
	





#bed_to_gff()

def extract_outlier_antisense_row():
    file = args.count
    new_file_name = file + "_outliers"
    new_file = open(new_file_name, 'w+')
    label = ""
    with open(file, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            for row in csv_reader:
                if int(row[3]) > int(row[1]):
		    print_string = ""
                    if int(row[3]) > 2*int(row[1]):
		        label = "outlier"
			print_string = row[0] + "\t" + row[1] + "\t" + row[3] + "\t" + label + "\n"
			#print row
		    else:
			label = "-"
                        print_string = row[0] + "\t" + row[1] + "\t" + row[3] + "\t" + label + "\n"
		    new_file.write(print_string)


extract_outlier_antisense_row()
