__author__ = 'alipirani'


import os
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description='Convert Genbank to Fasta')
parser.add_argument('-dir', action='store', dest="dir", help='Directory where Genbank files are located')
parser.add_argument('-filenames', action='store', dest="filenames", help='Genbank filenames: One name per line')
parser.add_argument('-out', action='store', dest="out", help='Output directory')
args = parser.parse_args()
directory = args.dir
filenames = args.filenames

filenames_array = []

with open(filenames) as fp:
    for line in fp:
        line = line.strip()
        line = directory + "/" + line
        filenames_array.append(line)

for file in filenames_array:
    print "Converting %s ......" % file
    input_handle = open(file, "rU")
    #print os.path.basename(file.rsplit( ".", 1 )[ 0 ])
    #out_file = args.out + "/" + os.path.basename(file).replace('.gb', '.fasta')
    out_file = args.out + "/" + os.path.basename(file.rsplit( ".", 1 )[ 0 ]) + ".fasta"
    output_handle = open(out_file, "w")
    sequences = SeqIO.parse(input_handle, "genbank")
    count = SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()
    input_handle.close()

