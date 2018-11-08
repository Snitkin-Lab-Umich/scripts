_author__ = 'alipirani'
import os
import argparse
import re
import subprocess
import statistics
from collections import defaultdict
from collections import OrderedDict
from Bio import SeqIO
parser = argparse.ArgumentParser(description='Extract Specific Contigs')
parser.add_argument('-filename', action='store', dest="filename", help='fasta file names')
parser.add_argument('-out', action='store', dest="out", help='out_directory')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fasta files')
parser.add_argument('-contigs', action='store', dest="contigs", help='Contig names; File    ContigID')
args = parser.parse_args()
filenames_array = []
with open(args.filename) as fp:
    for line in fp:
        line = line.strip()
        line = args.dir + "/" + line
        filenames_array.append(line)

contigid_dict = defaultdict(list)
with open(args.contigs) as fp:
    for line in fp:
        line = line.strip()
        linesplit = line.split('\t')
        contigid_dict[linesplit[0]].append(linesplit[1])




for file in filenames_array:
    handle = open(file, "rU")
    out_file = args.out + "/" + os.path.basename(file)
    with open(out_file, 'w+') as out:
        for record in SeqIO.parse(handle, "fasta"):
            contig_name_split = record.id.split('_')
            contig_name = str(contig_name_split[6]) + "_" + str(contig_name_split[7])
            #print file
            temp_str = str(contigid_dict[os.path.basename(file)])
            #print contig_name
            #print str(temp_str)
            if contig_name in temp_str:
                write_string = ">" + str(record.id) + "\n" + str(record.seq) + "\n"
                out.write(write_string)
        handle.close()



# for i in Rush_KPC_*.fasta; do cat $i | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | cat <( echo ">$i") - > All_pseudomolecule/$i; done