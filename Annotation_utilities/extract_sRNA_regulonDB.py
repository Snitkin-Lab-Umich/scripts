__author__ = 'alipirani'
# This scripts requires Biopython.
# Use this script to extract small RNA sequences a regulonDb gene sequences file using the key text file.


from Bio import SeqIO
import sys
import argparse
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbitblastnCommandline

parser = argparse.ArgumentParser(description='Extract small RNA sequences from RegulonDB sequences File')
parser.add_argument('-gene_sequences', action='store', dest="gene_sequences_file", help='File with gene sequences downloaded from RegulonDB')
parser.add_argument('-map', action='store', dest="map_file", help='Map File with small RNA sequences and ID downloaded fron RegulonDB')
parser.add_argument('-out_file', action='store', dest="output", help='Output file')
parser.add_argument('-genome', action='store', dest="genome", help='Genome fasta file')
args = parser.parse_args()
regulonDB_out_file = "./regulonDB.faa"
genome_file =args.genome
genome_filename_split = genome_file.split('.')
first_part = genome_filename_split[0]
sequence_length = OrderedDict()
map_file_dict = OrderedDict()
map_prod_name = OrderedDict()
def create_genome_db():
    makeblastdb_cmd = "makeblastdb -in %s -out %s -dbtype nucl" % (args.genome, args.genome)
    os.system(makeblastdb_cmd)


def annotate_sRNA_regulonDB():

    with open(args.map_file) as fp:
        for line in fp:
            if not line.startswith('#'):
                line = line.strip()
                line_Split = line.split('\t')
                print len(line_Split)
                if len(line_Split) > 6:
                    prod_name = line_Split[6]
                else:
                    prod_name = ""
                map_file_dict[line_Split[0]] = line_Split[1]
                map_prod_name[line_Split[0]] = prod_name
        fp.close()

    f1 = open(regulonDB_out_file, 'w+')
    with open(args.gene_sequences_file) as fp:
        for line in fp:
            if not line.startswith('#'):
                line = line.strip()
                line_Split = line.split('\t')
                if line_Split[0] in map_file_dict.keys():
                    if line_Split[1] == "()":
                        prod_name = "small RNA"
                    else:
                        prod_name = line_Split[6]
                    gene_name = map_file_dict[line_Split[0]]
                    seq = line_Split[9]
                    coding_dna = Seq(seq, IUPAC.unambiguous_dna)
                    translated_seq = coding_dna.translate(table="Bacterial")
                    sequence_length[line_Split[0]] = len(translated_seq)
                    print_string_header = ">%s ~~~%s~~~%s\n" % (line_Split[0], gene_name, prod_name)
                    print_string_seq = str(translated_seq)
                    final_print_string = print_string_header + print_string_seq + "\n"
                    f1.write(final_print_string)
                    #print sequence_length


def blast_and_parse():
    regulon_db_blast_output = "%s_regulaon_db_blast" % first_part
    tblastn_cmd_line = "tblastn -query %s -db %s -outfmt 6 -evalue 1e-10 | sort -k1,1 -k12,12nr -k11,11n | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > %s" % (regulonDB_out_file, genome_file, regulon_db_blast_output)
    os.system(tblastn_cmd_line)
    f2 = open(args.output, 'w+')
    gff_header = "##gff-version 3\n"
    f2.write(gff_header)
    with open(regulon_db_blast_output) as fp:
        for line in fp:
            line = line.strip()
            line_Split = line.split('\t')
            if float(line_Split[2]) > 95.00:
                query_length = int(sequence_length[line_Split[0]])
                if line_Split[3] > (0.95*query_length):
                    if int(line_Split[9] > line_Split[8]):
                        strand = "+"
                        start = line_Split[8]
                        end = line_Split[9]
                    else:
                        strand = "-"
                        start = line_Split[9]
                        end = line_Split[8]
                    gff_string = "%s\tRegulonDB\tsRNA\t%s\t%s\t.\t%s\t0\tID=%s;inference=RegulonDB;locus_tag=%s;gene=%s;product=%s\n" % (line_Split[1], start, end, strand, line_Split[0], line_Split[0], map_file_dict[line_Split[0]], map_prod_name[line_Split[0]])
                    #print gff_string
                    f2.write(gff_string)



#def prokka_annotation():
    
create_genome_db()
annotate_sRNA_regulonDB()
blast_and_parse()