# This scripts requires Biopython.
# Use this script to extract all the CDS/rRNA genes from a genbank file or to extract user-specified genes
# Each genes will be extracted in fasta format in output directory.
# All CDS genes in fasta format can be used for LS-BSR analysis

from Bio import SeqIO
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Extract Genes from Genbank File: Either all CDS genes or specific genes')
parser.add_argument('-extract', action='store', dest="extract", help='A/S: A for extracting all CDS genes. S for specific genes mentioned in file provided with argument -genes')
parser.add_argument('-genes', action='store', dest="genenames_file", help='Specific genes to extract from Genbank file. One gene name per line.')
parser.add_argument('-gb', action='store', dest="genbank_file", help='Genbank File')
parser.add_argument('-out', action='store', dest="output", help='Output directory')
args = parser.parse_args()

filename = args.genbank_file
out_directory = args.output
filebase = os.path.splitext(filename)[0]
handle = open(filename, 'rU')
log_handle = open(out_directory + '/' + 'extract_genes.log', 'w+')
def extract_cds_genes():
    features_extracted = 0
    for record in SeqIO.parse(handle, 'genbank') :
        seq = str(record.seq)
        for feature in record.features:
            log_handle.write(str(feature.type) + '\n')
            log_handle.write(str(feature.location) + '\n')
            if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'misc_RNA' or feature.type == 'repeat_region' or feature.type == 'rRNA' or feature.type == 'tmRNA' or feature.type == 'tRNA':
                if 'gene' in feature.qualifiers:
                    geneName = feature.qualifiers['gene'][0]
                elif 'product' in feature.qualifiers:
                    geneName = feature.qualifiers['product'][0]
                else:
		    if 'note' in feature.qualifiers:
		        geneName = (feature.qualifiers['note'][0]).replace(' ', '')
		    else:
			geneName = "UNKNOWN"
                    print 'Gene/Product not found. ERROR when parsing feature:'
                    print feature.qualifiers
                    #geneName = "UNKNOWN"
                    #quit()
                geneName = geneName.replace(' ', '-')
                geneName = geneName.replace('/', '-')
                log_handle.write(geneName + '\n')
                geneFile = open(out_directory + '/' + geneName + '.fa', 'w+')
                geneFile.write('>')
                geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                geneFile.write("\n\n")
                features_extracted += 1
            log_handle.write('----------------------------' + '\n')
    log_handle.write('\n\nTotal features extracted from Genbank file: %s' % features_extracted + '\n')

def extract_individual_gene():
    genome_name = os.path.basename(filebase)
    absent_genomes = []
    for record in SeqIO.parse(handle, 'genbank') :
        seq = str(record.seq)
        for feature in record.features:
            if feature.type == 'CDS':
                #print "CDS"
                #print feature
                if 'gene' not in feature.qualifiers:
                    if filename not in absent_genomes:
                        absent_genomes.append(filename)
                    continue
                    #print "No gene"
                if feature.qualifiers['gene'][0] in gene_array:
                    #print feature
                    geneName = feature.qualifiers['gene'][0]
                    print geneName
                    log_handle.write(genome_name + '\n')
                    log_handle.write(str(feature.type) + ':\t' + geneName + '\t' + str(feature.location) + '\n')
                    geneName = geneName.replace(' ', '-')
                    geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".fna"
                    geneFile = open(geneFilename, 'w+')
                    geneFile.write('>')
                    geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                    geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                    geneFile.write("\n\n")
                    translation_seq = feature.qualifiers['translation'][0]
                    proteinFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                    proteinFile = open(proteinFilename, 'w+')
                    proteinFile.write('>')
                    proteinFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                    proteinFile.write(translation_seq)
                    proteinFile.write("\n\n")


if args.extract == "A":
    print "\nExtracting all CDS genes from %s into directory: %s\n" % (filename, out_directory)
    extract_cds_genes()
elif args.extract == "S":

    if args.genenames_file:
        print "\nExtracting Speciifc genes from %s into directory: %s\n" % (filename, out_directory)
        print "Reading gene names..."
        gene_array = []
        with open(args.genenames_file) as fp:
            for line in fp:
                line = line.strip()
                gene_array.append(line)
            fp.close()
        print gene_array
        extract_individual_gene()

    else:
        print "Asked for specific genes in -extract argument. Please provide a file listing name of specific genes to extract with -genes argument\n"
        exit()






