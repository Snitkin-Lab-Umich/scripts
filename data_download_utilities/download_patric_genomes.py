"""
Description:

Genome name and Genome ID can be obtained in the following manner: https://www.patricbrc.org/view/Taxonomy/2#view_tab=genomes
Click Download on this website to download information in csv format.

Features that can be extracted:

fna                 : Genome Assembly in Fasta format
PATRIC.cds.tab      : cds annotations in tab format
PATRIC.faa          : Amino acid sequences of all the features
PATRIC.features.tab : Annotation features in tab format
PATRIC.gff          : GFF annotations
PATRIC.pathway.tab  : Pathway annotation for each feature in tab format
PATRIC.rna.tab      : rna annotations in tab format
PATRIC.spgene.tab   : Special genes annotations such as AR genes/Drug target in tab format
RefSeq.cds.tab      : Refseq cds annotations in tab format
RefSeq.faa          : Refseq Amino acid sequences of all the genes
RefSeq.gbf          : Refseq genbank file
RefSeq.gff          : Refseq GFF file
RefSeq.rna.tab      : Refseq rna annotation in tab format

"""

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import glob
import readline
import errno
from datetime import datetime

parser = argparse.ArgumentParser(description='Download Genome data from PATRIC Database. '
                                             'Note: The script requires Genome name and Genome Id for the genomes that you want it to download from Patric. '
                                             'Both of this field are included in Patric metadata csv file that can be downloaded for your desired genomes from this link:'
                                             'https://www.patricbrc.org/view/Taxonomy/2#view_tab=genomes')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-csv', action='store', dest="csv", help='Patric csv metadata file containing Genome ID and Genome name resp. in first and second column')
required.add_argument('-features', action='store', dest="features", help='comma seperated list of features to download. ex: fna,PATRIC.faa,RefSeq.gbf')
required.add_argument('-out', action='store', dest="output", help='Output directory')
required.add_argument('-parallel', action='store', dest="parallel", help='Download genome data in parallel: yes/no')
args = parser.parse_args()

def read_csv_metadata(csv_metadata):
    genome_id_name_dict = {}
    with open(args.csv, 'rU') as csv_file:
        print "Reading Patric CSV Metadata file: %s\n" % args.csv
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader, None)
        for row in csv_reader:
            genome_id_name_dict[row[0]] = row[1]

    return genome_id_name_dict


def run_command(i):
    print "Running: %s" % i
    os.system(i)
    #subprocess.call(i, shell=True)
    done = "Downloaded: %s" % i
    return done


def download_data(genome_id_names, features_list, csv_metadata, log_dir, out_dir, parallel):
    print "Downloading genome data from PATRIC...\n"
    link_list = []
    for genome_id in genome_id_names.keys():
        for feature in features_list:
            genome_name = re.sub(r'\W+', '_', genome_id_names[genome_id])
            link = "wget ftp://ftp.patricbrc.org/patric2/genomes/%s/%s.%s --output-file=%s/%s.%s.log --output-document=%s/%s.%s" % (genome_id, genome_id, feature, log_dir, genome_name, feature, out_dir, genome_name, feature)
            link_list.append(link)
            print link
    if parallel == "yes":
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in link_list)
        # Check if the file was downloaded
        for genome_id in genome_id_names.keys():
            genome_name = re.sub(r'\W+', '_', genome_id_names[genome_id])
            if not os.path.isfile("%s/%s.%s" % (out_dir, genome_name, feature)):
                print "Fail to download:    %s   %s" % (genome_name, genome_id)
    else:
        for i in link_list:
            os.system(i)


def make_sure_path_exists(out_path):
    """
    Make sure the output folder exists or create at given path
    :param out_path:
    :return:
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()

if __name__ == '__main__':

    """Start Timer"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()

    print "The Script started at: %s\n" % start_time

    print "Start: Help Message:\nThis script will download genome data from PATRIC Database.\n\n" \
          "Note: The script requires Genome name and Genome Id for the genomes that you want it to download from Patric." \
          "Both of this field are included in Patric metadata csv file that can be downloaded for your desired genomes from this link: " \
          "https://www.patricbrc.org/view/Taxonomy/2#view_tab=genomes\n\n" \
          "\n" \
          "Genome name and Genome ID in csv format can be obtained in the following manner: \n" \
          "1. Go to https://www.patricbrc.org/view/Taxonomy/2#view_tab=genomes\n" \
          "2. Enter the name of organism in search bar.\n" \
          "3. Click Download on this website to download this organism's metadata information in csv format.\n" \
          "\nFeatures that can be extracted from PATRIC using this script:\n" \
          "fna                 : Genome Assembly in Fasta format\n" \
          "PATRIC.cds.tab      : cds annotations in tab format\n" \
          "PATRIC.faa          : Amino acid sequences of all the features\n" \
          "PATRIC.features.tab : Annotation features in tab format\n" \
          "PATRIC.gff          : GFF annotations\n" \
          "PATRIC.pathway.tab  : Pathway annotation for each feature in tab format\n" \
          "PATRIC.rna.tab      : rna annotations in tab format\n" \
          "PATRIC.spgene.tab   : Special genes annotations such as AR genes/Drug target in tab format\n" \
          "RefSeq.cds.tab      : Refseq cds annotations in tab format\n" \
          "RefSeq.faa          : Refseq Amino acid sequences of all the genes\n" \
          "RefSeq.gbf          : Refseq genbank file\n" \
          "RefSeq.gff          : Refseq GFF file\n" \
          "RefSeq.rna.tab      : Refseq rna annotation in tab format\n" \
          "\nEnd: Help Message:\n"

    print "Generating an output/log folders\n"
    make_sure_path_exists(args.output)
    log_directory = "%s/log_files" % args.output
    make_sure_path_exists(log_directory)
    genome_id_names = read_csv_metadata(args.csv)
    print genome_id_names
    features_list = args.features.split(',')
    download_data(genome_id_names, features_list, args.csv, log_directory, args.output, args.parallel)
