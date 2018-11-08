from __future__ import division
import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import thread
import glob
import readline
import pandas as pd
import errno
from pyfasta import Fasta
from datetime import datetime
import threading


parser = argparse.ArgumentParser(description='Extract Only reference and variant positions and generate a fasta file out of it.')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
required.add_argument('-filter2_only_snp_vcf_filename', action='store', dest="filter2_only_snp_vcf_filename",
                    help='Name of filter2 only SNP vcf file')
required.add_argument('-reference', action='store', dest="reference",
                    help='Path to Reference Fasta File')
args = parser.parse_args()


def extract_only_ref_variant_fasta():
    f = Fasta(args.reference)
    if len(f.keys()) == 1:
        ref_id = str(f.keys())
        #print ref_id

    ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir).readlines()
    core_vcf_file = args.filter2_only_snp_vcf_filename.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_core.vcf.gz')
    print core_vcf_file
    #fasta_string = ">%s\n" % core_vcf_file.replace('_filter2_final.vcf_core.vcf.gz', '')
    fasta_string = ""
    firstLine = ffp.pop(0)
    print len(ffp)
    count = 0
    for lines in ffp:
        #next(ffp)
        lines = lines.strip()
        #grep_position = "zcat %s | grep -w \'%s\' | awk -F\'\\t\' \'{print $5}\' | wc -l" % (core_vcf_file, lines)
        grep_position = "zcat %s | grep -v \'#\' | awk -F\'\\t\' \'{ if ($2 == %s) print $0 }\' | awk -F\'\\t\' \'{print $5}\'" % (core_vcf_file, lines)
        #print grep_position
        proc = subprocess.Popen([grep_position], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        if out and "," not in out:
            #print out
            fasta_string = fasta_string + out
            count += 1
            #print fasta_string
            #test = "1"
        else:
            #extract_base = "tr -d \'\\n\' < %s | cut -b%s" % (args.reference, lines)
            extract_base = "grep -v \'>\' %s | tr -d \'\\n\'| cut -b%s" % (args.reference, lines)
            #print extract_base
            # proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
            # (out, err) = proc.communicate()
            # out = out.strip()
            # fasta_string = fasta_string + out
            # if not out:
            #     print extract_base
            #     print "Error extracting reference allele"
            #out = str(f.sequence({'chr': str(f.keys()[0]), 'start': int(lines), 'stop': int(lines)}))
            #print str(f.sequence({'chr': str(f.keys()[0]), 'start': int(lines), 'stop': int(lines)}))
            count += 1
            fasta_string = fasta_string + str(f.sequence({'chr': str(f.keys()[0]), 'start': int(lines), 'stop': int(lines)}))
            #print fasta_string


    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(core_vcf_file.replace('_filter2_final.vcf_core.vcf.gz', '')) + fasta_string
    fp = open("%s/%s_variants.fa" % (args.filter2_only_snp_vcf_dir, os.path.basename(core_vcf_file.replace('_filter2_final.vcf_core.vcf.gz', ''))), 'w+')
    print final_fasta_string
    print "Count: %s " % count
    print "Length: %s " % len(fasta_string)
    fp.write(final_fasta_string + '\n')
    fp.close()
    #ffp.close()
    #sequence_lgth_cmd = "for i in %s/*_variants.fa; do bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % args.filter2_only_snp_vcf_dir
    #os.system(sequence_lgth_cmd)















extract_only_ref_variant_fasta()
