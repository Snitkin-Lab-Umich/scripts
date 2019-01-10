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
import errno
parser = argparse.ArgumentParser(description='Parse nucmer coordinates and delta files')
parser.add_argument('-folder', action='store', dest="folder", help='filename of pseudomolecule assembly fasta file')
parser.add_argument('-out', action='store', dest="out", help='out_directory')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fasta files')
args = parser.parse_args()


def parse_coord_aggregate(folder, folder_dir):
    folder_base = os.path.basename(folder)
    list_command = "ls %s/%s/*.coords" % (args.out, folder_base.replace('.fasta', ''))
    list_of_coords_files = subprocess.check_output(list_command, shell=True)
    comparison_file_anno = args.out + (os.path.basename(folder)) + "_anno.score"
    f_out_anno_score = open(comparison_file_anno, 'w+')
    for lines in list_of_coords_files.splitlines():
        print lines
        file_prefix = lines.replace('.coords', '')
        aligned_bases = 0
        anno_aligned_bases = 0
        f = open(lines,'r')
        r_line = f.readlines()[4:]
        # query_name = ((os.path.basename(lines).replace(contig_name, '')).replace('.coords', '')).replace('_final_ordered_', '')
        if len(r_line) > 0:
            cmd1 = "tail -n1 %s | awk -F\'\t\' \'{print $8}\'" % lines
            r_length = subprocess.check_output(cmd1, shell=True)
            cmd2 = "tail -n1 %s | awk -F\'\t\' \'{print $9}\'" % lines
            q_length = subprocess.check_output(cmd2, shell=True)
            r_name = "tail -n1 %s | awk -F\'\t\' \'{print $12}\'" % lines
            q_name = "tail -n1 %s | awk -F\'\t\' \'{print $13}\'" % lines
            reference_name = subprocess.check_output(r_name, shell=True)
            query_name = subprocess.check_output(q_name, shell=True)
            for i in r_line:
                i = i.strip()
                i = i.split('\t')
                gaps = ""
                line_name = os.path.basename(lines)
                # changed to 95. Hard coded
                if float(i[6]) >= 95.00 and int(i[4]) > 2000:
                    if float(i[4]) == float(i[5]):
                        gaps = "No gaps"
                        anno_aligned_bases = anno_aligned_bases + int(i[4])
                    elif float(i[4]) > float(i[5]):
                        gaps = float(i[4]) - float(i[5])
                        anno_aligned_bases = anno_aligned_bases + int(i[4])
                    elif float(i[5]) > float(i[4]):
                        gaps = float(i[5]) - float(i[4])
                        anno_aligned_bases = anno_aligned_bases + int(i[5])



            fragment_length = int(q_length.strip())
            score_anno = float(anno_aligned_bases / fragment_length)
            if score_anno > 0.9999:
                score_anno = 1.0
            elif score_anno < 0:
                score_anno = 0.0
            f_out_anno_score.write(str(score_anno) + '\n')

        else:
            #score = "0"
            aligned_bases = "0"
            anno_aligned_bases = "0"
            score_anno = "0"
            f_out_anno_score.write(str(score_anno) + '\n')

        #if anno_aligned_bases > 0.99:
        #    if r_length > q_length:
        #        print "%s" % (reference_name.strip())
        #    elif q_length > r_length:
        #        print "%s" % (query_name.strip())
        #    else:
        #        print "%s" % (query_name.strip())




def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

folder = args.folder
file_name = os.path.basename(folder)


#folder_dir = args.out + "/facilities/%s/%s/%s/" % (facility, facility_id, patient_id)
folder_dir = args.out + "/parse_nucmer_db"
make_sure_path_exists(folder_dir)


parse_coord_aggregate(folder, folder_dir)
