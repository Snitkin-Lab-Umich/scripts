__author__ = 'alipirani'

import argparse
import re
import os
from joblib import Parallel, delayed
import multiprocessing
# Command line argument parsing
# Command line argument parsing
parser = argparse.ArgumentParser(description='Run Raxml Parallel')
parser.add_argument('-f1', action='store', dest="file_1", help='phylip file names')
parser.add_argument('-d', action='store', dest="directory", help='directory where phylip files are located')
parser.add_argument('-o', action='store', dest="output", help='output directory where the raxml results will be saved')
args = parser.parse_args()
alignment_file = args.file_1
#out_dir = args.output + "/"
arr = []
muscle_cmd = []
directory = args.directory + "/"
phylip_filenames = []
#parameter = "-f a -x 12345 -p 12345 -N 50 -k -m GTRCAT"
parameter = "-f a -x 12345 -p 12345 -N 50 -k -m PROTCATGTR"
with open(alignment_file) as fp:
    for line in fp:
        line = line.strip()
        line = args.directory + line
        phylip_filenames.append(line)


def create_directory(phylip_filenames):
    for file in phylip_filenames:
        cmd = "mkdir %s" % file + "_raxml"
        os.system(cmd)

#create_directory(phylip_filenames)
num_cores = multiprocessing.cpu_count()

def run_raxml(i):
    dir = i + "_raxml"
    change_dir_cmd = "cd %s" % dir
    os.system(change_dir_cmd)
    filename = os.path.basename(i)
    cmdstring = "/home2/apirani/bin/raxmlHPC %s -s %s -n %s_final.tre" % (parameter, i, filename)
    os.system(cmdstring)



results = Parallel(n_jobs=num_cores)(delayed(run_raxml)(i) for i in phylip_filenames)


