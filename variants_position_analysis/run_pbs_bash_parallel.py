__author__ = 'alipirani'


import os
import readline
import argparse
from joblib import Parallel, delayed
import multiprocessing
parser = argparse.ArgumentParser(description='Alignment using Muscle')
parser.add_argument('-command', action='store', dest="command", help='command file containing commands for aligning each fasta file')
args = parser.parse_args()

#print args.filenames
command_file = args.command
command_array = []
with open(command_file) as fp:
    for line in fp:
        line = line.strip()
        command_array.append(line)
num_cores = multiprocessing.cpu_count()

def run_command(i):
    os.system(i)
    done = "done"
    return done
results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in command_array)
