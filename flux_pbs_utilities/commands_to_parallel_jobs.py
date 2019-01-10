__author__ = 'alipirani'

import os
import re
import readline
import argparse
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser(description='This script takes a file containing commands(-command command_file) and generates a PBS job for each commands')
parser.add_argument('-command', action='store', dest="command", help='command file containing one command on each line')
parser.add_argument('-pbs_memory_arg', action='store', type=str, dest="pbs_memory_arg", help='pbs memory requirement such as nodes=1:ppn=4,pmem=4000mb,walltime:72:00:00')
parser.add_argument('-email', action='store', dest="email", help='pbs email argument for notifications')
parser.add_argument('-dir', action='store', dest="dir", help='Directory to save job files and PBS working directory')
parser.add_argument('-prefix', action='store', dest="prefix", help='prefix for PBS job file')
parser.add_argument('-acc', action='store', dest="acc", help='Account for PBS job file: fluxod or flux')
args = parser.parse_args()

command_array = []
with open(args.command) as fp:
    for line in fp:
        line = line.strip()
        command_array.append(line)
    fp.close()

count = 1
for i in command_array:
    job_name = "command_" + str(count)
    job_file_name = args.dir + "/" + job_name + ".pbs"
    job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m a\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A esnitkin_%s\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.email, args.pbs_memory_arg, args.acc, args.acc, args.dir, i)
    job_file_name = "%s_%s.pbs" % (args.prefix, job_name)
    f1=open(job_file_name, 'w+')
    f1.write(job_print_string)
    count += 1
    f1.close()
