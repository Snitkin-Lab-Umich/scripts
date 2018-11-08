from __future__ import division
import os
import argparse
import re
import subprocess
import errno
import readline
from joblib import Parallel, delayed
import multiprocessing
import glob
import csv
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict
import time


parser = argparse.ArgumentParser(description='PBS Batch Scheduler: Manages and Submits batch jobs')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-pbs_job_scripts', action='store', dest="pbs_job_scripts", help='file with PBS job scripts to run', required=True)
args = parser.parse_args()

# GENERATE PBS job script file ARRAY
filenames_array = []
with open(args.pbs_job_scripts) as fp:
    for line in fp:
        line = line.strip()
        line = "qsub " + line
        filenames_array.append(line)



while True:
    get_no_of_jobs = "showq -w acct=\"esnitkin_fluxod\" -r | grep \"^Total jobs:\" | cut -d' ' -f4"
    proc = subprocess.Popen([get_no_of_jobs], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    n = 10
    if int(out) < 1000:
        print "No of Jobs is < 1000\n"
        print len(filenames_array)
        newlist = filenames_array[:n]
        print newlist
        del filenames_array[:n]
        print len(filenames_array)
        time.sleep(20)