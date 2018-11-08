from __future__ import division
__author__ = 'alipirani'
import statistics
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np

parser = argparse.ArgumentParser(description='GATK Depth of Coverage Analysis')
parser.add_argument('-coverage_file', action='store', dest="coverage_file", help='coverage_file')
parser.add_argument('-bin', action='store', dest="bin", help='Number of Bins. It represents how many values should be there in a bin.')


args = parser.parse_args()
n = int(args.bin)
read_counts = []
out_file = args.coverage_file + "_depth_bins"
with open(args.coverage_file) as fp:
    next(fp)
    for line in fp:
        line = line.strip()
        line_split = line.split('\t')
        counts = int(line_split[1])
        read_counts.append(int(counts))







def generate_mean_stats():
    # Partition/divide the read count array into bins of size 'n'.
    # The below command sets the value of i by dividing the read counts array each time starting with 0 and increasing to step 'n'.
    # So if n=1000 then the first bin will contain values from 0 to 1000 index of read counts array.
    partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
    # This below line creates a seperate array that was left out from partition due to lack of range
    partitions = partitions if len(partitions[-1]) == n else partitions[:-1]

    mean_array = []

    for x in partitions:
        mean = statistics.mean(x)
        mean_array.append(mean)

    with open(out_file, 'w') as out:
        header = "bin\tcount\n"
        count = 0
        out.write(header)
        for i in mean_array:
            count += 1
            out.write(str(count)+'\t'+str(i)+'\n')

generate_mean_stats()