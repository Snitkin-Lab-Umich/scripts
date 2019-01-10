from __future__ import division
__author__ = 'alipirani'
import statistics
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np
import time
from datetime import datetime
start = time.time()

#Argument Parser
parser = argparse.ArgumentParser(description='Read BED file and generate PTR dataframe')
parser.add_argument('-bed', action='store', dest="bed", help='bedfile')
parser.add_argument('-out', action='store', dest="out", help='output file')
parser.add_argument('-bin', action='store', dest="bin", help='Number of Bins. It represents how many values should be there in a bin.')
parser.add_argument('-total_aligned_reads', action='store', dest="total", help='Total Aligned Reads')
parser.add_argument('-window_bin', action='store', dest="window_bin", help='Size of Window bin')
args = parser.parse_args()

# Initiate variables and array
window = int(args.window_bin)
perc_of_window = 0.6*window
#n = int(args.bin)
read_counts = []
print "[%s] Processing: %s" % (time.ctime(), args.bed)

#Read Bed file
with open(args.bed) as fp:
    for line in fp:
        line = line.strip()
        line_split = line.split('\t')
        counts = int(line_split[2])
        read_counts.append(int(counts))

# Methods
def generate_moving_sum_results(moving_sum_array):
    out_file = args.out + "_moving_sum_bins"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in moving_sum_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')
    peak = max(moving_sum_array)
    through = min([x for x in moving_sum_array if x !=0])
    PTR_moving = peak/through
    #print "The peak and trough values for moving_sum_array is: %s; %s" % (peak, through)
    out_file_ptr = args.out + "_PTR"
    with open(out_file_ptr, 'w+') as out:
        out.write(str(args.out)+' moving_sum_array :\t'+str(PTR_moving)+'\n')

    stats_file = args.bed.replace('_coverage.bed', '_alignment_stats')
    out_file = args.out + "_moving_sum_bins_normalize"
    cmd = "grep \'mapped (\' " + stats_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    out_split = out.split(' ')
    mapped_reads = out_split[0]
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in moving_sum_array:
            count += 1
            ratio = int(i)/int(mapped_reads)
            out.write(str(count)+','+str(ratio)+'\n')
    out.close()




def calculate_per_of_reads(median_sliding_window_array):
    stats_file = args.bed.replace('_coverage.bed', '_alignment_stats')
    cmd = "grep \'mapped (\' " + stats_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    out_split = out.split(' ')
    mapped_reads = out_split[0]
    print "The no of mapped reads are: %s" % mapped_reads
    out_file = args.bed + "_perc_bins"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                perc = (i*100)/int(mapped_reads)
                out.write(str(count)+','+str(perc)+'\n')

# Main Method
def smoothing():
    #moving_array_length = []
    # Create an array of Raw Counts from Bed file with a sliding window of 100 and 10000 bins
    raw_count_sliding_window_array = []
    raw_count_sliding_window_array_file = args.out + "raw_count_sliding_window_arrays"
    f1=open(raw_count_sliding_window_array_file, 'w+')
    raw_count_sliding_window_array_bins_with_0 = 0
    for i in xrange(0, len(read_counts), 500):
        start = i
        end = i + window
        f1.write(str(read_counts[start:end])+'\n')
        if read_counts[start:end].count(0) > perc_of_window:
            raw_count_sliding_window_array_bins_with_0 += 1
        else:
            raw_count_sliding_window_array.append(read_counts[start:end])
    #print "The length of raw_count_sliding_window_array is %s" % len(raw_count_sliding_window_array)
    print "The number of arrays in Raw count Sliding window array with a count of 0 more than 60 percent is: %s" % raw_count_sliding_window_array_bins_with_0
    # Get the Sum of each counts bins in raw_count_sliding_window_array
    moving_sum_array = []
    for i in raw_count_sliding_window_array:
        moving_sum_array.append(sum(i))
    generate_moving_sum_results(moving_sum_array)

    # Calculate Median of each 10000 moving_sum_array bins by sliding 100 bins
    median_sliding_window_array = []
    diff = 0
    for i in xrange(0, len(read_counts), 500):
        start = i
        end = i + window
        #moving_array_length.append(len(moving_sum_array[start:end]))
        if len(moving_sum_array[start:end]) > 500:
            median_sliding_window_array.append(statistics.median(moving_sum_array[start:end]))

    # Calculate PTR
    peak = max(median_sliding_window_array)
    through = min([x for x in median_sliding_window_array if x !=0])
    PTR_median = peak/through
    print "The peak and trough values for median_sliding_window_array is: %s; %s" % (peak, through)
    #print "#####%s , %s" % (median_sliding_window_array[393], median_sliding_window_array[105])
    print "PTR=%s" % PTR_median
    print "The genomic locations for peak and trough values: %s, %s" % (median_sliding_window_array.index(peak), median_sliding_window_array.index(through))
    out_file_ptr = args.out + "_PTR"
    with open(out_file_ptr, 'a') as out:
        out.write(str(args.out)+' median_sliding_window_array :\t'+str(PTR_median)+'\n')
    with open(args.out, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')

    stats_file = args.bed.replace('_coverage.bed', '_alignment_stats')
    out_file = args.out + "_median_bins_normalize"
    cmd = "grep \'mapped (\' " + stats_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    out_split = out.split(' ')
    mapped_reads = out_split[0]
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
            count += 1
            ratio = int(i)/int(mapped_reads)
            out.write(str(count)+','+str(ratio)+'\n')
    out.close()





    calculate_per_of_reads(median_sliding_window_array)
    #print moving_array_length

def generate_coverage_normalized():
    partitions = [read_counts[i:i+100000] for i in xrange(0, len(read_counts), 100000)]
    partitions = partitions if len(partitions[-1]) == 100000 else partitions[:-1]
    stats_file = args.bed.replace('_coverage.bed', '_alignment_stats')
    cmd = "grep \'mapped (\' " + stats_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    out_split = out.split(' ')
    mapped_reads = out_split[0]
    normalize_array = []
    #print mapped_reads
    for x in partitions:
        sum_of_bin = sum(x)
        ratio = sum_of_bin/int(mapped_reads)
        normalize_array.append(ratio)
    coverage_out_file = args.bed + "_coverage_normalize"
    with open(coverage_out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in normalize_array:
            count += 1
            out.write(str(count)+','+str(i)+'\n')

#smoothing()
generate_coverage_normalized()

total_time = time.time() - start
print "Total time: %s seconds" % total_time

















########################################################################################################################
##################################################Deprecated############################################################
# def generate_median_stats():
#     # Partition/divide the read count array into bins of size 'n'.
#     # The below command sets the value of i by dividing the read counts array each time starting with 0 and increasing to step 'n'.
#     # So if n=1000 then the first bin will contain values from 0 to 1000 index of read counts array.
#     partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
#     # This below line creates a seperate array that was left out from partition due to lack of range
#     partitions = partitions if len(partitions[-1]) == n else partitions[:-1]
#     median_array = []
#     for x in partitions:
#         median = statistics.median(x)
#         median_array.append(median)
#
#     with open(args.out, 'w') as out:
#         header = "bin\tcount\n"
#         count = 0
#         out.write(header)
#         for i in median_array:
#             count += 1
#             out.write(str(count)+'\t'+str(i)+'\n')
# #generate_median_stats()
#
# def generate_mean_stats():
#     partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
#     partitions = partitions if len(partitions[-1]) == n else partitions[:-1]
#     mean_array = []
#     for x in partitions:
#         mean = statistics.mean(x)
#         mean_array.append(mean)
#     with open(args.out, 'w') as out:
#         header = "bin\tcount\n"
#         count = 0
#         out.write(header)
#         for i in mean_array:
#             count += 1
#             out.write(str(count)+'\t'+str(i)+'\n')
# #generate_mean_stats()
#
# def generate_percentage_stats():
#     partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
#     partitions = partitions if len(partitions[-1]) == n else partitions[:-1]
#     total = 0
#     for x in partitions:
#         total = total + sum(x)
#     #print total
# def generate_coverage():
#     partitions = [read_counts[i:i+13079] for i in xrange(0, len(read_counts), 13079)]
#     partitions = partitions if len(partitions[-1]) == 13079 else partitions[:-1]
#
#     mean_array = []
#
#     for x in partitions:
#         mean = statistics.mean(x)
#         mean_array.append(mean)
#     coverage_out_file = args.bed + "_coverage"
#     with open(coverage_out_file, 'w') as out:
#         header = "bin,count\n"
#         count = 0
#         out.write(header)
#         for i in mean_array:
#             count += 1
#             out.write(str(count)+','+str(i)+'\n')
# #generate_coverage()