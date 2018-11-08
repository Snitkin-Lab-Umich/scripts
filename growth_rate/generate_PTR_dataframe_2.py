from __future__ import division
__author__ = 'alipirani'
import statistics
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np

parser = argparse.ArgumentParser(description='Read BED file and generate PTR dataframe')
parser.add_argument('-bed', action='store', dest="bed", help='bedfile')
parser.add_argument('-out', action='store', dest="out", help='output file')
parser.add_argument('-bin', action='store', dest="bin", help='Number of Bins. It represents how many values should be there in a bin.')
parser.add_argument('-total_aligned_reads', action='store', dest="total", help='Total Aligned Reads')
parser.add_argument('-window_bin', action='store', dest="window_bin", help='Size of Window bin')
args = parser.parse_args()
window = int(args.window_bin)
#n = int(args.bin)
read_counts = []
print args.bed
with open(args.bed) as fp:
    for line in fp:
        line = line.strip()
        line_split = line.split('\t')
        counts = int(line_split[2])
        read_counts.append(int(counts))





def generate_median_stats():
    # Partition/divide the read count array into bins of size 'n'.
    # The below command sets the value of i by dividing the read counts array each time starting with 0 and increasing to step 'n'.
    # So if n=1000 then the first bin will contain values from 0 to 1000 index of read counts array.
    partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
    # This below line creates a seperate array that was left out from partition due to lack of range
    partitions = partitions if len(partitions[-1]) == n else partitions[:-1]

    median_array = []

    for x in partitions:
        median = statistics.median(x)
        median_array.append(median)

    with open(args.out, 'w') as out:
        header = "bin\tcount\n"
        count = 0
        out.write(header)
        for i in median_array:
            count += 1
            out.write(str(count)+'\t'+str(i)+'\n')

#generate_median_stats()

def generate_mean_stats():
    partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
    partitions = partitions if len(partitions[-1]) == n else partitions[:-1]

    mean_array = []

    for x in partitions:
        mean = statistics.mean(x)
        mean_array.append(mean)

    with open(args.out, 'w') as out:
        header = "bin\tcount\n"
        count = 0
        out.write(header)
        for i in mean_array:
            count += 1
            out.write(str(count)+'\t'+str(i)+'\n')


#generate_mean_stats()



def generate_percentage_stats():
    partitions = [read_counts[i:i+n] for i in xrange(0, len(read_counts), n)]
    partitions = partitions if len(partitions[-1]) == n else partitions[:-1]
    total = 0
    for x in partitions:
        total = total + sum(x)

    #print total

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
    #through = np.min(moving_sum_array[np.nonzero(moving_sum_array)])
    through = min([x for x in moving_sum_array if x !=0])
    PTR_moving = peak/through
    print "The peak and trough values for moving_sum_array is: %s; %s" % (peak, through)
    #print through
    out_file_ptr = args.out + "_PTR"
    with open(out_file_ptr, 'w+') as out:
        out.write(str(args.out)+' moving_sum_array :\t'+str(PTR_moving)+'\n')





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








def smoothing_1():
    raw_count_sliding_window_array = []
    zero_bins = 0
    less_median_coverage = 0
    for i in xrange(0, len(read_counts), 100):
        start = i
        end = i + window
        #print str(start) + "\t" + str(end)
        #if len(read_counts[start:end]) == 10000:
        #print read_counts[start:end].count(0)
        # print statistics.median(read_counts[start:end])
        if read_counts[start:end].count(0) == window:
            zero_bins += 1
        # if float(statistics.median(read_counts[start:end])) < 0.05:
        #     less_median_coverage += 1
            #print statistics.median(read_counts[start:end])
        raw_count_sliding_window_array.append(read_counts[start:end])
    print "The length of raw_count_sliding_window_array is %s" % len(raw_count_sliding_window_array)
    print "The number of bins with no mapped reads: %s" % str(zero_bins)
    # print "less median coverage: %s" % str(less_median_coverage)
    moving_sum_array = []
    for i in raw_count_sliding_window_array:
        moving_sum_array.append(sum(i))
    #print len(moving_sum_array)
    generate_moving_sum_results(moving_sum_array)





    median_sliding_window_array = []
    diff = 0
    for i in xrange(0, len(read_counts), 100):
        start = i
        end = i + window
        #print str(start) + "\t" + str(end) + "\n"
        #print len(moving_sum_array[start:end])
        #print moving_sum_array[start:end]
        #print len(moving_sum_array[start:end])
        if len(moving_sum_array[start:end]) > 5000:
            #print len(moving_sum_array[start:end])

            #value = statistics.median(moving_sum_array[start:end]) - diff
            median_sliding_window_array.append(statistics.median(moving_sum_array[start:end]))
            #diff = statistics.median(moving_sum_array[start:end])


    print len(median_sliding_window_array)
    #print median_sliding_window_array
    peak = max(median_sliding_window_array)
    #through = np.min(median_sliding_window_array[np.nonzero(median_sliding_window_array)])
    through = min([x for x in median_sliding_window_array if x !=0])
    PTR_median = peak/through
    print "The peak and trough values for median_sliding_window_array is: %s; %s" % (peak, through)
    print "The genomic location for peak values: %s" % (median_sliding_window_array.index(peak))
    print "The genomic location for trough values: %s" % (median_sliding_window_array.index(through))
    #print through
    #print median_sliding_window_array[117]
    out_file_ptr = args.out + "_PTR"
    with open(out_file_ptr, 'a') as out:
        out.write(str(args.out)+' median_sliding_window_array :\t'+str(PTR_median)+'\n')
        out.write(str(args.out)+' Peak and trough location:\t' + str(median_sliding_window_array.index(peak)) + '\t' + str(median_sliding_window_array.index(through)) + '\n')





    with open(args.out, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')

    calculate_per_of_reads(median_sliding_window_array)
smoothing_1()






#generate_percentage_stats()






def generate_coverage():
    partitions = [read_counts[i:i+13079] for i in xrange(0, len(read_counts), 13079)]
    partitions = partitions if len(partitions[-1]) == 13079 else partitions[:-1]

    mean_array = []

    for x in partitions:
        mean = statistics.mean(x)
        mean_array.append(mean)
    coverage_out_file = args.bed + "_coverage"
    with open(coverage_out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in mean_array:
            count += 1
            out.write(str(count)+','+str(i)+'\n')




generate_coverage()





