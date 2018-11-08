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
import glob
import readline
import pandas as pd
import errno
import gzip
import random
from datetime import datetime

parser = argparse.ArgumentParser(description='Downsample the existing raw sequencing data into specified coverage(genome size given) or number of reads specified')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-fastq_dir', action='store', dest="fastq_dir",
                    help='Directory where all the original fastq files are located.')
required.add_argument('-fastq_files', action='store', dest="fastq_files",
                    help='Name of fastq files(one per line) in fastq directory.')
required.add_argument('-type', action='store', dest="type",
                    help='Type of reads. PE/SE')
required.add_argument('-out_folder', action='store', dest="out_folder",
                    help='Output folder to save output')
optional.add_argument('-coverage', action='store', dest="coverage",
                    help='Downsampling the original data to this coverage')
optional.add_argument('-genome_size', action='store', dest="genome_size",
                    help='Genome size')
optional.add_argument('-coverage_reads', action='store', dest="coverage_reads",
                    help='Downsampling the original data to this many reads')
optional.add_argument('-jobrun', action='store', dest="jobrun",
                    help='Running a job on Cluster, Running Parallel jobs, Run jobs/commands locally (default): cluster, local, parallel-local, parallel-single-cluster')
optional.add_argument('-cluster_type', action='store', dest="cluster_type",
                    help='Type of Cluster: torque, pbs, sgd')
optional.add_argument('-cluster_resources', action='store', dest="cluster_resources",
                    help='Cluster Resources to use. for example nodes,core. Ex: 1,4')
optional.add_argument('-numcores', action='store', dest="numcores",
                    help='Number of cores to use on local system for parallel-local parameter')
optional.add_argument('-remove_temp', action='store', dest="remove_temp",
                    help='Remove Temporary files generated during the run')
args = parser.parse_args()


def coverage(filenames_array, awk, output_folder, type, samples, size):
    temp_forward_coverage = output_folder + "/temp_forward.txt"
    temp_reverse_coverage = output_folder + "/temp_reverse.txt"
    f1=open(temp_forward_coverage, 'w+')
    f2=open(temp_reverse_coverage, 'w+')
    file_coverage = {}
    if type == "PE":
        for file in filenames_array:
            if file.endswith('.gz'):
                coverage_cmd_forward = "zcat %s | %s" % (file, awk)
                coverage_cmd_reverse = "zcat %s | %s" % (file.replace('_R1_', '_R2_'), awk)
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                coverage_msg_reverse = "Calculating coverage for file: %s\n" % file.replace('_R1_', '_R2_')
            else:
                coverage_cmd_forward = "cat %s | %s" % (file, awk)
                coverage_cmd_reverse = "cat %s | %s" % (file.replace('_R1_', '_R2_'), awk)
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                coverage_msg_reverse = "Calculating coverage for file: %s\n" % file.replace('_R1_', '_R2_')

            proc = subprocess.Popen([coverage_cmd_forward], stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            f1.write(out2)
            proc2 = subprocess.Popen([coverage_cmd_reverse], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc2.communicate()
            f2.write(out)
        temp_final_file = "%s/temp_final.txt" % output_folder
        sample_file = "%s/%s" % (output_folder, samples)
        with open(temp_final_file, 'w') as res, open(sample_file) as f1, open(temp_forward_coverage) as f2, open(temp_reverse_coverage) as f3:
            for line1, line2, line3 in zip(f1, f2, f3):
                res.write("{}\t{}\t{}\n".format(line1.rstrip(), line2.rstrip(), line3.rstrip()))
        final_coverage_file = "%s/Final_Coverage.txt" % output_folder
        f3=open(final_coverage_file, 'w+')
        header = "Sample_name\tForward_reads\tForward_unique_reads\tPerc_Forward_unique_reads\tMost_abundant_Forward_read(adapter)\tMost_abundant_Forward_read_count\tPerc_Most_abundant_Forward_read\tAvg_Forward_readLength\tReverse_reads\tReverse_unique_reads\tPerc_Reverse_unique_reads\tMost_abundant_Reverse_read(adapter)\tMost_abundant_Reverse_read_count\tPerc_Most_abundant_Reverse_read\tAvg_Reverse_readLength\tCoverage\n"
        f3.write(header)
        with open(temp_final_file, 'a+') as fp:
            for line in fp:
                line = line.strip()
                line_split = line.split('\t')
                avg_read_length = statistics.mean([float(line_split[7]), float(line_split[14])])
                final_coverage = (float(line_split[1]) * 2 * avg_read_length) / float(size)
                print "The coverage for %s: %s" % (line_split[0], final_coverage)
                file_coverage[line_split[0]] = float(final_coverage)
                print_string = line + "\t" + str(final_coverage) + "\n"
                f3.write(print_string)
        f3.close()

        return file_coverage

    elif type == "SE":
        for file in filenames_array:
            coverage_msg_forward = "Calculating coverage for file: %s\n" % file
            if file.endswith('.gz'):
                coverage_cmd_forward = "zcat %s | %s" % (file, awk)
            else:
                coverage_cmd_forward = "cat %s | %s" % (file, awk)

            proc = subprocess.Popen([coverage_cmd_forward], stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            f1.write(out2)
        temp_final_file = "%s/temp_final.txt" % output_folder
        sample_file = "%s/%s" % (output_folder, os.path.basename(samples))
        with open(temp_final_file, 'w') as res, open(sample_file) as f1, open(temp_forward_coverage) as f2:
            for line1, line2 in zip(f1, f2):
                res.write("{}\t{}\n".format(line1.rstrip(), line2.rstrip()))
        final_coverage_file = "%s/Final_Coverage.txt" % output_folder
        f3=open(final_coverage_file, 'w+')
        header = "Sample_name\tForward_reads\tForward_unique_reads\tPerc_Forward_unique_reads\tMost_abundant_Forward_read(adapter)\tMost_abundant_Forward_read_count\tPerc_Most_abundant_Forward_read\tAvg_Forward_readLength\tCoverage\n"
        f3.write(header)
        with open(temp_final_file, 'a+') as fp:
            for line in fp:
                line = line.strip()
                line_split = line.split('\t')
                avg_read_length = float(line_split[7])
                final_coverage = (float(line_split[1]) * avg_read_length) / float(size)
                print "The coverage for %s: %s" % (line_split[0], float(final_coverage))
                file_coverage[line_split[0]] = float(final_coverage)
                print_string = line + "\t" + str(final_coverage) + "\n"
                f3.write(print_string)
        f3.close()

        return file_coverage

def downsample_coverage(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage):
    # proc = subprocess.Popen(["zcat %s | wc -l" % files], stdout=subprocess.PIPE, shell=True)
    # (out, err) = proc.communicate()
    # records = int(out) / 4
    def random_coverage_sampling(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage):
        if type == "SE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | %s" % (files, awk)], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                out_split = out.split('\t')
                records = float(out_split[0].strip())
                length = 52
                ori_coverage = float(file_coverage[os.path.basename(files)])
                #downsampled_reads = (records * coverage) / ori_coverage
                downsampled_reads = float((coverage * size) / length)
                rand_records = []
                c = float(downsampled_reads % 4)
                if c != 0:
                    downsampled_reads = downsampled_reads + c
                print "The downsampled reads for %s: %s" % (files, float(downsampled_reads))
                while len(rand_records) < float(downsampled_reads): #int(downsampled_reads):
                    r = random.randint(0,records-1)
                    if r not in rand_records:
                        rand_records.append(r)
                rand_records = sorted(rand_records)
                # outfile = output_folder + os.path.basename(files).replace('.gz', '') + ".subset"
                # fha = gzip.open(files)
                # suba = open(output_folder + os.path.basename(files).replace('.gz', '') + ".subset", "w")

                outfile = output_folder + os.path.basename(files) + ".subset"
                fha = open(files)
                suba = open(output_folder + os.path.basename(files) + ".subset", "w")
                rec_no = -1
                written = 0
                print len(rand_records)
                print len(set(rand_records))
                for rr in rand_records:
                    while rec_no < rr:
                        for i in range(4):
                            fha.readline()
                        rec_no += 1
                    for i in range(4):
                        suba.write(fha.readline())
                    rec_no += 1
                    written += 1
                print written
                assert written == float(downsampled_reads) #int(downsampled_reads)






                gzip_cmd = "gzip -f %s" % outfile
                os.system(gzip_cmd)
        ###Pending
        elif type == "PE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | wc -l" % files], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                records = int(out) / 4
                rand_records = []
                while len(rand_records) < 1000:
                    r = random.randint(0,records-1)
                    if r not in rand_records:
                        rand_records.append(r)
                rand_records = sorted(rand_records)
                #rand_records = sorted([random.randint(0, int(records) - 1) for _ in xrange(10000)])
                #print rand_records
                outfile = args.out_folder + files.replace('.gz', '')
                fha = gzip.open(files)
                suba = open(files.replace('.gz', '') + "_random_subset.fastq", "w")
                rec_no = -1
                written = 0
                for rr in rand_records:
                    while rec_no < rr:
                        for i in range(4):
                            fha.readline()
                        rec_no += 1
                    for i in range(4):
                        suba.write(fha.readline())
                    rec_no += 1
                    written += 1
                assert written == 1000
                gzip_cmd = "gzip -f %s" % files.replace('.gz', '') + "_random_subset.fastq"
                os.system(gzip_cmd)

    def get_head(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage):
        if type == "SE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | %s" % (files, awk)], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                out_split = out.split('\t')
                records = float(out_split[0].strip())
                length = 52
                ori_coverage = float(file_coverage[os.path.basename(files)])
                #downsampled_reads = (records * coverage) / ori_coverage
                downsampled_reads = float((coverage * size) / length)
                c = int(downsampled_reads % 4)
                print c
                if c != 0:
                    downsampled_reads = downsampled_reads + c
                print downsampled_reads

                downsampled_reads = int(downsampled_reads) * 4
                print "The downsampled reads for %s: %s" % (files, int(downsampled_reads))
                out_file = output_folder + os.path.basename(files).replace('_mapped.fastq', '_head_subset.fastq')
                head_cmd = "head -n%s %s > %s" %(int(downsampled_reads), files, out_file)
                print head_cmd
                os.system(head_cmd)
                gzip_cmd = "gzip -f %s" % out_file
                os.system(gzip_cmd)

        ###Pending
        elif type == "PE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | wc -l" % files], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()


    def seqtk(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage):
        if type == "SE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | %s" % (files, awk)], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                out_split = out.split('\t')
                records = float(out_split[0].strip())
                length = 52
                ori_coverage = float(file_coverage[os.path.basename(files)])
                #downsampled_reads = (records * coverage) / ori_coverage
                downsampled_reads = float((coverage * size) / length)
                c = float(downsampled_reads % 4)
                if c != 0:
                    downsampled_reads = downsampled_reads + c
                print "The downsampled reads for %s: %s" % (files, float(downsampled_reads))
                out_file = output_folder + os.path.basename(files).replace('_mapped.fastq', '_seqtk_subset.fastq')
                seqtk_cmd = "/nfs/esnitkin/bin_group/seqtk/seqtk sample -s100 %s %s > %s" %(files, float(downsampled_reads), out_file)
                print seqtk_cmd
                os.system(seqtk_cmd)
                gzip_cmd = "gzip -f %s" % out_file
                os.system(gzip_cmd)

        ###Pending
        elif type == "PE":
            for files in filenames_array:
                proc = subprocess.Popen(["cat %s | wc -l" % files], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()


    #random_coverage_sampling(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage)
    get_head(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage)
    #seqtk(filenames_array, output_folder, type, samples, size, coverage, awk, file_coverage)








# Make sure the output folder exists or create at given path
def make_sure_path_exists(out_path):
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

    print "\nThe Script started at: %s\n" % start_time

    print "\nThe Script: coverage_sampling.py will down-sample the existing raw sequencing data into specified coverage(genome size given) or number of reads specified:\n\n" \
          "1. Calculate the coverage of original data and Downsample it to specified coverage. requires genome size parameter\n" \
          "2. Downsample the original data to specified number of reads. requires coverage reads parameter\n" \

    """ Gather fastq filenames from the directory """
    global fastq_filenames
    fastq_filenames = []
    with open(args.fastq_files) as fp:
        for line in fp:
            line = line.strip()
            line = args.fastq_dir + line
            fastq_filenames.append(line)
        fp.close()
    if args.out_folder != '':
        args.out_folder += '/'
        make_sure_path_exists(args.out_folder)

    cp_cmd = "cp %s %s" % (args.fastq_files, args.out_folder)
    os.system(cp_cmd)

    """ Calculate Original file coverage """
    awk = "awk \'BEGIN{OFS=\"\t\"};((NR-2)%4==0){read=$1;total++;count[read]++;len+=length(read)}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total,len/total}\'"
    if args.genome_size:
        file_coverage = coverage(fastq_filenames, awk, args.out_folder, args.type, args.fastq_files, args.genome_size)
        print file_coverage
    if args.coverage:
        out_coverage = float(args.coverage)
        downsample_coverage(fastq_filenames, args.out_folder, args.type, args.fastq_files, int(args.genome_size), float(args.coverage), awk, file_coverage)
    # elif args.coverage_reads:
    #     downsample_reads()

    time_taken = datetime.now() - start_time_2
    if args.remove_temp:
        del_command = "rm -r %s" % temp_dir
        os.system(del_command)
