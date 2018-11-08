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
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict

parser = argparse.ArgumentParser(description='Recombination analysis pipeline.')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filename', action='store', dest="filename", help='filename of pseudomolecule assembly fasta file', required=True)
optional.add_argument('-filename_db', action='store', dest="filename_db", help='filenames of nucmer db', required=False)
required.add_argument('-out', action='store', dest="out", help='out_directory', required=True)
required.add_argument('-dir', action='store', dest="dir", help='directory of fasta files', required=True)
optional.add_argument('-dir_db', action='store', dest="dir_db", help='directory of nucmer db fasta files', required=False)
required.add_argument('-prokka_dir', action='store', dest="prokka_dir", help='directory of prokka annotations', required=True)
required.add_argument('-jobrun', action='store', dest="jobrun",
                    help='Running a job on Cluster, Running Parallel jobs, Run jobs/commands locally (default): cluster, local, parallel-local, parallel-single-cluster')
required.add_argument('-steps', action='store', dest="steps",
                    help='Analysis Steps to be performed.'
                         '   Step 1: Run Nucmer on all samples. All-vs-All.'
                         '   Step 2: Parse the nucmer generated coords files and generate annotated results of aligned fragments'
                         '   Step 3: ')
#required.add_argument('-db', action='store', dest="db", help='Nucmer reference genome to be used as pseudomolecule.', required=True)
args = parser.parse_args()


def create_job(jobrun, commands_list):

    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        job_directory = args.out + "/" + "temp_jobs"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = job_directory + "/nucmer_job_command_*.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_file = commands_list
        print len(command_file)
        # if args.numcores:
        #     num_cores = int(num_cores)
        # else:
        #     num_cores = multiprocessing.cpu_count()
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)

    elif jobrun == "parallel-single-cluster":
        print "  "
    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            os.system(i)

def create_job_parse(jobrun, commands_list):

    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        job_directory = args.out + "/" + "temp_jobs"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_coordinates_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = job_directory + "/nucmer_coordinates_job_command_*.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_file = commands_list
        print len(command_file)
        # if args.numcores:
        #     num_cores = int(num_cores)
        # else:
        #     num_cores = multiprocessing.cpu_count()
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)

    elif jobrun == "parallel-single-cluster":
        print "  "
    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            os.system(i)

def make_sure_path_exists(out_path):
    """
    Make sure the output folder exists or create at given path
    :param out_path:
    :return:
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()

def run_nucmer(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
        os.system(mkdir_command)
        for file_again in filenames_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=%s %s %s && show-coords -r -c -l -T %s.delta > %s.coords && mv %s.delta %s.coords %s/%s\n" % (prefix, file, file_again, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
            command_array.append(command)
            f.write(command)
    f.close()
    return command_array

def run_nucmer_db(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
        os.system(mkdir_command)
        for file_again in filenames_db_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + ":" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=%s %s %s && show-coords -r -c -l -T -I 99 -L 2000 %s.delta > %s.coords && mv %s.delta %s.coords %s/%s\n" % (prefix, file, file_again, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
            command_array.append(command)
            f.write(command)
    f.close()
    return command_array


# def run_nucmer_db(temp_cmd):
#     f = open(temp_cmd, 'w+')
#     command_array = []
#     for file in filenames_array:
#         filebase = os.path.basename(file)
#         mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
#         os.system(mkdir_command)
#         prefix = filebase.replace('.fasta', '')
#         command = "nucmer --maxmatch --prefix=%s %s %s && show-coords -r -c -l -T -I 99 -L 2000 %s.delta > %s.coords && mv %s.delta %s.coords %s/%s\n" % (prefix, file, args.db, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
#         command_array.append(command)
#         f.write(command)
#     f.close()
#     return command_array



def run_command(i):
    os.system(i)
    done = "done"
    return done


def generate_parse_coord_aggregate_jobs(jobrun, filenames_array, temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = args.out + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates.py -folder %s -out %s -dir %s -prokka_dir %s\n" % (folder, args.out, args.dir, args.prokka_dir)
        #print cmd
        command_array.append(cmd)
        f.write(cmd)
        #job_name = os.path.basename(folder)
        #job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
        #job_file_name = "%s/%s_run_nucmer_coordinates.pbs" % (job_directory, job_name)
        #f1=open(job_file_name, 'w+')
        #f1.write(job_print_string)
    #f1.close()
    f.close()
    print len(command_array)
    return command_array



def generate_parse_coord_db_aggregate_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = args.out + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates_db.py -folder %s -out %s -dir %s -prokka_dir %s\n" % (folder, args.out, args.dir, args.prokka_dir)
        #print cmd
        command_array.append(cmd)
        f.write(cmd)
        #job_name = os.path.basename(folder)
        #job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
        #job_file_name = "%s/%s_run_nucmer_coordinates.pbs" % (job_directory, job_name)
        #f1=open(job_file_name, 'w+')
        #f1.write(job_print_string)
    #f1.close()
    f.close()
    print len(command_array)
    return command_array
#generate_parse_coord_aggregate_jobs()


#Main Steps
if __name__ == '__main__':

    """Start Timer"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()

    print "\nThe Script started at: %s\n" % start_time

    print "\nThe Script: recombination_analysis.py will run the following steps:\n\n" \
          "1. Runs nucmer to align each genome against each other provided with filename argument\n" \
          "2. Parse the coordinates file generated by nucmer and generates an annotation file containing annotations of aligned region.\n" \
          "3. Calculates the alignment score which is defined as number of bases aligned per the genome length.\n"

    # GENERATE FASTA FILENAMES ARRAY
    filenames_array = []
    with open(args.filename) as fp:
        for line in fp:
            line = line.strip()
            line = args.dir + "/" + line
            filenames_array.append(line)

    # Run pipeline steps
    if "1" in args.steps:
        make_sure_path_exists(args.out)
        temp_cmd = "%s/temp_commands" % args.out
        command_array = run_nucmer(temp_cmd)

        if args.jobrun:
            create_job(args.jobrun, command_array)

    if "2" in args.steps:
        temp_cmd = "%s/temp_commands_parse_coord" % args.out
        command_array = generate_parse_coord_aggregate_jobs(args.jobrun, filenames_array, temp_cmd)
        if args.jobrun:
            create_job_parse(args.jobrun, command_array)


    if "3" in args.steps:
        print "The Results output directory containing species-pair folders is %s" % args.out
        list_command = "ls -d %s/*_* | rev | cut -d'/' -f1 | rev" % args.out
        species_pair_uniq_array = []
        species_array = []
        #print list_command
        list_of_species_pai_dir = []
        list_of_species_pai_dir = subprocess.check_output(list_command, shell=True)
        for i in list_of_species_pai_dir.split('\n'):
            species_name = (i.split('_'))[0]
            if species_name not in species_array and species_name != "":
                species_array.append(species_name)
        os.system("mkdir %s/redundant_species_pairs" % args.out)
        f1=open("%s/move_redundant.sh" % args.out, 'w+')
        for file in species_array:
            for file_again in species_array:
                if file != file_again:
                    generate_string = str(file) + "_" + str(file_again)
                    reverse_string = str(file_again) + "_" + str(file)
                    if reverse_string not in species_pair_uniq_array:
                        species_pair_uniq_array.append(generate_string)
                        move = "mv %s/%s %s/redundant_species_pairs/" % (args.out, reverse_string, args.out)
                        f1.write(move + '\n')
        f1.close()
    if "4" in args.steps:
        print "Comparing each fasta file with the nucmer reference database generated from extracted aligned fragments"
        make_sure_path_exists(args.out)
        temp_cmd = "%s/temp_commands" % args.out
        # GENERATE FASTA FILENAMES ARRAY
        filenames_db_array = []
        with open(args.filename_db) as fp:
            for line in fp:
                line = line.strip()
                line = args.dir_db + "/" + line
                filenames_db_array.append(line)



        command_array = run_nucmer_db(temp_cmd)
        # print command_array[20]
        # print len(command_array)
        if args.jobrun:
            create_job(args.jobrun, command_array)

    if "5" in args.steps:
        temp_cmd = "%s/temp_commands_parse_coord_db" % args.out
        filenames_db_array = []
        with open(args.filename_db) as fp:
            for line in fp:
                line = line.strip()
                line = args.dir_db + "/" + line
                filenames_db_array.append(line)
        command_array = generate_parse_coord_db_aggregate_jobs(args.jobrun, filenames_array, filenames_db_array, temp_cmd)
        if args.jobrun:
            create_job_parse(args.jobrun, command_array)