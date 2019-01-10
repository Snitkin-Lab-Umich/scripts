__author__ = 'alipirani'

import os
import argparse
import re
import glob
parser = argparse.ArgumentParser(description='Generate PBS scripts for Different Pipelines for multiple samples: varcall/assembly/new_assembly/rna/ptr')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fastq files')
parser.add_argument('-filenames', action='store', dest="filenames", help='These file should contain name of forward fastq files. \
One file per line. \
These can be obtained by running: ls *_R1_001.fastq.gz > filenames \
Make Sure your Forward reads extension ends with \"_R1_001.fastq.gz\" ')
parser.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
parser.add_argument('-pipeline', action='store', dest="pipeline", help='Generating Jobs for which pipeline? varcall/assembly/new_assembly/new_assembly_SE/rna/ptr/ariba')
parser.add_argument('-reference', action='store', dest="reference", help='Reference Genome to be used for pipeline')
parser.add_argument('-type', action='store', dest="type", help='Type of Fastq files: PE or SE')
parser.add_argument('-ariba_db', action='store', dest="ariba_db", help='Path to Ariba Database')
parser.add_argument('-reference_key', action='store', dest="reference_key", help='Reference Genome to be used for each sample. \
Provide a tab-seperated file with fastq filename in first column and reference genome name in corresponding second column')
parser.add_argument('-depthofcoverage', action='store', dest="depthofcoverage", help='Varcall pipeline until GATK Depth of Coverage')
parser.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All : This will run all the steps starting from cleaning the reads to variant calling;\n'
                                                                     '2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. \nYou can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.\nThe order is required to be sequential. Also, while skipping any of the step make sure you have results already present in your output folder.\n'
                                                                     '3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps')
parser.add_argument('-ariba', action='store', dest="ariba", help='Run Ariba to find AMR / MLST / both: possible options are AMR, MLST, both')
parser.add_argument('-email', action='store', dest="email", help='email for PBS job notifications')
parser.add_argument('-resources', action='store', dest="resources", help='PBS resources for jobs')
parser.add_argument('-assembly', action='store', dest="assembly", help='assembly: wga/plasmid/both. Default is both ')
parser.add_argument('-suffix', action='store', dest="suffix", help='suffix for fastq read files ')
args = parser.parse_args()

# # read filenames file
# filenames_array = []
# with open(args.filenames) as fp:
#     for line in fp:
#         line = line.strip()
#         line = args.dir + "/" + line
#         filenames_array.append(line)

def get_filenames(dir, type, filenames, suffix):
    if not filenames:
        if not suffix:
            suffix = ".fastq.gz"
        try:
            list_of_files = glob.glob("%s/*%s" % (dir, suffix))
            if len(list_of_files) < 1:
                print "No fastq files with suffix %s found in reads directory %s" % (suffix, dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                keep_logging('Error while listing files in reads directory.', 'Error while listing files in reads directory.', logger, 'exception')
                exit()
    else:
        list_of_files = []
        with open(filenames) as fp:
            for line in fp:
                line = line.strip()
                line = dir + "/" + line
                list_of_files.append(line)
    return list_of_files

# Function to create RNA Seq pipeline pbs scripts for each fastq files
def create_RNA_seq_jobs():
    # read the reference key file: Filename -> reference genome
    reference_key_dict = {}
    fp = open(args.reference_key, 'r')
    for line in fp:
        line = line.strip()
        line_split = line.split('\t')
        line_split_file = args.dir + "/" + line_split[0]
        reference_key_dict[line_split_file] = line_split[1]

    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
    "#PBS -m a\n" \
    "#PBS -V\n" \
    "#PBS -l %s\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in reference_key_dict.keys():
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('_')
        first_part = str(first_part_split[0]) + "_"

        #first_file = file
        #second_part = filename_base.replace("_1_", "_2_")
        #second_part = filename_base.replace("_R1", "_R2")
        #second_part = filename_base.replace("forward", "reverse")
        #second_part = filename_base.replace("1_sequence", "2_sequence")
        #second_file = args.dir + "/" + second_part

        # Change these directory path to where your pipeline is located
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/RNA_seq_pipeline/"
        command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -s yes" % (file, args.out_dir, first_part, first_part, reference_key_dict[file])
        job_name = "./" + first_part + ".pbs"

        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create (old)assembly pipeline pbs scripts for each fastq files
def create_assembly_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        # Forward reads file name and get analysis name from its name
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('.')
        first_part = first_part_split[0]
        first_file = file

        second_part = filename_base.replace("_R1_", "_R2_")
        #second_part = filename_base.replace("forward", "reverse")
        #second_part = filename_base.replace("1_sequence", "2_sequence")
        second_file = args.dir + "/" + second_part

        # Change these directory path to where your pipeline code is located
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/assembly_umich/"
        command = "/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s" % (first_file, second_file, args.out_dir, first_part, args.reference)
        job_name = "./" + first_part + ".pbs"

        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create Variant Calling pipeline pbs scripts for each fastq files
def create_varcall_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
    "#PBS -m a\n" \
    "#PBS -V\n" \
    "#PBS -l %s\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        # Forward reads file name and get analysis name from its name
        filename_base = os.path.basename(file)
        #first_part_split = filename_base.split('_R1_001.fastq.gz')
        # first_part_split = filename_base.split('_R1.fastq.gz')
        # first_part = first_part_split[0].replace('_L001', '')
        first_file = file

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        #second_part = filename_base.replace("_R1", "_R2")
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
        elif ".1.fastq.gz" in filename_base:
            second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
            first_part_split = filename_base.split('.1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
	    first_part = re.sub("_S.*_", "", first_part)
	elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        else:
                print "Using Standard second file naming convention"
                second_part = filename_base.replace("_R1_", "_R2_")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
		first_part = re.sub("_S.*_", "", first_part)
        second_file = args.dir + "/" + second_part

        # Change these directory path to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/varcall_umich"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/varcall_umich/config"
        job_name = "./" + first_part + ".pbs"

        if args.depthofcoverage:
            command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -coverage_depth_stats yes" % (first_file, second_file, args.out_dir, first_part, first_part, args.reference, config)
        elif args.type == "SE":
	    command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s" % (first_file, args.out_dir, first_part, first_part, args.reference, config, args.steps)
        else:    
            command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s" % (first_file, second_file, args.out_dir, first_part, first_part, args.reference, config, args.steps)

        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create PTR pipeline pbs scripts for each fastq files
def create_PTR_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
    "#PBS -m a\n" \
    "#PBS -V\n" \
    "#PBS -l %s\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        # Forward reads file name and get analysis name from its name
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('.')
        first_part = first_part_split[0]
        first_file = file

        # Change these directory path to where your pipeline code is located.
        #cd_command = "cd /nfs/esnitkin/bin_group/pipeline/PTR_analysis/"
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/Growth-rate-analysis/"
        #command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -SE %s -o %s/%s -analysis %s -index CFT073" % (first_file, args.out_dir, first_part, first_part)
        command = "/nfs/esnitkin/bin_group/anaconda2/bin/python pipeline.py -type SE -PE1 %s -o %s/%s -analysis %s -index %s" % (first_file, args.out_dir, first_part, first_part, args.reference)
        job_name = "./" + first_part + ".pbs"

        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_jobs(list_of_files):
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in list_of_files:
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/config"
        job_name = "%s" % args.out_dir + first_part + ".pbs"
        if args.assembly:
            assembly_para = "-assembly " + args.assembly
        else:
            assembly_para = "-assembly both"
        if args.reference:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, args.reference, first_part, config, assembly_para)
        else:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, first_part, config, assembly_para)
        if args.ariba:
            ariba_para = "-ariba " + args.ariba
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s %s" % (first_file, second_file, args.out_dir, first_part, args.reference, first_part, config, assembly_para, ariba_para)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_SE_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_unclassified.fastq.gz" in filename_base:
            second_part = filename_base.replace("_unclassified.fastq.gz", "_unclassified.fastq.gz")
            first_part_split = filename_base.split('_unclassified.fastq.gz')
            first_part = first_part_split[0].replace('_L00', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/config"
        job_name = "./" + first_part + ".pbs"

        if args.reference:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type SE -analysis %s -config %s" % (first_file, args.out_dir, first_part, args.reference, first_part, config)
        else:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type SE -analysis %s -config %s" % (first_file, args.out_dir, first_part, first_part, config)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')



# Function to create assembly pbs scripts for each fastq files
def create_a5_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        mkdir_cmd = "mkdir %s/%s" % (args.out_dir, first_part)
        cd_command = "cd %s/%s" % (args.out_dir, first_part)
        job_name = "./" + first_part + ".pbs"

        command = "module load perl-modules\na5_pipeline.pl %s %s %s" % (first_file, second_file, first_part)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                out.write(mkdir_cmd+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')


# Function to create assembly pbs scripts for each fastq files
def create_ariba_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=60:00:00\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        #mkdir_cmd = "mkdir %s/%s" % (args.out_dir, first_part)
        cd_command = "cd %s" % (args.out_dir)
        job_name = "./" + first_part + ".pbs"

        command = "module load cd-hit\n/nfs/esnitkin/bin_group/anaconda3/bin/ariba run %s %s %s %s" % (args.ariba_db, first_file, second_file, first_part)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                #out.write(mkdir_cmd+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')


# Function to create srst2 pbs scripts for each fastq files
def create_srst2_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=60:00:00\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        #mkdir_cmd = "mkdir %s/%s" % (args.out_dir, first_part)
        cd_command = "cd %s" % (args.out_dir)
        job_name = "./" + first_part + ".pbs"

        command = "/nfs/esnitkin/bin_group/anaconda2/bin//srst2 --input_pe %s %s --output %s --log --gene_db %s --save_scores --report_all_consensus --forward _1_sequence --reverse _1_sequence" % (first_file, second_file, first_part, args.ariba_db)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                #out.write(mkdir_cmd+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')


def create_recombination_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        cd_command = "cd %s" % (args.out_dir)
        job_name = file + ".pbs"
        command = "/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/recombination_analysis.py -filename %s -out /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/generate_matrix/same_patient_different_species_same_facility/ -dir /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/fasta_files/ -steps 4 -prokka_dir /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations/ -jobrun parallel-local -filename_db same_patient_different_species_same_facility_nucmer_db_filenames -dir_db /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/Results/Results_parse_3/same_patient_different_species_same_facility/ariba_database/dedup_cluster/extracted_fasta/fasta_headers/" % file
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % os.path.basename(file)
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            #out.write(mkdir_cmd+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')














# Main Function; Check which Pipeline function to use based on command-line argument: pipeline
if args.pipeline == "varcall":
    print "Generating Variant Calling PBS scripts for samples in args.filenames"
    create_varcall_jobs()
elif args.pipeline == "old_assembly":
    print "Generating Assembly PBS scripts for samples in args.filenames"
    create_assembly_jobs()
elif args.pipeline == "assembly":
    print "Generating new Assembly PBS scripts for samples in args.filenames"
    list_of_files = get_filenames(args.dir, args.type, args.filenames, args.suffix)
    create_new_assembly_jobs(list_of_files)
elif args.pipeline == "a5_assembly":
    print "Generating a5 Assembly PBS scripts for samples in args.filenames"
    create_a5_jobs()
elif args.pipeline == "rna":
    print "Generating RNA Seq PBS scripts for samples in args.filenames"
    create_RNA_seq_jobs()
elif args.pipeline == "ptr":
    print "Generating PTR pipeline PBS scripts for samples in args.filenames"
    create_PTR_jobs()
elif args.pipeline == "new_assembly_SE":
    print "Generating new Assembly PBS scripts for samples in args.filenames"
    create_new_assembly_SE_jobs()
elif args.pipeline == "ariba":
    print "Generating new Ariba PBS scripts for samples in args.filenames"
    create_ariba_jobs()
elif args.pipeline == "srst2_plasmid":
    print "Generating new SRST2 Plasmid PBS scripts for samples in args.filenames"
    create_srst2_jobs()
elif args.pipeline == "recombination":
    print "Generating recombination PBS scripts for samples in args.filenames"
    create_recombination_jobs()
else:
    print "Please provide args.pipeline command-line argument"
