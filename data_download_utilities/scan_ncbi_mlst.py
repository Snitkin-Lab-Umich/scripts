__author__ = 'alipirani'

import sys
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
import select
import os
import argparse
import errno
import ConfigParser
import glob
import re
import multiprocessing
import logging
import subprocess as subprocess
from datetime import datetime
from joblib import Parallel, delayed
import xml.etree.ElementTree as ET
from argparse import RawTextHelpFormatter

""" Command Line Argument Parsing """
def parser():
    parser = argparse.ArgumentParser(description='\nParse Biosample XML file, Download fastq, Run Ariba\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-accession', action='store', dest="accession", help='Biosample accession list text file', required=True)
    required.add_argument('-dir', action='store', dest="dir", help='Output Directory Path', required=True)
    required.add_argument('-analysis', action='store', dest="analysis",
                          help='Unique analysis name',
                          required=True)
    optional.add_argument('-mlst_db', action='store', dest="mlst_db", help='Ariba Database to use for MLST scan', required=False)
    optional.add_argument('-coverage', action='store', dest="coverage", help='Minimum raw data download Coverage', required=False)
    optional.add_argument('-sra_file', action='store', dest="sra_file", help='List of SRA accession numbers',
                          required=False)

    return parser


def keep_logging(pmessage, lmessage, logger, mode):
    print(pmessage)
    if mode == 'warning':
        logger.warning(lmessage)
    elif mode == 'info':
        logger.info(lmessage)
    elif mode == 'exception':
        logger.exception(lmessage)
    elif mode == 'debug':
        logger.debug(lmessage)
    else:
        logger.error(lmessage)


def generate_logger(output_folder, analysis_name, log_unique_time):
    # Create a new logger; create a file handler; create a logging format; add the handlers to the logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #handler = logging.FileHandler('{}/{}_{}_.log.txt'.format(args.output_folder, args.analysis_name, datetime.now().strftime('%Y_%m_%d_%H_%M_%S')))
    handler = logging.FileHandler('{}/{}_{}.log.txt'.format(output_folder, log_unique_time, analysis_name))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def call(popenargs, logger, stderr_log_level=logging.ERROR, stdout_log_level=logging.DEBUG, **kwargs):
    """
    Variant of subprocess.call that accepts a logger instead of stdout/stderr,
    and logs stdout messages via logger.debug and stderr messages via
    logger.error.
    """

    try:
        child = sp.Popen(popenargs, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, **kwargs)
        log_level = {child.stdout: stdout_log_level, child.stderr: stderr_log_level}
        def check_io():
            ready_to_read = select.select([child.stdout, child.stderr], [], [], 1000)[0]
            for io in ready_to_read:
                line = io.readline()
                if line:
                    logger.log(log_level[io], line[:-1])
        # keep checking stdout/stderr until the child exits
        while child.poll() is None:
            check_io()
        check_io()  # check again to catch anything after the process exits
    except Exception as e:
        raise(e)
    return child.wait()

def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

# def parse_xml(xml, dir, analysis, logger):
#     keep_logging('Parsing XML file: %s' % xml,
#                  'Parsing XML file: %s' % xml, logger, 'info')
#     tree = ET.parse(xml)
#     root = tree.getroot()
#
#     print ET.tostring(root, encoding='utf8').decode('utf8')
#
#     # for movie in root.findall("./genre/decade/movie/[year='1992']"):
#     #     print(movie.attrib)
#
#     # for country in root.findall('geographic location'):
#     #     name = country.get('USA')
#     #     print name




def parse_accesion_file(accession, dir, analysis, logger):
    biosample_accession = []
    with open(accession) as fp:
        for line in fp:
            line = line.strip()
            if line != "":
                biosample_accession.append(line)
        fp.close()
    keep_logging('START: Parse Biosample accession file and generate a Biosample list to download.',
                 'START: Parse Biosample accession file and generate a Biosample list to download.', logger, 'info')


    print "No. of Biosample SRA Accession to download: %s" % len(biosample_accession)

    return biosample_accession

def run_command(i):
    """
    Function to run each command and is run as a part of python Parallel mutiprocessing method.
    :param: command
    :return: Done status
    """
    keep_logging('Running: %s' % i,
                 'Running: %s' % i, logger, 'info')
    call("%s" % i, logger)
    done = "Completed: %s" % i
    return done


def download_biosample(biosample_accession, dir, analysis, logger):
    keep_logging('START: Generate Fastq Dump commands and download fastq files.',
                 'START: Generate Fastq Dump commands and download fastq files.', logger, 'info')
    fastq_dump_cmd_list = []
    f1 = open("%s/fastq_dump_command_list.sh" % dir, 'w+')
    sra_file = open("%s/SRA_acc_list.txt" % dir, 'w+')
    #sra_acc = []
    num_cores = multiprocessing.cpu_count()
    for acc in biosample_accession:
        cmd2 = "/nfs/esnitkin/bin_group/EDirect/edirect/esearch -db sra -query %s </dev/null | /nfs/esnitkin/bin_group/EDirect/edirect/efetch -format docsum | xtract -pattern Runs -element Run@acc" % acc
        proc = subprocess.Popen([cmd2], stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc.communicate()
        if out2 and out2.startswith('SR'):
            fastq_dump_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/parallel-fastq-dump --sra-id %s --threads 4 --outdir %s --split-files --gzip -X %s" % (out2.strip(), dir, max_cov_spots)
            print out2.strip()
            f1.write(fastq_dump_cmd + '\n')
            sra_file.write("%s\n" % out2.strip())
            #sra_acc.append(out2.strip())
            fastq_dump_cmd_list.append(fastq_dump_cmd)
    f1.close()
    #results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in fastq_dump_cmd_list)
    #print len(sra_acc)
    #for sra in sra_acc:






""" Start of Main Method/Pipeline """
if __name__ == '__main__':

    # Set up logging modules and config file
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()

    args = parser().parse_args()

    global log_unique_time
    if args.dir != '':
        args.dir += '/'
    make_sure_path_exists(args.dir)
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    global logger
    logs_folder = args.dir + "/Logs"
    make_sure_path_exists(logs_folder)
    # logger = generate_logger(logs_folder, args.analysis_name, log_unique_time)

    logger = generate_logger(args.dir, args.analysis, log_unique_time)

    keep_logging('START: Download Biosamples From NCBI and Scan against Ariba MLST Db', 'START: Download Biosamples From NCBI and Scan against Ariba MLST Db', logger, 'info')
    #parse_xml(args.xml, args.dir, args.analysis, logger)
    biosample_accession= parse_accesion_file(args.accession, args.dir, args.analysis, logger)


    if args.coverage:
        max_cov_spots = args.coverage
    else:
        max_cov_spots = 1000000

    if args.sra_file:
        sra_acc = []
        with open(args.sra_file) as fp:
            for line in fp:
                line = line.strip()
                sra_acc.append(line)
        wget_cmd_list = []
        fastq_dump_cmd_list = []
        ariba_cmd_list = []
        if args.mlst_db:
            make_sure_path_exists("%s/ariba_results/")
        wget = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/"
        for acc in sra_acc:
            code = str(acc[:3])
            acc_code = str(acc[:6])
            wget_command = "%s/%s/%s/%s/%s.sra" % (wget, code, acc_code, acc, acc)
            fastq_dump_cmd = "fastq-dump --split-files --gzip -X %s %s.sra" % (max_cov_spots, acc)
            #print wget_command
            #print fastq_dump_cmd
            wget_cmd_list.append(wget_command)
            fastq_dump_cmd_list.append(fastq_dump_cmd)
            if args.mlst_db:
                ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/ariba run --verbose --force %s/ %s/%s_1.fastq.gz %s/%s_2.fastq.gz %s/ariba_results/%s" % (args.mlst_db, args.dir, acc, args.dir, acc, args.dir, acc)
                print ariba_cmd
                ariba_cmd_list.append(ariba_cmd)
        num_cores = multiprocessing.cpu_count()

        #results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in wget_cmd_list)
        #results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in fastq_dump_cmd_list)

        pending_wget_sra = []
        for acc in sra_acc:
            if not os.path.isfile("%s/%s.sra" % (args.dir, acc)):
                print "%s was not downloaded" % acc
                code = str(acc[:3])
                acc_code = str(acc[:6])
                wget_command = "%s/%s/%s/%s/%s.sra -O %s/%s" % (wget, code, acc_code, acc, acc, args.dir, acc)
                pending_wget_sra.append(wget_command)
                #os.system(wget_command)

        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in pending_wget_sra)

        pending_fastq_dump = []
        for acc in sra_acc:
            if not os.path.isfile("%s/%s_1.fastq.gz" % (args.dir, acc)):
                print "%s was not downloaded" % acc
                code = str(acc[:3])
                acc_code = str(acc[:6])
                fastq_dump_cmd = "fastq-dump --split-files --gzip -X %s %s.sra" % (max_cov_spots, acc)
                pending_fastq_dump.append(fastq_dump_cmd)
                print fastq_dump_cmd
                #os.system(fastq_dump_cmd)
        print "Pending fastq-dump commands: %s" % len(pending_fastq_dump)
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in pending_fastq_dump)

        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in ariba_cmd_list)

    else:
        #keep_logging('Exit: Download Biosamples From NCBI and Scan against Ariba MLST Db',
                     #'Exit: Download Biosamples From NCBI and Scan against Ariba MLST Db', logger, 'info')
        download_biosample(biosample_accession, args.dir, args.analysis, logger)
