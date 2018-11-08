from __future__ import division
import os
import argparse
import re
import subprocess
import readline
from joblib import Parallel, delayed
import multiprocessing
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser(description='Run Nucmer on Sample 2008 All-vs-All')
parser.add_argument('-filename', action='store', dest="filename", help='filename of pseudomolecule assembly fasta file')
parser.add_argument('-out', action='store', dest="out", help='out_directory')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fasta files')
parser.add_argument('-prokka_dir', action='store', dest="prokka_dir", help='directory of prokka annotations')
args = parser.parse_args()

filenames_array = []

with open(args.filename) as fp:
    for line in fp:
        line = line.strip()
        line = args.dir + "/" + line
        filenames_array.append(line)

temp_cmd = "./temp_commands"
def run_nucmer():
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
            #os.system(command)
	    f.write(command)
            command_array.append(command)
    #print command_array
    return command_array
    # with open(temp_cmd) as fp:
    #     for line in fp:
    #         line = line.strip()
    #         command_array.append(line)

command_array = run_nucmer()

num_cores = multiprocessing.cpu_count()

def run_command(i):
    os.system(i)
    done = "done"
    return done
results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in command_array)




def parse_coord():
    for folder in filenames_array:
        folder_base = os.path.basename(folder)
        list_command = "ls %s/%s/*.coords" % (args.out, folder_base.replace('.fasta', ''))
        list_of_coords_files = subprocess.check_output(list_command, shell=True)
        comparison_file = args.out + (os.path.basename(folder)) + ".score"
        #print comparison_file
        f_out = open(comparison_file, 'w+')
        for lines in list_of_coords_files.splitlines():
            aligned_bases = 0
            f = open(lines,'r')
            r_line = f.readlines()[4:]
            if len(r_line) > 0:
                bed_out_file = lines + ".bed"
                f1 = open(bed_out_file, 'w+')
                cmd1 = "tail -n1 %s | awk -F\'\t\' \'{print $8}\'" % lines
                r_length = subprocess.check_output(cmd1, shell=True)
                cmd2 = "tail -n1 %s | awk -F\'\t\' \'{print $9}\'" % lines
                q_length = subprocess.check_output(cmd2, shell=True)
                for i in r_line:
                    i = i.strip()
                    i = i.split('\t')
                    name = os.path.basename(lines)
                    write_string = name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                    #r_length = i[7]
                    #q_length = i[8]
                    f1.write(write_string)
                f1.close()
                bed_cmd = "bedtools merge -i %s" % bed_out_file
                merged_bed_positions = subprocess.check_output(bed_cmd, shell=True)
                for rows in merged_bed_positions.splitlines():
                    rows_split = rows.split('\t')
                    sum = int(rows_split[2]) - int(rows_split[1]) + 1
                    aligned_bases = aligned_bases + sum
                    uniq_ref = int(r_length.strip()) - aligned_bases
                uniq_query = int(q_length.strip()) - aligned_bases
                deno = aligned_bases + uniq_ref + uniq_query
                score = float(aligned_bases / deno)
                #score = aligned_bases/(aligned_bases + uniq_ref + uniq_query)
                #print score
                f_out.write(str(score) + '\n')
            else:
                score = "0"
                #print score
                f_out.write(str(score) + '\n')



#parse_coord()

# for i in *.fa*; do cat $i | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | cat <( echo ">$i") - > All_pseudomolecule/$i
# for i in */*.delta; do show-coords -r -c -l -T $i > $i.coords; done

# Create pseudomolecule from contigs and save it in a directory
# for i in *.fa*; do cat $i | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | cat <( echo ">$i") - > All_pseudomolecule/$i
# Run this script from the directory where pseudomolecule is saved
# It will create score files in the output directory.
# example: ~/anaconda/bin/python run_nucmer.py -filename filenames -out alignment_delta/ -dir ./All_pseudomolecule

def parse_coord_aggregate():
    for folder in filenames_array:
        folder_base = os.path.basename(folder)
        list_command = "ls %s/%s/*.coords" % (args.out, folder_base.replace('.fasta', ''))
        list_of_coords_files = subprocess.check_output(list_command, shell=True)
        comparison_file = args.out + (os.path.basename(folder)) + ".score"
        aligned_bases_file = args.out + (os.path.basename(folder)) + ".aligned"
        comparison_file_anno = args.out + (os.path.basename(folder)) + "_anno.score"
        aligned_bases_file_anno = args.out + (os.path.basename(folder)) + "_anno.aligned"
        #print comparison_file
        f_out = open(comparison_file, 'w+')
        f_out_aligned = open(aligned_bases_file, 'w+')
        f_out_anno_score = open(comparison_file_anno, 'w+')
        f_out_aligned_anno = open(aligned_bases_file_anno, 'w+')

        anno_folder = os.path.basename((folder.replace('_l500_contigs_final_ordered.fasta', '')).replace('_final_ordered.fasta',''))
        list_bed_cmd = "ls %s/%s*/*.bed" % (args.prokka_dir, anno_folder)
        #print list_bed_cmd
        bed_file = (subprocess.check_output(list_bed_cmd, shell=True)).strip()
        #print bed_file
        print "here"
        get_contig_name_cmd = "head -n1 %s | awk -F'\t' '{print $1}'" % bed_file
        contig_name = (subprocess.check_output(get_contig_name_cmd, shell=True)).strip()
        #contig_name = "619-G008-1-R-E-coli-CIP-R_l500_contigs"
        #print contig_name
        for lines in list_of_coords_files.splitlines():
            aligned_bases = 0
            anno_aligned_bases = 0
            f = open(lines,'r')
            r_line = f.readlines()[4:]
            print lines
            if len(r_line) > 0:

                bed_out_file = lines + ".bed"
                annotation_file = lines + ".annotations.txt"
                #bed_out_file_thr = lines + "_thr.bed"
                f1 = open(bed_out_file, 'w+')
                f_out_anno = open(annotation_file, 'w+')
                #f2 = open(bed_out_file_thr, 'w+')
                cmd1 = "tail -n1 %s | awk -F\'\t\' \'{print $8}\'" % lines
                r_length = subprocess.check_output(cmd1, shell=True)
                cmd2 = "tail -n1 %s | awk -F\'\t\' \'{print $9}\'" % lines
                q_length = subprocess.check_output(cmd2, shell=True)
                for i in r_line:
                    i = i.strip()
                    i = i.split('\t')
                    gaps = ""
                    line_name = os.path.basename(lines)
                    #print "%s\t%s" % (float(i[6]), int(i[4]))
                    if float(i[6]) >= 95.00 and int(i[4]) > 2000:
                        #print "%s\t%s" % (float(i[6]), int(i[4]))
                        if float(i[4]) == float(i[5]):
                            gaps = "No gaps"
                        elif float(i[4]) > float(i[5]):
                            gaps = float(i[4]) - float(i[5])
                        elif float(i[5]) > float(i[4]):
                            gaps = float(i[5]) - float(i[4])
                        if float(i[7]) - float(i[4]) > 10:
                            bedop_cmd = "echo -e \"%s\\t%s\\t%s\" | bedmap --echo --echo-map-id - %s" % (contig_name, str(i[0]), str(i[1]), bed_file)
                            #print bedop_cmd
                            bedop_results = (subprocess.check_output(bedop_cmd, shell=True)).split('|')
                            prokka_id_list = str(bedop_results[1]).split(';')
                            gene_list = ""
                            gene_type = []
                            #print prokka_id_list
                            rna_genes_len = 0
                            for id in prokka_id_list:
                                if id.startswith('PROKKA_'):

                                    get_gene_details_cmd = "grep -w '%s' %s" % (id.strip(), bed_file)
                                    get_gene_details = subprocess.check_output(get_gene_details_cmd, shell=True)
                                    get_gene_details_list = get_gene_details.split('\t')
                                    get_gene_name = "grep -w '%s' %s | awk -F'\\t' '{print $10}' | sed 's/ID=.*;product=//g' | tr -d '\\n'" % (id.strip(), bed_file)
                                    name = subprocess.check_output(get_gene_name, shell=True)
                                    gene_list = gene_list + ";" + (str(name)).strip()
                                    gene_types = get_gene_details_list[7].strip()
                                    gene_type.append(gene_types.strip())
                                    gene_start = get_gene_details_list[1]
                                    gene_end = get_gene_details_list[2]

                                    #get_gene_name = "grep -w '%s' %s | awk -F'\\t' '{print $10}' | sed 's/ID=.*;product=//g' | tr -d '\\n'" % (id.strip(), bed_file)
                                    #name = subprocess.check_output(get_gene_name, shell=True)
                                    #gene_list = gene_list + ";" + (str(name)).strip()
                                    #get_gene_type = "grep -w '%s' %s | awk -F'\\t' '{print $8}'" % (id.strip(), bed_file)
                                    #gene_types = subprocess.check_output(get_gene_type, shell=True)
                                    #gene_type.append(gene_types.strip())
                                    #gene_start_cmd = "grep -w '%s' %s | awk -F'\\t' '{print $2}'" % (id.strip(), bed_file)
                                    #gene_end_cmd = "grep -w '%s' %s | awk -F'\\t' '{print $3}'" % (id.strip(), bed_file)
                                    #gene_start = subprocess.check_output(gene_start_cmd, shell=True)
                                    #gene_end = subprocess.check_output(gene_end_cmd, shell=True)
                                    #print "%s; %s" % (str(gene_start).strip(), str(gene_end).strip())
                                    #print gene_start_cmd
                                    #print id.strip()
                                    if "RNA" in gene_types:
                                        #print "RNA"
                                        gene_length = int(gene_end) - int(gene_start) + 1
                                        #print gene_length
                                        rna_genes_len = rna_genes_len + gene_length
                                else:
                                    print "Problem: %s" % prokka_id_list
                            #print "RNA" + str(rna_genes_len)
                            #print "region" + str(i[4])
                            rna_excluded_aligned_region = int(i[4]) - rna_genes_len
                            if rna_excluded_aligned_region < 2000:
                                rna_excluded_aligned_region = 0
                            anno_aligned_bases = int(anno_aligned_bases + rna_excluded_aligned_region)
                            #print rna_excluded_aligned_region
                            #print str(i[4])
                            no_of_features = len(prokka_id_list)
                            #print set(gene_type)
                            write_string = line_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                            anno_string = contig_name + "\t" + query_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]) + "\t" + str(i[6]) + "\t" + str(gaps) + "\t" + str(set(gene_type)) + "\t" + str(no_of_features) + "\t" + gene_list + "\t" + str(rna_excluded_aligned_region) + "\n"
                            #print (str(anno_string)).strip()
                            f1.write(write_string)
                            #f1.close()
                            f_out_anno.write(anno_string)
                            #f_out.close()
                        else:
                            #print "All genes"
                            write_string = line_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                            anno_string = contig_name + "\t" + query_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5])+ "\t" + str(i[6]) + "\t" + str(gaps) + "\t" + "All" + "\t" + "All" + "\t" + "All" + "\t" + str(i[4]) + "\n"
                            f1.write(write_string)
                            anno_aligned_bases = int(i[4])
                            #f1.close()
                            f_out_anno.write(anno_string)
                            #f_out.close()
                f_out_anno.close()
                f1.close()
                bed_cmd = "bedtools merge -i %s" % bed_out_file
                merged_bed_positions = subprocess.check_output(bed_cmd, shell=True)
                for rows in merged_bed_positions.splitlines():
                    rows_split = rows.split('\t')
                    sum = int(rows_split[2]) - int(rows_split[1]) + 1
                    aligned_bases = aligned_bases + sum
                #print aligned_bases
                #print anno_aligned_bases
                uniq_ref = int(r_length.strip()) - aligned_bases
                uniq_query = int(q_length.strip()) - aligned_bases
                deno = aligned_bases + uniq_ref + uniq_query
                score = float(aligned_bases / deno)
                uniq_ref_anno = int(r_length.strip()) - anno_aligned_bases
                uniq_query_anno = int(q_length.strip()) - anno_aligned_bases
                deno_anno = anno_aligned_bases + uniq_ref_anno + uniq_query_anno
                score_anno = float(anno_aligned_bases / deno_anno)
                if score_anno > 0.9999:
                    score_anno = 1
                elif score_anno < 0:
                    score_anno = 0
                #score = aligned_bases/(aligned_bases + uniq_ref + uniq_query)
                #print score

                aligned_bases_kb = int(aligned_bases) / 1000
                aligned_bases_kb_anno = int(anno_aligned_bases) / 1000
                #print aligned_bases_kb
                f_out.write(str(score) + '\n')
                f_out_aligned.write(str(aligned_bases_kb) + '\n')
                f_out_anno_score.write(str(score_anno) + '\n')
                f_out_aligned_anno.write(str(aligned_bases_kb_anno) + '\n')
            else:
                score = "0"
                aligned_bases = "0"
                f_out.write(str(score) + '\n')
                f_out_aligned.write(str(aligned_bases) + '\n')
                f_out_anno_score.write(str(score_anno) + '\n')
                f_out_aligned_anno.write(str(aligned_bases_kb_anno) + '\n')
        f_out.close()
        f_out_aligned.close()
        f_out_anno.close()
        f_out_aligned_anno.close()

#parse_coord_aggregate()


def generate_parse_coord_aggregate_jobs():
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/plasmid_analysis/run_nucmer_coordinates.py -folder %s -out %s -dir %s -prokka_dir %s" % (folder, args.out, args.dir, args.prokka_dir)
        print cmd
        job_name = os.path.basename(folder)
        job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
        job_file_name = "%s.pbs" % (job_name)
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()






#generate_parse_coord_aggregate_jobs()