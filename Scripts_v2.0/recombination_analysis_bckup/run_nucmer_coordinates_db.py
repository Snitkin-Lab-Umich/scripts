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
import errno
parser = argparse.ArgumentParser(description='Parse nucmer coordinates and delta files')
parser.add_argument('-folder', action='store', dest="folder", help='filename of pseudomolecule assembly fasta file')
parser.add_argument('-out', action='store', dest="out", help='out_directory')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fasta files')
parser.add_argument('-prokka_dir', action='store', dest="prokka_dir", help='directory of prokka annotations')
args = parser.parse_args()


def parse_coord_aggregate(folder, folder_dir):
    folder_base = os.path.basename(folder)
    list_command = "ls %s/%s_temp/*.coords" % (args.out, folder_base.replace('.fasta', ''))
    #print list_command
    list_of_coords_files = subprocess.check_output(list_command, shell=True)
    #print list_of_coords_files
    #comparison_file = args.out + (os.path.basename(folder)) + ".score"
    #f_out = open(comparison_file, 'w+')
    comparison_file_anno = args.out + (os.path.basename(folder)) + "_anno.score"
    f_out_anno_score = open(comparison_file_anno, 'w+')

    anno_folder = os.path.basename((folder.replace('_l500_contigs_final_ordered.fasta', '')).replace('_final_ordered.fasta',''))
    list_bed_cmd = "ls %s/%s*/*.bed" % (args.prokka_dir, anno_folder)
    bed_file = (subprocess.check_output(list_bed_cmd, shell=True)).strip()
    #print bed_file
    get_contig_name_cmd = "head -n1 %s | awk -F'\t' '{print $1}'" % bed_file
    contig_name = (subprocess.check_output(get_contig_name_cmd, shell=True)).strip()
    #print contig_name
    for lines in list_of_coords_files.splitlines():
        # condition = ""
        # sample = (os.path.basename(lines)).replace('.coords', '')
        # sample_split = sample.split('-')
        # sample_1 = sample_split[0:7]
        # sample_2 = sample_split[7:]
        # if sample_1[5] in sample_2:
        #     condition = condition + "_same_species"
        # else:
        #     condition = condition + "_different_species"
        #
        # #Set condition parameter
        # if sample_1[1] == sample_2[1]:
        #     condition = condition + "same_patient"
        # else:
        #     condition = condition + "different_patient"
        # if sample_1[5] == sample_2[5]:
        #     condition = condition + "_same_species"
        # else:
        #     condition = condition + "_different_species"
        # if sample_1[1][0] == sample_2[1][0]:
        #     condition = condition + "_same_facility"
        # else:
        #     condition = condition + "_different_facility"










        file_prefix = lines.replace('.coords', '')
        bed_out_file = file_prefix + ".bed"
        #f1 = open(bed_out_file, 'w+')
        aligned_bases = 0
        anno_aligned_bases = 0
        f = open(lines,'r')
        r_line = f.readlines()[4:]
        #print lines
        query_name = ((os.path.basename(lines).replace(contig_name, '')).replace('.coords', '')).replace('_final_ordered_', '')
        print query_name
        if len(r_line) > 0:
            cmd1 = "tail -n1 %s | awk -F\'\t\' \'{print $8}\'" % lines
            r_length = subprocess.check_output(cmd1, shell=True)
            cmd2 = "tail -n1 %s | awk -F\'\t\' \'{print $9}\'" % lines
            q_length = subprocess.check_output(cmd2, shell=True)
            for i in r_line:
                i = i.strip()
                i = i.split('\t')
                gaps = ""
                line_name = os.path.basename(lines)
                if float(i[6]) >= 95.00 and int(i[4]) > 2000:
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
                                # rna_genes_len = 0
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

                                if "RNA" in gene_types:
                                    gene_length = int(gene_end) - int(gene_start) + 1
                                    rna_genes_len = rna_genes_len + gene_length
                            else:
                                gene_list = gene_list + ""

                        rna_excluded_aligned_region = int(i[4]) - rna_genes_len
                        if rna_excluded_aligned_region < 2000:
                            rna_excluded_aligned_region = 0
                        if "CDS" in gene_type:
                            #fasta_header = ">%s:%s:%s-%s:%s" % ((contig_name.replace('_l500_contigs', '')).replace('_final_ordered', ''), (query_name.replace('_l500_contigs', '')).replace('_final_ordered', ''), str(i[0]), str(i[1]), str(i[4]))
                            #f_out_anno_aligned.write(fasta_header + '\n' + str(out) + '\n')
                            if gene_list == "":
                                gene_list = "."
                            #meta_line = "%s:%s:%s-%s:%s\t0\t1\t.\t.\t%s\t%s\t%s\n" % ((contig_name.replace('_l500_contigs', '')).replace('_final_ordered', ''), (query_name.replace('_l500_contigs', '')).replace('_final_ordered', ''), str(i[0]), str(i[1]), str(i[4]), gene_list, gene_list_flanking, gene_list_flanking_down)
                            #print "here"
                            #print meta_line
                            #f_out_anno_aligned_meta.write(meta_line)
                        anno_aligned_bases = int(anno_aligned_bases + rna_excluded_aligned_region)
                        #no_of_features = len(prokka_id_list)
                        #write_string = line_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                        #anno_string = contig_name + "\t" + query_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]) + "\t" + str(i[6]) + "\t" + str(gaps) + "\t" + str(set(gene_type)) + "\t" + str(no_of_features) + "\t" + gene_list + "\t" + str(rna_excluded_aligned_region) + "\n"
                        #f1.write(write_string)
                        #f_out_anno.write(anno_string)
                    else:
                        # extract_aligned_region_cmd = "tr -d \'\\n\' < %s | cut -b%s-%s" % (folder, str(i[0]), str(i[1]))
                        # proc = subprocess.Popen([extract_aligned_region_cmd], stdout=subprocess.PIPE, shell=True)
                        # (out, err) = proc.communicate()
                        # out = out.strip()
                        # #print extract_aligned_region_cmd
                        #write_string = line_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n"
                        #anno_string = contig_name + "\t" + query_name + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5])+ "\t" + str(i[6]) + "\t" + str(gaps) + "\t" + "All" + "\t" + "All" + "\t" + "All" + "\t" + str(i[4]) + "\t" + "All" + "\n"
                        #f1.write(write_string)
                        anno_aligned_bases = int(i[4])
                        #f_out_anno.write(anno_string)
                        #f_out_anno.close()
            #f_out_anno.close()
            #f1.close()
            # bed_cmd = "bedtools merge -i %s" % bed_out_file
            # merged_bed_positions = subprocess.check_output(bed_cmd, shell=True)
            # for rows in merged_bed_positions.splitlines():
            #     rows_split = rows.split('\t')
            #     sum = int(rows_split[2]) - int(rows_split[1]) + 1
            #     aligned_bases = aligned_bases + sum
            # uniq_ref = int(r_length.strip()) - aligned_bases
            # uniq_query = int(q_length.strip()) - aligned_bases
            # deno = aligned_bases + uniq_ref + uniq_query
            # score = float(aligned_bases / deno)
            #uniq_ref_anno = int(r_length.strip()) - anno_aligned_bases
            #uniq_query_anno = int(q_length.strip()) - anno_aligned_bases
            #deno_anno = anno_aligned_bases + uniq_ref_anno + uniq_query_anno

            fragment_length = int(q_length.strip())
            score_anno = float(anno_aligned_bases / fragment_length)
            if score_anno > 0.9999:
                score_anno = 1.0
            elif score_anno < 0:
                score_anno = 0.0
            #print anno_aligned_bases
            #print q_length.strip()
            #print score_anno
            # aligned_bases_kb = int(aligned_bases) / 1000
            # aligned_bases_kb_anno = int(anno_aligned_bases) / 1000
            # f_out.write(str(score) + '\n')
            #f_out_aligned.write(str(aligned_bases_kb) + '\n')
            f_out_anno_score.write(str(score_anno) + '\n')
            #f_out_aligned_anno.write(str(aligned_bases_kb_anno) + '\n')
        else:
            #score = "0"
            aligned_bases = "0"
            score_anno = "0"
            #f_out.write(str(score) + '\n')
            #f_out_aligned.write(str(aligned_bases) + '\n')
            f_out_anno_score.write(str(score_anno) + '\n')
            #f_out_aligned_anno.write(str(aligned_bases_kb_anno) + '\n')
    #f_out.close()
    f_out_anno_score
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

folder = args.folder
file_name = os.path.basename(folder)


#folder_dir = args.out + "/facilities/%s/%s/%s/" % (facility, facility_id, patient_id)
folder_dir = args.out + "/parse_nucmer_db"
make_sure_path_exists(folder_dir)


parse_coord_aggregate(folder, folder_dir)
