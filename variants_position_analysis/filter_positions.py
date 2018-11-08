__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict


parser = argparse.ArgumentParser(description='Parsing All position with label file and investigating positions to determine the reason why it was filtered out from the final list')
#All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-positions_file_dir', action='store', dest="positions_file_dir", help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-positions_filenames', action='store', dest="positions_filenames", help='Names of filter2 only SNP vcf files with name per line.')
parser.add_argument('-unique_positions', action='store', dest="unique_positions", help='Names of unique_positions_file')

args = parser.parse_args()

All_position_file = args.positions_filenames
position_label = OrderedDict()


with open(All_position_file, 'rU') as csv_file:
        #reader = csv.reader(current_fp, delimiter="\t")
        # file_name = file + "_positions"
        # f1=open(file_name, 'w+')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #position = row[0]
            #print row
            for row in csv_reader:
                position_label[row[0]] = row[1:]

##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def generate_heatmap_position():
    f1=open("Only_ref_variant_positions", 'w+')
    for value in position_label:
        # key_array = value
        # key_array = []
        #lll = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        lll = ['0', '2', '3', '4', '5', '6', '7']
        ref_var = ['1', '1']
        # print cmp(lll, position_label[value])
        # print position_label[value]
        #if x in lll in position_label[value]:
        if set(ref_var) & set(position_label[value]):
            if set(lll) & set(position_label[value]):
                print value + "\t" + str(position_label[value])
         	#print "bakwass"   
	    #print value + "\t" + str(position_label[value])
            else:
                 strr = value + "\n"
                 f1.write(strr)
            # data = ""
            # for i in position_label[value]:
            #     data = data + str(i) + "\t"
            #
            # string_print = value + "\t" + data + "\n"
            # f1.write(string_print)


    # f1=open("All_label_file_proximate_suffix_raw_added_final_highproximate_only_position", 'w+')
    # for value in position_label:
    #     lll = ['0', '2', '3', '4', '5', '6']
    #     high_prox = ['7']
    #     prox_but_filtered_out = 0
    #     if set(high_prox) & set(position_label[value]):
    #         #Check if any positions contains highFQ proximate snps
    #         if set(lll) & set(position_label[value]):
    #             #check if that position contains any of the filtered variants
    #             prox_but_filtered_out += 1
    #         else:
    #             data = ""
    #             for i in position_label[value]:
    #                 data = data + str(i) + "\t"
    #
    #             string_print = value + "\t" + data + "\n"
    #             f1.write(string_print)

generate_heatmap_position()




####generate list of positions that dont include those positions which were filtered out
#### Even if it is filtered out in one it wont be included
def generate_position_array():
    for value in position_label:
        # key_array = value
        # key_array = []
        lll = ['LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ']
        # print cmp(lll, position_label[value])
        # print position_label[value]
        #if x in lll in position_label[value]:
        if not set(lll) & set(position_label[value]):
            print value
            #print value + "\t" + str(position_label[value])





#generate_position_array()




def generate_label_report():
    f1=open("./label_report_fom_unique", 'w+')
    unmapped_array = ['0']
    unmapped_count = 0
    #proximate_count = 0
    #proximate_filter_array = ['LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
    #filtered_array = ['LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ']
    reference_variant_array = ['1']
    ref_var_count = 0
    #filtered_count = 0
    filtered_out_low_FQ = ['2']
    filtered_out_low_FQ_count = 0
    filtered_out_high_FQ_but_low_DP_QUAL = ['3']
    filtered_out_high_FQ_but_low_DP_QUAL_count = 0
    filtered_out_high_FQ_but_low_DP_QUAL_proximate = ['4']
    filtered_out_high_FQ_but_low_DP_QUAL_proximate_count = 0
    lowFQ = ['5']
    lowFQ_count = 0
    highFQ = ['6']
    highFQ_count = 0
    highfq_proximate = ['7']
    highfq_proximate_count = 0
    filtered_out_positions = ['0', '2', '3', '4', '5', '6', '7']
    filtered_out_positions_count = 0
    for value in position_label:
    #     if len(set(reference_variant_array).intersection(position_label[value])) == 1:
    #         ref_var_count += 1
    #     #print ref_var_count
    #     if len(set(unmapped_array).intersection(position_label[value])) == 1:
    #         unmapped_count += 1
    # #print unmapped_count
    #     if len(set(filtered_out_low_FQ).intersection(position_label[value])) == 1:
    #         filtered_out_low_FQ_count += 1
    #     if len(set(filtered_out_high_FQ_but_low_DP_QUAL).intersection(position_label[value])) == 1:
    #         filtered_out_high_FQ_but_low_DP_QUAL_count += 1
    #     if len(set(filtered_out_high_FQ_but_low_DP_QUAL_proximate).intersection(position_label[value])) == 1:
    #         filtered_out_high_FQ_but_low_DP_QUAL_proximate_count += 1
    #     if len(set(lowFQ).intersection(position_label[value])) == 1:
    #         lowFQ_count += 1
    #     if len(set(highFQ).intersection(position_label[value])) == 1:
    #         highFQ_count += 1
    #     if len(set(highfq_proximate).intersection(position_label[value])) == 1:
    #         highfq_proximate_count += 1
    #     if len(set(filtered_out_positions).intersection(position_label[value])) == 1:
    #         filtered_out_positions_count += 1

        if set(reference_variant_array) & set(position_label[value]):
            ref_var_count += 1
        if set(unmapped_array) & set(position_label[value]):
            unmapped_count += 1
        if set(filtered_out_low_FQ) & set(position_label[value]):
            filtered_out_low_FQ_count += 1
        if set(filtered_out_high_FQ_but_low_DP_QUAL) & set(position_label[value]):
            filtered_out_high_FQ_but_low_DP_QUAL_count += 1
        if set(filtered_out_high_FQ_but_low_DP_QUAL_proximate) & set(position_label[value]):
            filtered_out_high_FQ_but_low_DP_QUAL_proximate_count += 1
        if set(lowFQ) & set(position_label[value]):
            lowFQ_count += 1
        if set(highFQ) & set(position_label[value]):
            highFQ_count += 1
        if set(highfq_proximate) & set(position_label[value]):
            highfq_proximate_count += 1
        if set(filtered_out_positions) & set(position_label[value]):
            filtered_out_positions_count += 1


    f1.write("The number of unmapped positions are: %s\n" % unmapped_count)
    f1.write("The number of Ref_Variant positions are: %s\n" % ref_var_count)
    f1.write("The number of filtered_out_positions_count positions are: %s\n\n" % filtered_out_positions_count)
    f1.write("The number of filtered_out_low_FQ positions are: %s\n" % filtered_out_low_FQ_count)
    f1.write("The number of filtered_out_high_FQ_but_low_DP_QUAL positions are: %s\n" % filtered_out_high_FQ_but_low_DP_QUAL_count)
    f1.write("The number of filtered_out_high_FQ_but_low_DP_QUAL_proximate positions are: %s\n" % filtered_out_high_FQ_but_low_DP_QUAL_proximate_count)
    f1.write("The number of lowFQ_count positions are: %s\n" % lowFQ_count)
    f1.write("The number of highFQ_count positions are: %s\n" % highFQ_count)
    f1.write("The number of highfq_proximate_count positions are: %s\n" % highfq_proximate_count)


    # f1=open("./label_report", 'w+')
    # unmapped_array = ['reference_unmapped_position']
    # unmapped_count = 0
    # proximate_count = 0
    # proximate_filter_array = ['LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
    # filtered_array = ['LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ']
    # reference_variant_array = ['reference_allele', 'VARIANT']
    # ref_var_count = 0
    # filtered_count = 0
    # variant_array = ['VARIANT']
    # variant_count = 0
    # only_variant = 0
    # for value in position_label:
    #     if set(reference_variant_array) & set(position_label[value]):
    #         ref_var_count += 1
    #     if set(filtered_array) & set(position_label[value]):
    #         filtered_count += 1
    #     if set(unmapped_array) & set(position_label[value]):
    #         unmapped_count += 1
    #     if set(proximate_filter_array) & set(position_label[value]):
    #         proximate_count += 1
    #     if set(variant_array) & set(position_label[value]):
    #         variant_count += 1
    #     if not set(unmapped_array) & set(position_label[value]):
    #         if not set(filtered_array) & set(position_label[value]):
    #             if not set(proximate_filter_array) & set(position_label[value]):
    #                 only_variant +=1
    #                 print value
    # f1.write("The number of unmapped positions are: %s\n" % unmapped_count)
    # f1.write("The number of proximate positions are: %s\n" % proximate_count)
    # f1.write("The number of filtered positions are: %s\n" % filtered_count)
    # f1.write("The number of Ref_Variant positions are: %s\n" % ref_var_count)
    # f1.write("The number of Variant positions are: %s\n" % variant_count)
    # f1.write("The number of only Variant positions are: %s\n" % only_variant)


#generate_label_report()



def keep_unique_positions_only():
    filename = args.unique_positions
    unique_position = []
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            unique_position.append(line)
    f1=open("All_label_file_unique_positions_only", 'w+')
    for value in position_label:
        if value in unique_position:
            data = ""
            for i in position_label[value]:
                data = data + str(i) + "\t"
            string_print = value + "\t" + data + "\n"
            f1.write(string_print)
    f1.close()

#keep_unique_positions_only()




