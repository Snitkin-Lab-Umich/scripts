
filtered_out_vcf_file = "filtered_out_vcf_files"

filtered_out_vcf_file_array = []

import os

with open(filtered_out_vcf_file) as fp:
    for line in fp:
            line = line.strip()
            filtered_out_vcf_file_array.append(line)


filename = "consensus.sh"


for file in filtered_out_vcf_file_array:
    f1 = open(filename, 'a+')
    bgzip_cmd = "bgzip %s\n" % file
    f1.write(bgzip_cmd)
    tabix_cmd = "tabix -p vcf %s.gz\n" % file
    f1.write(tabix_cmd)
    fasta_cmd = "cat /home2/apirani/bin/reference/MRSA_USA_300/MRSA_USA_300.fasta | /home/apirani/bin/vcftools_0.1.12b/bin/vcf-consensus %s.gz > %s_filter2_consensus.fa\n" % (file, file)
    f1.write(fasta_cmd)
    base = os.path.basename(file)
    sed_command = "sed -i 's/>.*/>%s/g' %s_filter2_consensus.fa\n" % (base, file)
    f1.write(sed_command)



