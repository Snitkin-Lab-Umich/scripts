# input:
# arg1 = reference fasta file
# arg2 = path to bam files to run pilon on IN QUOTES 
# ex. generate_pilon_commands.sh ./kpnih1.fa "/path/to/files/*/*.bam"

for i in $(ls $2); do 

pref=$(echo $i | rev | cut -d/ -f1 | rev | cut -d_ -f1-2); 

echo "java -Xmx16G -jar /sw/med/centos7/pilon/1.22/pilon-1.22.jar --genome $1 --frags $i  --output $pref --outdir ${pref}_pilon --variant > ${pref}_pilon.out"; 

done > pilon_commands.sh
