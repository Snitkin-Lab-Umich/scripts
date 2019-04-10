# Scripts
This github repository contain frequently used scripts being used by lab members.

The scripts are arranged in different directories based on their specific utilities.

## Alignment_utilities

- convertAlignment.pl
Converts an alignment to another alignment format. The script will make a guess if the input format is one of these formats: fasta, clustalw, phylip, xmfa

## Annotation_utilities

- ecoli_assign_phylogroup.py
This script assigns phylogroup based on the combination of presence and absence of E. coli specific primer target sequences such as: arpA,chuA,tspe4,yjaA,C_specific,E_specific.

- extract_cds.py

- extract_genes_from_genbank.py
- extract_sRNA_regulonDB.py
- gb2ptt.pl
- gbtofasta.py
- gbtogff.pl
- gff2bed
- gff_to_genbank.py
- gtf2bed
- run_roary.pbs

variant_parser_functions.R 
To parse the SNP matrix output from the variant calling pipeline

render_variant_matrix_qc.R
To perform QC on variant matrices (SNV and Indel)
