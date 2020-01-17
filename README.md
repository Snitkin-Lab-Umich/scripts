# Scripts

This repository contain scripts that were written or are actively being developed by Snitkin Lab members. The purpose of this repository is to make tracking of collaborative and personal projects easy and hold scripts required for certain type of analyses under one umbrella repository. This will help members to search for a script that the lab member has already developed for a task and avoid reinventing the wheel.

Lab members are highly encourage to add their frequently used scripts and analysis scripts to this repo. They can be arranged in different directories with directory name being the name of analysis. For example: Alignment_utilities folder will hol scripts that are useful in working with alignment format files such as converting one type of alignment format to another. 

## SOP for creating a new repository:

- Before creating a new directory, please go through the existing directories and decide if your analysis can a part of an existing analysis directory. If it requires a new directory, go ahead and make one with a README file.
- It is mandatory to Initialise each analysis repo with a README file - explaining what scripts/functions does it contain and how this scripts can be used by other lab members or integrated into their workflow.
- Since this scripts will be used by other lab members, make sure that you have thoroughly tested it and can work without breaking in a clean environment. SOP for testing the code will follow in the next iteration of this document update.
- Since this is a public repository, Before committing the code, make sure you dont have any sensitive patient data hardcoded. 
- Having a test data that can be run with the scripts or a test fuction that can run on its own can speed up the debugging and testing process. This will also help members to format their input data accordingly.


<!---
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


--->
