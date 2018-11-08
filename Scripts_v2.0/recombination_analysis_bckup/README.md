Synopsis:

The pipeline takes fasta files and their annotations as inputs, compares each fasta file to one another and extracts aligned regions that are at least 99% and 2000bp long. 

Since this pipeline is made to compare patients from same as well as different facilities, it requires the fasta file names to be in certain naming specifications.

example: 1388-G009-14-G-P-mirabilis-Aim2-CIP_R_final_ordered 

The filename needs to be divided into 8 parts, out of which only second (patient id: G009) and sixth part (species: mirabilis) will be used to generate different type of condition directories:

- different_patient_different_species_different_facility
- different_patient_different_species_same_facility
- different_patient_same_species_different_facility
- different_patient_same_species_same_facility
- same_patient_different_species_same_facility
- same_patient_same_species_same_facility

## Preparing input

1. Making a pseudomolecule of each fasta files

Run abacas on each fasta file against each respective reference genome (species: ex: order all mirabilis fasta files against a single mirabilis reference genome) to order it against a reference genome. 
(If you decide not to order you can also join all contigs into one pseudomolecule)

NH: these inputs are already generated and placed here: /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/fasta_files



2. Annotations:

The pipeline requires prokka annotation, 






The pipeline needs fasta files and their prokka annotations. Since we are comparing ea