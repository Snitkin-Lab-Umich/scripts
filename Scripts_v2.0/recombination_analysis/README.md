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

The pipeline requires prokka annotation. Annotate pseudomolecules using prokka and place the results in one folder. 

NH: these annotations are placed in /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations

## Run pipeline:

step 1: Run nucmer to compare each sample against each and generate coordinates file.

- generate a file containing fasta filenames, one line per line.

cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis
mkdir Results
cd Results

ls /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/fasta_files/*.fasta | awk -F'/' '{print $(NF)}' > filenames


Put the below command inside a pbs script and run it. Try giving it morfe than 8 cores so that it will run nucmer commands in parallel.
```
cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results
~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/recombination_analysis.py -filename filenames -out /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results/Results_parse_3 -dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/fasta_files/ -steps 1 -prokka_dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations/ -jobrun parallel-local
```

step 2: parse \*.coords file to generate condition directories and extract aligned fasta sequences for each species pairs

```
cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results
~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/recombination_analysis.py -filename filenames -out /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results/Results_parse_3 -dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/fasta_files/ -steps 2 -prokka_dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations/ -jobrun parallel-local
```

step 3: remove redundant species-pairs folders from each condition directories. we will use same_patient_different_species_same_facility as an example for condition directory 

```
cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results
~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/recombination_analysis.py -filename filenames -out /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results/Results_parse_3/same_patient_different_species_same_facility/ -dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/fasta_files/ -steps 3 -prokka_dir /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations/ -jobrun parallel-local
```

No go to same_patient_different_species_same_facility and make a new directory

```
mkdir ariba_database
cat */*.aligned.fasta > same_patient_different_species_same_facility.fasta
cat */*.fasta_meta.tsv > same_patient_different_species_same_facility.fasta.tsv

Run dedupe from BBtools to remove duplicate sequences and exact sequences.
pending


```

step 4: compare each aligned fragments to each fasta file and generate a score.



```
cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/
mkdir generate_matrix

cd generate_matrix
```

Put the below command in a seperate job:

```
cd /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/generate_matrix
~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/recombination_analysis.py -filename filenames -out /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/generate_matrix/same_patient_different_species_same_facility/ -dir /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/fasta_files/ -steps 4 -prokka_dir /scratch/esnitkin_fluxod/apirani/Project_Nursing_Home/Analysis/2017_10_28_HGT_Mummer_analysis/2017_10_28_HGT_Mummer_analysis_annotations/ -jobrun parallel-local -filename_db nucmer_db_filenames -dir_db /nfs/esnitkin/Project_Nursing_Home/Analysis/Regional_MDRO_populations/2017_10_28_HGT_Mummer_analysis/Results/Results_parse_3/same_patient_different_species_same_facility/ariba_database/dedup_cluster/extracted_fasta/fasta_headers/
```