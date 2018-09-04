#!/bin/bash

# $1 = path to fasta file
# $2 = gubbins (1) or no gubbins (0)
#
# Usage:
# Gubbins:
# gubbins_iqtree_raxml.sh /path/to/whole/genome/alignment/fasta/file 1
# No gubbins:
# gubbins_iqtree_raxml.sh /path/to/snps/only/fasta/file 0

# get prefix for output files
pref=$(echo $1 | cut -d. -f1 | rev | cut -d/ -f1 | rev)

# modules to load (some of these might not be necessary)
modules=$(echo python-anaconda2/201607 biopython fasttree dendropy reportlab RAxML raxml bioperl fastml/gub gubbins openmpi/1.10.2/gcc/4.8.5 gcc/4.8.5)

# if gubbins
if [ $2 = 1 ]; then
  echo Will run gubbins.
  # gubbins command
  gub=$(echo run_gubbins.py --prefix $pref --threads 12 $1)
  echo $gub > gubbins_command.sh

  # raxml command
  raxml=$(echo mpirun -np 2 raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -N autoMRE -m ASC_GTRGAMMA --asc-corr=lewis -s $pref.filtered_polymorphic_sites.fasta -n $pref.filtered_polymorphic_sites_raxML-fa-x12345-p12345-NautoMRE-mGTRGAM-lew -T 6)
echo $raxml > raxml_command.sh

  # iqtree command
  iqtree=$(echo /nfs/esnitkin/bin_group/anaconda3/bin/iqtree -s $pref.filtered_polymorphic_sites.fasta -nt AUTO -bb 1000 -m MFP -pre $pref)
  echo $iqtree > iqtree_command.sh
  
  # generate pbs scripts for gubbins, raxml, iqtree
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c gubbins_command.sh -o ${pref}_gubbins.pbs -M "$modules"
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c iqtree_command.sh -o ${pref}_iqtree.pbs -M "$modules"
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c raxml_command.sh -o ${pref}_raxml.pbs -M "$modules"  

  # start gubbins, iqtree, raxml jobs
  echo qsub ${pref}_gubbins.pbs
  gub_qsub=$(qsub ${pref}_gubbins.pbs)
  gub_id=$(echo $gub_qsub | cut -d. -f1)
  echo gubbins job id: $gub_id
  echo qsub -W depend=afterok:$gub_id ${pref}_iqtree.pbs
  echo qsub -W depend=afterok:$gub_id ${pref}_raxml.pbs
  iqtree_qsub=$(qsub -W depend=afterok:$gub_id ${pref}_iqtree.pbs)
  raxml_qsub=$(qsub -W depend=afterok:$gub_id ${pref}_raxml.pbs)
else
  echo Will not run gubbins.
  # raxml command
  raxml=$(echo mpirun -np 2 raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -N autoMRE -m ASC_GTRGAMMA --asc-corr=lewis -s $1 -n $pref.filtered_polymorphic_sites_raxML-fa-x12345-p12345-NautoMRE-mGTRGAM-lew -T 6)
echo $raxml > raxml_command.sh

  # iqtree command
  iqtree=$(echo /nfs/esnitkin/bin_group/anaconda3/bin/iqtree -s $1 -nt AUTO -bb 1000 -m MFP -pre $pref)
  echo $iqtree > iqtree_command.sh

  # generate pbs scripts for gubbins, iqtree, raxml
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c gubbins_command.sh -o ${pref}_gubbins.pbs -M "$modules"
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c iqtree_command.sh -o ${pref}_iqtree.pbs -M "$modules"
  /nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/pipeline/Github/scripts/pbs_script_maker.py -c raxml_command.sh -o ${pref}_raxml.pbs -M "$modules"

  # start iqtree, raxml jobs
  echo qsub ${pref}_iqtree.pbs
  echo qsub ${pref}_raxml.pbs
  iqtree_qsub=$(qsub ${pref}_iqtree.pbs)
  raxml_qsub=$(qsub ${pref}_raxml.pbs)

fi

