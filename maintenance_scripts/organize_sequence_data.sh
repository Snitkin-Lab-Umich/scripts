#!/usr/bin/env bash
# This script generates the Sequence_data and its sub-directories.

echo "Generating Sequence_data and its sub-directories for $USER"
PWD=`pwd`
echo "Generating Sequence_data directory in $PWD/"
mkdir $PWD/fastq $PWD/metadata $PWD/Reports $PWD/assembly  $PWD/assembly_annotation  $PWD/output_files  $PWD/fastq
