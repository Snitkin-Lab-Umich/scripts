find . -type f -name "*.sam" -exec rm -rf {} \;
find . -type f -name "*_aln.bam" -exec rm -rf {} \;
find . -type f -name "*_marked.bam" -exec rm -rf {} \;
find . -type f -name "*.bedcov" -exec rm -rf {} \;
find . -type f -name "*.fq.gz" -exec rm -rf {} \;
