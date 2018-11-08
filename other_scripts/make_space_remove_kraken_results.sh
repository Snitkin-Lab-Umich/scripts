find . -type f -name "*.gz_unclassified.txt" -exec rm -f {} \;
find . -type f -name "*.gz_kraken" -exec rm -f {} \;
find . -type f -name "*_krona.input" -exec rm -f {} \;
find . -type f -name "*_krona.out.html.files" -exec rm -f {} \;
