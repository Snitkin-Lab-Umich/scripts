package glimmer;

#sub parse_glimmer_fasta - Takes as input: 1) glimmer_file      - A fasta file produced by my glimmer shell script
#                       - Returns a hash where the keys are gene names and values are: start, end
sub parse_glimmer_fasta
{
        my ($glimmer_file) = @_;

        my %gene_hash;


        #GO THROUGH EACH LINE IN FILE 
        open INPUT, $glimmer_file;

        foreach my $line (<INPUT>)
        {

                if($line =~ />/)
                {

                        chomp $line;

                        my ($gene_name, $coords, $length) = split /\s\s/, $line;

                        $gene_name =~ s/>//;

			#IF GENE NAME HAS PIPE IN IT, THEN ASSUME THAT FIRST PART IS GENE NAME AND SECOND PART IS CONTIG THAT IT RESIDES ON
			#NOTE, THIS WILL MEAN THAT THE COORDINATES ARE RELATIVE TO THE CONTIG START!!!!!!!!!!!!!!!
			my @gene_name = split /\|/, $gene_name;
		

                        @coords = split /\s/, $coords;

                        $gene_hash{$gene_name[0]}{'START'} = $coords[0];
                        $gene_hash{$gene_name[0]}{'END'} = $coords[1];
                        $gene_hash{$gene_name[0]}{'TYPE'} = "protein";
                        $gene_hash{$gene_name[0]}{'CONTIG'} = $gene_name[1];
			
                }#end if

        }#end foreach


        return \%gene_hash;

}#end parse_glimmer_fasta


#sub parse_ncbi_fasta - Takes as input: 1) ncbi_file      - A cds fasta file produced by ncbi
#                     - Returns a hash where the keys are gene names and values are: start, end
sub parse_ncbi_fasta
{
        my ($ncbi_file) = @_;

        my %gene_hash;
	my $gene_count = 1;


        #GO THROUGH EACH LINE IN FILE 
        open INPUT, $ncbi_file;

        foreach my $line (<INPUT>)
        {

                if($line =~ />/)
                {

                        chomp $line;

                        my ($ref, $gene_name, $coords) = split /\|/, $line;

			$coords =~ s/^:c?//;
                        @coords = split /-/, $coords;

			$gene_name = $gene_name . "_" . $gene_count;

                        $gene_hash{$gene_name}{'START'} = $coords[0];
                        $gene_hash{$gene_name}{'END'} = $coords[1];
                        $gene_hash{$gene_name}{'TYPE'} = "protein";

			$gene_count ++;

                }#end if

        }#end foreach


        return \%gene_hash;

}#end parse_ncbi_fasta



1;
