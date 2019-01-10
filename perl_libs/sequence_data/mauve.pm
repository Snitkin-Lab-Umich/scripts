package mauve;

#sub parse_mauve_contigs - Takes as input: contig_file - The output of the mauve contig mover process
#			 - Returns a list of contigs in order and  hash where the keys are contigs and
#			   the values are 'forward' or 'coplement' depending on how the sequence is oriented
sub parse_mauve_contigs
{
	my ($contig_file) = @_;

	my @contigs;
	my %contig2orient;


	#GO THROUGH FILE AND RETRIVE ORDERED CONTIGS
	my $contig_flag = 0;

	open MAUVE, $contig_file;

	foreach my $line (<MAUVE>)
	{

		#CHECK IF IN CONTIG SECTION
		if( ($contig_flag == 1)  && $line =~ /^contig/)
		{
			
			chomp $line;

			my @line = split /\t/, $line;
		
			push @contigs, $line[1];
			$contig2orient{$line[1]} = $line[3];

		}elsif($contig_flag == 1 && $line !~ /type/){

			last;
	
		}#end if


		#CHECK FOR BEGINNING OF ORDERED CONTIG SECTION
		if($line =~ /Ordered Contigs/)
		{
			
			$contig_flag = 1;

		}#end if

	}#end foreach


	return (\@contigs, \%contig2orient);

}#end parse_mauve_contigs




1;
