package assembly_file_parsing;

#sub parse_A5_asm_stats - Takes as input: 1) A5_stats_file - A file produced by A5 with a summary of the assembly
#			- Returns a hash with the column headers as keys and entries as values
sub parse_A5_asm_stats
{
	my ($A5_stats_file) = @_;

	my %stats_hash;

	#GO THROUGH FILE AND EXTRACT INFORMATION
	open A5, $A5_stats_file;

	my $header = <A5>;
	chomp $header;
	my @header = split /\t/, $header;

	my $values = <A5>;
	chomp $values;
	my @values = split /\t/, $values;

	close A5;

	@stats_hash{@header} = @values;


	return \%stats_hash;

}#end parse_A5_asm_stats


#sub parse_mira_asm_stats - Takes as input: 1) mira_stats_file - A file produced by mira with a summary of the assembly
#			  - Returns a hash with the column headers as keys and entries as values
sub parse_mira_asm_stats
{
	my ($mira_stats_file) = @_;

	my %stats_hash;

	#GO THROUGH FILE AND EXTRACT INFORMATION
	open MIRA, $mira_stats_file;
	my $start = 0;

	foreach my $line (<MIRA>)
	{

		chomp $line;
		
		#CHECK IF IN LARGE CONTIGS SECTION
		if($line =~ /^Num. reads assembled/)
		{
	
			$start = 1;

		}elsif($line =~ /^All contigs/){

			last;

		}#end if

		#CHECK FOR STATS OF INTEREST
		if($line =~ /Number of contigs|Largest contig|N50 contig size|Num\. reads assembled|Total consensus/)
		{

			my $key = $line;
			$key =~ s/^\s*([^:]+):.*$/$1/;

			my $value = $line;	
			$value =~ s/^.*:\s+(\d+)\s*$/$1/;

			$stats_hash{$key} = $value;

		}#end if

	}#end foreach

	return \%stats_hash;

}#end parse_mira_asm_stats

1;
