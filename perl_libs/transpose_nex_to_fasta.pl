#Converts transpose nexus matrix to fasta format
#
#USAGE
if(@ARGV == 0)
{

	print "\n\nprint transpose_nex_to_fasta.pl nex_file\n\n";
	
	exit;

}#end if


#GET COMMAND LINE ARGUMENTS
my $nex_file = shift @ARGV;


#GO THROUGH FILE AND CONSTRUCT SEQUENCE BY STRAIN
open NEXUS, $nex_file;

my @taxa;
my $junk;
my $seq = 0;

foreach my $line (<NEXUS>)
{

	chomp $line;

	#CHECK IF END OF FILE
	if($line =~ /^;/)
	{
		
		last;

	}#end if
	
	#CHECK IF LINE HAS TAXA LABELS
	if($line =~ /^taxlabels/)
	{

		$line =~ s/;$//;
		($junk, @taxa) = split / /, $line;

	}#end if


	#IF SEQUENCE HAS STARTED THEN PARSE
	if($seq == 1)
	{
	
		my ($label, @alleles) = split / /, $line;

		for(my $i = 0; $i < @taxa; $i ++)
		{

			$seq_hash{$taxa[$i]} .= $alleles[$i];

		}#end foreach
		
	}#end if

	#CHECK IF SEQUENCE IS STARTING
	if($line =~ /^matrix/)
	{

		$seq = 1;

	}#end if

	

}#end foreach


#PRINT OUT FASTA FILE
my $fasta_file = $nex_file;
$fasta_file =~ s/.nex/.fasta/;

open FASTA, ">$fasta_file";

foreach my $taxa (@taxa)
{

	print FASTA ">$taxa\n$seq_hash{$taxa}\n"

}#end foreach
