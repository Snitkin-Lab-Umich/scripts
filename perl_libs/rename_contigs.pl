#CHANGES NAMES OF CONTIGS (ASSUMES THE ID NUMBER IS THE FINAL DIGITS AFTER A _)

#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nrename_contigs.pl contig_file renamed_contig_file prefix\n\n";

	exit;

}#end if

#READ IN COMMAND LINE VARIABLES
my $contig_file = shift @ARGV;
my $renamed_contig_file = shift @ARGV;
my $prefix = shift @ARGV;


#GO THROUGH FILE AND RENAME FASTA HEADERS
open IN_CONTIGS, $contig_file;
open OUT_CONTIGS, ">$renamed_contig_file";

foreach my $line (<IN_CONTIGS>)
{

	if($line =~ /^>/)
	{

		#HEADER LINE, SO REPLACE ACCORDINGLY
		$line =~ s/>.*(_\d+)$/>$prefix$1/;
		print OUT_CONTIGS $line;
		

	}else{
		#NOT HEADER, SO PRINT AS IS
		print OUT_CONTIGS $line;

	}#end if

}#end foreach

close IN_CONTIGS;
close OUT_CONTIGS;
