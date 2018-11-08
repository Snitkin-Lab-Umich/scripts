#Counts kmers in a file of sequences
#
#Takes as input:
#		1) seq_file	- A file of sequences, one per line
#		2) kmer_length	- THe length of kmers to counts
#
#Returns kmer counts in the sequences

#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl count_kmers.pl seq_file kmer_length\n\n";

	exit;

}#end if

use seq_stats;


#GET COMMAND LINE ARGUMENTS
my $seq_file = shift @ARGV;
my $kmer_length = shift @ARGV;

my %temp_kmer_counts;
my $r_temp_kmer_counts;
my %kmer_counts;


#GO THROUGH FILE COUNTING KMERS IN EACH SEQUENCE
open SEQS, $seq_file;

foreach my $seq (<SEQS>)
{

	chomp $seq;
	$seq =~ tr/[a-z]/[A-Z]/;

	$r_temp_kmer_counts = &seq_stats::kmer_count($seq, $kmer_length);
	%temp_kmer_counts = %$r_temp_kmer_counts;

	foreach my $kmer (keys %temp_kmer_counts)
	{

		$kmer_counts{$kmer} = $kmer_counts{$kmer} + $temp_kmer_counts{$kmer};

	}#end foreach


}#end foreach



#PRINT OUT KMER COUNTS
foreach my $kmer (sort keys %kmer_counts)
{

	if($kmer !~ /N/)
	{

		print "$kmer\t$kmer_counts{$kmer}\n";


	}#end if
}#end foreach


close SEQS;
