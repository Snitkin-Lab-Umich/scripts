#PARSES BED COVERAGAE AND PRODUCES BINARY FASTA INDICATING PRESENCE/ABSENCE OF GENOMIC REGION

use List::Util qw(sum);


#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl bedcov_to_var_fasta.pl bed_file fasta_file\n\n";

	exit;

}#end if


#GET COMMAND LINE ARGUMENTS
my $bed_file = shift @ARGV;
my $fasta_file = shift @ARGV;

#GO THOROUGH BED FILE AND CREATE LIST OF BINARY STRINGS
open BED, $bed_file;

my $header = <BED>;
chomp $header;
my ($tab1, $tab2, $tab3,  @taxa) = split /\t/, $header;

my $prev_cov = -1;
my %cov_taxa;

foreach my $line (<BED>)
{

	chomp $line;

	#$line =~ s/\t0\.9[5-9][^\t]+/\t1/g;
	$line =~ s/\t0\.[5-9][^\t]+/\t1/g;
	$line =~ s/\t0\.[0-4][^\t]+/\t0/g;
	$line =~ s/\t1\.0[^\t]+/\t1/g;

	my ($chr, $start, $end,  @cov) = split /\t/, $line;

	my $sum = sum(@cov);
	my $cur_cov = join("", @cov);

	if($cur_cov != $prev_cov && $sum != 0 && $sum != @cov)
	{
		print "$start\n";
		for(my $i = 0; $i < @taxa; $i ++)
                {

                        $seq_hash{$taxa[$i]} .= $cov[$i];

                }#end foreach


	}#end if

	$prev_cov = $cur_cov;

}#end foreach

#PRINT OUT FASTA FILE
open FASTA, ">$fasta_file";

foreach my $taxa (@taxa)
{

        print FASTA ">$taxa\n$seq_hash{$taxa}\n"

}#end foreach
