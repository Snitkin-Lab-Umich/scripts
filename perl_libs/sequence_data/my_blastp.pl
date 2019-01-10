#Takes as input a fasta file and a database and blasts the sequences against the database
#
#Takes as input
#		1) gene_fasta	- Fasta file of gene sequences to be blasted
#		2) blast_db	- Database tp blast against
#		3) blast_output	- An prefix for output files with blast results
#

#PRINT OUT USAGE
if(@ARGV == 0)
{

	print "\n\nmy_blastp.pl gene_fasta blast_db blast_output\n\n";
	exit;

}#end if


use blast;
use seq_io;


#GET COMMAND LINE ARGUMENTS
my $gene_fasta = shift @ARGV;
my $blast_db = shift @ARGV;
my $blast_output = shift @ARGV;


#CREATE OUTPUT FILES
my $bo_file = $blast_output . ".bo";
my $pbo_file = $blast_output . ".pbo";

#READ IN SEQUENCE FILE
my ($r_gene_hash, $r_gene_list) = &seq_io::fasta2hash($gene_fasta);


#BLAST
my ($r_blast_results_hash, $r_blast_results_list) = &blast::blastp_vs_userDB($blast_db, $bo_file, $r_gene_hash);
my %blast_results = %$r_blast_results_list;


#OUTPUT PARSED RESULTS
open PBO, ">$pbo_file";

foreach my $query(sort keys %blast_results)
{

	foreach my $hit (@{$blast_results{$query}})
	{

		my ($hit_name, $significance, $q_length, $h_length, $h_frac_al, $q_frac_al, $q_start, $q_end, $h_start, $h_end, $hit_desc) = split /\t/, $hit;

		print PBO "$query\t$hit_name\t$hit_desc\t$significance\t$q_length\t$h_length\t$q_frac_al\t$h_frac_al\n";
		

	}#end foreach
	
	print PBO "\n";


}#end foreach
