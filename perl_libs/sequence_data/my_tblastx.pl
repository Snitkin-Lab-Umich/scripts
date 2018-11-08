#Takes as input a fasta file and a database and blasts the sequences against the database
#
#Takes as input
#		1) gene_fasta	- Fasta file of gene sequences to be blasted
#		2) blast_db	- Database tp blast against
#		3) blast_output	- An output file for blast results
#

#PRINT OUT USAGE
if(@ARGV == 0)
{

	print "\n\nmy_tblastx.pl gene_fasta blast_db blast_output\n\n";
	exit;

}#end if


use blast;
use seq_io;


#GET COMMAND LINE ARGUMENTS
my $gene_fasta = shift @ARGV;
my $blast_db = shift @ARGV;
my $blast_output = shift @ARGV;


#READ IN SEQUENCE FILE
my ($r_gene_hash, $r_gene_list) = &seq_io::fasta2hash($gene_fasta);


#BLAST
&blast::tblastx_vs_userDB($blast_db, $blast_output, $r_gene_hash);
