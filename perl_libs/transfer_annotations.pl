#TRANSFERS ANNOTATIONS TO LS-BSR CLUSTERS
#
#Takes as input:
#		1) annot_blast_db	- A nucleotide database of gene sequences
#		2) annot_file		- A file of annotations for blasta db sequences
#		3) ls_bsr_mat		- A matrix produced by ls bsr
#		4) ls_bsr_fasta		- A fasta file of consensous sequences produced by ls-bsr
#
#Produces a modifed ls-bsr mat where annotations have been added based on best hist blast

use blast;
use seq_io;

#USAGE
if(@ARGV == 0)
{

	print "\n\nperl transfer_annotations.pl annot_blast_db annot_file ls_bsr_mat ls_bsr_fasta\n\n";

	exit;

}#end if

#GET COMMAND LINE ARGUMENTS
my $annot_blast_db = shift @ARGV;
my $annot_file = shift @ARGV;
my $ls_bsr_mat = shift @ARGV;
my $ls_bsr_fasta = shift @ARGV;


#READ IN LS-BSR SEQUENCES
my ($r_gene_hash, $r_gene_list) = &seq_io::fasta2hash($ls_bsr_fasta);


#READ IN ANNOTATIONS
open ANNOT, $annot_file;
my $annot_hash;

foreach my $line (<ANNOT>)
{

	my @line = split /\t/, $line;
	$annot_hash{$line[1]} = $line[3];

}#end foreach


#BLAST QUERY FASTA AGAINST DATABASE
my $r_blast_results = &blast::get_top_blast_hit($annot_blast_db, "temp_blast_results", $r_gene_hash);
my %blast_results = %$r_blast_results;


#GO THROUGH LS-BSR MATRIX, TRANSFERING ANNOTATION
my $ls_bsr_annot_mat = $ls_bsr_mat;
$ls_bsr_annot_mat =~ s/.txt/_annot.txt/;

open IN_BSR, $ls_bsr_mat;
my $header = <IN_BSR>;

open OUT_BSR, ">$ls_bsr_annot_mat";
print OUT_BSR $header;

foreach my $line (<IN_BSR>)
{

	my @line = split /\t/, $line;

	if(exists($blast_results{$line[0]}))
	{

		my @hit = keys(%{$blast_results{$line[0]}});
		$hit[0] =~ s/^(.*)\|.*$/$1/;

		my $annotation = $annot_hash{$hit[0]};
		$annotation =~ s/^([^\[]+) \[.*$/$1/;
		
		print  "*" . $hit[0] . "*\n";
		print "*" . $annotation . "*\n";
		$line[0] .= " ($annotation)";

	}

	$line = join "\t", @line;
	print OUT_BSR $line;

}#end foreach

close IN_BSR;
close OUT_BSR;
