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

	print "\n\nmy_blastn.pl gene_fasta blast_db blast_output\n\n";
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
my $tab_file = $blast_output . ".tab";


#READ IN SEQUENCE FILE
my ($r_gene_hash, $r_gene_list) = &seq_io::fasta2hash($gene_fasta);


#BLAST
my $r_blast_results = &blast::blastn_vs_userDB($blast_db, $bo_file, $r_gene_hash);
my %blast_results = %$r_blast_results;


#OUTPUT PARSED BLAST RESULTS, AND KEEP TRACK OF HITS TO MAKE TABLE
open PBO, ">$pbo_file";

my %hits;
my %hit_qaln;

foreach my $query(sort keys %blast_results)
{

        foreach my $hit (keys %{$blast_results{$query}})
        {

                my ($significance, $q_length, $h_length, $h_frac_al, $q_frac_al, $q_start, $q_end, $h_start, $h_end, $score, $hit_desc) = @{$blast_results{$query}{$hit}};

                print PBO "$query\t$hit\t$hit_desc\t$significance\t$q_length\t$h_length\t$q_frac_al\t$h_frac_al\n";

		$hit_qaln{$query}{$hit} = $q_frac_al;
		$hits{$hit} = 1;

        }#end foreach

        print PBO "\n";

}#end foreach


#OUTPUT TABLE OUTPUT
open TAB, ">$tab_file";

#PRINT HEADER
print TAB "\t" . join("\t", sort(keys(%hits))) . "\n";


#GO THROUGH EACH QUERY AND PRINT OUT FRACTION ALIGNED TO DIFFERENT HITS
foreach my $query(sort keys %blast_results)
{

	print TAB $query;

	foreach my $hit (sort keys %hits)
	{

		if(exists($blast_results{$query}{$hit}))
		{

			print TAB "\t$hit_qaln{$query}{$hit}";	

		}else{

			print TAB "\t0";	

		}#end if

	}#end foreach

	print TAB "\n";
}#end foreach
