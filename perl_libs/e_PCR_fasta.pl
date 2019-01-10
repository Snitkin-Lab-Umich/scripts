#Performs ePCR on genomes in database
#
#Takes as input:
#		1) genome_dir	- Directory with genome files
#		2) primers	- A file of primers (4-column file -> Primer set name|forward seq|reverse seq|amplicon size)
#		3) output_file	- A file to output the results of the e-PCR
#
#Returns a genome by primer matrix indicating which primers produced amplicons on the correct size
#

BEGIN {push @INC, '/nfs/esnitkin/bin_group/perl_libs/util';push @INC, '/nfs/esnitkin/bin_group/perl_libs/sequence_data'; push @INC, '/nfs/esnitkin/bin_group/perl_libs/math'}

use hash_lib;
use seq_io;


#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl e_PCR.pl genome_dir primers output_file\n\n";
	
	exit;

}#end if

#GET COMMAND LINE ARGUMENTS
my $genome_dir = shift @ARGV;
my $primers = shift @ARGV;
my $output_file = shift @ARGV;


#CONCATENATE CONTIGS INTO SINGLE SEQUENCE AND OUTPUT SEQUENCE OF CHROMOSOME TO NEW FILE
my @genomes = `ls $genome_dir`;

foreach my $genome (@genomes)
{

	chomp $genome;

	#READ IN FILE
	my ($r_genome_hash, $r_contig_names) = &seq_io::fasta2hash("$genome_dir/$genome");
	my %genome_hash = %$r_genome_hash;
	my @contig_names = @$r_contig_names;


	#CONCATENATE SEQUENCES
	my $concat_seq = join "", @genome_hash{@contig_names};


	#APPEND TO FILE
	&seq_io::append_fasta_file($genome, $concat_seq, 'temp_chr.fa');


}#end foreach


#PERFORM e-PCR
`e-PCR -w9 -f 1 -m100 -u+ $primers D=100-400 temp_chr.fa V=+ O=temp_pcr N=1 G=1 T=3`;


#GET THE NAMES OF THE PRIMERS
my $r_primer_hash = &hash_lib::dict_to_hash($primers);
my @primers = sort keys(%$r_primer_hash);


#PARSE OUTPUT FILE
my %pcr_results;

open PCR, "temp_pcr";

foreach my $line(<PCR>)
{

         #SPLIT UP LINE                                                                                                                                                                                                    
         chomp $line;

         my ($chr, $primer, $strand, $start, $end, $size, @junk) = split /\t/, $line;
         my @size = split /[-\/]/, $size;


         #MAKE SURE AMPLICON MATCHES EXPECTED SIZE                                                                                                                                                                         
         if(abs($size[1] - $size[2]) < 10)
         {

                 $pcr_results{$chr}{$primer} = "$start-$end";

         }else{

                 print "$chr/$primer SIZE DISCREPANCY ($size[1] vs. $size[2])!!!\n\n";
                 $pcr_results{$chr}{$primer} = "$start-$end";

         }#end if                                                                                                                                                                                                          

 }#end foreach                                                                                                                                                                                                             
 close PCR;


#CREATE OUTPUT MATRIX
open OUTPUT, ">$output_file";

print OUTPUT "\t" . join("\t", @primers) . "\n";

foreach my $chr(sort @genomes)
{

	chomp $chr;

	print $chr . "\n";
	print OUTPUT "$chr";

	foreach my $primer (@primers)
	{

		if(exists($pcr_results{$chr}{$primer}))
		{

			print OUTPUT "\t1";

		}else{

			print OUTPUT "\t0";

		}#end if

	}#end foreach

	print OUTPUT "\n";

}#end foreach

close OUTPUT;
