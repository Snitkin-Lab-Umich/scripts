#PRODUCES AN OUTPUT FILE WHICH CONTAINS THE FIRST N COLUMNS OF THE INPUT FILE

#Takes as input:
#		1) in_file	- A tab-delimited input file
#		2) start_col	- The first column to print out (zero-based)
#		3) end_col	- The last column to print out (zerp-based)
#		4) out_file	- An output file 
#
#Writes to output file the desired column range

#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl awk.pl in_file num_cols out_file\n\n";

	exit;

}#end if


my $in_file = shift @ARGV;
my $start_col = shift @ARGV;
my $end_col = shift @ARGV;
my $out_file = shift @ARGV;

open IN_FILE, "$in_file";
open OUT_FILE, ">$out_file";

my @file = <IN_FILE>;

foreach my $line (@file)
{

	chomp $line;

	my @line = split "\t", $line;

	my @line_subset = @line[$start_col..($end_col)];

	my $line_subset = join "\t", @line_subset;

	print OUT_FILE $line_subset . "\n";

}#end foreach

close IN_FILE;
close OUT_FILE;
