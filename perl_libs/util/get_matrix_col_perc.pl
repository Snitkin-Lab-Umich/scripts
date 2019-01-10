#Takes in a matrix and returns the same matrix with column percentages in parens for each element
#
#Takes as input:
#		1) mat_file		- The matrix file
#		2) num_col_headers	- The number of rows that have descriptive data, that should not be analyzed
#		3) num_row_headers	- The number of columns that have descriptive data, that should not be analyzed
#		4) output_file		- An output file to print the matrix to
#
#Returns the input matrix, where the values in each data column now have the percent of the column total in parens


#PRINT USAGE STATEMENT
if(@ARGV == 0)
{

	print "\n\nperl get_matrix_col_perc.pl mat_file num_col_headers num_row_headers\n\n";

	exit;

}#end if

use list_ops;


#GET COMMAND LINE ARGUMENTS
my $mat_file = shift @ARGV;
my $num_col_headers = shift @ARGV;
my $num_row_headers = shift @ARGV;
my $output_file = shift @ARGV;


#OPEN FILE
open MAT, $mat_file;
my @rows = <MAT>;
close MAT;


#REMOVE HEADER LINES
my @col_headers;

for(my $i = 0; $i < $num_col_headers; $i ++)
{

	push @col_headers, $rows[$i];

}#end for


#GO THROUGH EACH ROW AND BUILD COLUMN LISTS
my $num_cols = 0;
my $num_rows = @rows;

for(my $r = $num_col_headers; $r < $num_rows; $r ++)
{


	chomp $rows[$r];
	my @row = split /\t/, $rows[$r];

	#GET THE NUMBER OF COLUMNS
	if($r == $num_col_headers)
	{

		$num_cols = @row;

	}#end if

	#GO THROUGH COLUMNS AND PLACE IN LIST
	for(my $c = 0; $c < @row; $c ++)
	{
	
		push @{$col_hash{$c}}, $row[$c];

	}#end for

}#end for


#GET THE PERCENTAGES FOR EACH COLUMN 
my %col_perc_hash;

for(my $c = $num_row_headers; $c < $num_cols; $c ++)
{

	$col_perc_hash{$c} = &list_ops::rnd(&list_ops::list_perc_total($col_hash{$c}), 2);

}#end foreach


#PRINT OUT NEW MATRIX
open OUTPUT, ">$output_file";

#PRINT OUT COLUMN HEADERS
foreach my $col_head (@col_headers)
{

	print OUTPUT "$col_head";

}#end foreach


#PRINT OUT EACH ROW
for(my $r = 0; $r < ($num_rows - $num_col_headers); $r ++)
{

	for(my $c = 0; $c < $num_cols; $c ++)
	{

		if( ($c == 0) && (exists($col_perc_hash{$c})) )
		{

			print OUTPUT @{$col_hash{$c}}[$r] . " (" . @{$col_perc_hash{$c}}[$r] . "%)";

		}elsif($c == 0){

			print OUTPUT @{$col_hash{$c}}[$r];

		}elsif(exists($col_perc_hash{$c})){

			print OUTPUT "\t" . @{$col_hash{$c}}[$r] . " (" . @{$col_perc_hash{$c}}[$r] . "%)";

		}else{

			print OUTPUT "\t" .  @{$col_hash{$c}}[$r];

		}#end if


	}#end for	

	print OUTPUT "\n";

}#end for

close OUTPUT;
