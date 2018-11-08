package file_lib;

#sub list2file - Takes as input : 1) file 	- The name of a file
#				  2) list	- A reference to a list
#		- Outputs the contents of teh list to the file
sub list2file
{
	my ($file , $r_list) = @_;
	my @list = @$r_list;


	open OUTPUT , ">$file";


	foreach my $entry (@list)
	{

		print OUTPUT $entry . "\n";

	}#end foeach	


	close OUTPUT;


}#end list2file	



#sub matrix2file - Takes as input : 1) mat_hash	- A hash where the first value is the row value and the second is the column value
#				    2) row_list	- A list of values for the rows
#				    3) col_list	- A list of values for the columns
#				    4) mat_file	- The name of the file to output the matrix to
#		 - Outputs matrix specified by input to file
sub matrix2file
{
	my ($r_mat_hash , $r_row_list , $r_col_list , $mat_file) = @_;

	my %mat_hash = %$r_mat_hash;
	my @row_list = @$r_row_list;
	my @col_list = @$r_col_list;


	#OUTPUT MATRIX TO FILE
	open MAT_FILE , ">$mat_file";


	foreach my $row (@row_list)
	{


		foreach my $col (@col_list)
		{


			if( exists( $mat_hash{$row} -> {$col} ) )
			{

				print MAT_FILE $mat_hash{$row} -> {$col} . "\t"; 	

			}else{

				print MAT_FILE "0\t";
	
			}#end if

		}#end foreach


		print MAT_FILE "\n";

	}#end foreach 


	close MAT_FILE;


}#end matrix2file


#sub file_to_list - Takes as input: 1) file - The name of a file
#		  - Returns a reference to a list with each entry containing a line in the file without a newline
sub file_to_list
{
	my ($file) = @_;
	my @list;


	open FILE , $file;

	foreach my $line (<FILE>)
	{

		chomp $line;

		push @list , $line;

	}#end foreach


	close FILE;

	return \@list;

}#end file_to_list


#sub matrix_file_to_hash - Takes as input: 1) file - THe name of the file
#			 - Returns a 2D hash where keys are indicies in matrix and values are entries from
#			   the LOWER TRIANGLE OF THE MATRIX. Note that indicies start with 0, not 1.
sub matrix_file_to_hash
{
	my ($file) = @_;

	my $line_num = 0;
	my %matrix_hash;

	#PARSE THE FILE AND STORE IN THE HASH
	open FILE , $file;

	
	foreach my $line (<FILE>)
	{

		chomp $line;


		my @line = split /\t/ , $line;

		for(my $col = 0; $col < (@line); $col ++)
		{


			$matrix_hash{$line_num}{$col} = $line[$col];
			$matrix_hash{$col}{$line_num} = $line[$col];

		}#end for

	
		$line_num ++;

	}#end foreach


	close FILE;


	return \%matrix_hash;

}#end matrix_file_to_hash


#sub matrix_file_wLabels_to_hash - Takes as input: 1) file - THe name of the file
#			 	 - Returns a 2D hash where keys are row/column names in matrix and values are entries from matrix
sub matrix_file_wLabels_to_hash
{
	my ($file) = @_;

	my %matrix_hash;

	#PARSE THE FILE AND STORE IN THE HASH
	open FILE , $file;

	my $cols = <FILE>;
	chomp $cols;
	my @cols = split /\t/, $cols;

	
	foreach my $line (<FILE>)
	{

		chomp $line;


		my ($row, @line) = split /\t/ , $line;

		for(my $col = 0; $col < (@line); $col ++)
		{

			$matrix_hash{$row}{$cols[$col]} = $line[$col];
			$matrix_hash{$cols[$col]}{$row} = $line[$col];

		}#end for


	}#end foreach


	close FILE;


	return \%matrix_hash;

}#end matrix_file_to_hash



#sub LT_matrix_file_to_hash - Takes as input: 1) file - THe name of the file
#			    - Returns a 2D hash where keys are indicies in matrix and values are entries from
#			      the LOWER TRIANGLE OF THE MATRIX. Note that indicies start with 0, not 1.
sub LT_matrix_file_to_hash
{
	my ($file) = @_;

	my $line_num = 0;
	my %matrix_hash;

	#PARSE THE FILE AND STORE IN THE HASH
	open FILE , $file;

	
	foreach my $line (<FILE>)
	{

		chomp $line;


		my @line = split /\t/ , $line;

		for(my $col = 0; $col < $line_num; $col ++)
		{


			$matrix_hash{$line_num}{$col} = $line[$col];
			$matrix_hash{$col}{$line_num} = $line[$col];

		}#end for

	
		$line_num ++;

	}#end foreach


	close FILE;


	return \%matrix_hash;

}#end LT_matrix_file_to_hash


1;
