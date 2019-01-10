package hash_lib;

BEGIN {push @INC, '/nfs/esnitkin/bin_group/perl_libs/math'}

use arith_ops;

#sub dict_to_hash - Takes as input the name of a dictionary file which is in two
#                   column tab delimited format
#                 - Returns a hash with column one as the key and column
#                   two as the value
sub dict_to_hash
{
        my ($dict_file) = @_;
        my %dict_hash;

        open INPUT , $dict_file;

        foreach my $entry (<INPUT>)
        {
                my($key , $value) = split /\t/ , $entry;
		chomp($value);
                $dict_hash{$key} =$value;
        }#end foreach

        close INPUT;

        return \%dict_hash;

}#end dict_to_hash


#sub dict_to_2D_hash - Takes as input a 3 column file
#	             - Returns a reference to a hash with the first 2 cols being keys and
#	               the third being the value
sub dict_to_2Dhash
{
	my ($dict_file) = @_;
        my %dict_hash;

        open INPUT , $dict_file;

        foreach my $entry (<INPUT>)
        {
                my($key1 , $key2, $value) = split /\t/ , $entry;
                chomp($value);
                $dict_hash{$key1}{$key2} =$value;
        }#end foreach

        close INPUT;

        return \%dict_hash;

}#end sub dict_to_2Dhash


#sub dict_to_hash_valC1 - Takes as input the name of a dictionary file which is in two
#                         column tab delimited format
#                       - Returns a hash with column one as the value and the part of column
#                         two before the first space as the key
sub dict_to_hash_valC1
{
        my ($dict_file) = @_;
        my %dict_hash;

        open INPUT , $dict_file;

        foreach my $entry (<INPUT>)
        {
                my($value,$key) = split /\t/ , $entry;
                my @key = split /\s+/ , $key;
                $dict_hash{$key[0]} =$value;
        }#end foreach

        close INPUT;

        return \%dict_hash;

}#end dict_to_hash_valC1


#sub dict_to_hash_keyC1 - Takes as input the name of a dictionary file which is in two
#                         column tab delimited format
#                       - Returns a hash with column one as the key and the part of column
#                         two before the first space as the value
sub dict_to_hash_keyC1
{
        my ($dict_file) = @_;
        my %dict_hash;

        open INPUT , $dict_file;

        foreach my $entry (<INPUT>)
        {
                my($key , $value) = split /\t/ , $entry;
                my @value = split /\s+/ , $value;
                $dict_hash{$key} =$value[0];
        }#end foreach

        close INPUT;

        return \%dict_hash;

}#end dict_to_hash_keyC1


#sub hash_of_lists_to_file - Takes as input a hash containing array references as values
#                 	   - Appends the result results in the hash to an output file
sub hash_of_lists_to_file
{
        my($output_file , $r_result_hash) = @_;
        %result_hash = %$r_result_hash;

        open OUTPUT , ">>$output_file";

        foreach $key(sort keys %result_hash)
        {
                if((@{$result_hash{$key}}) != 0)
                {
                        foreach my $result ( @{$result_hash{$key}} )
                        {
                                print OUTPUT $result;
                        }
                }
        }#end foreach

        close OUTPUT;

}#end hash_of_lists_to_file


#sub twoD_hash_to_file - Takes as input: 1) r_hash 	- A reference to a 2D hash
#					 2) output_file	- A file to write the hash to
#		       - Writes the 2D hash to a file in the form of a matrix
sub twoD_hash_to_file
{
	my ($r_hash, $output_file) = @_;

	my %hash = %$r_hash;


	#PRINT OUT MATRIX
	open OUTPUT, ">$output_file";

	#PRINT HEADER
	my $cols = join "\t", (sort keys %hash);
	print OUTPUT "\t$cols\n";

	#PRINT MATRIC
	foreach my $key1 (sort keys %hash)
	{

	        print OUTPUT "$key1";

	        foreach my $key2 (sort keys %hash)
	        {

	                print OUTPUT "\t" . $hash{$key1}{$key2};

	        }#end foreach

	        print OUTPUT "\n";

	}#end foreach


	close OUTPUT;

}#end twoD_hash_to_file


#sub file_to_hash - Takes as input a single column file
#		  - Returns a reference to a hash with the file entries as keys and a value of 1
sub file_to_hash
{

	my ($file) = @_;
	my %hash;

	open INPUT , $file;

	foreach $line (<INPUT>)
	{

		chomp($line);

		$hash{$line} = 1;
	
	}#end foreach

	return \%hash;

}#end file_to_hash


#sub init_hash - Takes as input: 1) r_keys	- A list of keys
#	       - Returns a reference to a hash whose keys are keys and values are 0
sub init_hash
{
	my ($r_keys) = @_;

	my %hash;

	foreach my $key (@$r_keys)
	{

		$hash{$key} = 0;

	}#end foreach


	return \%hash;

}#end init_hash


#sub hash_sum - Takes as input: 1) r_hash	- A reference to a hash
#	      - Returns a sum of all the elements in the hash
sub hash_sum
{
	my ($r_hash) = @_;

	my %hash = %$r_hash;
	my $hash_sum = 0;


	foreach my $key (keys %hash)
	{

		$hash_sum += $hash{$key};

	}#end foreach


	return $hash_sum;

}#end hash_sum


#sub hash_add - Takes as input - 1) r_hash1	 - A reference to a hash
#	      			 2) r_hash2	 - A reference to a second hash
#	      - Returns a hash that represents the sum of elements in the two hashes that share the same key
sub hash_add
{
	my ($r_hash1, $r_hash2) = @_;

	my %hash1 = %$r_hash1;
	my %hash2 = %$r_hash2;

	my @all_keys = (keys(%hash1), keys(%hash2));

	my $r_add_hash = &init_hash(\@all_keys);
	my %add_hash = %$r_add_hash;


	foreach my $key (keys %add_hash)
	{

		$add_hash{$key} = $hash1{$key} + $hash2{$key};

	}#end foreach
	

	return \%add_hash;

}#end hash_add


1;
