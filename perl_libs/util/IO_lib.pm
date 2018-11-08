package IO_lib;


#sub get_params_from_file - Takes as input 1) The name of a file containing parameter assignments in
#					      which the first column is the name of a parameter and \
#					      the second column is its desired value
#			  - Returns a hash reference in which the keys is the parameter name and 
#			    the values its assigned value in the file
sub get_params_from_file
{
	my ($param_file) = @_;
	my %param_hash;


	open PARAMS , $param_file;


	foreach my $param_pair (<PARAMS>)
	{

		chomp $param_pair;

		my ($param_name , $param_value) = split /\t/ , $param_pair;

		$param_hash{$param_name} = $param_value;	
		

	}#end forach



	return \%param_hash;


}#end get_params_from_file



1;

