package db_lib;


use DBI;


#sub DB_connect - Establishes a database connection and returns a handle
sub DB_connect
{

	my ($data_source , $user, $password) = @_;

        #CONNECT TO DB
        my $dbh = DBI -> connect($data_source, $user, $password)
                         or die "Can't connect to $data_source: $DBI::errstr";

        return $dbh;

}#end DB_connect


#sub add_place_holders - Takes as input a partailly completed query and the number of place holders t0 ass
#		       - Returns the same query string with the desired number of ? added on the the end
sub add_place_holders
{
	my ($query , $num_ph) = @_;

	if($num_ph > 0)
	{

		$query = $query . "? ";
	
	}else{

		return $query;

	}#end if


	for(my $i = 1; $i < $num_ph; $i++)
	{

		$query = $query . ", ?";

	}#end for

	return $query;

}#end add_place_holders

1;
