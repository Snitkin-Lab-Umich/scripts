package visML_coords_db_lib;


use DBI;
use db_lib;


###################################################### DB CONNECTION #######################################################

#sub visML_coords_DB_connect - Establishes a database connection and returns a handle
sub visML_coords_DB_connect
{
        $data_source = "dbi:mysql:visML_coords";

        my $user = "root";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end visML_coords_DB_connect


################################################### BASIC DB INSERTION ROUTINES ###########################################


#sub visML_coords_to_DB - Takes as input: 1) dbh	- A database handle
#					  2) coord_hash	- A hash whose key is the node name and value is X/Y coordinates in 
#							  the network layout
#					  3) model	- The model for which this layout applies to
#					  4) layout	- The name of the layout being inserted
#			- Places the data into the appropriate database tables
sub visML_coords_to_DB
{
	my ($dbh , $r_coord_hash, $model , $layout) = @_;

	my %coord_hash = %$r_coord_hash;

	my $lid;
	my $mid;


	#PLACES THE MODEL NAME INTO THE MODEL TABLE, IF IT DOES NOT ALREADY EXIST
	my $select_model = "SELECT (m.model_id)
                            FROM model m
                            WHERE m.model_name = ?";


        my $sth = $dbh -> prepare($select_model)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model)
                         or die "Can't execute statment: $DBI::errstr";

        $mid =  $sth -> fetchrow_array;


	#INSERT THE MODEL NAME INTO THE DATABASE IF IT IS NOT ALREADY PRESENT
	if( $mid eq "" )
	{

	        my $insert_model = "INSERT INTO model (model_name)
                                    VALUES (?)";


                my $sth = $dbh -> prepare($insert_model)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($model)
                                         or die "Can't execute statment: $DBI::errstr";

		$mid = $dbh->{ q{mysql_insertid}};

	}#end if


	#PLACES THE LAYOUT NAME INTO THE LAYOUT TABLE IF IT DOES NOT ALREADY EXIST
        my $select_layout = "SELECT (l.layout_id)
                             FROM layout l
                             WHERE l.layout_name = ?";


        my $sth = $dbh -> prepare($select_layout)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($layout)
                         or die "Can't execute statment: $DBI::errstr";

        $lid =  $sth -> fetchrow_array;


        #INSERT THE LAYOUT NAME INTO THE DATABASE IF IT IS NOT ALREADY PRESENT
        if( $lid eq "" )
        {

                my $insert_layout = "INSERT INTO layout (layout_name)
                                     VALUES (?)";


                my $sth = $dbh -> prepare($insert_layout)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($layout)
                                         or die "Can't execute statment: $DBI::errstr";

		$lid = $dbh->{ q{mysql_insertid}};

        }#end if



	#PLACES THE COORDINATES INTO THE COORDINATES TABLE
	foreach my $node_name (keys %coord_hash)
	{

		my $X_coord = $coord_hash{$node_name}->{"X"};
		my $Y_coord = $coord_hash{$node_name}->{"Y"};


                my $insert_coords = "INSERT INTO coordinates (model_id , layout_id, node_name, X_coord, Y_coord)
                                     VALUES (?,?,?,?,?)";


                my $sth = $dbh -> prepare($insert_coords)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($mid, $lid, $node_name,$X_coord, $Y_coord)
                                         or die "Can't execute statment: $DBI::errstr";



	}#end foreach


}#end visML_coords_to_DB



################################################### BASIC DB SELECTION ROUTINES ###########################################

#sub get_visML_coords - Takes as input: 1) dbh		- A database handle
#					2) layout	- The name of a layout in the database
#					3) model	- THe name of a model for which layouts have been stored
#		      - Returns a hash in which the keys are node names and the values are X and Y coordinates
sub get_visML_coords
{
	my ($dbh , $layout, $model) = @_;

	my %node_coords;

	
	#SELECT ALL NODE POSITIONS FOR THE DESIRED LAYOUT AND MODEL
	my $select_coords = "SELECT c.node_name, c.X_coord, c.Y_coord
                             FROM coordinates c, layout l, model m
                             WHERE c.model_id = m.model_id AND
				   c.layout_id = l.layout_id AND 
				   l.layout_name = ? AND
				   m.model_name = ?";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($layout , $model)
                         or die "Can't execute statment: $DBI::errstr";


	#GET THE COORDINATES FOR EACH NODE
	while(my ($node, $X_coord,  $Y_coord) =  $sth -> fetchrow_array)
        {

                $node_coords{$node} -> {"X"} = $X_coord;;
                $node_coords{$node} -> {"Y"} = $Y_coord;;

        }#end while



	return \%node_coords;

}#end get_visML_coords


#sub get_layouts
#
sub get_layouts
{
	my ($dbh) = @_;

        my @layouts;


        #SELECT ALL NODE POSITIONS FOR THE DESIRED LAYOUT AND MODEL
        my $select_layouts = "SELECT l.layout_name
                              FROM layout l";


        my $sth = $dbh -> prepare($select_layouts)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET THE COORDINATES FOR EACH NODE
        while(my ($layout) =  $sth -> fetchrow_array)
        {

		push @layouts , $layout;

        }#end while


        return \@layouts


}#end get_layouts




1;
