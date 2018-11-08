package KO_ANALYSIS_db_lib;



use DBI;
use db_lib;


###################################################### DB CONNECTION #######################################################

#sub DB_connect - Establishes a database connection and returns a handle
sub DB_connect
{
        my $data_source = "dbi:mysql:SC_KO_ANALYSIS";
        #my ($data_source) = @_;

        my $user = "root";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end DB_connect


################################################### BASIC DB INSERTION ROUTINES ###########################################

#sub SGD_info_to_gene - Takes as input : 1) name 	- Common name of gene
#					 2) ORF_size 	- The nucleotide size of the ORF
#				         3) GO_BP	- GO biological process of protein
#				         4) GO_MF	- GO molecular function of protein
#				         5) GO_CC	- GO cellular component of protein
#                                        6) ORF_id 	- The orf id of the gene which is already in thr table
#                                        7) SGD_id 	- THe SGD ID of the gene
#                                        8) dbh 	- A database handle
#
#                     - Puts the gene information into the gene table
sub SGD_info_to_gene
{
        my ($ORF_id , $name , $ORF_size , $GO_BP , $GO_MF , $GO_CC , $SGD_id, $dbh) = @_;


        #INSERT GENE INFO INTO DATABASE
        my $insert_gene = "INSERT INTO gene (ORF_id , name ,  ORF_size , GO_BP , GO_MF , GO_CC , SGD_id)
			   VALUES (? , ? ,? ,? ,? ,? , ?)";	


        my $sth = $dbh -> prepare($insert_gene)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($ORF_id , $name , $ORF_size , $GO_BP , $GO_MF , $GO_CC , $SGD_id)
                         or die "Can't execute statment: $DBI::errstr";


}#end SGD_info_to_gene


#sub insert_KO_data - Takes as intput: 1) exp_vals	- A reference to a list of experimental values
#				       2) ORF		- The ORF ID of the gene KO
#				       3) conds		- A list of conditions which are columns headers in the table
#							  and also in the same order as the exp_vals
#				       4) table_name	- The name of the table to insert the data into
#				       5) dbh		- A database handle
#
#		    - Puts KO data into appropriate table
sub insert_KO_data
{
	my ( $r_exp_vals , $ORF , $r_conds , $table_name , $dbh) = @_;
	my $gid;
	my $conds = join " , " , @$r_conds;

	
	#GET THE GID FOR THE CURRENT ORF
	my $select_gid = "SELECT (g.gid)
                          FROM gene g
                          WHERE g.ORF_id = ?";


        my $sth = $dbh -> prepare($select_gid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ORF)
                         or die "Can't execute statment: $DBI::errstr";

        $gid =  $sth -> fetchrow_array;


	#INSERT THE EXPERIMETNAL DATA INTO THE DATABASE
	if ($gid)
	{


		my $insert_KO = "INSERT INTO $table_name (gid , $conds )
			         VALUES (? ,"; 

		$insert_KO = &db_lib::add_place_holders($insert_KO , scalar(@$r_conds));

		$insert_KO = $insert_KO . ")";
	

		my $sth = $dbh -> prepare($insert_KO)
        	                  or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($gid , @$r_exp_vals)
        	                 or die "Can't execute statment: $DBI::errstr";

	}#end if


}#end insert_KO_data



################################################### BASIC DB SELECTION ROUTINES ###########################################

#sub get_KO_data - Takes as input: 1) table	- THe name of a table in the database
#				   2) condition	- The name of a condition field in the table
#				   3) dbh	- A database handle
#		 - Returns a hash where the keys are gene names and the values are deletion phenotypes
sub get_KO_data
{
	my ($table, $condition, $dbh) = @_;

	my %KO_data;

       #SELECT KO DATA FROM RELEVENT TABLE UNDER RELEVENT CONDITION
        my $select_KO = "SELECT g.ORF_id, ko.$condition
                           FROM gene g, $table ko
                           WHERE g.gid = ko.gid";

        my $sth = $dbh -> prepare($select_KO)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($ORF_id, $KO_phen) =  $sth -> fetchrow_array)
        {

                $KO_data{$ORF_id} = $KO_phen;

        }#end while

	return \%KO_data;

}#end get_KO_data


#sub protName2gid - Takes as input: 1) prot_names	- A list of protein names
#				    2) dbh		- A database handle
#		    Returns gids corresponding to the protein names inputted 
sub protName2gid
{
	my ($r_prot_names , $dbh) = @_;
	my @gids;
	
	#SELECT ALL GIDS CORRESPONDING TO THE LIST OF INPUTTED PROTEINS
        my $select_gids = "SELECT DISTINCT g.gid
                           FROM gene g
                           WHERE g.name IN (";

        #ADD PLACE HOLDERS FOR GENES
        $select_gids = &db_lib::add_place_holders($select_gids , scalar(@$r_prot_names));

        $select_gids = $select_gids . ") OR
			        g.ORF_id IN (";

	$select_gids = &db_lib::add_place_holders($select_gids , scalar(@$r_prot_names));

        $select_gids = $select_gids . ")";


        my $sth = $dbh -> prepare($select_gids)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_prot_names , @$r_prot_names)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid) =  $sth -> fetchrow_array)
        {

                push @gids , $gid;

        }#end while

	if( scalar(@$r_prot_names) != scalar(@gids) )
	{

		print "Prot anmes = @$r_prot_names\n"

	}

	
	return \@gids;

}#end protName2gid


1;
