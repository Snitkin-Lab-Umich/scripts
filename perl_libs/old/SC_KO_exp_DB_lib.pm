package SC_KO_exp_DB_lib;



use DBI;
use db_lib;


###################################################### DB CONNECTION #######################################################

#sub SC_exp_DB_connect - Establishes a database connection and returns a handle
sub SC_KO_exp_DB_connect
{
        my $data_source = "dbi:mysql:SC_KO_exp";
        #my ($data_source) = @_;

        my $user = "root";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end DB_connect



#sub get_sauer_GR2005_data - Takes as input 1) ORFS - A list of ORD ids to retrieve KO data for
#			   - Returns a hash with the keys being ORF ids and the values being tab
#			     delimited KO data
sub get_sauer_GR2005_data
{
	my ($r_ORFS , $dbh) = @_;
	my @ORFS = @$r_ORFS;

	my %KO_data;

	
	#IF LSIT OF ORFS PROVIDED IS EMPTY THEN SELECT ALL
	my $select_ORFS;
	my $sth;
	my $rc;

	if( (@ORFS) == 0 )
	{


                #GET KO DATA FOR EACH ORF IN THE LIST
                $select_ORFS = "SELECT ORF_id, glucose_aerob , glucose_anaerob , galactose , glycerol , ethanol
                                   FROM sauer_GR2005";


                $sth = $dbh -> prepare($select_ORFS)
                                 or die "Can't prepare statment: $DBI::errstr";

                $rc = $sth -> execute()
                                 or die "Can't execute statment: $DBI::errstr";


	}else{


		#GET KO DATA FOR EACH ORF IN THE LIST
	        $select_ORFS = "SELECT ORF_id, glucose_aerob , glucose_anaerob , galactose , glycerol , ethanol
	                           FROM sauer_GR2005
	                           WHERE ORF_id IN (";


	        #ADD PLACE HOLDERS FOR GENES
	        $select_ORFS = &db_lib::add_place_holders($select_ORFS , scalar(@ORFS));

	        $select_ORFS = $select_ORFS . ")";



        	$sth = $dbh -> prepare($select_ORFS)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	$rc = $sth -> execute(@ORFS)
        	                 or die "Can't execute statment: $DBI::errstr";


	}#end if



        while(my ($ORF_id , $gluc_aer, $gluc_anaer, $galac , $glyc , $etha) =  $sth -> fetchrow_array)
        {


		#REPLACE EMPTY ENTRIES WITH 1'S (VIABLE)
		$gluc_aer =~ s/^$/1/;
		$gluc_anaer =~ s/^$/1/; 
		$galac =~ s/^$/1/;
		$glyc =~ s/^$/1/;
		$etha =~ s/^$/1/;                


		#REPLACE M'S (MUTATIONAL ADJUSTMENT) WITH 0'S (LETHAL)
		$gluc_aer =~ s/m/0/;
		$gluc_anaer =~ s/m/0/; 
		$galac =~ s/m/0/;
		$glyc =~ s/m/0/;
		$etha =~ s/m/0/;


		#MAKE A TAB DELIMITED LIST OF THE DATA
		my $curr_KO_data = join "\t" , ($gluc_aer, $gluc_anaer, $galac , $glyc , $etha);
		$KO_data{$ORF_id} = $curr_KO_data;
           

        }#end while



	return \%KO_data;


}#end get_sauer_GR2005_data



1;
