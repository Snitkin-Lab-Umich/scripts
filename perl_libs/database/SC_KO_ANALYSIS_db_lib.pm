package SC_KO_ANALYSIS_db_lib;



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



#sub insert_mips_cid2complex - Takes as input: 1) cid		- The MIPS ID of the complex
#					       2) complex name	- The name of the MIPS complex
#					       3) dbh		- A database handle
#			     - Inserts data into MIPS_cid2complex table
sub insert_mips_cid2complex
{
	my ($cid , $cname , $dbh) = @_;


	#INSERT DATA
	my $insert_mips = "INSERT INTO MIPS_cid2complex (mips_code , mips_cname )
                                 VALUES (? , ?)";


        my $sth = $dbh -> prepare($insert_mips)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($cid , $cname)
                          or die "Can't execute statment: $DBI::errstr";



}#end insert_mips_cid2complex



#sub insert_mips_cid2gid - Takes as input: 1) cid	- The ID of the MIPS complex
#					   2) ORF_ID	- The ORF ID of the gene part of the complex
#				           3) dbh	- A database handle
#			 - Inserts the data into the MIPS_gid2cid table
sub insert_mips_cid2gid
{
	my ($cid , $ORF_ID , $dbh)= @_;


	#GET THE DATABASE GENE ID FOR THE CURRENT ORF
        my $select_gid = "SELECT (g.gid)
                          FROM gene g
                          WHERE g.ORF_id = ?";


        my $sth = $dbh -> prepare($select_gid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ORF_ID)
                         or die "Can't execute statment: $DBI::errstr";

        $gid =  $sth -> fetchrow_array;


	#GET THE DATABASE CID FOR THE CURRENT COMPLEX
        my $select_cid = "SELECT (m.mips_cid)
                          FROM MIPS_cid2complex m
                          WHERE m.mips_code = ?";


        my $sth = $dbh -> prepare($select_cid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($cid)
                         or die "Can't execute statment: $DBI::errstr";

        $db_cid =  $sth -> fetchrow_array;

	print "gid = $gid  cid = $db_cid\n";


	#INSERT THE DATA IF THE GENE IS IN THE DATABASE
	if($gid)
	{

	        my $insert_mips = "INSERT INTO MIPS_gid2cid (mips_cid , gid)
	                                 VALUES (? , ?)";


	        my $sth = $dbh -> prepare($insert_mips)
	                          or die "Can't prepare statment: $DBI::errstr";


	        my $rc = $sth -> execute($db_cid , $gid)
	                          or die "Can't execute statment: $DBI::errstr";


	}#end if



}#end insert_mips_cid2gid


#sub insert_gavin06_module - Takes as input : 1) mod_id		- The gavin ID of the module
#					      2) gene_names	- A | delimited list of gene ids in the module
#					      3) dbh		- A database handle
#			   - Inserts data into the GAVIN06_module table
sub insert_gavin06_module
{
	my ($mod_id , $gene_names , $dbh) = @_;


        #INSERT DATA
        my $insert_mod = "INSERT INTO GAVIN06_module (gavin_mid , mod_prots)
                                 VALUES (? , ?)";


        my $sth = $dbh -> prepare($insert_mod)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($mod_id , $gene_names)
                          or die "Can't execute statment: $DBI::errstr";


}#end insert_gavin06_module


#sub insert_gavin06_complex - Takes as input : 1) complex_id	- Gavin complex ID
#					       2) complex_name	- Gavin complex name
#					       3) core_prots	- A | delimited list of gids
#					       4) mods		- A | delimited list of module ids
#					       5) att_prots	- A bar delimited list of gids
#					       6) loc		- Cellular localization
#					       7) dbh		- Database handle
#                           - Inserts data into the GAVIN06_complex table
sub insert_gavin06_complex
{
        my ($complex_id , $complex_name , $core_prots , $mods , $att_prots , $loc ,  $dbh) = @_;


        #INSERT DATA
        my $insert_com = "INSERT INTO GAVIN06_complex (gavin_cid , gavin_cname , core_prots , mod_prots , att_prots , loc)
                          VALUES (? , ? , ? , ? , ? , ?)";


        my $sth = $dbh -> prepare($insert_com)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($complex_id , $complex_name , $core_prots , $mods , $att_prots , $loc)
                          or die "Can't execute statment: $DBI::errstr";


}#end insert_gavin06_module


################################################### BASIC DB SELECTION ROUTINES ###########################################

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



#sub get_gavin_cid2complex - Takes as input: 1) dbh	- A database handle 
#
#			   - Returns a reference to a hash where the key is the gavin complex ID and the 
#			     values are the corresponding complex descriptions
sub get_gavin_cid2complex
{
	my ($dbh) = @_;
	my %cid2complex;


	#SELECT ALL CID TO COMPLEX PAIRS
        my $select_cid2complex = "SELECT DISTINCT gc.gavin_cid , gc.gavin_cname
                           	  FROM GAVIN06_complex gc";


        my $sth = $dbh -> prepare($select_cid2complex)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($cid , $cname) =  $sth -> fetchrow_array)
        {

                $cid2complex{$cid} = $cname;

        }#end while


        return \%cid2complex;


}#end get_gavin_cid2complex


#sub get_gavin_cid2orfs - Takes as input: 1) dbh 	- A database handle
#
#			- Returns a reference to a hash where the keys are gavin complex ids and the values
#			  are references to lists of ORFS which make up the core complexes  
sub get_gavin_cid2orfs
{
	my ($dbh) = @_;

        my %cid2core;
        my %cid2orfs;


        #SELECT ALL CIDS AND ASSOCIATED CORE COMPLEX COMPONENTS
        my $select_cid2complex = "SELECT DISTINCT gc.gavin_cid , gc.core_prots
                                  FROM GAVIN06_complex gc";


        my $sth = $dbh -> prepare($select_cid2complex)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($cid , $core_prots) =  $sth -> fetchrow_array)
        {

                $cid2core{$cid} = $core_prots;

        }#end while



	#FOREACH COMPLEX GET ORFS
	foreach my $cid (keys %cid2core)
	{
		

		$cid2core{$cid} =~ s/^\|(.*)\|$/$1/;
		my @ORFS = split /\|/ ,  $cid2core{$cid}; 


	        #SELECT ALL GIDS CORRESPONDING TO THE LIST OF INPUTTED PROTEINS
	        my $select_orfs = "SELECT DISTINCT g.ORF_id
	                           FROM gene g
	                           WHERE g.gid IN (";

	        #ADD PLACE HOLDERS FOR GENES
	        $select_orfs = &db_lib::add_place_holders($select_orfs , scalar(@ORFS));

	        $select_orfs = $select_orfs . ")";


	        my $sth = $dbh -> prepare($select_orfs)
	                         or die "Can't prepare statment: $DBI::errstr";

	        my $rc = $sth -> execute(@ORFS)
	                         or die "Can't execute statment: $DBI::errstr";



        	while(my ($ORF_id) =  $sth -> fetchrow_array)
        	{

                	push @{$cid2orfs{$cid}} , $ORF_id;

        	}#end while


	}#end foreach


        return \%cid2orfs;


}#end get_gavin_cid2orfs


#sub get_gavin_cid2atts - Takes as input: 1) dbh 	- A database handle
#
#			- Returns a reference to a hash where the keys are gavin complex ids and the values
#			  are references to lists of ORFS which make up the complex attatchments 
sub get_gavin_cid2atts
{
	my ($dbh) = @_;

        my %cid2atts;
        my %cid2orfs;


        #SELECT ALL CIDS AND ASSOCIATED CORE COMPLEX COMPONENTS
        my $select_cid2atts = "SELECT DISTINCT gc.gavin_cid , gc.att_prots
                                  FROM GAVIN06_complex gc";


        my $sth = $dbh -> prepare($select_cid2atts)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($cid , $att_prots) =  $sth -> fetchrow_array)
        {

                $cid2atts{$cid} = $att_prots;

        }#end while



	#FOREACH COMPLEX GET ORFS
	foreach my $cid (keys %cid2atts)
	{
		

		$cid2atts{$cid} =~ s/^\|(.*)\|$/$1/;
		my @ORFS = split /\|/ ,  $cid2atts{$cid}; 


	        #SELECT ALL GIDS CORRESPONDING TO THE LIST OF INPUTTED PROTEINS
	        my $select_orfs = "SELECT DISTINCT g.ORF_id
	                           FROM gene g
	                           WHERE g.gid IN (";

	        #ADD PLACE HOLDERS FOR GENES
	        $select_orfs = &db_lib::add_place_holders($select_orfs , scalar(@ORFS));

	        $select_orfs = $select_orfs . ")";


	        my $sth = $dbh -> prepare($select_orfs)
	                         or die "Can't prepare statment: $DBI::errstr";

	        my $rc = $sth -> execute(@ORFS)
	                         or die "Can't execute statment: $DBI::errstr";



        	while(my ($ORF_id) =  $sth -> fetchrow_array)
        	{

                	push @{$cid2orfs{$cid}} , $ORF_id;

        	}#end while


	}#end foreach


        return \%cid2orfs;


}#end get_gavin_cid2atts

#sub get_gavin_mid2orfs - Takes as input: 1) dbh 	- A database handle
#
#			- Returns a reference to a hash where the keys are gavin module ids and the values
#			  are references to lists of ORFS which make up the modules 
sub get_gavin_mid2orfs
{
	my ($dbh) = @_;

        my %mid2prots;
        my %mid2orfs;


        #SELECT ALL CIDS AND ASSOCIATED CORE COMPLEX COMPONENTS
        my $select_mid2orfs = "SELECT DISTINCT gm.gavin_mid , gm.mod_prots
                                  FROM GAVIN06_module gm";


        my $sth = $dbh -> prepare($select_mid2orfs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($mid , $mod_prots) =  $sth -> fetchrow_array)
        {

                $mid2prots{$mid} = $mod_prots;

        }#end while



	#FOREACH COMPLEX GET ORFS
	foreach my $mid (keys %mid2prots)
	{
		

		my @ORFS = split /\|/ ,  $mid2prots{$mid}; 


	        #SELECT ALL GIDS CORRESPONDING TO THE LIST OF INPUTTED PROTEINS
	        my $select_orfs = "SELECT DISTINCT g.ORF_id
	                           FROM gene g
	                           WHERE g.gid IN (";

	        #ADD PLACE HOLDERS FOR GENES
	        $select_orfs = &db_lib::add_place_holders($select_orfs , scalar(@ORFS));

	        $select_orfs = $select_orfs . ")";


	        my $sth = $dbh -> prepare($select_orfs)
	                         or die "Can't prepare statment: $DBI::errstr";

	        my $rc = $sth -> execute(@ORFS)
	                         or die "Can't execute statment: $DBI::errstr";



        	while(my ($ORF_id) =  $sth -> fetchrow_array)
        	{

                	push @{$mid2orfs{$mid}} , $ORF_id;

        	}#end while


	}#end foreach


        return \%mid2orfs;


}#end get_gavin_mid2orfs




#sub get_gavin_cid2mids - Takes as input: 1) dbh 	- A database handle
#
#			- Returns a reference to a hash where the keys are gavin complex ids and the values
#			  are references to lists of gavin module IDs
sub get_gavin_cid2mids
{
	my ($dbh) = @_;

        my %cid2mids;


        #SELECT ALL CIDS AND ASSOCIATED CORE COMPLEX COMPONENTS
        my $select_cid2module = "SELECT DISTINCT gc.gavin_cid , gc.mod_prots
                                  FROM GAVIN06_complex gc";


        my $sth = $dbh -> prepare($select_cid2module)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($cid , $mod_prots) =  $sth -> fetchrow_array)
        {

		$mod_prots =~ s/^\|(.*)\|$/$1/;

                @{$cid2mids{$cid}} = split /\|/ , $mod_prots;

        }#end while


	return \%cid2mids;

}#end get_gavin_cid2mids



1;
