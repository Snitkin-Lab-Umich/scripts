package gene_annotation_db_lib;

use DBI;
use db_lib;
use list_ops;


###################################################### DB CONNECTION #######################################################

#sub gene_annotation_db_connect - Establishes a database connection and returns a handle
sub gene_annotation_db_connect
{
        my ($data_source) = @_;

	$data_source = "dbi:mysql:snitkines_gene_annotation;host=biobase";
        my $user = "snitkines";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end FBA_model_DB_connect


################################################### BASIC DB INSERTION ROUTINES ###########################################

#sub populate_COG_org - Takes as input: 1) cog_org_hash - A reference to a hash whose key is an org_id and values are other
#							  table columns
#		      - Populates table
sub populate_COG_org
{
	my ($r_org_hash, $dbh) = @_;

	my %org_hash = %$r_org_hash;


        #INSERT ALL ORGS INTO  COG_org TABLE 
        foreach my $org_id (sort keys %org_hash)
        {


        	#PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_org = "INSERT INTO COG_org (COG_org_id , org_class, tax_id , species)
                                  VALUES (? , ?, ?, ?)";


                my $sth = $dbh -> prepare($insert_org)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($org_id, $org_hash{$org_id}{'CLASS'}, $org_hash{$org_id}{'TAX ID'}, $org_hash{$org_id}{'SPECIES'})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end populate_COG_org


#sub populate_COG_func - Takes as input: 1) cog_func_hash - A reference to a hash whose key is an func code and values are other
#                                                           table columns
#                      - Populates table
sub populate_COG_func
{
	my ($r_func_hash, $dbh) = @_;

	my %func_hash = %$r_func_hash;


       #INSERT ALL FUNCTIONS INTO COG_func TABLE 
        foreach my $func_code (sort keys %func_hash)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_func = "INSERT INTO COG_func (COG_func_code , class, description)
                                  VALUES (? , ?, ?)";


                my $sth = $dbh -> prepare($insert_func)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($func_code, $func_hash{$func_code}{'CLASS'}, $func_hash{$func_code}{'DESC'})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach



}#end populate_COG_func


#sub populate_COG_cog2gi - Takes as input: 1) cog2gi_hash - A reference to a hash whose key is cog gene id  and values are NCBI gi
#                        - Populates table
sub populate_COG_cog2gi
{
	my ($r_cog2gi_hash, $dbh) = @_;

	my %cog2gi_hash = %$r_cog2gi_hash;


	#INSERT ALL RELATIONSHIPS INTO COG_prot_ID_to_gi TABLE 
        foreach my $prot_id (sort keys %cog2gi_hash)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_cog2gi = "INSERT INTO COG_prot_id_to_gi (COG_prot_id , gi_id)
                                     VALUES (? , ?)";


                my $sth = $dbh -> prepare($insert_cog2gi)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($prot_id, $cog2gi_hash{$prot_id})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach



}#end populate_COG_cog2gi


#sub populate_COG_info - Takes as input: 1) cog_info_hash - A reference to a hash whose key is a COG ID and values are other
#                                                           table columns
#                      - Populates table
sub populate_COG_info
{
	my ($r_cog_info_hash, $dbh) = @_;

	my %cog_info_hash = %$r_cog_info_hash;


	#INSERT ALL FUNCTIONS INTO COG_info TABLE 
        foreach my $cog_id (sort keys %cog_info_hash)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_cog_info = "INSERT INTO COG_info (COG_id, description, COG_func_code)
                                       VALUES (? , ?, ?)";


                my $sth = $dbh -> prepare($insert_cog_info)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($cog_id, $cog_info_hash{$cog_id}{'DESC'}, $cog_info_hash{$cog_id}{'FUNC'})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end populate_COG_info


#sub populate_COG_annot - Takes as input: 1) cog_annot_hash - A reference to a hash whose key is a COG gene id and values are other
#                                                             table collumns reevent to its COG annotation
#                       - Populates table
sub populate_COG_annot
{
	my ($r_cog_annot_hash, $dbh) = @_;

	my %cog_annot_hash = %$r_cog_annot_hash;


	#INSERT ALL COG ASSIGNMENTS INTO COG_annotations TABLE 
        foreach my $prot_id (sort keys %cog_annot_hash)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_cog = "INSERT INTO COG_annotation (COG_id , COG_prot_id, COG_org_id)
                                  VALUES (? , ?, ?)";


                my $sth = $dbh -> prepare($insert_cog)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($cog_annot_hash{$prot_id}{'COG'}, $prot_id, $cog_annot_hash{$prot_id}{'SPECIES'})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end populate_COG_annot


###############################################################BASIC SELECTION ROUTINES###################################################

#sub get_prot2cog - Takes as input: 1) dbh	- A database handle
#		  - Returns a hash where the keys are proteins IDs and the values are COGs
sub get_prot2cog
{
	my ($dbh) = @_;

	my %prot2cog_hash;


       #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cogs = "SELECT COG_prot_id, COG_id
                            FROM COG_annotation";


        my $sth = $dbh -> prepare($select_cogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($prot_id, $COG_id)  =  $sth -> fetchrow_array)
        {

                $prot2cog_hash{$prot_id} = $COG_id;

        }#end while


	return \%prot2cog_hash;

}#end get_prot2cog


#sub get_prot2org - Takes as input: 1) dbh      - A database handle
#                 - Returns a hash where the keys are proteins IDs and the values are org IDs
sub get_prot2org
{
        my ($dbh) = @_;

        my %prot2org_hash;


       #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cogs = "SELECT COG_prot_id, COG_org_id
                            FROM COG_annotation";


        my $sth = $dbh -> prepare($select_cogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($prot_id, $COG_org_id)  =  $sth -> fetchrow_array)
        {

                $prot2org_hash{$prot_id} = $COG_org_id;

        }#end while


        return \%prot2org_hash;

}#end get_prot2cog

#sub get_cog_info - Takes as input: 1) dbh      - A database handle
#                 - Returns a hash where the keys are COG ids and the values are codes/descriptions
sub get_cog_info
{
        my ($dbh) = @_;

        my %cog_info_hash;


       #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cogs = "SELECT COG_id, description, COG_func_code
                            FROM COG_info";


        my $sth = $dbh -> prepare($select_cogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($COG_id, $desc, $func_code)  =  $sth -> fetchrow_array)
        {

                $cog_info_hash{$COG_id}{'FUNC'} = $func_code;
                $cog_info_hash{$COG_id}{'DESC'} = $desc;

        }#end while


        return \%cog_info_hash;

}#end get_cog_info


#sub get_cog_codes - Takes as input: 1) dbh      - A database handle
#                   - Returns a hash where the keys are COG function codes and the values are descriptions
sub get_cog_codes
{
        my ($dbh) = @_;

        my %cog_code_hash;


       #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cogs = "SELECT description, COG_func_code, class
                            FROM COG_func";


        my $sth = $dbh -> prepare($select_cogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($desc, $func_code, $class)  =  $sth -> fetchrow_array)
        {

                $cog_code_hash{$func_code}{'DESC'} = $desc;
                $cog_code_hash{$func_code}{'CLASS'} = $class;

        }#end while


        return \%cog_code_hash;

}#end get_cog_codes


#sub get_pc_info - Takes as input: 1) dbh- A database handle
#		 - Returns a hash where the keys are protein cluster asseccions and the values are relevent info
sub get_pc_info
{
	my ($dbh) = @_;


        my %pc_info_hash;


       #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_pc = "SELECT ID ,accession ,proteins,conserved_in ,genera ,organisms,paralogs,publications,description,COG
                         FROM PC_info";


        my $sth = $dbh -> prepare($select_pc)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my (@pc_info)  =  $sth -> fetchrow_array)
        {

                $pc_info_hash{$pc_info[1]}{'ID'} = $pc_info[0];
                $pc_info_hash{$pc_info[1]}{'PROTS'} = $pc_info[2];
                $pc_info_hash{$pc_info[1]}{'CONS_IN'} = $pc_info[3];
                $pc_info_hash{$pc_info[1]}{'GENERA'} = $pc_info[4];
                $pc_info_hash{$pc_info[1]}{'ORG'} = $pc_info[5];
                $pc_info_hash{$pc_info[1]}{'PARALOG'} = $pc_info[6];
                $pc_info_hash{$pc_info[1]}{'PUB'} = $pc_info[7];
                $pc_info_hash{$pc_info[1]}{'DESC'} = $pc_info[8];
                $pc_info_hash{$pc_info[1]}{'COG'} = $pc_info[9];

        }#end while


        return \%pc_info_hash;


}#end get_pc_info

1;
