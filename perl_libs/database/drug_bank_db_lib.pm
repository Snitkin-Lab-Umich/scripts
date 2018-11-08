package drug_bank_db_lib;


use DBI;
use db_lib;
use list_ops;


###################################################### DB CONNECTION #######################################################

#sub drug_bank_db_connect - Establishes a database connection and returns a handle
sub drug_bank_db_connect
{
        #my $data_source = "dbi:mysql:test";
        my ($data_source) = @_;

        my $user = "root";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end drug_bank_db_connect


################################################### CREATE TABLE STATEMENTS ###########################################

#sub create_drug_cards_table - Takes as input a hash linking column IDs to descriptions of the clumns
#			    - Creates a table called drug_cards with the appropriate column names
sub create_drug_cards_table
{
	my ($dbh, $r_drug_card_key_hash) = @_;

	my %drug_card_key_hash = %$r_drug_card_key_hash;


       my $create_table = "CREATE TABLE IF NOT EXISTS drug_cards (\n";

       #ADD REACTION FIELDS
       foreach my $col (sort keys %drug_card_key_hash)
       {

		my $col2 = $col;
		$drug_card_key_hash{$col} =~ s/\'//g;

               $create_table = $create_table . "\t$col2  TEXT COMMENT \'" . $drug_card_key_hash{$col} . "\',\n";

       }#end foeach     

       $create_table = $create_table . "\tdrug_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY
                                                       )";

	print $create_table . "\n\n";

       #CREATE THE TABLE
       my $cth = $dbh -> prepare($create_table)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $cth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


}#end create_drug_cards_table

################################################### BASIC DB INSERTION ROUTINES ###########################################

#sub insert_entry_into_drug_card_table - Takes as input: 1) $drug_card_cols	- A comma delimited list of the columns on the drug card table
#							 2) $r_drug_card	- A list of values to insert in the table, in the same order as 
#										  columns in $drug_card_cols
sub insert_entry_into_drug_card_table
{
	my ($drug_card_cols, $r_drug_card, $dbh) = @_;

	my @drug_card = @$r_drug_card;


	#PREPARE INSERT STATEMENT
        my $insert_drug = "INSERT INTO drug_cards ($drug_card_cols)
                                       VALUES (";

	$insert_drug = &db_lib::add_place_holders($insert_drug , scalar(@drug_card));

	$insert_drug = $insert_drug . ")";


	print "$insert_drug\n";

	#EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($insert_drug)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@drug_card)
                         or die "Can't execute statment: $DBI::errstr";


}#end insert_entry_into_drug_bank_table

#sub insert_entry_into_KO_table - Takes as input: 1) 
sub insert_entry_into_KO_table
{
	my ($KO_id, $name, $definition, $class, $db_links, $genes, $sc_genes, $hs_genes, $dbh) = @_;

	#PREPARE INSERT STATEMENT
        my $insert_KO = "INSERT INTO KEGG_orthology (KO_id,name,definition,class,db_links,genes,sc_genes,hs_genes)
                                       VALUES (?,?,?,?,?,?,?,?)";


	#EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($insert_KO)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($KO_id, $name, $definition, $class, $db_links, $genes, $sc_genes, $hs_genes)
                         or die "Can't execute statment: $DBI::errstr";

}#end insert_entry_into_KO_table


#sub insert_entry_into_OMIM_table
sub insert_entry_into_OMIM_table
{
	my ($OMIM_id, $genes, $chr_loc, $desc, $dbh) = @_;

	#PREPARE INSERT STATEMENT
        my $insert_OMIM = "INSERT INTO morbid_map (OMIM_id, gene_symbols, chr_loc, description)
                                       VALUES (?,?,?,?)";


        #EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($insert_OMIM)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($OMIM_id, $genes, $chr_loc, $desc)
                         or die "Can't execute statment: $DBI::errstr";


}#end insert_entry_into_OMIM_table

################################################DB SELEXTION ROUTINES###########################################

#sub get_ortholog_y2h - Takes as input: 1) dbh		- A database handle
#					2) yeast_ORF 	- A yeast ORF
#		      - Returns a reference to a list of human genes in the same orthology group
sub get_ortholog_y2h
{
	my ($yeast_ORF, $dbh) = @_;

        my @human_genes;


        #GET HUMAN GENES THAT ARE IN SAME KO GROUP AS YEAST ORF
        my $select_genes = "SELECT hs_genes
                           FROM KEGG_orthology
			   WHERE sc_genes LIKE \'%$yeast_ORF%\'";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


	my $human_genes =  $sth -> fetchrow_array;
	@human_genes = split /\|/ , $human_genes;

        return \@human_genes;

}#end get_ortholog_y2h


#sub get_drug_given_target - Takes as input: 1) gene	- A gene name
#					     2) dbh	- A database handle
#			   - Returns the names of all drug that target the given gene
sub get_drug_given_target
{
	my ($gene, $dbh) = @_;

	my %drugs;

	#GET HUMAN GENES THAT ARE IN SAME KO GROUP AS YEAST ORF
        my $select_drugs = "SELECT Drug_ID, Generic_Name
                           FROM drug_cards
                           WHERE (Drug_Target_1_Gene_Name LIKE \'%\n$gene\n%\' OR Drug_Target_1_Gene_Name LIKE \'%\n$gene\' OR Drug_Target_1_Gene_Name LIKE \'$gene\n%\')
				 AND Organisms_Affected LIKE \'%Humans%\'
				 AND Drug_Category NOT LIKE \'%Dietary%\'";

        my $sth = $dbh -> prepare($select_drugs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($drug_id , $name) =  $sth -> fetchrow_array)
        {

                $drugs{$drug_id} = $name;

        }#end while

        return \%drugs;


}# end get_drug_given_target


#sub get_number_of_targets - Takes as input: 1) drug_id	- The id of the drug
#					     2) dbh	- A database handle
#			   - Returns the number of targets associated with the gene
sub get_number_of_targets
{
	my ($drug_id, $dbh) = @_;


        #GET THE TARGET INFO
	my $select_targets = "SELECT Drug_Target_1_Gene_Name 
                              FROM drug_cards
                              WHERE drug_id = ?";

        my $sth = $dbh -> prepare($select_targets)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($drug_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $drug_targets = $sth -> fetchrow_array;

	my @drug_targets = split /\n=\n/, $drug_targets;

	return \@drug_targets;

}#end get_number_of_targets


#sub get_number_of_interactions - Takes as input: 1) drug_id 	- The id of the drug
#						  2) dbh	- A database handle
#				- Returns the number of interactions the drug has
sub get_number_of_interactions
{
	my ($drug_id, $dbh) = @_;

	my $num_interactions = 0;

	#GET THE INTERACTION INFO
	my $select_interactions = "SELECT Drug_Interactions 
                              	   FROM drug_cards
                              	   WHERE drug_id = ?";

        my $sth = $dbh -> prepare($select_interactions)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($drug_id)
                         or die "Can't execute statment: $DBI::errstr";



        $drug_interactions =$sth -> fetchrow_array;

	if($drug_interactions !~ /Not Available/)
	{

		$num_interactions = $drug_interactions =~ tr/\t/\t/;

	}#end if

	
	return $num_interactions;


}#end get_number_of_interactions


#sub get_num_phenotypes - Takes as input: 1) gene	- A gene symbol
#					  2) dbh	- A database handle
#			- Returns the number of occurences of the gene symbol in the database
sub get_num_phenotypes
{
	my ($gene, $dbh) = @_;

	my $num_phenotypes = 0;


	#GET THE PHENOTYPE INFO
        my $select_phenotypes = "SELECT COUNT(*) 
                                 FROM morbid_map
                                 WHERE gene_symbols LIKE \'\"$gene\"\' OR
					gene_symbols LIKE \'% $gene\"\' OR	
					gene_symbols LIKE \'% $gene,%\' OR
					gene_symbols LIKE \'\"$gene,%\' OR
					gene_symbols LIKE \'$gene\'";

        my $sth = $dbh -> prepare($select_phenotypes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


	$num_phenotypes = $sth -> fetchrow_array;

	return $num_phenotypes;

}#end get_num_phenotypes


1;
