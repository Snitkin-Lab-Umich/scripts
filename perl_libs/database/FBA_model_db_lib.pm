package FBA_model_db_lib;


use DBI;
use db_lib;
use list_ops;


###################################################### DB CONNECTION #######################################################

#sub DB_connect - Establishes a database connection and returns a handle
sub FBA_model_DB_connect
{
	my ($data_source) = @_;

        my $user = "root";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end FBA_model_DB_connect


################################################### BASIC DB INSERTION ROUTINES ###########################################

#sub genes_to_FBA_DB - Takes as input: 1) genes	- A list of genes
#				       2) dbh	- A database handle
#		     - Puts genes in gene table where AI gid will be assigned. Only puts gene in database if 
#		       it does not already exist.
sub genes_to_FBA_DB
{
	my($r_genes , $dbh) = @_;


	#GET ALL GENES CURRENTLY IN DATABASE
	my $r_existing_genes = &get_all_genes_from_FBA_DB($dbh);	
	my %existing_genes = %$r_existing_genes;


	#INSERT ALL GENES INTO GENE TABLE WHERE AN AUTO INC ID WILL BE ASSIGNED
	foreach my $gene (sort keys %$r_genes)
	{


		#CHECK IF GENE EXISTS IN DB, OTHERWISE INSERT IT
		if( !exists($existing_genes{$gene}) )
		{

			#PREPARE INSERT STATEMENT AND EXECUTE
	        	my $insert_gene = "INSERT INTO gene (gid , ORF_id)
	        	                   VALUES (gid , ?)";


	        	my $sth = $dbh -> prepare($insert_gene)
	        	                  or die "Can't prepare statment: $DBI::errstr";


	        	my $rc = $sth -> execute($gene)
	                                 or die "Can't execute statment: $DBI::errstr";


		}#end if

	}#end foreach


}#end genes_to_FBA_DB


#sub SGD_info_to_gene - Takes as input : 1) name - Common name of gene
#					 2) alais - An alias if one exists
#					 3) description - A short description of the gene
#					 4) phenotype - Phenotype of KO
#					 5) gene_prod - 
#					 6) ORF_id - The orf id of the gene which is already in thr table
#					 7) SGD_id - THe SGD ID of the gene
#					 8) dbh - A database handle
#
#		      - Puts the gene information into the gene table
sub SGD_info_to_gene
{
	my ($name , $alias , $description , $phenotype , $gene_prod , $ORF_id , $SGD_id , $dbh) = @_;


        #UPDATE ENTRY WITH APPROPRIATE ID
        my $update_gene = "UPDATE gene
			   SET name = ? ,
			       alias = ? ,
			       SGD_id = ? ,
			       description = ? ,
			       phenotype = ? ,
			       gene_prod = ?
                           WHERE ORF_id = ?";


        my $sth = $dbh -> prepare($update_gene)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($name , $alias , $SGD_id , $description , $phenotype , $gene_prod , $ORF_id)
                         or die "Can't execute statment: $DBI::errstr";


}#end SGD_info_to_gene


#sub proteins_to_FBA_DB - Takes as input: 1) proteins	- A list of proteins
#					  2) dbh	- A datbase handle
#			 - Puts proteins in protein table where AI pid will be assigned. Only puts gene in database if 
#                          it does not already exist.
sub proteins_to_FBA_DB
{
        my($r_proteins , $dbh) = @_;


	#GET ALL PROTEINS CURRENTLY IN DATABASE
	my $r_existing_proteins = &get_all_proteins_from_FBA_DB($dbh);	
	my %existing_proteins = %$r_existing_proteins;


        #INSERT ALL PROTEINS INTO PROTEIN TABLE WHERE AN AUTO INC ID WILL BE ASSIGNED
        foreach my $protein (sort keys %$r_proteins)
        {

		if( !exists($existing_proteins{$protein}) )
		{

                	my $insert_protein = "INSERT INTO protein (protein_name)
                	                       VALUES (?)";


                	my $sth = $dbh -> prepare($insert_protein)
        	                          or die "Can't prepare statment: $DBI::errstr";
	

                	my $rc = $sth -> execute($protein)
                                         or die "Can't execute statment: $DBI::errstr";

		}#end if

        }#end foreach


}#end proteins_to_FBA_DB



#sub genes2protein_to_FBA_DB - Takes as input: 1) protein	-  A protein 
#					       2) genes		- A lst of genes forming this protein
#					       3) model_id	- An ID from the model table
#					       4) dbh		- A database handle
#			     - Puts pid/gid pairs in gene2protein table
sub genes2protein_to_FBA_DB
{
	my ($protein , $r_genes , $model_id,  $dbh) = @_;
	my $pid;
	my %gids = ();

	
	#GET THE PID FOR THE GIVEN COMPLEX
        my $select_pid = "SELECT (p.pid)
                          FROM protein p
                          WHERE p.protein_name = ?";


        my $sth = $dbh -> prepare($select_pid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($protein)
                         or die "Can't execute statment: $DBI::errstr";

        $pid =  $sth -> fetchrow_array;



	#GET THE GID FOR ALL THE GENES
	foreach my $gene (@$r_genes)
	{

		#PREPARE SELECT STATEMENT AND EXECUTE
        	my $select_gid = "SELECT (g.gid)
        	                  FROM gene g
        	                  WHERE g.ORF_id = ?";


        	my $sth = $dbh -> prepare($select_gid)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($gene)
        	                 or die "Can't execute statment: $DBI::errstr";
	
        	my $gid =  $sth -> fetchrow_array;

		$gids{$gid} = 1;


	}#end foreach



	#FOREACH GID INSERT pID/GID PAIR
        foreach my $gid (keys %gids)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_g2p = "INSERT INTO gene2protein (gid, pid, model_id)
                                       VALUES (? , ?, ?)";


                my $sth = $dbh -> prepare($insert_g2p)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($gid , $pid, $model_id)
                                         or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end genes2protein_to_FBA_DB



#sub rxn2expr_to_FBA_DB - Takes as input: 1) rxn		- A reaction code 
#					  2) expressions	- A list of logical expressions indicating which protein 
#								  combinations(complexes) are sufficient to catalyze the reaction
#					  3) model_id		- An ID from the mode table
#					  4) dbh		- A database handle
#			 - Puts rid and expression pairs in the rxn2complex_expr table
sub rxn2expr_to_FBA_DB
{
	my ($rxn , $r_expressions , $model_id, $dbh) = @_;
	my $rid;


	#GET THE REACTION ID OF THE GIVEN REACTIOM

        my $select_rid = "SELECT (r.rid)
                          FROM rxn r
                          WHERE r.rxn_code = ? AND
				r.model_id = ?";


        my $sth = $dbh -> prepare($select_rid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($rxn, $model_id)
                         or die "Can't execute statment: $DBI::errstr";

        $rid =  $sth -> fetchrow_array;



	#FOR EACH EXPRESSION MAKE ENTRY IN TABLE WITH RID/EXPRESSION PAIR
	my %expressions;

	foreach my $expr (@$r_expressions)
	{

		if(!exists($expressions{$expr}))
		{	
	
			#PREPARE INSERT STATEMENT AND EXECUTE
                	my $insert_r2e = "INSERT INTO rxn2complex_expr (rid, complex_expr)
                	                       VALUES (? , ?)";


               		my $sth = $dbh -> prepare($insert_r2e)
                        	          or die "Can't prepare statment: $DBI::errstr";


                	my $rc = $sth -> execute($rid , $expr)
                        	                 or die "Can't execute statment: $DBI::errstr";

			$expressions{$expr} = 1;

		}#end if



	}#end foreach


}#end rxn2expr_to_FBA_DB



#sub rxn2protein_to_FBA_DB - Takes as input: 1) rxn		- A reaction code 
#					     2) proteins	- A list of proteins participating in the reaction
#					     3) model_id	- An ID from the model table
#					     4) dbh		- A database handle
#			    - Puts rid/pid pairs in the rxn2protein table (NOTE- protein does not have to be sufficient
#			      or neccasary to catalyze reaction, just can be part)
sub rxn2protein_to_FBA_DB
{
	my ($rxn , $r_proteins , $model_id, $dbh) = @_;
	my $rid;
	my %pids;


        #GET THE RID FOR THE GIVEN REACTION
        my $select_rid = "SELECT (r.rid)
                          FROM rxn r
                          WHERE r.rxn_code = ? AND
				r.model_id = ?";


        my $sth = $dbh -> prepare($select_rid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($rxn, $model_id)
                         or die "Can't execute statment: $DBI::errstr";

        $rid =  $sth -> fetchrow_array;


        #GET THE PID FOR ALL THE PROTEINS PARTICIPATING IN THE REACTION
        foreach my $protein (@$r_proteins)
        {

                my $select_pid = "SELECT (p.pid)
                                  FROM protein p
                                  WHERE p.protein_name = ?";


                my $sth = $dbh -> prepare($select_pid)
                                 or die "Can't prepare statment: $DBI::errstr";

                my $rc = $sth -> execute($protein)
                                 or die "Can't execute statment: $DBI::errstr";

                my $pid =  $sth -> fetchrow_array;

                $pids{$pid} = 1;


        }#end foreach


        #FOREACH GID INSERT PID/GID PAIR
        foreach my $pid (keys %pids)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_r2p = "INSERT INTO rxn2protein (rid,pid)
                                  VALUES (? , ?)";


                my $sth = $dbh -> prepare($insert_r2p)
                                  or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($rid , $pid)
                                         or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end rxn2protein_to_FBA_DB



#sub rxn_to_FBA_DB - Takes as input: 1) rxn_code	- A reaction code
#				     2) rxn_name	- Name of a reaction
#				     3) rxn		- The stoichiometry
#				     4) subsystem	- The process the reaction participates in
#				     5) protein_class	- EC information
#				     6) model_id	- An ID from the model table
#				     7) dbh		- A database handle	
#		    - Generates an AI rid and puts associated data in the rxn table
sub rxn_to_FBA_DB
{
	my ($rxn_code , $rxn_name , $rxn,  $subsystem , $protein_class , $model_id, $dbh) = @_;

	#INSERT REACTION AND ASSOCIATED INFORMATION INTO THE DATABASE

        #PREPARE INSERT STATEMENT AND EXECUTE
        my $insert_rxn = "INSERT INTO rxn (rxn_code , rxn_name , rxn , subsystem , prot_class, model_id)
                          VALUES (? , ? , ? , ?, ?, ?)";


        my $sth = $dbh -> prepare($insert_rxn)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($rxn_code , $rxn_name, $rxn , $subsystem, $protein_class, $model_id)
                          or die "Can't execute statment: $DBI::errstr";



}#end rxn_to_FBA_DB


#sub compartments_to_FBA_DB - Takes as inputL 1) comp_name	- The name of the compartment
#					     2) dbh		- A database handle
#			   - Enters compartment into database
sub compartments_to_FBA_DB
{
	my ($r_comp_names, $dbh) = @_;

	my %comp_names = %$r_comp_names;

	
	#GET EXISTING COMPARTMENTS
	my $r_existing_comps = &get_all_compartments_from_FBA_DB($dbh);
	my %existing_comps = %$r_existing_comps;

	
	#INSERT COMPARTMENTS INTO DATABASE IF THEY DO NOT ALREADY EXIST
	foreach my $comp_id (keys %comp_names)
	{

		my $comp_name = $comp_names{$comp_id};

 		if( !exists($existing_comps{$comp_name}) )
		{
	
			#PREPARE INSERT STATEMENT AND EXECUTE
        		my $insert_comp = "INSERT INTO compartment (comp_id , comp_name)
        	        	           VALUES (comp_id, ?)";
                

        		my $sth = $dbh -> prepare($insert_comp)
        	        	          or die "Can't prepare statment: $DBI::errstr";

        
       			my $rc = $sth -> execute($comp_name)
        	        	          or die "Can't execute statment: $DBI::errstr";

		}#end if


	}#end foreach


}#end compartment_to_FBA_DB


#sub metab_to_FBA_DB - Takes as input:	1) metab	- The code for the metabolite in a reaction
#					2) desc		- A description of the metabolite
#					3) dbh		- A database handle
#					4) comp_name	- The name of the compartment where the metabolite resides
#                     - Inserts metabolite data into database
sub metab_to_FBA_DB
{
        my ($metab , $desc , $dbh, $comp_name, $formula) = @_;


	my %metabs;

	#GET ALL METABOLITES CURRENTLY IN THE DATABASE

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_metabs = "SELECT m.metab_symbol
                             FROM metabolite m";

        
        my $sth = $dbh -> prepare($select_metabs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";
        
                         

        while(my ($metab_symbol) =  $sth -> fetchrow_array)
        {                

                $metabs{$metab_symbol} = 1;

        }#end while


	#INSERT METABOLITE INTO TBALE IF IT IS NOT ALREADY IN IT	
	if( !exists($metabs{$metab}) )
	{

		#GET THE ID FOR THE CORRESPONDING COMPARTMENT NAME
		my $r_comps = &get_all_compartments_from_FBA_DB($dbh);		
		my %comps = %$r_comps;


        	#PREPARE INSERT STATEMENT AND EXECUTE
        	my $insert_metab = "INSERT INTO metabolite (metab_id, metab_symbol , metab_desc, formula, comp_id)
              	 	             VALUES (metab_id, ? , ?, ?, ?)";


        	my $sth = $dbh -> prepare($insert_metab)
                	          or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($metab , $desc, $formula, $comps{$comp_name})
                	          or die "Can't execute statment: $DBI::errstr";

	}#end if


}#end metab_to_FBA_DB


#sub cond_to_FBA_DB - Takes as input : 1) cond_name 	- The name of the condition to be inserted
#                                       2) flux_Lbound	- A tab delimited list of flux lower bounds in the alphabetical
#							  order by reaction codes 
#                                       3) flux_Ubound	- A tab delimited list of flux upper bounds in the alphabetical
#							  order by reaction codes 
#					4) model_id	- An ID in the model table
#					5) dbh		- A database handle
#
#                    - Inserts condition name into condition table and bounds into condition bounds table
sub cond_to_FBA_DB
{
	my ($cond_name , $flux_Lbound , $flux_Ubound , $model_id, $dbh) = @_;

	my $cond_id;
	my %conds;

	#GET CONDITION NAME FROM TABLE SO CURRENT CONDITION CAN BE CHECKED FOR

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_con = "SELECT c.condition_name, c.con_id
                            FROM `condition` c";


        my $sth = $dbh -> prepare($select_con)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($cond_name, $cond_id) =  $sth -> fetchrow_array)
        {

                $conds{$cond_name} = $cond_id;

        }#end while




        #INSERT ENTRY INTO CONDITION TABLE AND RETRIEVE ID IF DOESN't ALREADY EXIST
	if( !exists($conds{$cond_name}))
	{

	        my $insert_cond = "INSERT INTO `condition` (condition_name)
        	                   VALUES (?)";


        	my $sth = $dbh -> prepare($insert_cond)
        	                  or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($cond_name)
                	          or die "Can't execute statment: $DBI::errstr";

       		$cond_id = $dbh->{ q{mysql_insertid}};

	}else{

	
		$cond_id = $conds{$cond_name};

	}

	#INSERT INTO CONDITION BOUNDS TABLE THE FLUX BOUNDARIES FOR THE CONDITION
	my $insert_condB = "INSERT INTO condition_bounds (con_id , flux_LB , flux_UB, model_id)
                            VALUES (? , ? , ?, ?)";


        $sth = $dbh -> prepare($insert_condB)
                          or die "Can't prepare statment: $DBI::errstr";


        $rc = $sth -> execute($cond_id , $flux_Lbound , $flux_Ubound, $model_id)
                          or die "Can't execute statment: $DBI::errstr";


}#end cond_to_FBA_DB


################################################### BASIC SELECTION ROUTINES ###########################################

#sub get_all_genes_from_FBA_DB - Takes as input: 1) dbh             - A database handle
#                          - Queries the FBA model DB and returns all the ORF_ids in the gene table
#
sub get_all_genes_from_FBA_DB
{
        my ($dbh) = @_;

        my %genes;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_genes = "SELECT gid, g.ORF_id
                            FROM gene g";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid,$ORF_id) =  $sth -> fetchrow_array)
        {

                $genes{$ORF_id} = $gid;

        }#end while


        return \%genes;

}#end get_all_genes_from_FBA_DB


#sub get_genes_from_FBA_DB - Takes as input: 1) model_id	- An ID from the model table
#					     2) dbh		- A database handle
#			   - Queries the FBA model DB and returns all the ORF_ids in the gene table
#
sub get_genes_from_FBA_DB
{
	my ($model_id, $dbh) = @_;
	my @genes;
	my @sorted_genes;


	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_genes = "SELECT DISTINCT g.ORF_id
                            FROM gene g, rxn r, rxn2protein rp,  gene2protein gp
                            WHERE r.model_id = ? AND
                                  rp.rid = r.rid AND
                                  rp.pid = gp.pid AND
                                  gp.model_id = r.model_id AND
                                  g.gid = gp.gid";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my $ORF_id =  $sth -> fetchrow_array)
	{

		push @genes , $ORF_id;

	}#end while


	#SORT GENES BEFORE RETURNING
	@sorted_genes = sort @genes;
	

	return \@sorted_genes;


}#end get_genes_from_FBA_DB


#sub get_genes_and_ids_from_FBA_DB - Takes as input: 1) model_id	- An ID from the model table
#						     2) dbh		- A database handle
#				   - Queries the FBA model DB and returns the ORF_ids and ids in the gene table
#
sub get_genes_and_ids_from_FBA_DB
{
        my ($model_id, $dbh) = @_;
        my %genes;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_genes = "SELECT g.gid , g.ORF_id
                            FROM gene g, rxn r, rxn2protein rp,  gene2protein gp
                            WHERE r.model_id = ? AND
                                  rp.rid = r.rid AND
                                  rp.pid = gp.pid AND
                                  gp.model_id = r.model_id AND
                                  g.gid = gp.gid";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid , $ORF_id) =  $sth -> fetchrow_array)
        {

                $genes{$ORF_id} = $gid;

        }#end while


        return \%genes;


}#end get_genes_and_ids_from_FBA_DB


#sub get_ids_and_genes_from_FBA_DB - Takes as input: 1) model_id	- An ID from the model database
#						     2) dbh		- A database handle
#				   - Queries the FBA model DB and returns the ORF_ids and ids in the gene table
#
sub get_ids_and_genes_from_FBA_DB
{
        my ($model_id, $dbh) = @_;
        my %genes;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_genes = "SELECT g.gid , g.ORF_id
                            FROM gene g, rxn r, rxn2protein rp,  gene2protein gp
			    WHERE r.model_id = ? AND
				  rp.rid = r.rid AND
				  rp.pid = gp.pid AND
				  gp.model_id = r.model_id AND
				  g.gid = gp.gid";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid , $ORF_id) =  $sth -> fetchrow_array)
        {

                $genes{$gid} = $ORF_id;

        }#end while


        return \%genes;


}#end get_ids_and_genes_from_FBA_DB


#sub ORF_to_geneName - Takes as input : 1) ORF_id 	- An ORF ID present in the database
#					2) dbh		- A datbase handle
#		     - Returns the value in gene_name field if one exists, if not returns the ORF ID
sub ORF_to_geneName
{
	my ($ORF_id , $dbh) = @_;

	my $gene_name;


	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_gene = "SELECT g.name
                           FROM gene g
			   WHERE g.ORF_id = ?";


        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ORF_id)
                         or die "Can't execute statment: $DBI::errstr";



        $gene_name  =  $sth -> fetchrow_array;


	if($gene_name eq "")
	{

		return $ORF_id;

	}else{

		return $gene_name;

	}#end if


}#end ORF_to_geneName


#sub geneName_to_ORF - Takes as input : 1) geneName 	- An gene name  present in the database
#					2) dbh		- A database handle
#                    - Returns the value in ORF_ID field if one exists, if not returns an empty string
sub geneName_to_ORF
{
        my ($geneName , $dbh) = @_;

        my $ORF_id = "";
	my $geneName_like =  $geneName . "%";


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_gene = "SELECT g.ORF_id
                           FROM gene g
                           WHERE (g.name = ? OR
				  g.alias LIKE ? OR
				  g.ORF_id = ?)";


        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($geneName , $geneName_like , $geneName)
                         or die "Can't execute statment: $DBI::errstr";



        $ORF_id  =  $sth -> fetchrow_array;


        return $ORF_id;


}#end geneName_to_ORF


#sub get_genes_and_proteins_from_FBA_DB - Takes as input:
#
#					- Returns hash relating proteins to genes
sub get_genes_and_proteins_from_FBA_DB
{
	my ($model_id , $dbh) = @_;

	my %prots2genes;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_prots = "SELECT p.protein_name, g.ORF_id
                            FROM protein p, gene g, gene2protein gp
			    WHERE gp.model_id = ? AND
				  g.gid = gp.gid AND
				  p.pid = gp.pid";

        
        my $sth = $dbh -> prepare($select_prots)
                         or die "Can't prepare statment: $DBI::errstr";
        
        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";


        
        while(my ($prot_name, $gene_name) =  $sth -> fetchrow_array)
        {
            
                $prots2genes{$prot_name} = $gene_name;
        
        }#end while



	return \%prots2genes;

}#end get_genes_and_proteins_from_FBA_DB  


#sub get_all_proteins_from_FBA_DB - Takes as input: 1) dbh 	- A database handle
#                          	  - Queries the FBA model DB and returns all the protein_names in the gene table
#
sub get_all_proteins_from_FBA_DB
{
        my ($dbh) = @_;

        my %proteins;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_prots = "SELECT pid, protein_name
                            FROM protein p";

        
        my $sth = $dbh -> prepare($select_prots)
                         or die "Can't prepare statment: $DBI::errstr";
        
        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        
        while(my ($pid, $prot_name) =  $sth -> fetchrow_array)
        {
                
                $proteins{$prot_name} = $pid;
        
        }#end while

        
        return \%proteins; 

}#end get_all_genes_from_FBA_DB


#sub get_all_compartments_from_FBA_DB - Takes as input: 1) dbh             - A database handle
#                         	      - Queries the FBA model DB and returns all the comp_names in the compartment table
#
sub get_all_compartments_from_FBA_DB
{
        my ($dbh) = @_;

        my %comps;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_comps = "SELECT comp_id, comp_name
                            FROM compartment";


        my $sth = $dbh -> prepare($select_comps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($comp_id,$comp_name) =  $sth -> fetchrow_array)
        {

                $comps{$comp_name} = $comp_id;

        }#end while


        return \%comps;

}#end get_all_compartments_from_FBA_DB



#sub get_condition - Takes as input the name of a condition and returns the id from the condition table
#
sub get_condition
{
	my ($condition , $dbh) = @_;
	my $cond_id;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_condid = "SELECT (c.con_id)
                            FROM `condition` c
			    WHERE c.condition_name = ?";


        my $sth = $dbh -> prepare($select_condid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition)
                         or die "Can't execute statment: $DBI::errstr";



        $cond_id =  $sth -> fetchrow_array;



        return $cond_id;


}#end get_condition


#sub get_conditions - Takes as input a dataase handle and returns the names of all conditions in the database
#
sub get_conditions
{
        my ($dbh) = @_;
	my @conditions;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_con = "SELECT c.condition_name
                            FROM `condition` c";


        my $sth = $dbh -> prepare($select_con)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



	 while(my $cond_name =  $sth -> fetchrow_array)
        {

                push @conditions , $cond_name;

        }#end while


        return \@conditions;


}#end get_conditions


#sub get_conditions_and_ids - Takes as input a dataase handle and returns the names of all conditions in the database
#
sub get_conditions
{
        my ($dbh) = @_;

        my %conditions;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_con = "SELECT c.con_id, c.condition_name
                            FROM `condition` c";


        my $sth = $dbh -> prepare($select_con)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



         while(my ($cond_id, $cond_name) =  $sth -> fetchrow_array)
        {

                $conditions{$cond_name} = $cond_id;

        }#end while


        return \%conditions;


}#end get_conditions


#sub get_condition_bounds - Takes as input: 1) cond	- The name of a condition in the condition table
#					    2) model_id	- The id of the model for which the bounds apply
#					    3) dbh	- A database handle
#			  - Returns a tab delimited list for the lower bounds and for the upper bounds
sub get_condition_bounds
{
	my ($cond , $model_id,  $dbh) = @_;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_bounds = "SELECT cb.flux_LB , cb.flux_UB
                             FROM condition_bounds cb , `condition` c
			     WHERE c.condition_name = ? AND
				   c.con_id = cb.con_id AND
				   cb.model_id = ?";


        my $sth = $dbh -> prepare($select_bounds)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($cond, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        my ($flux_LB , $flux_UB) =  $sth -> fetchrow_array;


        return ($flux_LB , $flux_UB);


}#end get_condition_bounds


#sub get_rxns_from_FBA_DB - Queries the FBA model DB and returns are the rxn_codes in the rxn table
#
sub get_rxns_from_FBA_DB
{
        my ($model_id, $dbh) = @_;
        my (@rxns , @sorted_rxns);


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_rxns = "SELECT (r.rxn_code)
                           FROM rxn r
			   WHERE r.model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my $rxn_code =  $sth -> fetchrow_array)
        {

                push @rxns , $rxn_code;

        }#end while


	@sorted_rxns = sort(@rxns);

        return \@sorted_rxns;


}#end get_rxns_from_FBA_DB


#sub get_rxns_and_ids_from_FBA_DB - Queries the FBA model DB and returns the rxn_codes and ids in the gene table
#
sub get_rxns_and_ids_from_FBA_DB
{
        my ($model_id, $dbh) = @_;
        my %rxns;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_rxns = "SELECT r.rid , r.rxn_code
                           FROM rxn r
			   WHERE model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rid , $rxn_code) =  $sth -> fetchrow_array)
        {

                $rxns{$rxn_code} = $rid;

        }#end while


        return \%rxns;


}#end get_rxns_and_ids_from_FBA_DB


#sub get_rxn_equations_from_FBA_DB - Takes as input : 1) model_id	- An ID of a model in the database
#						       2) dbh		- A database handle
#				    - Returns a hash with the key being the reaction code and the value being
#				      the reaction itself (e.g. A + B -> C + D)
sub get_rxn_equations_from_FBA_DB
{
	my ($model_id, $dbh) = @_;

        my %rxns;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_rxns = "SELECT r.rxn_code , r.rxn
                           FROM rxn r
			   WHERE model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $rxn) =  $sth -> fetchrow_array)
        {

                $rxns{$rxn_code} = $rxn;

        }#end while


        return \%rxns;


}#end get_rxn_equations_from_FBA_DB



#sub get_biomass_reaction_from_FBA_DB - Takes as input : 1) model_id       - An ID of a model in the database
#                                                      2) dbh           - A database handle
#				      - Returns the name of the biomass reaction in the given model
sub get_biomass_reaction_from_FBA_DB
{
	my ($model_id, $dbh) = @_;

	my $biomass_rxn;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_rxns = "SELECT r.rxn_code
                           FROM rxn r
                           WHERE model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";
                         
        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";
                         


        while(my ($rxn_code) =  $sth -> fetchrow_array)
        {
        
		if( $rxn_code =~ /biomass/i)
		{

			$biomass_rxn = $rxn_code;

		}#end if

        }#end while



	return $biomass_rxn

}#end 

#sub get_gene_subsystem - Takes as input : 1) gid	- gid from gene table
#					   2) model_id	- ID from the model table
#					   3) dbh	- database handle
sub get_gene_subsystem
{
	my ($gid , $model_id, $dbh) = @_;

	my @sub_sys;


	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_sub_sys = "SELECT DISTINCT r.subsystem
                              FROM rxn r , gene g , gene2protein g2p ,rxn2protein r2p
			      WHERE g.gid = ? AND 
				    g.gid = g2p.gid AND
				    g2p.pid = r2p.pid AND
			            r2p.rid = r.rid AND
				    r.model_id = ?";


        my $sth = $dbh -> prepare($select_sub_sys)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($gid, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($sub_sys) =  $sth -> fetchrow_array)
        {

		if( $sub_sys =~ /[a-zA-Z]+/ )
		{

			push @sub_sys , $sub_sys;

		}#end if

        }#end while


        return \@sub_sys;


}#end get_gene_subsystem



#sub get_processes - Takes as input :1) model_id 	- A ID from the model table
#				     2) dbh 		- Database handle
#		   - Returns a list of processes present in the reaction table
sub get_processes
{
	my ($model_id , $dbh) = @_;

	my @processes;

	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_proc = "SELECT DISTINCT r.subsystem
                           FROM rxn r
			   WHERE r.model_id = ?";


        my $sth = $dbh -> prepare($select_proc)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($proc) =  $sth -> fetchrow_array)
        {

		if( $proc =~ /[a-zA-Z]+/ )
		{

                	push @processes , $proc;

		}#end if

        }#end while


        return \@processes;


}#end get_processes


#sub get_rxns_in_processes - Takes as input: 1) r_processes 	- A reference to a list of processes
#					     2) model_id	- An ID from the model table
#					     3) dbh 		- A database handle
#			   - Returns a list of rxns which are annotated as participating in the given process
sub get_rxns_in_processes
{
	my ($r_processes , $model_id, $dbh) = @_;
	my @rxns;

        #GET ALL REACTIONS WHICH PARTICIPATE IN ONE OF THE PROCESSES OF INTEREST
        my $select_rxns = "SELECT DISTINCT r.rxn_code
                           FROM rxn r
                           WHERE r.model_id = ? AND 
				 r.subsystem IN (";

        #ADD PLACE HOLDERS FOR GENES
        $select_rxns = &db_lib::add_place_holders($select_rxns , scalar(@$r_processes));

        $select_rxns = $select_rxns . ")";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id, @$r_processes)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code) =  $sth -> fetchrow_array)
        {

		push @rxns , $rxn_code;

	}#end while


	return \@rxns;


}#end get_rxns_in_processes


#sub get_rxns_and_processes - Takes as input:  1) model_id	  - An ID from the model table
#					       2) dbh             - A database handle
#
#                          - Returns a list of rxns which are annotated as participating in the given process
sub get_rxns_and_processes
{
        my ($model_id, $dbh) = @_;
        my %rxn2proc;

        #GET ALL REACTIONS WHICH PARTICIPATE IN ONE OF THE PROCESSES OF INTEREST
        my $select_rxns = "SELECT DISTINCT r.rxn_code , r.subsystem
                           FROM rxn r
			   WHERE r.model_id = ?";

        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $proc) =  $sth -> fetchrow_array)
        {

                $rxn2proc{$rxn_code} = $proc;

        }#end while


        return \%rxn2proc;


}#end get_rxns_and_processes


#sub get_genes_and_proccess - Takes as input: 1) model_id	- An ID from the model table
#					      2) dbh 		- A database handle
#			    - Returns a reference to a list containing the name of processes in the reaction table 
#			      in the database and a reference to a hash which links gene names to indicies in the 
#			      process list which correspond to the processes the reactions which they catalyze belong
#			      to
sub get_genes_and_processes
{
	my ($model_id, $dbh) = @_; 

	my %proc2ind;
	my $process_ind = 1;
	my @processes;

	my %gene2proc;


	#GET A LIST OF PROCESSES IN THE DATABASE
        my $select_procs = "SELECT DISTINCT r.subsystem
                            FROM rxn r
			    WHERE r.subsystem != 'NULL' AND
				  r.model_id = ?
                            ORDER BY r.subsystem ASC";


        my $sth = $dbh -> prepare($select_procs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($proc) =  $sth -> fetchrow_array)
        {

		push @processes , $proc;		

                $proc2ind{$proc} = $process_ind;

		$process_ind ++;

        }#end while
	

	#GET ALL PAIRS OF GENE/PROCESS ASSOCIATIONS
        my $select_gp_ints = "SELECT DISTINCT g.ORF_id , r.subsystem 
                              FROM rxn r, gene g , gene2protein gp , rxn2protein rp
			      WHERE r.rid = rp.rid AND
			            g.gid = gp.gid AND
				    gp.pid = rp.pid AND
				    r.model_id = ? AND
				    gp.model_id = r.model_id"; 


        my $sth = $dbh -> prepare($select_gp_ints)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene , $proc) =  $sth -> fetchrow_array)
        {

		push @{$gene2proc{$gene}} , $proc2ind{$proc};

        }#end while


	return (\%gene2proc , \@processes);


}#end get_genes_and_processes


#sub get_rxns_and_ECnums - Takes as input: 1) model_id	- An ID from the model table
#					   2) dbh	- A database handle
#			 - Returns a reference to a hash where the keys are rxn_codes and the values are 
#			   the EC numbers of the corresponding reactions
sub get_rxns_and_ECnums
{
	my ($model_id , $dbh) = @_;
	my %rxn2EC;


	#GET ALL REACTION CODES AND EC NUMBERS IF THE EC NUMBER EXISTS
        my $select_rxns = "SELECT DISTINCT r.rxn_code , r.prot_class
                           FROM rxn r
                           WHERE r.prot_class LIKE '%.%' AND
				 r.model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $EC) =  $sth -> fetchrow_array)
        {

		$EC =~ s/^([^;]+);.*$/$1/;

        	$rxn2EC{$rxn_code} = $EC;

        }#end while


	return \%rxn2EC;


}#end get_rxns_and_ECnums


#sub get_rxns_procs_and_EC - Takes as input: 1) model_id	- An ID from the model table
#					     2) dbh 		- A database handle
#                          - Returns a reference to a hash where the keys are rxn_codes and the values are
#                            the EC numbers of the corresponding reactions
sub get_rxns_procs_and_EC
{
        my ($model_id, $dbh) = @_;
        my %rxn2ECproc;


        #GET ALL REACTION CODES AND EC NUMBERS IF THE EC NUMBER EXISTS
        my $select_rxns = "SELECT DISTINCT r.rxn , r.prot_class, r.subsystem
                           FROM rxn r
			   WHERE r.model = ?";

        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn , $EC , $proc) =  $sth -> fetchrow_array)
        {


                $rxn2ECproc{$rxn} -> {'EC'} = $EC;
                $rxn2ECproc{$rxn} -> {'PROC'} = $proc;

        }#end while


        return \%rxn2ECproc;


}#end get_rxns_procs_and_EC


#sub get_rxns_and_ORFS	 - Takes as input: 1) process	- The name of the process for which the ORFs corresponding to
#							  the reactions should be retrieved
#					   2) model_id	- An ID from the model table
#					   3) dbh       - A database handle
#                        - Returns a reference to a hash where the keys are rxn_codes and the values are
#                          the ORFS catalyzing the corresponding reactions
sub get_rxns_and_ORFs
{
        my ($process , $model_id, $dbh) = @_;
        my %rxn2ORF;


        #GET ALL REACTION CODES AND EC NUMBERS IF THE EC NUMBER EXISTS
        my $select_rxns = "SELECT DISTINCT r.rxn_code , g.ORF_id
                           FROM gene g, protein p , gene2protein gp, rxn r , rxn2protein rp
                           WHERE g.gid = gp.gid AND
				 gp.pid = rp.pid AND
				 rp.rid = r.rid AND
				 r.subsystem = ? AND
				 r.model_id = ? AND
				 gp.model_id = r.model_id";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($process, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $ORF) =  $sth -> fetchrow_array)
        {


                $rxn2ORF{$rxn_code} = $ORF;

        }#end while


        return \%rxn2ORF;


}#end get_rxns_and_ORFs


#sub get_rxns_and_ORFS_GL - Takes as input: 1) gene_list	- The name of the ORFS for which catalyzed reactions should be 
#							   	  found
#					    2) model_id		- An ID from the model table
#					    3) dbh       	- A database handle
#                        - Returns a reference to a hash where the keys are rxn_codes and the values are
#                          the ORFS catalyzing the corresponding reactions
sub get_rxns_and_ORFs_GL
{
        my ($r_gene_list , $model_id, $dbh) = @_;
        my %ORF2rxn;


        #GET ALL REACTION CODES AND EC NUMBERS IF THE EC NUMBER EXISTS
        my $select_rxns = "SELECT DISTINCT r.rxn_code , g.ORF_id
                           FROM gene g, protein p , gene2protein gp, rxn r , rxn2protein rp
                           WHERE g.gid = gp.gid AND
				 gp.pid = rp.pid AND
				 rp.rid = r.rid AND
				 r.model_id = ? AND
				 gp.model_id = r.model_id AND
				 g.ORF_id IN (";

	#ADD PLACE HOLDERS FOR GENES
        $select_rxns = &db_lib::add_place_holders($select_rxns , scalar(@$r_gene_list));

        $select_rxns = $select_rxns . ")";

        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id, @$r_gene_list)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $ORF) =  $sth -> fetchrow_array)
        {


                push @{$ORF2rxn{$ORF}} , $rxn_code;

        }#end while


        return \%ORF2rxn;


}#end get_rxns_and_ORFs


#sub get_rxns_and_ORFS_ALL - Takes as input: 1) process - The name of the process for which the ORFs corresponding to
#                                                         the reactions should be retrieved
#                                          2) dbh       - A database handle
#                        - Returns a reference to a hash where the keys are rxn_codes and the values are
#                          the ORFS catalyzing the corresponding reactions
sub get_rxns_and_ORFs_ALL
{
        my ($model_id, $dbh) = @_;

        my %rxn2ORF;
        my %rxn2gene;


        #GET ALL REACTION CODES AND EC NUMBERS IF THE EC NUMBER EXISTS
        my $select_rxns = "SELECT DISTINCT r.rxn_code , g.SGD_id , g.name
                           FROM gene g, gene2protein gp, rxn r , rxn2protein rp
                           WHERE g.gid = gp.gid AND
                                 gp.pid = rp.pid AND
                                 rp.rid = r.rid";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rxn_code , $ORF , $name) =  $sth -> fetchrow_array)
        {


                $rxn2ORF{$rxn_code} = $ORF;

                $rxn2gene{$rxn_code} = $name;

        }#end while


        return (\%rxn2ORF , \%rxn2gene);


}#end get_rxns_and_ORFs_ALL


#sub get_model_id_from_FBA_DB - Takes as input: 1) model	- The name of a model in the model table
#						2) dbh		- A database handle
#			      - Returns the id corresponding to the model
sub get_model_id_from_FBA_DB
{
	my ($model, $dbh) = @_;

	my $model_id = 0;


        #GET MODEL ID FROM THE DATABASE
        my $select_model = "SELECT m.model_id
                           FROM model m
                           WHERE m.model_name = ?";

        my $sth = $dbh -> prepare($select_model)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model)
                         or die "Can't execute statment: $DBI::errstr";


        $model_id =  $sth -> fetchrow_array;


	return $model_id;

}#end get_model_id_from_FBA_DB

################################################### COMPLEX DB INSERTION ROUTINES ###########################################


#sub insert_singleKO - Takes as input : 1) gid 			-gid from gene table
#					2) r_flux_vector	- a reference to a flux vector
#					3) r_rxns 		- a refernce to a reaction vector in the same order as the fluc vector
#					4) con_id 		- con_id from the condition table
#					5) opt_type 		- the type of optimization performed (eg. FBA or MOMA)
#					6) status 		- the exit status for the glpk FBA
#					7) model_id		- An ID from the model table
#					8) dbh 			- database handle
#		     - Inserts data in appropriate format into singleKO_rxn and singleKO_flux tables
sub insert_singleKO
{
	my ($gid , $r_flux_vector , $con_id , $opt_type , $status , $model_id,  $dbh) = @_;

        #REMOVE LEADING NEGATIVE SIGNS FROM 0 AND CONDENSE SCIENTIFIC NOTATION (-0 -> 0 , 1E-001 -> 1E-1)
        my $flux_vector = join "\t" , @$r_flux_vector;

        $flux_vector =~ s/-0\t/0\t/g;
        $flux_vector =~ s/e-0+([0-9])/e-$1/g;


	#MAKE INSERTION INTO THE SINGLE_KO_RXN TABLE
        my $insert_fluxv = "INSERT INTO single_KO_flux (gid , con_id , flux_vector, opt_type, status, model_id)
                            VALUES (? , ? , ? , ? , ?, ?)";


        my $sth = $dbh -> prepare($insert_fluxv)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($gid , $con_id , $flux_vector, $opt_type, $status, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



}#end insert_singleKO



#sub insert_doubleKO - Takes as input : 1) gid1 		- A gid from gene table
#					2) gid2 		- A second gid from gene table
#                                       3) r_flux_vector 	- a reference to a flux vector
#                                       4) ccon_id 		- on_id from the condition table
#					5) status 		- the exit status for the glpk FBA
#					6) opt_type 		- the optimization type (eg. FBA, MOMA)
#					7) model_id		- An ID from the model table
#                                       8) dbh 			- database handle
#                    - Inserts data in appropriate format into the doubleKO_flux table
sub insert_doubleKO
{
        my ($gid1 , $gid2 , $r_flux_vector , $con_id , $status, $opt_type , $model_id, $dbh) = @_;

	#REMOVE LEADING NEGATIVE SIGNS FROM 0 AND CONDENSE SCIENTIFIC NOTATION (-0 -> 0 , 1E-001 -> 1E-1)
        my $flux_vector = join "\t" , @$r_flux_vector;

	$flux_vector =~ s/-0\t/0\t/g;
	$flux_vector =~ s/e-0+([0-9])/e-$1/g;



        #MAKE INSERTION INTO THE SINGLE_KO_RXN TABLE
        my $insert_fluxv = "INSERT INTO double_KO_flux (gid1 , gid2 , con_id , flux_vector, status , opt_type, model_id)
                            VALUES (? , ? , ? , ? , ? , ?, ?)";


        my $sth = $dbh -> prepare($insert_fluxv)
                          or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($gid1 , $gid2 , $con_id , $flux_vector, $status , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";


}#end insert_doubleKO



#sub insert_wildtype - Takes as input: 1) r_flux_vector	- a reference to a flux vector
#                                      2) r_rxns 	- a refernce to a reaction vector in the same order as the fluc vector
#                                      3) con_id 	- con_id from the condition table
#				       4) status 	- the exit status for the glpk FBA
#				       5) model_id	- An ID from the model table
#                                      6) dbh 		- database handle
#                    - Inserts data in appropriate format into wildtype table
#
sub insert_wildtype
{
        my ($r_flux_vector , $con_id ,$status, $opt_type, $model_id, $dbh) = @_;

	my @rxns = @$r_rxns;
    
 
	#REMOVE LEADING NEGATIVE SIGNS FROM 0 AND CONDENSE SCIENTIFIC NOTATION (-0 -> 0 , 1E-001 -> 1E-1)
        my $flux_vector = join "\t" , @$r_flux_vector;

        $flux_vector =~ s/-0\t/0\t/g;
        $flux_vector =~ s/e-0+([0-9])/e-$1/g;


        my @fluxes = split /\t/ , $flux_vector;


	#MAKE INSERTION INTO WILDTYPE_FLUX TABLE
        my $insert_flux_vector = "INSERT INTO wildtype_flux ( con_id , flux_vector , status , opt_type, model_id)
                                  VALUES ( ? , ? , ? , ?, ?)";


        my $sth = $dbh -> prepare($insert_flux_vector)
                                  or die "Can't prepare statment: $DBI::errstr";


        my $rc = $sth -> execute($con_id ,  $flux_vector , $status, $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



}#end insert_wildtype


################################################### COMPLEX SELECTION ROUTINES ###########################################


#sub get_dependent_rxns - Takes as input: 1) genes	- A list of gene names
#					  2) model_id	- An ID from the model table
#					  3) dbh	- A database handle
#			- Returns a list of all reactions which this these genes are neccasary for
sub get_dependent_rxns
{
	my ($r_genes , $model_id, $dbh) = @_;
	my @pids;
	my %protein_names;
	my %rxns;
	my @dep_rxns;
	my $protein_expr;
	


	#GET ALL PROTIENS WHICH CONTAIN ONE OF THE GENES OF INTEREST
	my $select_proteins = "SELECT DISTINCT gp.pid , p.protein_name
                               FROM gene2protein gp , gene g , protein p
			       WHERE gp.gid = g.gid AND
				     g.ORF_id IN ("; 

	#ADD PLACE HOLDERS FOR GENES
 	$select_proteins = &db_lib::add_place_holders($select_proteins , scalar(@$r_genes));

	$select_proteins = $select_proteins . ") AND
				     gp.pid = p.pid AND
				     gp.model_id = ?";


        my $sth = $dbh -> prepare($select_proteins)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_genes, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($pid , $protein_name) =  $sth -> fetchrow_array)
        {

                push @pids , $pid;

		$protein_names{$protein_name} = 1;

        }#end while



	#GET ALL REACTIONS IN WHICH EACH PROTEIN PARTICIPATES IN 
        my $select_rxns = "SELECT DISTINCT rp.rid , r.rxn_code 
                           FROM rxn2protein rp ,  rxn r
                           WHERE r.rid = rp.rid AND 
				 rp.pid IN (";

 	#ADD PLACE HOLDERS FOR PIDS
        $select_rxns = &db_lib::add_place_holders($select_rxns , scalar(@pids));

        $select_rxns = $select_rxns . " ) AND
			        r.model_id = ?";


        my $sth = $dbh -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@pids, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rid , $rcode) =  $sth -> fetchrow_array)
        {


		$rxns{$rid} = $rcode;


        }#end while



	#DETERMINE FOR EACH REACTION IF THERE IS A COMPLEX WHICH CAN CATALYZE IT WHICH DOES NOT
	#REQUIRE THE GENE OF INTEREST
	RXN : foreach my $rid ( keys %rxns )
	{

		#SET PARAMETER TO KEEP TRACK OF WETHER GENE IS REQUIRED FOR REACTION
		my $is_dependent = 1;


	        #GET ALL REACTIONS IN WHICH EACH PROTEIN PARTICIPATES IN
	        my $select_expr = "SELECT rc.complex_expr
	                           FROM rxn2complex_expr rc
	      	                   WHERE rc.rid = ?";
	
	        my $sth = $dbh -> prepare($select_expr)
	                         or die "Can't prepare statment: $DBI::errstr";
	
	        my $rc = $sth -> execute($rid)
	                         or die "Can't execute statment: $DBI::errstr";



		#GO THROUGH EACH EXPRESSION, IF ONE DOES NOT REQUIRE ANY OF THE PROTEINS WHICH THE 
		#GENE IS A PART OF THEN GO TO NEXT REACTION
	        EXPR : while(my ($expr) =  $sth -> fetchrow_array)
	        {
	
			my @expr = split / and / , $expr;

		
			#CHECK IF CURRENT EXPRESSION CONTAINS A PROTEIN REQUIRING THE GENE OF INTEREST
			foreach my $rxn_prot (@expr)
			{

				if( exists($protein_names{$rxn_prot}) )
				{
										
					next EXPR;

				}#end if
	

			}#end foreach
						
			
			#IF GET TO HERE IT MEANS THAT THE CURRENT EXPRESSION DOES NOT REQUIRE THE GENE
			#OF INTEREST AND THEREFORE THE REACTION IS NOT KO'D
			$is_dependent = 0;
			last EXPR;
			

	        }#end while

		
		#CHECK IF FLAG HAS BEEN SET TO ZERO INDICATING THAT THE GENE IS NOT REQUIRED FOR REACTION	
		if( $is_dependent == 1 )
		{

			push @dep_rxns , $rxns{$rid};

		}#end if


	}#end foreach



	return \@dep_rxns;


}#end get_dependent_rxns



#sub get_essential_genes - Takes as input: 1) WT_cond_id 	- The condition wildtype fluxes should be retirieved for
#					   2) KO_cond_id 	- The condition KO fluxes should be retrieved for
#					   3) biomass_rxn 	- The name of the biomass reaction
#					   4) growth_cutoff 	- The percent reduction to be considered essential
#					   5) opt_type 		- The optimization method used (e.g. FBA , MOMA)
#					   6) r_rxn_list 	- A list of reactions in the same order as the fluc vector
#					   7) model_id		- An ID from the model table
#					   8) dbh 		- A database handle
#			 - Returns a hash containing as keys all essential genes as defined by having
#			   a decrease in biomass flux of at least the percentage of 'growth_cutoff'
sub get_essential_genes
{
	my ($WT_cond_id , $KO_cond_id , $biomass_rxn , $growth_cutoff , $WT_opt_type ,  $KO_opt_type, $r_rxn_list, $model_id,  $dbh) = @_;

	my @biomass_rxn = ($biomass_rxn);

	my %dbid2ORFid;
	my $r_dbid2ORFid;

	my $r_KO_biomass;
	my %KO_biomass;
	my %essential_genes;

	my $r_WT_biomass;
	my @WT_biomass;
	my $WT_biomass;
	my $biomass_cutoff;

	
	#REALTE THE ORF IDS OF THE GENES TO THE DATABASE IDS
	$r_dbid2ORFid = &get_ids_and_genes_from_FBA_DB($model_id,$dbh);
	%dbid2ORFid = %$r_dbid2ORFid;

	#GET THE WILDTYPE VALUE FOR THE BIOMASS FLUX UNDER THE GIVEN CONDITION
	$r_WT_biomass = &get_wildtype_fluxes($WT_cond_id , \@biomass_rxn , $r_rxn_list ,  $dbh , $WT_opt_type , $model_id);
	@WT_biomass = @$r_WT_biomass;
	$WT_biomass = $WT_biomass[0];

	print "Biomass rxn = $biomass_rxn   WT_cond_id = $WT_cond_id  WT biomass = $WT_biomass Growth cutoff = $growth_cutoff\n";



	#GET ALL GENES WHOSE KO UNDER THE GIVEN CONDITION RESULTS IN A DECREASE IN THE
	#BIOMASS FLUX OF AT LEAST 'GROWTH_CUTOFF' PERCENT

	$r_KO_biomass = &get_single_KO_fluxes($KO_cond_id , \@biomass_rxn , $KO_opt_type , $r_rxn_list , $model_id, $dbh);
	%KO_biomass = %$r_KO_biomass;


	foreach my $gene(keys %KO_biomass)
	{

		my @KO_biomass = @{$KO_biomass{$gene}};
		my $KO_biomass = $KO_biomass[0];

		my $norm_biomass = $KO_biomass / $WT_biomass;
		print "gene = $dbid2ORFid{$gene}  norm biomass = $norm_biomass\n";


		if( ($KO_biomass / $WT_biomass) < $growth_cutoff)
		{

			$essential_genes{$dbid2ORFid{$gene}} = 1;


		}#end if


	}#end foreach

	

	return \%essential_genes;


}#end get_essential_genes



#sub get_influential_genes - Takes as input: 1) cond 		- The condition essentiality is being measured under
#					     2) model_id	- An ID from the model table
#					     3) dbh		- A database handle
#                          - Returns a hash containing as keys all influential genes as defined by having
#                            a flux distribution different from wildtype when knocked out
sub get_influential_genes
{
        my ($cond, $model_id, $dbh) = @_;
        my %influential_genes;


        my $select_infl_genes = "SELECT DISTINCT g.ORF_id
                                 FROM gene g, single_KO_flux ko, `condition` c, wildtype_flux wt
                                 WHERE g.gid = ko.gid AND
				       wt.con_id = ko.con_id AND
				       c.condition_name = ? AND
				       c.con_id = ko.con_id AND
				       ko.flux_vector != wt.flux_vector AND
				       wt.model_id = ko.model_id AND
				       ko.model_id = ?";


        my $sth = $dbh -> prepare($select_infl_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($cond, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my $gene =  $sth -> fetchrow_array)
        {

                $influential_genes{$gene} = 1;

        }#end while



        return \%influential_genes;


}#end get_influential_genes


#sub get_phenotypic_altering_genes - Takes as input: 1) cond		- The condition essentiality is being measured under
#						     2) model_id	- An ID from the model table
#						     3) dbh		- A database handle
#                          	   - Returns a hash containing as keys all phenotypic altering genes and pairs as 
#				     defined by having a double KO flux distribution different than both of the 
#				     single KOs
sub get_phenotypic_altering_genes
{
        my ($cond, $model_id,  $dbh) = @_;
        my %phenotypic_altering_genes;


	#SELECT GENE ID 2'S WHICH AFFECT AT LEAST ONE SINGLE KO FLUX DISTRIBUTION
        my $select_pa_genes = "(SELECT DISTINCT g.ORF_id
				FROM single_KO_flux ko1 , double_KO_flux ko2, `condition` c, gene g
				WHERE ko1.gid = ko2.gid1 AND
       				      ko2.flux_vector != ko1.flux_vector AND
       				      ko1.con_id = ko2.con_id AND
       				      ko2.con_id = c.con_id AND
				      g.gid = ko2.gid2 AND
				      c.condition_name = ? AND
				      ko1.model_id = ko2.model_id AND
			              ko2.model_id = ?)"; 


        my $sth = $dbh -> prepare($select_pa_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($cond, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my $gene =  $sth -> fetchrow_array)
        {

                $phenotypic_altering_genes{$gene} = 1;

        }#end while



	#SELECT GENE ID 1'S WHICH AFFECT AT LEAST ONE SINGLE KO FLUX DISTRIBUTION
	my $select_pa_genes = "	(SELECT DISTINCT g.ORF_id
				FROM single_KO_flux ko1 , double_KO_flux ko2, `condition` c, gene g
				WHERE ko1.gid = ko2.gid2 AND
      				      ko2.flux_vector != ko1.flux_vector AND
       				      ko1.con_id = ko2.con_id AND
                                      ko2.con_id = c.con_id AND
				      g.gid = ko2.gid1 AND
                                      c.condition_name = ? AND
				      ko1.model_id = ko2.model_id AND
                                      ko2.model_id = ?)";


        my $sth = $dbh -> prepare($select_pa_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($cond, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my $gene =  $sth -> fetchrow_array)
        {

                $phenotypic_altering_genes{$genes} = 1;

        }#end while



        return (\%phenotypic_altering_genes);

}#end get_phenotypic_altering_genes


#sub get_complex_rep_gene_list - Takes as input: 1) model_id	- An ID from the model table
#						 2) dbh		- A database handle
#			       - Returns a list of genes with only a single gene from a complex being used
#				 to represent all the genes in that complex
sub get_complex_rep_gene_list
{
	my ($model_id, $dbh) = @_;

	my $r_prots2genes;
	my %prots2genes;

	my %genes;

	my %complex_sets;


	#GET A LIST OF ALL GENES IN THE DATABASE
	$r_prots2genes = &get_genes_and_proteins_from_FBA_DB($model_id, $dbh);
	%prots2genes = %$r_prots2genes;


	#FOR EACH GENE GET THE REACTIONS WHICH ARE DEPENDANT ON IT
	foreach my $prot (keys %prots2genes)
	{

		#GET ALL REACTION EXPRESSION THAT A GENE IS PART OF
        	my $select_expr = "SELECT DISTINCT rc.complex_expr
                 	           FROM rxn2complex_expr rc, rxn r
                        	   WHERE (rc.complex_expr LIKE '$prot \%' OR
					 rc.complex_expr LIKE '\% $prot \%' OR
					 rc.complex_expr LIKE '\% $prot' OR
					 rc.complex_expr = '$prot') AND
					 r.rid = rc.rid AND
					 r.model_id = ?
				   ORDER BY rc.complex_expr";


        	my $sth = $dbh -> prepare($select_expr)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($model_id)
        	                 or die "Can't execute statment: $DBI::errstr";


		my @complexes;
        	while(my ($complex) =  $sth -> fetchrow_array)
        	{

			push @complexes , $complex;

		}#end while

		my $complexes = join "_" , @complexes;

		#ADD GENE TO LIST IF THE EXACT SAME EXPRESSION SET DOES NOT ALREADY EXIST
		if( !exists( $complex_sets{$complexes} ) )
		{

			$genes{$prots2genes{$prot}} = 1;

			$complex_sets{$complexes} = 1;

		}#end if


	}#end foreach


	return \%genes;

}#end get_complex_rep_gene_list


#sub get_wildtype_fluxes - Takes as input : 1) con_id 		- The ID for the condition of interest
#                                           2) r_rxn_codes 	- The name of the flux whose value should be retrieved
#					    3) r_rxn_list 	- A list of reactions in order of readtions in the flux vector	
#                                           4) dbh 		- A database handle
#					    5) opt_type 	- The optimization type used to generate the flux data
#					    6) model_id		- An ID from the model table
#                        - Returns the value of the desired flux
sub get_wildtype_fluxes
{
        my($con_id , $r_rxn_codes , $r_rxn_list , $dbh , $opt_type, $model_id) = @_;

	my @wt_fluxes;

	my $r_rxn_code2rid;
	my %rxn_code2rid;
	
	my %rxn_indexes;
	my $rxn_index = 0;

	my @rxn_codes = @$r_rxn_codes;
	my @rxn_list = @$r_rxn_list;


        #IF OPTIMIZATION TYPE NOT DEFINED, THEN DEFAULT TO FBA
        if( $opt_type eq "")
        {

                $opt_type = "FBA";

        }#end if


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



        #GET THE VALUE OF THE WILDTYPE FLUXES FOR THE DESIRED CONDITIONS
        my $select_wt_flux = "SELECT wt.flux_vector
                                FROM wildtype_flux wt 
                                WHERE wt.con_id = ? AND
				      wt.opt_type = ? AND
				      wt.model_id = ?";


        my $sth = $dbh -> prepare($select_wt_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($con_id, $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUXES
        while(my ($flux_vector) =  $sth -> fetchrow_array)
        {


                chomp $flux_vector;
		#$flux_vector =~ s/ //g;
		$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;

                my @flux_vector = split /\t/ , $flux_vector;

                foreach my $rxn_code (@rxn_codes)
                {

                        push @wt_fluxes , $flux_vector[ $rxn_indexes{$rxn_code} ];

			#print "rxn code = $rxn_code  rxn ind = $rxn_indexes{$rxn_code}  flux = $wt_fluxes[0]\n";
                }


        }#end while


        return \@wt_fluxes;


}#end get_wildtype_fluxes



#sub get_wildtype_fluxes_RND - Takes as input : 1) con_id 	- The ID for the condition of interest
#                                               2) r_rxn_codes 	- The name of the flux whose value should be retrieved
#				  	        3) r_rxn_list 	- A list of reactions in order of readtions in the flux vector	
#						4) precision 	- The number of places after the decimal places to round the fluxes 
#							          to
#                                               5) dbh 		- A database handle
#					        6) opt_type 	- The optimization type used to generate the flux data
#						7) model_id	- An ID from the model table
#                        - Returns the value of the desired flux
sub get_wildtype_fluxes_RND
{
        my($con_id , $r_rxn_codes , $r_rxn_list , $precision , $dbh , $opt_type, $model_id) = @_;

	my @wt_fluxes;

	my $r_rxn_code2rid;
	my %rxn_code2rid;
	
	my %rxn_indexes;
	my $rxn_index = 0;

	my @rxn_codes = @$r_rxn_codes;
	my @rxn_list = @$r_rxn_list;


        #IF OPTIMIZATION TYPE NOT DEFINED, THEN DEFAULT TO FBA
        if( $opt_type eq "")
        {

                $opt_type = "FBA";

        }#end if


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



        #GET THE VALUE OF THE WILDTYPE FLUXES FOR THE DESIRED CONDITIONS
        my $select_wt_flux = "SELECT wt.flux_vector
                                FROM wildtype_flux wt 
                                WHERE wt.con_id = ? AND
				      wt.opt_type = ? AND
			              wt.model_id = ?";


        my $sth = $dbh -> prepare($select_wt_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($con_id, $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUXES
        while(my ($flux_vector) =  $sth -> fetchrow_array)
        {


                chomp $flux_vector;
		#$flux_vector =~ s/ //g;
		$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;

                my @flux_vector = split /\t/ , $flux_vector;
		
		my $r_flux_vector_rnd = &list_ops::rnd(\@flux_vector , $precision);
		my @flux_vector_rnd = @$r_flux_vector_rnd;

                foreach my $rxn_code (@rxn_codes)
                {

                        push @wt_fluxes , $flux_vector_rnd[ $rxn_indexes{$rxn_code} ];

			#print "rxn code = $rxn_code  rxn ind = $rxn_indexes{$rxn_code}  flux = $wt_fluxes[0]\n";
                }


        }#end while


        return \@wt_fluxes;


}#end get_wildtype_fluxes_RND



#sub get_wildtype_fluxes_hash - Takes as input: 1) con_id	- THe ID for the condition under which wildtype fluxes are desired
#						2) opt_type	- The opt_type by which wildtype fluxes were computed
#						3) model_id	- An ID from the model table
#						4) dbh		- A database handle
#			      - Returns a reference to a hash where the keys are rxn codes and values are corresponding wildtype
#				flux values
sub get_wildtype_fluxes_hash
{
	my($condition, $opt_type, $model_id, $dbh) = @_;
	my %WT_rxn2flux;


	#GET ALL REACTIONS FROM THE DATABASE IN THE SAME ORDER AS IN THE FLUX VECTOR
	my $r_rxns = &get_rxns_from_SBML_DB($model_id, $dbh);
	my @rxns = @$r_rxns;


	#GET THE WILDTYPE FLUXES
	my $WT_fluxV = &get_wildtype_flux_vector($condition , $dbh , $opt_type, $model_id);
	my @WT_fluxes = split /\t/ , $WT_fluxV;


	#PUT THE FLUXES IN A HASH KEYED BY THE REACTIONS
	for(my $r=0; $r < (@rxns); $r++)
	{

		#print "$rxns[$r] = $WT_fluxes[$r]\n";
		$WT_rxn2flux{ $rxns[$r] } = $WT_fluxes[$r];

	}#end for 	


	return \%WT_rxn2flux;

}#end  get_wildtype_fluxes_hash



#sub get_single_KO_fluxes - Takes as input : 1) con_id 		- The ID for the condition of interest
#                                            2) r_rxn_codes 	- The name of the fluxes whose value should be retrieved
#			  		     3) opt_type 	- The optimization method (e.g. FBA , MOMA)
#					     4) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#							          flux vector in the database
#					     5) model_id	- An ID from the model table
#                                            6) dbh 		- A database handle
#                         - Returns a reference to a hash with the key being gene ids and the
#                           values being the value of the flux in a given KO
sub get_single_KO_fluxes
{
        my($con_id , $r_rxn_codes , $opt_type , $r_rxn_list , $model_id,  $dbh) = @_;

        my %flux_KO;

        my $r_rxn_code2rid;
        my %rxn_code2rid;
        my $rxn_index = 0;
        my %rxn_indexes;


	my @rxn_codes = @$r_rxn_codes;
	my @rxn_list = @$r_rxn_list;


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach


        #GET ALL THE FLUX VECTORS FOR SINGLE KOS UNDER THE DESIRED CONDITION
        my $select_1KOflux = "SELECT ko1.gid , ko1.flux_vector
                                FROM single_KO_flux ko1
                                WHERE ko1.con_id = ? AND
				      ko1.opt_type = ? AND
				      ko1.model_id = ?";


        my $sth = $dbh -> prepare($select_1KOflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($con_id , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($gid , $flux_vector) =  $sth -> fetchrow_array)
        {

		
		chomp $flux_vector;
		$flux_vector =~ s/(\d)\s+([-\d])/$1\t$2/g;
		#$flux_vector =~ s/ //g;

                my @flux_vector = split /\t/ , $flux_vector;

		foreach my $rxn_code (@rxn_codes)
		{
                
			push @{$flux_KO{$gid}} , $flux_vector[ $rxn_indexes{$rxn_code} ];
	
		}


        }#end while



        return \%flux_KO;


}#end get_single_KO_fluxes



#sub get_single_KO_fluxes_RND - Takes as input : 1) con_id 	- The ID for the condition of interest
#                                                2) r_rxn_codes - The name of the fluxes whose value should be retrieved
#			  	   	         3) opt_type 	- The optimization method (e.g. FBA , MOMA)
#					         4) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#				  	 	 	          flux vector in the database
#						 5) precision 	- The number of places after the decimal point which the fluxes 
#							          should be rounded to
#						 6) model_id	- An ID from th model table
#                                                7) dbh 	- A database handle
#                         - Returns a reference to a hash with the key being gene ids and the
#                           values being the value of the flux in a given KO
sub get_single_KO_fluxes_RND
{
        my($con_id , $r_rxn_codes , $opt_type , $r_rxn_list , $precision , $model_id,  $dbh) = @_;

        my %flux_KO;

        my $r_rxn_code2rid;
        my %rxn_code2rid;
        my $rxn_index = 0;
        my %rxn_indexes;


	my @rxn_codes = @$r_rxn_codes;
	my @rxn_list = @$r_rxn_list;


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



        #GET ALL THE FLUX VECTORS FOR SINGLE KOS UNDER THE DESIRED CONDITION
        my $select_1KOflux = "SELECT ko1.gid , ko1.flux_vector
                                FROM single_KO_flux ko1
                                WHERE ko1.con_id = ? AND
				      ko1.opt_type = ? AND
				      ko1.model_id = ?";


        my $sth = $dbh -> prepare($select_1KOflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($con_id , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($gid , $flux_vector) =  $sth -> fetchrow_array)
        {

		
		chomp $flux_vector;
		$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;
		#$flux_vector =~ s/ //g;

                my @flux_vector = split /\t/ , $flux_vector;

		my $r_flux_vector_rnd = &list_ops::rnd(\@flux_vector , $precision);
		my @flux_vector_rnd = @$r_flux_vector_rnd;


		foreach my $rxn_code (@rxn_codes)
		{
                
			push @{$flux_KO{$gid}} , $flux_vector_rnd[ $rxn_indexes{$rxn_code} ];
	
		}


        }#end while



        return \%flux_KO;


}#end get_single_KO_fluxes_RND




#sub get_single_KO_flux_GL - Takes as input : 1) con_id 	- The ID for the condition of interest
#                                             2) r_rxn_codes 	- The name of the fluxes whose value should be retrieved
#			   		      3) opt_type 	- The optimization method (e.g. FBA , MOMA)
#					      4) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#							          flux vector in the database
#					      5)gene_list 	- A list of gene KOs to retrive flux data for
#					      6) model_id	- An ID from the model table
#                                             7) dbh 		- A database handle
#                          - Returns a reference to a hash with the key being gene ids and the
#                            values being the value of the flux in a given KO
sub get_single_KO_fluxes_GL
{
        my($con_id , $r_rxn_codes , $opt_type , $r_rxn_list , $r_gene_list , $model_id,  $dbh) = @_;

        my %flux_KO;

        my $r_rxn_code2rid;
        my %rxn_code2rid;
        my $rxn_index = 0;
        my %rxn_indexes;


	my @rxn_codes = @$r_rxn_codes;
	my @rxn_list = @$r_rxn_list;


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



        #GET ALL THE FLUX VECTORS FOR SINGLE KOS UNDER THE DESIRED CONDITION
        my $select_1KOflux = "SELECT g.ORF_id , ko1.flux_vector
                                FROM single_KO_flux ko1 , gene g
                                WHERE ko1.con_id = ? AND
				      ko1.opt_type = ? AND
				      g.gid = ko1.gid AND
				      ko1.model_id = ? AND
				      g.ORF_id IN (";

        #ADD PLACE HOLDERS FOR GENES
        $select_1KOflux = &db_lib::add_place_holders($select_1KOflux , scalar(@$r_gene_list));

        $select_1KOflux = $select_1KOflux . ")";


        my $sth = $dbh -> prepare($select_1KOflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($con_id , $opt_type , $model_id,  @$r_gene_list)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($gene , $flux_vector) =  $sth -> fetchrow_array)
        {

		
		chomp $flux_vector;
		$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;
		#$flux_vector =~ s/ //g;

                my @flux_vector = split /\t/ , $flux_vector;

		foreach my $rxn_code (@rxn_codes)
		{
                
			push @{$flux_KO{$gene}} , $flux_vector[ $rxn_indexes{$rxn_code} ];
	
		}


        }#end while



        return \%flux_KO;


}#end get_single_KO_flux_GL



#sub get_double_KO_fluxes - Takes as input :  1) con_id 	- The ID for the condition of interest
#                                             2) r_rxn_codes 	- The name of the fluxes whose value should be retrieved
#                                             3) opt_type 	- The optimization method (e.g. FBA , MOMA)
#                                             4) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#                                                                 flux vector in the database
#					      5) model_id	- An ID from the model table
#                                             6) dbh 		- A database handle

#			  - Returns a reference to a 2-dimensional hash with the keys being genes and the
#			    values being the value of the flux in a double KO
sub get_double_KO_fluxes
{
	my($con_id , $r_rxn_codes , $opt_type , $r_rxn_list , $r_gene_list , $model_id,  $dbh) = @_;

	my %flux_KO2;
	
	my $r_rxn_code2rid;
	my %rxn_code2rid;
	my $rxn_index = 0;
	my %rxn_indexes;

	my %gene2dbid;
	my $r_gene2dbid;

	my @gene_list = @$r_gene_list;

        my @rxn_codes = @$r_rxn_codes;
        my @rxn_list = @$r_rxn_list;


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



	#GET A LIST OF GIDS SO THAT FLUX CAN BE EXTRACTED ONE GENE AT A TIME TO AVOID RETURNING 
	#HUGE AMOUNTS OF DATA!!!!
	$r_gene2dbid = &get_model_genes_and_ids_from_FBA_DB($model_id, $dbh);
	%gene2dbid = %$r_gene2dbid;


	foreach my $gene (@gene_list)
	#foreach my $gene (sort keys %gene2dbid )
	{

		#for(my $g = 1; $g <=2; $g++)
		#{
	
			#GET ALL THE FLUX VECTORS FOR DOUBLE KOS UNDER THE DESIRED CONDITION
			my $select_2KOflux = "SELECT ko2.gid1 , ko2.gid2 , ko2.flux_vector
		                                FROM double_KO_flux ko2
	        	                        WHERE  ko2.con_id = ? AND
						       ko2.gid1 = ? AND
						       ko2.opt_type = ? AND
						       ko2.model_id = ?";


	        	my $sth = $dbh -> prepare($select_2KOflux)
	        	                 or die "Can't prepare statment: $DBI::errstr";

	        	my $rc = $sth -> execute($con_id , $gene2dbid{$gene} , $opt_type, $model_id)
	        	                 or die "Can't execute statment: $DBI::errstr";



			#SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
	        	while(my ($gid1 , $gid2 , $flux_vector) =  $sth -> fetchrow_array)
	        	{


				chomp($flux_vector);
				$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;

	               	 	my @flux_vector = split /\t/ , $flux_vector;


				for(my $r =0; $r < (@rxn_codes); $r++)
	                	{
	
	 	                       	@{$flux_KO2{$gid1}{$gid2}}[$r] = $flux_vector[ $rxn_indexes{$rxn_codes[$r]} ];
	 	                       	@{$flux_KO2{$gid2}{$gid1}}[$r] = $flux_vector[ $rxn_indexes{$rxn_codes[$r]} ];

        	        	}


        		}#end while


		#}#end for


	}#end foreach



	return \%flux_KO2;


}#end get_double_KO_fluxes


#sub get_double_KO_fluxes_RND - Takes as input :  1) con_id 		- The ID for the condition of interest
#                                                 2) r_rxn_codes 	- The name of the fluxes whose value should be retrieved
#                                                 3) opt_type 		- The optimization method (e.g. FBA , MOMA)
#                                                 4) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#                                                              	 	  flux vector in the database
#						  5) precision 		- THe number of places after the decimal which the fluxes should
#							         	  be rounded to
#						  6) model_id		- An ID from the model table
#                                                 7) dbh 		- A database handle

#			      - Returns a reference to a 2-dimensional hash with the keys being genes and the
#			        values being the value of the flux in a double KO
sub get_double_KO_fluxes_RND
{
	my($con_id , $r_rxn_codes , $opt_type , $r_rxn_list , $r_gene_list , $precision , $model_id,  $dbh) = @_;

	my %flux_KO2;
	
	my $r_rxn_code2rid;
	my %rxn_code2rid;
	my $rxn_index = 0;
	my %rxn_indexes;

	my %gene2dbid;
	my $r_gene2dbid;

	my @gene_list = @$r_gene_list;

        my @rxn_codes = @$r_rxn_codes;
        my @rxn_list = @$r_rxn_list;


        #GET POSITION OF EACH FLUX IN THE FLUX VECTOR
        foreach my $curr_rxn_code (@rxn_list)
        {

                $rxn_indexes{$curr_rxn_code} = $rxn_index;

                $rxn_index ++;

        }#end foreach



	#GET A LIST OF GIDS SO THAT FLUX CAN BE EXTRACTED ONE GENE AT A TIME TO AVOID RETURNING 
	#HUGE AMOUNTS OF DATA!!!!
	$r_gene2dbid = &get_model_genes_and_ids_from_FBA_DB($model_id, $dbh);
	%gene2dbid = %$r_gene2dbid;


	foreach my $gene (@gene_list)
	#foreach my $gene (sort keys %gene2dbid )
	{

		#for(my $g = 1; $g <=2; $g++)
		#{
	
			#GET ALL THE FLUX VECTORS FOR DOUBLE KOS UNDER THE DESIRED CONDITION
			my $select_2KOflux = "SELECT ko2.gid1 , ko2.gid2 , ko2.flux_vector
		                                FROM double_KO_flux ko2
	        	                        WHERE  ko2.con_id = ? AND
						       ko2.gid1 = ? AND
						       ko2.opt_type = ? AND
						       ko2.model_id = ?";


	        	my $sth = $dbh -> prepare($select_2KOflux)
	        	                 or die "Can't prepare statment: $DBI::errstr";

	        	my $rc = $sth -> execute($con_id , $gene2dbid{$gene} , $opt_type, $model_id)
	        	                 or die "Can't execute statment: $DBI::errstr";



			#SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
	        	while(my ($gid1 , $gid2 , $flux_vector) =  $sth -> fetchrow_array)
	        	{


				chomp($flux_vector);
				$flux_vector =~ s/(\d)\s+(\d)/$1\t$2/g;

	               	 	my @flux_vector = split /\t/ , $flux_vector;
				
				my $r_flux_vector_rnd = &list_ops::rnd(\@flux_vector , $precision);
				my @flux_vector_rnd = @$r_flux_vector_rnd;


				for(my $r =0; $r < (@rxn_codes); $r++)
	                	{
	
	 	                       	@{$flux_KO2{$gid1}{$gid2}}[$r] = $flux_vector_rnd[ $rxn_indexes{$rxn_codes[$r]} ];
	 	                       	@{$flux_KO2{$gid2}{$gid1}}[$r] = $flux_vector_rnd[ $rxn_indexes{$rxn_codes[$r]} ];

        	        	}


        		}#end while


		#}#end for


	}#end foreach



	return \%flux_KO2;


}#end get_double_KO_fluxes_RND


#sub get_KO_fitness - Takes as input : 1) r_genes 	- A list of ORF ids
#				       2) r_conds 	- A list of condition names
#				       3) flux_name 	- The name of the flux to compute fitness in reference to
#				       4) opt_type 	- Optimization type (eg. FBA, MOMA)
#				       5) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#					                  flux vector
#				       6) model_id	- An ID from the model table
#				       7) dbh 		- A database handle
#		    - Returns a 2 dimensional hash where the keys are condition/gene and the value is
#		      the fitness of the KO of the gene under the given condition
sub get_KO_fitness
{
	my ($r_genes , $r_conds , $flux_name , $opt_type, $r_rxn_list , $model_id, $dbh) = @_;

	my @genes = @$r_genes;
	my @conds = @$r_conds;

	my $r_rxn_code2rid;
	my %rxn_code2rid;
	
	my $flux_index;
	my $rxn_index = 0;

	my %WT_flux;
	my %KO_fitness;
	
	my @rxn_list = @$r_rxn_list;


        foreach my $curr_rxn_code (@rxn_list)
        {


                if( $curr_rxn_code eq $flux_name )
                {

               		$flux_index = $rxn_index;

                }#end if


                $rxn_index ++;

        }#end foreach


	#GET THE WILDTYPE FLUX FOR EACH CONDITION
        my $select_wtflux = "SELECT c.condition_name , wt.flux_vector
                             FROM wildtype_flux wt , `condition` c
                             WHERE wt.con_id = c.con_id AND
				   wt.model_id = ? AND
				   c.condition_name IN (";

        #ADD PLACE HOLDERS FOR CONDITIONS
        $select_wtflux = &db_lib::add_place_holders($select_wtflux , scalar(@conds));

        $select_wtflux = $select_wtflux . ")";



        my $sth = $dbh -> prepare($select_wtflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id, @conds)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($cond_name , $flux_vector) =  $sth -> fetchrow_array)
        {


                chomp $flux_vector;
                my @flux_vector = split /\t/ , $flux_vector;

		$WT_flux{$cond_name} = $flux_vector[$flux_index];


        }#end while



	#GET THE SINGLE KO FLUX FOREACH GENE, UNDER EACH CONDITION, FOR THE GIVEN OPTIMIZATION TYPE AND COMPUTE THE FITNESS
        my $select_KOflux = "SELECT c.condition_name , g.ORF_id , ko.flux_vector
                             FROM single_KO_flux ko , `condition` c, gene g
                             WHERE ko.con_id = c.con_id AND
				   ko.opt_type = ? AND
				   ko.model_id = ? AND
                                   c.condition_name IN (";

        #ADD PLACE HOLDERS FOR CONDITIONS
        $select_KOflux = &db_lib::add_place_holders($select_KOflux , scalar(@conds));

        $select_KOflux = $select_KOflux . ") AND
				   ko.gid = g.gid AND
				   g.ORF_id IN (";

	#ADD PLACE HOLDERS FOR GENES
        $select_KOflux = &db_lib::add_place_holders($select_KOflux , scalar(@genes));

        $select_KOflux = $select_KOflux . ")";



        my $sth = $dbh -> prepare($select_KOflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($opt_type , $model_id, @conds , @genes)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($cond_name , $orf_id , $flux_vector) =  $sth -> fetchrow_array)
        {


                chomp $flux_vector;
                my @flux_vector = split /\t/ , $flux_vector;

		if( $WT_flux{$cond_name} != 0 )
		{

	                $KO_fitness{$orf_id}{$cond_name} = $flux_vector[$flux_index] / $WT_flux{$cond_name};

		}else{


			 $KO_fitness{$orf_id}{$cond_name} = 1;

		}


        }#end while


	return \%KO_fitness;

}#end get_KO_fitness


#sub get_KO_fitness_cond_spec - Takes as input : 1) r_genes 	- A list of ORF ids
#				  	         2) r_conds 	- A list of condition names
#				       		 3) flux_name 	- The name of the flux to compute fitness in reference to
#				       		 4) opt_type 	- Optimization type (eg. FBA, MOMA)
#				       		 5) r_rxn_list 	- A list of reactions in the same order as the reactions in the
#					          	         flux vector
#				       		 6) model_id	- An ID from the model table
#				       		 7) dbh 		- A database handle
#		    	      - Returnse condition/gene and the value is
#		      		the fitness of the KO of the gene under the given condition
sub get_KO_fitness_cond_spec
{
	my ($cond , $flux_name , $opt_type, $r_rxn_list , $model_id, $dbh) = @_;

	my $r_rxn_code2rid;
	my %rxn_code2rid;
	
	my $flux_index;
	my $rxn_index = 0;

	my %WT_flux;
	my %KO_fitness;
	
	my @rxn_list = @$r_rxn_list;

	#GET THE INDEX OF THE FITNES FLUX IN THE REACTIO LIST
        foreach my $curr_rxn_code (@rxn_list)
        {


                if( $curr_rxn_code eq $flux_name )
                {

               		$flux_index = $rxn_index;

                }#end if


                $rxn_index ++;

        }#end foreach


	#GET THE WILDTYPE FLUX FOR EACH CONDITION
        my $select_wtflux = "SELECT wt.flux_vector
                             FROM wildtype_flux wt , `condition` c
                             WHERE wt.con_id = c.con_id AND
				   wt.model_id = ? AND
				   c.condition_name = ?";


        my $sth = $dbh -> prepare($select_wtflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id, $cond)
                         or die "Can't execute statment: $DBI::errstr";



        my $flux_vector =  $sth -> fetchrow_array;
        chomp $flux_vector;

        my @flux_vector = split /\t/ , $flux_vector;
	$WT_flux = $flux_vector[$flux_index];


	#GET THE SINGLE KO FLUX FOREACH GENE, UNDER EACH CONDITION, FOR THE GIVEN OPTIMIZATION TYPE AND COMPUTE THE FITNESS
        my $select_KOflux = "SELECT g.ORF_id , ko.flux_vector
                             FROM single_KO_flux ko , `condition` c, gene g
                             WHERE ko.con_id = c.con_id AND
				   ko.opt_type = ? AND
				   ko.model_id = ? AND
                                   c.condition_name = ? AND
				   ko.gid = g.gid";


        my $sth = $dbh -> prepare($select_KOflux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($opt_type , $model_id, $cond)
                         or die "Can't execute statment: $DBI::errstr";



        #SAVE THE DESIRED FLUX FOR EVERY DOUBLE KO IN THE DATABASE
        while(my ($orf_id , $flux_vector) =  $sth -> fetchrow_array)
        {
                
		chomp $flux_vector;
                my @flux_vector = split /\t/ , $flux_vector;

		if( $WT_flux != 0 )
		{

	                $KO_fitness{$orf_id} = $flux_vector[$flux_index] / $WT_flux;

		}else{


			 $KO_fitness{$orf_id} = 1;

		}


        }#end while


	return \%KO_fitness;

}#end get_KO_fitness_cond_spec


#sub get_wildtype_flux_vector - Takes as input : 1) condition 	-The name of the condition of interest
#                                                2) dbh 	- A database handle
#					 	 3) opt_type	- An optimization crteria
#						 4) model_id	- An ID from the model table	
#                      	      - Returns a scalar containing the tab delimited fluxes for the 
#				desired condtion 
#                       
sub get_wildtype_flux_vector
{
        my($condition  , $dbh , $opt_type, $model_id) = @_;

        my $wt_flux;


	#IF OPTIMIZATION TYPE NOT DEFINED, THEN DEFAULT TO FBA
	if( $opt_type eq "")
	{

		$opt_type = "FBA";

	}#end if


        #GET THE VALUE OF THE WILDTYPE FLUX FOR THE DESIRED CONDITION
        my $select_wt_flux = "SELECT wt.flux_vector
                                FROM wildtype_flux wt , `condition` c
                                WHERE wt.con_id = c.con_id AND
				      c.condition_name = ? AND
				      wt.opt_type = ? AND
				      wt.model_id = ?";


        my $sth = $dbh -> prepare($select_wt_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition, $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        $wt_flux =  $sth -> fetchrow_array;
        chomp $wt_flux;

        return $wt_flux;


}#end get_wildtype_flux_vector



#sub get_singleKO_flux_vector - Takes as input : 1) condition 		- The name of the condition of interest
#						 2) gene_name 		- The name of the gene deletion of interest
#						 3) opt_type 		- The optimzation method used for KO flux calculation
#						 4) model_id		- An ID from the model table
#                                                5) dbh 		- A database handle
#                             - Returns a scalar containing the tab delimited fluxes for the
#                               desired condtion and KO
#
sub get_singleKO_flux_vector
{
        my($condition  , $gene_name , $opt_type, $model_id, $dbh) = @_;

        my $ko_flux;


        #GET THE VALUE OF THE WILDTYPE FLUX FOR THE DESIRED CONDITION
        my $select_KO_flux = "SELECT ko.flux_vector
                                FROM single_KO_flux ko , `condition` c, gene g
                                WHERE ko.con_id = c.con_id AND
				      ko.gid = g.gid AND
				      c.condition_name = ? AND
				      g.ORF_id = ? AND
				      ko.opt_type = ? AND
				      ko.model_id = ?";


        my $sth = $dbh -> prepare($select_KO_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition , $gene_name, $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        $ko_flux =  $sth -> fetchrow_array;
        chomp $ko_flux;

        return $ko_flux;


}#end get_singleKO_flux_vector


#sub get_all_singleKO_flux_vectors - Takes as input : 1) condition 	- The name of the condition of interest
#				 	  	      2) opt_type 	- The optimzation method used for KO flux calculation
#						      3) model_id	- An ID from the model table
#                                                     4) dbh 		- A database handle
#                                  - Returns a scalar containing the tab delimited fluxes for the
#                                    desired condtion and KO
#
sub get_all_singleKO_flux_vectors
{
        my($condition , $opt_type, $model_id,  $dbh) = @_;

        my %ko_fluxes;


        #GET THE VALUE OF THE WILDTYPE FLUX FOR THE DESIRED CONDITION
        my $select_KO_flux = "SELECT g.orf_id  , ko.flux_vector
                                FROM single_KO_flux ko , `condition` c, gene g
                                WHERE ko.con_id = c.con_id AND
				      ko.gid = g.gid AND
				      c.condition_name = ? AND
				      ko.opt_type = ? AND
				      ko.model_id = ?";


        my $sth = $dbh -> prepare($select_KO_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



	while(my ($orf_id , $flux_vector) =  $sth -> fetchrow_array)
        {

		chomp $flux_vector;

		$ko_fluxes{$orf_id} = $flux_vector;

	}#end while


        return \%ko_fluxes;


}#end get_all_singleKO_flux_vectors


#sub get_all_doubleKO_flux_vectors - Takes as input : 1) condition 	- The name of the condition of interest
#                                    		      2) gene_name 	- The name of the gene deletion of interest
#			  	 		      3) opt_type 	- The optimzation method used for KO flux calculation
#						      4) model_id	- An ID from the model table
#                                          	      5) dbh 		- A database handle
#
#                             	    - Returns a hash where the key is an ORF ID and the value is a tab delimited flux
#				      vector of the double mutant with the inputted and keyed genes removed. If the key 
#				      is the inputted gene then the value is the single KO flux distribution.
#                             
#
sub get_all_doubleKO_flux_vectors
{
        my($condition  , $gene_name ,  $opt_type ,$model_id,  $dbh) = @_;

        my %ko_fluxes;


	#GET THE SINGLE KO FLUX DISTRIBUTION
	$ko_fluxes{$gene_name} = &get_singleKO_flux_vector($condition  , $gene_name , $opt_type, $model_id,  $dbh);


        #GET THE VALUE OF THE WILDTYPE FLUX FOR THE DESIRED CONDITION
        my $select_KO_flux = "SELECT ko.flux_vector , g2.ORF_id
                                FROM double_KO_flux ko , `condition` c, gene g1 , gene g2
                                WHERE ko.con_id = c.con_id AND
                                      ( (ko.gid1 = g1.gid AND g2.gid = ko.gid2) 
					       OR
					(ko.gid2 = g1.gid AND g2.gid = ko.gid1) ) AND
				      c.condition_name = ? AND
				      g1.ORF_id = ? AND
				      ko.opt_type = ? AND
				      ko.model_id = ?";
				      

        my $sth = $dbh -> prepare($select_KO_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition, $gene_name , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



	while(my ($flux_vector , $orf_id) =  $sth -> fetchrow_array)
        {

		chomp $flux_vector;

		$ko_fluxes{$orf_id} = $flux_vector;

	}#end while


        return \%ko_fluxes;


}#end get_all_doubleKO_flux_vectors


#sub get_doubleKO_flux_vector - Takes as input : 1) condition 		- The name of the condition of interest
#                                                2) gene_name1 		- The name of the first gene deletion of interest
#                                                3) gene_name2 		- The name of the second gene deletion of interest
#						 4) opt_type 		- The optimzation method used for KO flux calculation
#						 5) model_id		- An ID form the model table
#                                                6) dbh 		- A database handle
#
#                             - Returns a scalar containing the tab delimited fluxes for the
#                               desired condtion and KO
#
sub get_doubleKO_flux_vector
{
        my($condition  , $gene_name1 , $gene_name2 ,  $opt_type , $model_id,  $dbh) = @_;

        my $ko_flux;


        #GET THE VALUE OF THE WILDTYPE FLUX FOR THE DESIRED CONDITION
        my $select_KO_flux = "SELECT ko.flux_vector
                                FROM double_KO_flux ko , `condition` c, gene g1, gene g2
                                WHERE ko.con_id = c.con_id AND
                                      ( (ko.gid1 = g1.gid AND
				         ko.gid2 = g2.gid) 
					       OR
					(ko.gid1 = g2.gid AND
                                      	ko.gid2 = g1.gid) ) AND
				      c.condition_name = ? AND
				      g1.ORF_id = ? AND
				      g2.ORF_id = ? AND
				      ko.opt_type = ? AND
				      ko.model_id = ?";
				      

        my $sth = $dbh -> prepare($select_KO_flux)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($condition, $gene_name1 , $gene_name2 , $opt_type, $model_id)
                         or die "Can't execute statment: $DBI::errstr";



        $ko_flux =  $sth -> fetchrow_array;
        chomp $ko_flux;

        return $ko_flux;


}#end get_doubleKO_flux_vector


################################################## UPDATE SUBROUTINES ######################################

#sub update_rxn_procs - Takes as input a database handle
#		      - Updates all subsystems in the rxn table to have no leading white space
sub update_rxn_procs
{
	my ($dbh) = @_;

	my %rxn2sys;


	#GET ALL REACTION IDS AND THE CORRESPONDING SUBSYSTEMS
        my $select_rxn = "SELECT rid , subsystem
                          FROM rxn"; 


        my $sth = $dbh -> prepare($select_rxn)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rid , $subsystem) =  $sth -> fetchrow_array)
        {

                $rxn2sys{$rid} = $subsystem;

        }#end while


	#UPDATE EACH REACTION
	foreach my $rid (keys %rxn2sys)
	{
		
		my $subsys = $rxn2sys{$rid};

		$subsys =~ s/^\s+//g;
		$subsys =~ s/\s+$//g;

		my $update_sys = "UPDATE rxn
				  SET subsystem = ?
				  WHERE rid = ?";

	        my $sth = $dbh -> prepare($update_sys)
                         or die "Can't prepare statment: $DBI::errstr";

	        my $rc = $sth -> execute($subsys , $rid)
                         or die "Can't execute statment: $DBI::errstr";


	}#end foreach


}#end update_rxn_procs


#sub update_rxn_eq_procs - Takes as input: 1) rxn2proc	- A hash where the keys are rxn equations and the value 'PROC' is the
#							  new process assignment for that reaction
#					   2) dbh 	-a database handle
#                        - Updates all the appropriate subsystems in the rxn table 
sub update_rxn_eq_procs
{
        my ($r_rxn2proc , $dbh) = @_;

	my %rxn2proc = %$r_rxn2proc;


        #UPDATE EACH REACTION
        foreach my $rxn (keys %rxn2proc)
        {

                my $proc = $rxn2proc{$rxn}{'PROC'};


                my $update_proc = "UPDATE rxn
                                  SET subsystem = ?
                                  WHERE rxn = ?";

                my $sth = $dbh -> prepare($update_proc)
                         or die "Can't prepare statment: $DBI::errstr";

                my $rc = $sth -> execute($proc , $rxn)
                         or die "Can't execute statment: $DBI::errstr";


        }#end foreach


}#end update_rxn_procs


#sub update_exchange_fluxes - Takes as input: 1) model_id	- ID of model
#					      2) dbh		- A database handle
#			    - Updates all exchange flux reactions to have a * before the reaction arrow 
sub update_exchange_rxns
{
	my ($model_id, $dbh) = @_;


        #GET ALL REACTION IDS AND THE CORRESPONDING equations
        my $select_rxn = "SELECT rid , rxn
                          FROM rxn
			  WHERE model_id = ? AND 
			        rxn_code LIKE \'%EX_%\'";


        my $sth = $dbh -> prepare($select_rxn)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($model_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($rid , $rxn) =  $sth -> fetchrow_array)
        {

                $rid2rxn{$rid} = $rxn;

        }#end while


	#UPDATE REACTION EQUATIONS
	foreach my $rid ( keys %rid2rxn)
	{

		my $rxn = $rid2rxn{$rid};

		$rxn =~ s/^(.*) <->/* <->/;

		my $update_rxn = "UPDATE rxn
                                  SET rxn = ?
                                  WHERE rid = ?";

                my $sth = $dbh -> prepare($update_rxn)
                         or die "Can't prepare statment: $DBI::errstr";

                my $rc = $sth -> execute($rxn , $rid)
                         or die "Can't execute statment: $DBI::errstr";


	}#end foreach


}#end update_exchange_rxns


########################################## SUBROUTINES TO TRANSFER DATA FROM SBML DATABASE ###########################

#sub tranfer_SBML_to_FBA - Takes as input: 1) model_id	- An ID from the model table 
#					   2) dbh_FBA	- A handle to the FBA model database
#					   3) dbh_SBML	- A handle to an SBML database

sub transfer_SBML_to_FBA
{
	my ($model_id, $dbh_FBA, $dbh_SBML) = @_;

	my %SBMLgid_to_FBAgid;
	my %SBMLpid_to_FBApid;
	my %SBMLrid_to_FBArid;
	my %SBMLcid_to_FBAcid;


	##RETRIEVE AND INSERT METABOLITE DATA
	
        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_metabs = "SELECT metab_id, metab_symbol, metab_desc
                             FROM metabolite";

        
        my $sth = $dbh_SBML -> prepare($select_metabs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";
        
                         

        while(my ($metab_id, $metab_symbol, $metab_desc) =  $sth -> fetchrow_array)
        {                

		&metab_to_FBA_DB($metab_symbol, $metab_desc, $dbh_FBA);

        }#end while


	##RETRIEVE AND INSERT GENE DATA AND MAINTAIN GENE TO GENE ID MAPPING

	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_genes = "SELECT gid, ORF_id, name, alias, SGD_id, description, phenotype, gene_prod
                            FROM gene";


        my $sth = $dbh_SBML -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


	my $r_existing_genes = &get_all_genes_from_FBA_DB($dbh_FBA);
        my %existing_genes = %$r_existing_genes;

        while(my ($gid, $ORF_id, $name, $alias, $SGD_id, $description, $phenotype, $gene_prod) =  $sth -> fetchrow_array)
        {

		if( !exists($existing_genes{$ORF_id}) )
		{

	        	#PREPARE INSERT STATEMENT AND EXECUTE
	                my $insert_gene = "INSERT INTO gene (gid , ORF_id, name, alias, SGD_id, description, phenotype, gene_prod)
        	                           VALUES (gid , ?, ?, ?, ?, ?, ?, ?)";


                	my $sth = $dbh_FBA -> prepare($insert_gene)
                	                  or die "Can't prepare statment: $DBI::errstr";


                	my $rc = $sth -> execute($ORF_id, $name, $alias, $SGD_id, $description, $phenotype, $gene_prod)
                        	         or die "Can't execute statment: $DBI::errstr";



			my $new_gid = $dbh_FBA ->{ q{mysql_insertid}};
			$SBMLgid_to_FBAgid{$gid} = $new_gid;

		}else{

			$SBMLgid_to_FBAgid{$gid} = $existing_genes{$ORF_id};

		}#end if


        }#end while


	##RETIREVE AND INSERT PROTEIN DATA AND MAINTAIN PROTEIN TO PROTEIN ID MAPPING

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_prots = "SELECT pid, protein_name 
                            FROM protein";


        my $sth = $dbh_SBML -> prepare($select_prots)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


	my $r_existing_prots = &get_all_proteins_from_FBA_DB($dbh_FBA);
        my %existing_prots = %$r_existing_prots;

        while(my ($pid, $protein_name) =  $sth -> fetchrow_array)
        {

		if( !exists($existing_prots{$protein_name}))
		{

	                #PREPARE INSERT STATEMENT AND EXECUTE
        	        my $insert_prot = "INSERT INTO protein (pid , protein_name)
                	                   VALUES (pid , ?)";


                	my $sth = $dbh_FBA -> prepare($insert_prot)
                        	          or die "Can't prepare statment: $DBI::errstr";


                	my $rc = $sth -> execute($protein_name)
                        	         or die "Can't execute statment: $DBI::errstr";



                	my $new_pid = $dbh_FBA ->{ q{mysql_insertid}};
                	$SBMLpid_to_FBApid{$pid} = $new_pid;


		}else{

			$SBMLpid_to_FBApid{$pid} = $existing_prots{$protein_name};

		}#end if

        
        }#end while


	##RETIREVE AND INSERT RXN DATA AND MAINTAIN RXN TO RXN ID MAPPING

	#PREPARE SELECT STATEMENT AND EXECUTE
        my $select_rxns = "SELECT rid, rxn_code, rxn_name, rxn, subsystem, prot_class 
                           FROM rxn";


        my $sth = $dbh_SBML -> prepare($select_rxns)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($rid, $rxn_code, $rxn_name, $rxn, $subsystem, $prot_class) =  $sth -> fetchrow_array)
        {


	        #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_rxn = "INSERT INTO rxn (rid, rxn_code, rxn_name, rxn, subsystem, prot_class, model_id)
                                   VALUES (rid , ?, ?, ?, ?, ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_rxn)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($rxn_code, $rxn_name, $rxn, $subsystem, $prot_class, $model_id)
                                 or die "Can't execute statment: $DBI::errstr";



                my $new_rid = $dbh_FBA ->{ q{mysql_insertid}};
                $SBMLrid_to_FBArid{$rid} = $new_rid;


        }#end while


	##RETRIEVE AND INSERT CONDITION DATA
        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_conds = "SELECT con_id, condition_name, description
                            FROM `condition`";


        my $sth = $dbh_SBML -> prepare($select_conds)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        my $r_existing_conds = &get_conditions($dbh_FBA);
        my %existing_conds = %$r_existing_conds;

        while(my ($con_id, $condition_name, $desc) =  $sth -> fetchrow_array)
        {

                if( !exists($existing_conds{$condition_name}))
                {

                        #PREPARE INSERT STATEMENT AND EXECUTE
                        my $insert_cond = "INSERT INTO `condition` (con_id, condition_name, description)
                                           VALUES (con_id , ?, ?)";


                        my $sth = $dbh_FBA -> prepare($insert_cond)
                                          or die "Can't prepare statment: $DBI::errstr";


                        my $rc = $sth -> execute($condition_name, $desc)
                                         or die "Can't execute statment: $DBI::errstr";



                        my $new_cid = $dbh_FBA ->{ q{mysql_insertid}};
                        $SBMLcid_to_FBAcid{$con_id} = $new_cid;

		}else{


			$SBMLcid_to_FBAcid{$con_id} = $existing_conds{$condition_name};

                }#end if


        }#end while


	##RETRIEVE AND INSERT CONDITION BOUND DATA, BUT ADJUST CONDITION MAPPING
        
        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_con = "SELECT con_id, flux_UB, flux_LB 
                          FROM condition_bounds";
        
        
        my $sth = $dbh_SBML -> prepare($select_con)
                         or die "Can't prepare statment: $DBI::errstr";
                
        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($con_id, $flux_UB , $flux_LB) =  $sth -> fetchrow_array)
        {                             

                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_con = "INSERT INTO condition_bounds (con_id, flux_UB, flux_LB, model_id)
                                  VALUES ( ?, ?, ?, ?)";
        

                my $sth = $dbh_FBA -> prepare($insert_con)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($SBMLcid_to_FBAcid{$con_id}, $flux_UB, $flux_LB, $model_id)
                                 or die "Can't execute statment: $DBI::errstr";


        }#end while


	##RETRIEVE AND INSERT WILDTYPE FLUX DATA, BUT ADJUST CONDITION MAPPING

        #PREPARE SELECT STATEMENT AND EXECUTE
        #my $select_fluxes = "SELECT con_id, flux_vector, opt_type, status  
        #                     FROM wildtype_flux";
        my $select_fluxes = "SELECT con_id, flux_vector, status, opt_type 
                             FROM wildtype_flux";


        my $sth = $dbh_SBML -> prepare($select_fluxes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #while(my ($con_id, $flux_vector, $opt_type, $status) =  $sth -> fetchrow_array)
        while(my ($con_id, $flux_vector, $status, $opt_type) =  $sth -> fetchrow_array)
        {

		#my $opt_type = "FBA";

                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_fluxes = "INSERT INTO wildtype_flux (con_id, flux_vector, opt_type, status, model_id)
                                     VALUES ( ?, ?, ?, ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_fluxes)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($SBMLcid_to_FBAcid{$con_id}, $flux_vector, $opt_type, $status, $model_id)
                                 or die "Can't execute statment: $DBI::errstr";


        }#end while


        ##RETRIEVE AND INSERT SINGLE KO FLUX DATA, BUT ADJUST CONDITION AND GENE MAPPING

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_fluxes = "SELECT gid, con_id, flux_vector, opt_type, status  
                             FROM single_KO_flux";


        my $sth = $dbh_SBML -> prepare($select_fluxes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($gid, $con_id, $flux_vector, $opt_type, $status) =  $sth -> fetchrow_array)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_fluxes = "INSERT INTO single_KO_flux (gid, con_id, flux_vector, opt_type, status, model_id)
                                     VALUES ( ?, ?, ?, ?, ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_fluxes)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($SBMLgid_to_FBAgid{$gid} , $SBMLcid_to_FBAcid{$con_id}, $flux_vector, $opt_type, $status, $model_id)
                                 or die "Can't execute statment: $DBI::errstr";

        }#end while


	##RETRIEVE AND INSERT DOUBLE KO FLUX DATA, BUT ADJUST CONDITION AND GENE MAPPING

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_fluxes = "SELECT gid1, gid2, con_id, flux_vector, opt_type, status  
                             FROM double_KO_flux";


        #my $sth = $dbh_SBML -> prepare($select_fluxes)
        #                 or die "Can't prepare statment: $DBI::errstr";

        #my $rc = $sth -> execute()
        #                 or die "Can't execute statment: $DBI::errstr";


        #while(my ($gid1, $gid2, $con_id, $flux_vector, $opt_type, $status) =  $sth -> fetchrow_array)
        #{


                #PREPARE INSERT STATEMENT AND EXECUTE
        #        my $insert_fluxes = "INSERT INTO double_KO_flux (gid1, gid2, con_id, flux_vector, opt_type, status, model_id)
        #                             VALUES ( ?, ?, ?, ?, ?, ?, ?)";


        #        my $sth = $dbh_FBA -> prepare($insert_fluxes)
        #                              or die "Can't prepare statment: $DBI::errstr";


                #my $rc = $sth -> execute($SBMLgid_to_FBAgid{$gid1}, $SBMLgid_to_FBAgid{$gid2} , $SBMLcid_to_FBAcid{$con_id}, $flux_vector, $opt_type, $status, $model_id)
                #                 or die "Can't execute statment: $DBI::errstr";

        #}#end while


	##RETRIEVE AND INSERT GENE/PROTEIN, BUT ADJUST PROTEIN MAPPING 
	my $select_gid2pid = "SELECT gid, pid 
                              FROM gene2protein";


        my $sth = $dbh_SBML -> prepare($select_gid2pid)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($gid, $pid) =  $sth -> fetchrow_array)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_gid2pid = "INSERT INTO gene2protein (gid, pid, model_id)
                                     VALUES (?, ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_gid2pid)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($SBMLgid_to_FBAgid{$gid}, $SBMLpid_to_FBApid{$pid}, $model_id)
                                 or die "Can't execute statment: $DBI::errstr";

        }#end while


	##RETRIEVE AND INSERT RXN/COMPLEX, BUT ADJUST REACTION MAPPING 
        my $select_rxn2com = "SELECT rid, complex_expr
                              FROM rxn2complex_expr";
        

        my $sth = $dbh_SBML -> prepare($select_rxn2com)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($rid, $complex_expr) =  $sth -> fetchrow_array)
        {


                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_rxn2com = "INSERT INTO rxn2complex_expr (rid, complex_expr)
                                     VALUES ( ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_rxn2com)
                                      or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($SBMLrid_to_FBArid{$rid}, $complex_expr)
                                 or die "Can't execute statment: $DBI::errstr";

        }#end while


	##RETRIEVE AND INSERT RXN/PROTEIN, BUT ADJUST REACTION AND PROTEIN MAPPING
	my $select_rid2pid = "SELECT rid, pid 
                              FROM rxn2protein";


        my $sth = $dbh_SBML -> prepare($select_rid2pid)
                         or die "Can't prepare statment: $DBI::errstr";
                         
        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";
                         

        while(my ($rid, $pid) =  $sth -> fetchrow_array)
        {
        

                #PREPARE INSERT STATEMENT AND EXECUTE
                my $insert_rid2pid = "INSERT INTO rxn2protein (rid, pid)
                                     VALUES ( ?, ?)";


                my $sth = $dbh_FBA -> prepare($insert_rid2pid)
                                      or die "Can't prepare statment: $DBI::errstr";
                                      

                my $rc = $sth -> execute($SBMLrid_to_FBArid{$rid}, $SBMLpid_to_FBApid{$pid})
                                 or die "Can't execute statment: $DBI::errstr";

        }#end while
		
	
}#end transfer_SBML_to_FBA


################################################### FOR USE BY VISANT WEB APP ######################################

#sub get_wildtype_conditions - Takes as input a dataase handle and returns the names of all conditions for which fluxes exist 
#                               in the wildtype table
#
sub get_wildtype_conditions
{
        my ($dbh) = @_;
        my @conditions;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cond = "SELECT DISTINCT c.condition_name
                            FROM wildtype_flux wt , `condition` c
                            WHERE c.con_id = wt.con_id";


        my $sth = $dbh -> prepare($select_cond)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



         while(my $cond =  $sth -> fetchrow_array)
        {

                push @conditions , $cond;

        }#end while


        return \@conditions;


}#end get_wildtype_conditions


#sub get_singleKO_conditions - Takes as input a dataase handle and returns the names of all conditions for which fluxes exist 
#                              in the singleKO_flux table
#
sub get_singleKO_conditions
{
        my ($dbh) = @_;
        my @conditions;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_cond = "SELECT DISTINCT c.condition_name
                            FROM single_KO_flux ko , `condition` c
                            WHERE c.con_id = ko.con_id";


        my $sth = $dbh -> prepare($select_cond)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



         while(my $cond =  $sth -> fetchrow_array)
        {

                push @conditions , $cond;

        }#end while



        return \@conditions;


}#end get_singleKO_conditions


#sub get_opt_types - Takes as input a dataase handle and returns the names of all opt_types in the single KO table
#
sub get_opt_types
{
        my ($dbh) = @_;
        my @opt_types;


        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_opt = "SELECT DISTINCT opt_type
                            FROM single_KO_flux ko";


        my $sth = $dbh -> prepare($select_opt)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



         while(my $opt_type =  $sth -> fetchrow_array)
        {

                push @opt_types , $opt_type;

        }#end while



        return \@opt_types;


}#end get_opt_types



1;

