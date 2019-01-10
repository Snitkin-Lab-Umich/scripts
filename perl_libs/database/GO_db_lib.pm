package GO_db_lib;

use DBI;
use db_lib;


#sub DB_connect - Establishes a database connection and returns a handle
sub GO_DB_connect
{
        my $data_source = "dbi:mysql:cagt6.bu.edu:GO";
        my $user = "snitkines";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);
                

        return $dbh;

}#end DB_connect


#sub get_seqs - Takes as input an NCBI taxonomy ID and a database handle
#             - Returns a reference to an array of sequences in the desired organism
sub get_seqs
{
        my ($ncbi_tax_ID , $dbh) = @_;
	my %seqs;

	#PREPARE SELECT STATEMENT AND EXECUTE
	my $select_seqs = "SELECT gp.symbol , t.acc
			   FROM gene_product gp , term t, association ass, species sp
			   WHERE sp.ncbi_taxa_id = \'$ncbi_tax_ID\' AND
      			   sp.id = gp.species_id AND 
      			   gp.id = ass.gene_product_id AND 
      			   ass.term_ID = t.id AND
			   gp.type_id = 18864";
	

       my $sth = $dbh -> prepare($select_seqs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute
                         or die "Can't execute statment: $DBI::errstr";

	#GET SEQUENCE ID/GO ID PAIRS AND PUT IN RESULT HASH
	while(my ($seq_ID , $GO_ID) =  $sth -> fetchrow_array)
        {

		push @{$seqs{$seq_ID}} , $GO_ID;
	
	}#end while

        return \%seqs;

}#end get_seqs


#sub get_ancestors_dist - Takes as input the GO acc code of a query term and the term ID of a root
#		        - Returns a reference to two lists. One contains GO acc codes of ancestors and the
#			  other contains the anscestors distance from the chosen root. (Note that the lists
#			  are in the same order to ensure that th values are paired
sub get_ancestors_dist
{
	my ($child_symbol, $root_term_ID, $dbh) = @_;
	my (@GO_IDs , @distances);

        #PREPARE SELECT STATEMENT AND EXECUTE
	my $select_anc = "SELECT DISTINCT gp2.distance , t.acc
			  FROM graph_path gp1 , graph_path gp2 ,term t , gene_product gene , association ass
			  WHERE gene.symbol = \"$child_symbol\" AND
      				gene.id = ass.gene_product_id AND
      				ass.term_id = gp1.term2_id AND
      				t.id = gp1.term1_id AND	
      				gp2.term1_id = \'$root_term_ID\' AND
      				gp2.term2_id = t.id
			  ORDER BY gp2.distance";


       my $sth = $dbh -> prepare($select_anc)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute
                         or die "Can't execute statment: $DBI::errstr";

        #GET SEQUENCE ID/GO ID PAIRS AND PUT IN RESULT HASH
        while(my ($distance , $GO_ID) =  $sth -> fetchrow_array)
        {
                push @GO_IDs , $GO_ID;
		push @distances , $distance;

        }#end while


        return (\@GO_IDs, \@distances);

}#end get_ancestors_dist


#sub get_distance - Takes as input a GO_ID and the term ID of an anscestor
#                 - Returns teh distance between the two in the graph
sub get_distance
{
	my ($child_GO_ID , $anc_term_ID , $dbh) = @_;

        #PREPARE SELECT STATEMENT AND EXECUTE
	my $select_dist = "SELECT gp.distance
			   FROM term t, graph_path gp
			   WHERE t.acc = \'$child_GO_ID\' AND
      				 t.id = gp.term2_id AND
       				 gp.term1_id = $anc_term_ID";

#       print "$select_dist\n";

       my $sth = $dbh -> prepare($select_dist)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute
                         or die "Can't execute statment: $DBI::errstr";

	#GET THE DISTANCE
	$distance = $sth -> fetchrow_array;

	return $distance

}#end get_distance


#sub get_parents - Takes as input a GO ID
#                - Returns a reference to a,list containing the inputs nodes parents
sub get_parents
{
	my ($GO_ID , $dbh) = @_;

        #PREPARE SELECT STATEMENT AND EXECUTE
 	my $select_parent = "SELECT t2.acc
			     FROM term t1 , term t2, term2term tt
			     WHERE t1.acc= \'$GO_ID\' AND
		                   t1.id = tt.term2_id AND
				   t2.id = tt.term1_id";	

       my $sth = $dbh -> prepare($select_parent)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute
                         or die "Can't execute statment: $DBI::errstr";


        #GET ALL THE PARENTS
        while(my ($parent) =  $sth -> fetchrow_array)
        {

                push @parents , $parent;

        }#end while

	return \@parents;

}#end get_parents


#sub get_term_ID - Takes as input a GO acc and a DB handle
#		 - Returns the term ID corresponding to the inputted GO acc
sub get_term_ID
{
	my ($GO_acc , $dbh) = @_;

        #PREPARE SELECT STATEMENT AND EXECUTE
        my $select_term = "SELECT t.id
                             FROM term t
                             WHERE t.acc= ?";

       my $sth = $dbh -> prepare($select_term)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($GO_acc)
                         or die "Can't execute statment: $DBI::errstr";


        #GET THE TERM ID
        my $term_ID = $sth -> fetchrow_array;

        return $term_ID


}#end get_term_ID


#sub get_gene_annotation - Takes as inpot a DB handle and an ncbi_taxa_ID
#			 - Returns references to two lists. One is a list of gene_product
#			   IDs and the other is a list of the associated terms
sub get_gene_annotation
{
	my ($ncbi_taxID , $term_type, $r_gene_product_type, $r_evidence , $r_db_id, $dbh) = @_;
	my @gene_ids;
	my @term_ids;
	my %term_counts;
	my @gene_product_type = @$r_gene_product_type;
	my @evidence = @$r_evidence;
	my @db_id = @$r_db_id;


	#CONSTRUCT SELECT STATEMENT
	my $select_ass = "SELECT DISTINCT t.id , g.id
			  FROM term t, species s, gene_product g, evidence e , association a, db d
			  WHERE s.id = g.species_id AND
      				s.ncbi_taxa_id = ? AND
				a.is_not = 0 AND
      				g.type_id IN ("; 

	#ADD PLACEHOLDERS FOR GENE PRODUCT TYPES CODES
        $select_ass = &db_lib::add_place_holders($select_ass , scalar(@gene_product_type));


	$select_ass = $select_ass . " ) AND
      				g.id = a.gene_product_id AND
      				t.id = a.term_id AND
      				t.term_type = ? AND
      				e.association_id = a.id AND
      				e.code IN (";

	#ADD PLACEHOLDERS FOR EVIDENCE CODES
	$select_ass = &db_lib::add_place_holders($select_ass , scalar(@evidence));

	$select_ass = $select_ass . " ) AND
				d.id = a.source_db_id AND
				d.id IN (";

	#ADD PLACE HOLDERS FOR DB IDS
	$select_ass = &db_lib::add_place_holders($select_ass , scalar(@db_id));

	$select_ass = $select_ass . " )";


	#PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_ass)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ncbi_taxID , @gene_product_type , $term_type, @evidence, @db_id) 
                         or die "Can't execute statment: $DBI::errstr";


        #GET ALL THE TERM ID AND GENE PRODUCT ID PAIRS
        while(my ($term_id , $gene_id) =  $sth -> fetchrow_array)
        {

		#STORE GENE PRODUCT ID/TERM ID PAIRS
                push @term_ids , $term_id;
                push @gene_ids , $gene_id;

		#KEEP TRACK OF THE NUMBER OF TIMES EACH TERM IS USED
		if ( exists($term_counts{$term_id}) )
		{

			$term_counts{$term_id} ++;

		}else{

			$term_counts{$term_id} = 1;

		}#end if

        }#end while


	return (\@term_ids , \@gene_ids, \%term_counts);

}#end get_gene_annotation


#sub get_gene_annotation_org_spec - Takes as inpot a DB handle, a letter code for GO aspect, a refernce
#				    to a list of acceptable evidence codes, the name of the table holding
#				    the annotation for the desired genome 
#			 	  - Returns references to two lists. One is a list of gene IDs as given
#			   	    in the annotation table and the other is a list of the associated terms
#				    term IDs
sub get_gene_annotation_org_spec
{
	my ($aspect , $r_evidence , $table_name, $dbh) = @_;
	my @gene_ids;
	my @term_ids;
	my %term_counts;
	my @evidence = @$r_evidence;


	#CONSTRUCT SELECT STATEMENT
	my $select_ass = "SELECT DISTINCT t.id , o.DB_object_id
			  FROM $table_name o, term t
			  WHERE t.acc = o.GO_acc AND
      				o.GO_aspect = ? AND
				o.qualifier != 'NOT' AND
      				o.evidence_code IN (";

	#ADD PLACEHOLDERS FOR EVIDENCE CODES
	$select_ass = &db_lib::add_place_holders($select_ass , scalar(@evidence));

	$select_ass = $select_ass . " )"; 


	#PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_ass)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($aspect , @evidence) 
                         or die "Can't execute statment: $DBI::errstr";


        #GET ALL THE TERM ID AND GENE PRODUCT ID PAIRS
        while(my ($term_id , $gene_id) =  $sth -> fetchrow_array)
        {

		#STORE GENE PRODUCT ID/TERM ID PAIRS
                push @term_ids , $term_id;
                push @gene_ids , $gene_id;

		#KEEP TRACK OF THE NUMBER OF TIMES EACH TERM IS USED
		if ( exists($term_counts{$term_id}) )
		{

			$term_counts{$term_id} ++;

		}else{

                        $term_counts{$term_id} = 1;

                }#end if

        }#end while


        return (\@term_ids , \@gene_ids, \%term_counts);


}#end get_gene_annotation_org_spec


#sub get_term_counts - Takes as input an NCBI taxonomy ID and a DB handle
#                    - Returns a reference to a hash containing term ID as key and counts
#                      of the number of times each term occurs as the value. A term occurs
#		       when any of its children occur
sub get_term_counts
{
        my ($r_term_counts , $dbh) = @_;

	my %anc_term_counts;
	my %term_counts = %$r_term_counts;
	my @term_ids = keys %term_counts;


	#CONSTRUCT SELECT STATEMENT
	my $select_term_counts = "SELECT DISTINCT t1.id, t2.id 
				  FROM term t1, term t2, graph_path gp
				  WHERE t1.id = gp.term1_id AND
				        t1.term_type = 'biological_process' AND
			       	 	gp.term2_id = t2.id AND
			        	t2.id IN (";

	#ADD PLACE HOLDERS FOR TERMS
	$select_term_counts = &db_lib::add_place_holders($select_term_counts , scalar(@term_ids));

	$select_term_counts = $select_term_counts . ")";


        #PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_term_counts)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@term_ids)
                         or die "Can't execute statment: $DBI::errstr";


        #GET ALL THE TERM ID AND GENE PRODUCT ID PAIRS
        while(my ($term1_id , $term2_id) =  $sth -> fetchrow_array)
        {

		if( exists($anc_term_counts{$term1_id}) )
		{

			$anc_term_counts{$term1_id} += $term_counts{$term2_id};

		}else{


			$anc_term_counts{$term1_id} = $term_counts{$term2_id};

		}#end if

        }#end while


	return(\%anc_term_counts);

}#end get_term_counts


#sub get_ancestors - Takes as input a reference to a list of terms and a DB handle
#		   - Returns a reference to a list of terms which includes the terms on
#		     the original list and any ancestor of a term on the list
sub get_ancestors
{
	my ($r_term_ids , $term_type, $dbh) = @_;

	my @term_ids = @$r_term_ids;
	my @anc_term_ids;

        #CONSTRUCT SELECT STATEMENT
        my $select_anc = "SELECT DISTINCT t1.id
                          FROM term t1, term t2, graph_path gp
                          WHERE t1.id = gp.term1_id AND
                                t1.term_type = ? AND
                                gp.term2_id = t2.id AND
                                t2.id IN (";

        #ADD PLACE HOLDERS FOR TERMS
        $select_anc = &db_lib::add_place_holders($select_anc , scalar( @term_ids));

        $select_anc = $select_anc . ")";


        #PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_anc)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($term_type, @term_ids)
                         or die "Can't execute statment: $DBI::errstr";

        #GET ALL THE TERM ID AND GENE PRODUCT ID PAIRS
        while(my ($term1_id) =  $sth -> fetchrow_array)
        {

		push @anc_term_ids , $term1_id;

        }#end while

	
	return \@anc_term_ids;

}#end get_ancestors

#sub get_common_ancestor -Takes as input two term IDs and a DB handle
#			 - Returns a reference to a lsit containing the term IDs of all common
#			   ancestors of the two terms
sub get_common_ancestor
{
	my ($term_id1 , $term_id2, $term_type, $dbh) = @_;
	my @anc_term_ids;


        #CONSTRUCT SELECT STATEMENT
	my $select_common_anc = "SELECT DISTINCT gp1.term1_id  
			 	 FROM graph_path gp1 , graph_path gp2, term t1 , term t2
				 WHERE gp1.term1_id = gp2.term1_id AND
				       gp1.term1_id = t1.id AND
				       gp2.term1_id = t2.id AND
				       t1.term_type = ? AND
				       t2.term_type = ? AND
      				       gp1.term2_id = ? AND
      				       gp2.term2_id = ?";


        #PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_common_anc)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($term_type, $term_type, $term_id1 , $term_id2 )
                         or die "Can't execute statment: $DBI::errstr";


        #GET ALL THE TERM ID AND DISTANCE PAIRS
        while(my ($anc_term_id) =  $sth -> fetchrow_array)
        {

		push @anc_term_ids , $anc_term_id;

        }#end while


	return (\@anc_term_ids);

}#end get_common_ancestor


#sub get_gene_symbols - Takes as input a list of gene product IDs and a DB handle
#		      - Returns a reference to a hash linking gene product IDs to their
#			symbols
sub get_gene_symbols
{
	my ($r_gpids , $dbh) = @_;

	my @gpids = @$r_gpids;
	my %gp_hash;

	my $select_gp = "SELECT DISTINCT gp.id , gp.symbol
			 FROM gene_product gp
			 WHERE gp.id IN ( ";

        #ADD PLACE HOLDERS FOR GENE PRODUCT IDS
        $select_gp = &db_lib::add_place_holders($select_gp , scalar( @gpids));

        $select_gp = $select_gp . ")";


        #PREPARE AND EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($select_gp)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@gpids)
                         or die "Can't execute statment: $DBI::errstr";

        #GET ALL THE TERM ID AND GENE PRODUCT ID PAIRS
        while(my ($gpid , $symbol) =  $sth -> fetchrow_array)
        {

                $gp_hash{$gpid} = $symbol;

        }#end while


	return \%gp_hash;

}#end get_gene_symbols

1;
