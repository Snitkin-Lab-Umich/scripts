package biowulf_blast;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::Run::StandAloneNCBIBlast;


#sub plasmid_blast - Takes as input: 1) query_file	- A fasta formated file of sequences to BLAST
#							  against the plasmid database
#				     2) raw_blaat_file	- The name of a file to hold the BLAST resuts
#		   - Returns a reference to a 2D hash where the first dimension is the query name
#		     and the second dimnsion are properties of the query and its best hit in the
#		     plasmid database
sub plasmid_blast
{
	my ($query_file,$raw_blast_file) = @_;
	

	#BLAST CONTIGS AGAINST PLASMID DATABASE
	my@params = ( program  => 'blastn',
	              database => '/data/snitkines/Databases/Plasmid_DB/BLAST_DB/all_plasmids_genomes.fasta',
	              output   => $raw_blast_file);

	my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

	my $blast_report = $factory->blastall($query_file);


	#PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
	my %blast_results;

	while( my $result = $blast_report->next_result )
	{

	        my $hit = $result->next_hit;

	        my $hsp = $hit->next_hsp;

	        $blast_results{$result->query_name}{'query_length'} = $result->query_length;
	        $blast_results{$result->query_name}{'hit_name'} = $hit->name;
	        $blast_results{$result->query_name}{'hit_desc'} = $hit->description;        
		$blast_results{$result->query_name}{'hit_length'} = $hit->length;        
		$blast_results{$result->query_name}{'hsp_length'} = $hsp->length('total');
	        $blast_results{$result->query_name}{'hsp_eval'} = $hsp->evalue;
        	$blast_results{$result->query_name}{'hit_cov'} = &arith_ops::rnd($hsp->length('total') / $hit->length,2);
        	$blast_results{$result->query_name}{'query_cov'} = &arith_ops::rnd($hsp->length('total') / $result->query_length,2);

	}#end while


	return \%blast_results;

}#end plasmid_blast


#sub get_cog_assignments - Takes as input: 1) query_file      	- A fasta formated file of sequences to BLAST
#                                 	                        against the plasmid database
#                       	           2) raw_blaat_file	- The name of a file to hold the BLAST resuts
#			 - Returns a reference to a hash where the keys are protein IDs and values are lists of COGs
sub get_cog_assignments
{
        my ($query_file,$raw_blast_file, $r_prot2cog_hash, $r_prot2org_hash) = @_;

	my %COG_counts;
	my %COG_assignments;

	my $hit_cov_thresh = .5;
	my $min_COG_hits = 2;


	#GET COGS CORRESPONDING TO COG PROTEIN IDS
	my %prot2cog_hash = %$r_prot2cog_hash;

	my %prot2org_hash = %$r_prot2org_hash;


        #BLAST SEQUENCES AGAINST COG DATABASE
        my@params = ( program  => 'blastp',
                      database => '/data/snitkines/Databases/COG/BLAST_DB/cog',
		      e	       => 1,
		      b	       => 500);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS
        while( my $result = $blast_report->next_result )
        {

                my $query_name = $result->query_name;	
		my $query_length = $result->query_length;
		
		my %best_hit_spec_hash;
		

		#GO THROUGH EACH HIT FOR THE CURRENT QUERY
		while ( my $hit = $result->next_hit )
		{
	                	
			#GET INFO ABOUT HIT
			my $hit_name = $hit->name;
	                my $hit_length = $hit->length;
	                my $hit_eval = $hit->significance;

			my $org_id = $prot2org_hash{$hit_name};

			my $hit_cov = $hit->frac_aligned_hit;
			my $query_cov = $hit->frac_aligned_query;


			#GO THROUGH HSPs JUST AS A PRECAUTION (e.g. IF THERE ARE MULTIPLE HSPs WHOSE CUMULATIVE SIGNIFICANCE IS GREATER THAN SOME OTHER HIT)
			while( my $hsp = $hit->next_hsp )
			{

				#GET INFO ABOUT CURRENT HSP		
				my $num_hsps = $hit->num_hsps;
	                	my $hsp_length = $hsp->length('total');
	                	my $hsp_eval = $hsp->evalue;

				#print "$query_name, $hit_name , $num_hsps , $hit_cov, $hit_cov2, $query_cov, $query_cov2\n";


				#IF HIT MEETS CRITERIA (1) HIT ASSIGNED TO COG 2) ALIGNMENT COV 3) BEST HIT IN ORG), THEN MAKE NOT OF HIT
				if( (exists($prot2cog_hash{$hit_name})) && ($hit_cov > $hit_cov_thresh) && (!exists($best_hit_spec_hash{$org_id})) ) 
				{
	

					#KEEP TALLY OF NUMBER OF TIMES A GIVEN COG IS HIT AND THAT A BEST HIT TO THE CURRENT ORG HAS BEEN SEEN
					if( exists($COG_counts{$query_name}{$prot2cog_hash{$hit_name}}) )
					{

						$COG_counts{$query_name}{$prot2cog_hash{$hit_name}} ++;
						$best_hit_spec_hash{$org_id} = 1;

					}else{
	
						$COG_counts{$query_name}{$prot2cog_hash{$hit_name}} = 1;
						$best_hit_spec_hash{$org_id} = 1;

					}#end if

				}else{
	
					#print "$query_name, $hit_name , $prot2cog_hash{$hit_name} , $org_id , $hit_cov\n";


				}#end if


				#THROW FLAG IF HIT HAS MULTIPLE HSPs
				if( ($num_hsps) > 1 )
				{

					#print "$query_name, $hit_name, $num_hsps\n";

				}#end if


			}#end while

	
		}#end while


        }#end while

	
	#ASSIGN PROTEINS TO COGS IF THERE ARE ENOUGH HITS TO THE COG
	foreach my $protein ( keys %COG_counts )
	{
		
		#GET COGS COUNTS FOR CURRENT PROTEIN
		my %prot_COG_counts = %{$COG_counts{$protein}};

		foreach my $COG ( keys %prot_COG_counts)
		{

			#CHECK IF CURRENT COG HAS BEEN HIT ENOUGH TIMES FOR ASSIGNMENT
			if( $prot_COG_counts{$COG} >= $min_COG_hits )
			{

				if( exists($COG_assignments{$protein}) )
				{

					push @{$COG_assignments{$protein}}, $COG; 

				}else{

					@{$COG_assignments{$protein}} = ($COG);

				}#end if

				#print "$protein, $COG\n";
	
			}#end if


		}#end foreach

		#print "\n";

	}#end foreach

        return \%COG_assignments;

}#end get_cog_assignments


#sub nr_blast - Takes as input: 1) query_file 	- The name of a fasta file of sequences
#				  2) raw_blast_file	- The name of a file to hold the blast results
#		- Returns a hash where the keys are query sequence IDs and the values are 
#		  details of the top BLAST hit 
sub nr_blast
{
	my ($query_file, $raw_blast_file) = @_;

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastp',
                      -database => '/data/snitkines/Databases/nr/nr',
                      #-database => '/fdb/blastdb/nr',
		      -e	=> .001,
		      -v 	=> 500,
		      -b	=> 500);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                if(my $hit = $result->next_hit)
		{

                	$blast_results{$result->query_name}{'hit_name'} = $hit->name;
                	$blast_results{$result->query_name}{'hit_desc'} = $hit->description;
                	$blast_results{$result->query_name}{'eval'} = $hit->significance;
                        $blast_results{$result->query_name}{'frac_aligned_q'} = &arith_ops::rnd($hit->frac_aligned_query,2);
                        $blast_results{$result->query_name}{'frac_aligned_h'} = &arith_ops::rnd($hit->frac_aligned_hit,2);

		}#end if

        }#end while


        return \%blast_results;

}#end nrDB_blast


#sub kegg_blast - Takes as input: 1) query_file 	- The name of a fasta file of sequences
#				  2) raw_blast_file	- The name of a file to hold the blast results
#		- Returns a hash where the keys are query sequence IDs and the values are 
#		  details of the top BLAST hit 
sub kegg_blast
{
	my ($query_file, $raw_blast_file) = @_;

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastp',
                      -database => '/data/snitkines/Databases/KEGG/kegg_prok',
		      -e	=> .001,
		      -v 	=> 50000,
		      -b	=> 50000);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                if(my $hit = $result->next_hit)
		{

                	$blast_results{$result->query_name}{'hit_name'} = $hit->name;
                	$blast_results{$result->query_name}{'hit_desc'} = $hit->description;
                	$blast_results{$result->query_name}{'eval'} = $hit->significance;
                        $blast_results{$result->query_name}{'frac_aligned_q'} = &arith_ops::rnd($hit->frac_aligned_query,2);
                        $blast_results{$result->query_name}{'frac_aligned_h'} = &arith_ops::rnd($hit->frac_aligned_hit,2);

		}#end if

        }#end while


        return \%blast_results;

}#end kegg_blast


#sub kegg_blast_list - Takes as input: 1) query_file 	- The name of a fasta file of sequences
#				  2) raw_blast_file	- The name of a file to hold the blast results
#		     - Returns a hash where the keys are query sequence IDs and the values are 
#		       lists of the top BLAST hits in order of their singificance 
sub kegg_blast_list
{
	my ($query_file, $raw_blast_file) = @_;

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastp',
                      -database => '/data/snitkines/Databases/KEGG/kegg_prok',
		      -e	=> .001,
		      -v 	=> 50000,
		      -b	=> 50000);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;
        my %blast_results_list;

        while( my $result = $blast_report->next_result )
        {

		while(my $hit = $result->next_hit)
                {

                        my $hit_name = $hit->name;
                        my $hsp = $hit->next_hsp;

                        #my @hit_summary = ($hit->significance, $hit->query_length, $hit-> length, $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'));
                        my @hit_summary = ($hit->significance, $hit->query_length, $hit->length,  $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'), $hit->start('hit'), $hit->end('hit'), $hit->description);

                        $blast_results{$result->query_name}{$hit_name} = \@hit_summary;
                        push @{$blast_results_list{$result->query_name}}, join "\t", ($hit_name,@hit_summary);

                }#end while

        }#end while


        return \%blast_results_list;

}#end kegg_blast_list


#sub cdd_blast - Takes as input: 1) query_file	- The name of a fasta file of sequences
#                                2) raw_blast_file     - The name of a file to hold the blast results
#               - Performs rps blast against the cdd database and returns a hash where 
#		  the keys are query sequence IDs and the values are details of the top BLAST hit 
sub cdd_blast
{
	my ($query_file, $raw_blast_file) = @_;


	#BLAST PROTEINS AGAINST CDD DATABASE
        my @params = ( -p => 'T',
		      -e => .001,
		      database => '/data/snitkines/Databases/CDD/RPS_BLAST_DB/cog_pfam_smart_tigr_domains');

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->rpsblast($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

		my %temp_results;

                while( my $hit = $result->next_hit)
		{

			if( exists($temp_results{'hit_name'}) )
			{

				push @{$temp_results{'hit_name'}} ,$hit->name;
                                push @{$temp_results{'hit_desc'}} , $hit->description;
                                push @{$temp_results{'eval'}} , $hit->significance;
                                push @{$temp_results{'frac_aligned_q'}} , &arith_ops::rnd($hit->frac_aligned_query,2);
                                push @{$temp_results{'frac_aligned_h'}} , &arith_ops::rnd($hit->frac_aligned_hit,2);


			}else{

                		@{$temp_results{'hit_name'}} = ($hit->name);
                		@{$temp_results{'hit_desc'}} = ($hit->description);
                		@{$temp_results{'eval'}} = ($hit->significance);
                                @{$temp_results{'frac_aligned_q'}} = &arith_ops::rnd($hit->frac_aligned_query,2);
                                @{$temp_results{'frac_aligned_h'}} = &arith_ops::rnd($hit->frac_aligned_hit,2);


			}#end if

		}#end while


		if(exists($temp_results{'hit_name'}))
		{

			$blast_results{$result->query_name}{'hit_name'} = join " , ", @{$temp_results{'hit_name'}} ;
			$blast_results{$result->query_name}{'hit_desc'} = join " , ", @{$temp_results{'hit_desc'}} ;
			$blast_results{$result->query_name}{'eval'} = join " , ", @{$temp_results{'eval'}} ;
			$blast_results{$result->query_name}{'frac_aligned_q'} = join " , ", @{$temp_results{'frac_aligned_q'}} ;
			$blast_results{$result->query_name}{'frac_aligned_h'} = join " , ", @{$temp_results{'frac_aligned_h'}} ;

		}#end if

        }#end while 


        return \%blast_results;
	
}#end cdd_blast


#sub pc_blast - Takes as input: 1) query_file  - The name of a fasta file of sequences
#                               2) raw_blast_file     - The name of a file to hold the blast results
#				3) r_pc_info_hash 	- A hash linking protein cluster id to DESC and COG
#             - Performs rps blast against the protein clusters database and returns a hash where 
#               the keys are query sequence IDs and the values are details of the top BLAST hit 
sub pc_blast
{
        my ($query_file, $raw_blast_file, $r_pc_info_hash) = @_;


	#GET PROTEIN ClUSTER INFO
	my %pc_info_hash = %$r_pc_info_hash;


        #BLAST PROTEINS AGAINST CDD DATABASE
        my@params = ( -p => "T",
		      -e => 0.001,
                      -database => '/data/snitkines/Databases/protein_clusters/RPS_BLAST_DB/PRK_PC_DB');

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->rpsblast($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                if(my $hit = $result->next_hit)
		{

			#ADJUST PC NAME TO BE ALL UPPER CASE AND REMOVE UNWANTED CHARS
			my $hit_name = $hit->name;
			$hit_name = uc($hit_name);
			$hit_name =~ s/,//g;


                	$blast_results{$result->query_name}{'hit_name'} = $hit_name;
                	$blast_results{$result->query_name}{'hit_desc'} = $pc_info_hash{$hit_name}{'DESC'};
                	$blast_results{$result->query_name}{'hit_COG'} = $pc_info_hash{$hit_name}{'COG'};
                	$blast_results{$result->query_name}{'eval'} = $hit->significance;
                        $blast_results{$result->query_name}{'frac_aligned_q'} = &arith_ops::rnd($hit->frac_aligned_query,2);
                        $blast_results{$result->query_name}{'frac_aligned_h'} = &arith_ops::rnd($hit->frac_aligned_hit,2);

		}#end if

        }#end while 


        return \%blast_results;

}#end cdd_blast


#sub ardb_blast - Takes as input: 1) query_file         - The name of a fasta file of sequences
#                                 2) raw_blast_file     - The name of a file to hold the blast results
#               - Returns a hash where the keys are query sequence IDs and the values are 
#                 details of the top BLAST hit 
sub ardb_blast
{
        my ($query_file, $raw_blast_file) = @_;

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastp',
                      -database => '/data/snitkines/Databases/ardb/BLAST_DB/resisGenes.pfasta',
                      -e        => .001);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                if(my $hit = $result->next_hit)
                {

                        $blast_results{$result->query_name}{'hit_name'} = $hit->name;
                        $blast_results{$result->query_name}{'hit_desc'} = $hit->description;
                        $blast_results{$result->query_name}{'eval'} = $hit->significance;
                        $blast_results{$result->query_name}{'frac_aligned_q'} = &arith_ops::rnd($hit->frac_aligned_query,2);
                        $blast_results{$result->query_name}{'frac_aligned_h'} = &arith_ops::rnd($hit->frac_aligned_hit,2);

                }#end if

        }#end while


        return \%blast_results;


}#end ardb_blast



1;
