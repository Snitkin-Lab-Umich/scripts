package blast;

use gene_annotation_db_lib;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::Run::StandAloneNCBIBlast;

use arith_ops;


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
	              database => '/Users/snitkines/Documents/BLAST_databases/Plasmid_DB/BLAST_DB/all_plasmids_genomes.fasta',
	              output   => $raw_blast_file);

	my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

	my $blast_report = $factory->blastall($query_file);


	#PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
	my %blast_results;

	while( my $result = $blast_report->next_result )
	{

	        if( my $hit = $result->next_hit)
		{

	        	my $hsp = $hit->next_hsp;

	        	$blast_results{$result->query_name}{'hit_name'} = $hit->name;
	        	$blast_results{$result->query_name}{'hit_desc'} = $hit->description;        
	        	$blast_results{$result->query_name}{'hsp_eval'} = $hsp->evalue;
        		$blast_results{$result->query_name}{'hit_cov'} = &arith_ops::rnd($hit->frac_aligned_hit ,2);
        		$blast_results{$result->query_name}{'query_cov'} = &arith_ops::rnd($hit->frac_aligned_query,2);

		}#end if

	}#end while


	return \%blast_results;

}#end plasmid_blast


#sub get_cog_assignments - Takes as input: 1) query_file      	- A fasta formated file of sequences to BLAST
#                                 	                        against the plasmid database
#                       	           2) raw_blaat_file	- The name of a file to hold the BLAST resuts
#			 - Returns a reference to a hash where the keys are protein IDs and values are lists of COGs
sub get_cog_assignments
{
        my ($query_file,$raw_blast_file) = @_;

	my %COG_counts;
	my %COG_assignments;

	my $hit_cov_thresh = .5;
	my $min_COG_hits = 2;

	my $dbh = &gene_annotation_db_lib::gene_annotation_db_connect('dbi:mysql:gene_annotation');


	#GET COGS CORRESPONDING TO COG PROTEIN IDS
	my $r_prot2cog_hash = &gene_annotation_db_lib::get_prot2cog($dbh);
	my %prot2cog_hash = %$r_prot2cog_hash;

	my $r_prot2org_hash = &gene_annotation_db_lib::get_prot2org($dbh);
	my %prot2org_hash = %$r_prot2org_hash;


        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( program  => 'blastx',
                      database => '/Users/snitkines/Documents/BLAST_databases/COG/BLAST_DB/cog',
                      output   => $raw_blast_file,
		      e	       => 1);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        while( my $result = $blast_report->next_result )
        {

                my $query_name = $result->query_name;	
		my $query_length = $result->query_length;
		
		my %best_hit_spec_hash;
		

		while ( my $hit = $result->next_hit )
		{
	                	
			my $hit_name = $hit->name;
	                my $hit_length = $hit->length;
	                my $hit_eval = $hit->significance;

			my $org_id = $prot2org_hash{$hit_name};

			my $hit_cov = $hit->frac_aligned_hit;
			my $query_cov = $hit->frac_aligned_query;

			while( my $hsp = $hit->next_hsp )
			{

				#GET INFO ABOUT CURRENT HIT		
				my $num_hsps = $hit->num_hsps;
	                	my $hsp_length = $hsp->length('total');
	                	my $hsp_eval = $hsp->evalue;

	                	#my $hit_cov = &arith_ops::rnd( ($hsp_length / $hit_length) ,2);
        	        	#my $query_cov = &arith_ops::rnd( ($hsp_length / $query_length) ,2);
				
				#print "$query_name, $hit_name , $num_hsps , $hit_cov, $hit_cov2, $query_cov, $query_cov2\n";


				#IF HIT MEETGS CRITERIA THEN MAKE NOTE OF COG ASSIGNMENT (1) HIT ASSIGNED TO COG 2) ALIGNMENT COV 3) BEST HIT IN ORG)
				if( (exists($prot2cog_hash{$hit_name})) && ($hit_cov > $hit_cov_thresh) && (!exists($best_hit_spec_hash{$org_id})) ) 
				{
	
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

	
	#ASSIGN PROTEINS TO COGS IF THERE ARE ENOUGH HITS TO THE COGS
	foreach my $protein ( keys %COG_counts )
	{

		my %prot_COG_counts = %{$COG_counts{$protein}};

		foreach my $COG ( keys %prot_COG_counts)
		{

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


#sub nrDB_blast - Takes as input: 1) query_file - The name of a fasta file of sequences
#		- Returns a hash where the keys are query sequence IDs and the values are 
#		  details of the top BLAST hit 
sub nr_blast
{
	my ($query_file) = @_;

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastx',
                      -database => '/Users/snitkines/Documents/BLAST_databases/nrDB/BLAST_DB/nr',
		      -output   => 'test.blast',
		      -e	=> .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                my $hit = $result->next_hit;

                $blast_results{$result->query_name}{'hit_name'} = $hit->name;
                $blast_results{$result->query_name}{'hit_desc'} = $hit->description;
                $blast_results{$result->query_name}{'eval'} = $hit->significance;

        }#end while


        return \%blast_results;


}#end nrDB_blast


#sub cdd_blast - Takes as input: 1) query_file	- The name of a fasta file of sequences
#               - Performs rps blast against the cdd database and returns a hash where 
#		  the keys are query sequence IDs and the values are details of the top BLAST hit 
sub cdd_blast
{
	my ($query_file) = @_;


	#BLAST PROTEINS AGAINST CDD DATABASE
        my@params = ( -p => 'F',
		      -e => .001,
		      database => '/Users/snitkines/Documents/BLAST_databases/CDD/RPS_BLAST_DB/Pfam',
		      output   => 'test.blast');

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


			}else{

                		@{$temp_results{'hit_name'}} = ($hit->name);
                		@{$temp_results{'hit_desc'}} = ($hit->description);
                		@{$temp_results{'eval'}} = ($hit->significance);


			}#end if

		}#end while

		$blast_results{$result->query_name}{'hit_name'} = join " , ", @{$temp_results{'hit_name'}} ;
		$blast_results{$result->query_name}{'hit_desc'} = join " , ", @{$temp_results{'hit_desc'}} ;
		$blast_results{$result->query_name}{'eval'} = join " , ", @{$temp_results{'eval'}} ;

        }#end while 


        return \%blast_results;
	
}#end cdd_blast


#sub pc_blast - Takes as input: 1) query_file  - The name of a fasta file of sequences
#             - Performs rps blast against the protein clusters database and returns a hash where 
#               the keys are query sequence IDs and the values are details of the top BLAST hit 
sub pc_blast
{
        my ($query_file) = @_;

	#GET PROTEIN ClusTER INFO
	my $dbh = &gene_annotation_db_lib::gene_annotation_db_connect('dbi:mysql:gene_annotation');

	my $r_pc_info_hash = &gene_annotation_db_lib::get_pc_info($dbh);
	my %pc_info_hash = %$r_pc_info_hash;


        #BLAST PROTEINS AGAINST CDD DATABASE
        my@params = ( -p => "F",
                      -database => '/Users/snitkines/Documents/BLAST_databases/protein_clusters/RPS_BLAST_DB/PRK_PC_DB',
                      -output   => 'test.blast');

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->rpsblast($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                my $hit = $result->next_hit;

		#ADJUST PC NAME TO BE ALL UPPER CASE AND REMOVE UNWANTED CHARS
		my $hit_name = $hit->name;
		$hit_name = uc($hit_name);
		$hit_name =~ s/,//g;


                $blast_results{$result->query_name}{'hit_name'} = $hit_name;
                $blast_results{$result->query_name}{'hit_desc'} = $pc_info_hash{$hit_name}{'DESC'};
                $blast_results{$result->query_name}{'hit_COG'} = $pc_info_hash{$hit_name}{'COG'};
                $blast_results{$result->query_name}{'eval'} = $hit->significance;

        }#end while 


        return \%blast_results;

}#end cdd_blast

#sub blastp_vs_userDB - Takes as input: 1) blast_db            - Database to blast againast
#                                       2) $blast_output       - A file to hold blast output
#                                       3) r_genes             - A hash where the keys are genes and values are there sequences
#                      - Returns results of the blasting
sub blastp_vs_userDB
{

        my ($blast_db, $blast_output, $r_genes) = @_;

        my %genes = %$r_genes;
        my @seq_objs;


        #CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        foreach my $gene (keys %genes)
        {

                my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},
                                            -display_id => $gene,
                                            -alphabet => "dna" );

                push @seq_objs, $seq_obj;

        }#end foreach


        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastp',
                      -database => $blast_db,
                      -output   => $blast_output,
                      -e        => .01,
		      -v	=> 20,
		      -b	=> 20);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


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

                }#end if

        }#end while


        return (\%blast_results, \%blast_results_list);


}#end blastp_vs_userDB


#sub blastn_vs_userDB - Takes as input: 1) blast_db            - Database to blast againast
#                                       2) $blast_output       - A file to hold blast output
#                                       3) r_genes             - A hash where the keys are genes and values are there sequences
#                      - Returns results of the blasting
sub blastn_vs_userDB
{

        my ($blast_db, $blast_output, $r_genes) = @_;

        my %genes = %$r_genes;
        my @seq_objs;


        #CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        foreach my $gene (keys %genes)
        {

                my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},
                                            -display_id => $gene,
                                            -alphabet => "dna" );

                push @seq_objs, $seq_obj;

        }#end foreach

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastn',
                      -database => $blast_db,
                      -output   => $blast_output,
#		      -F		=> "F",
#                      -e        => 1000,
#		      -X	=> 10000,
#		      -y	=> 10000,
#		      -Z	=> 10000,
                      -e        => 1e-10,
#		      -W	=> 7,
		      -b	=> 500);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                while(my $hit = $result->next_hit)
                {

                        my $hit_name = $hit->name;
                        my $hsp = $hit->next_hsp;
                        
			my @hit_summary = ($hit->significance, $hit->query_length, $hit->length,  $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'), $hit->start('hit'), $hit->end('hit'), $hit->raw_score, $hit->description);

                        $blast_results{$result->query_name}{$hit_name} = \@hit_summary;

                }#end if

        }#end while

        return \%blast_results;

}#end blastn_vs_userDB


#sub tblastx_vs_userDB - Takes as input: 1) blast_db		- Database to blast againast
#					 2) $blast_output	- A file to hold blast output
#					 3) r_genes		- A hash where the keys are genes and values are there sequences
#		       - Returns results of the blasting
sub tblastx_vs_userDB
{

	my ($blast_db, $blast_output, $r_genes) = @_;

	my %genes = %$r_genes;
	my @seq_objs;
	

	#CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
	foreach my $gene (keys %genes)
	{

		my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},                        
                          		    -display_id => $gene,                        
                          		    -alphabet => "dna" );

		push @seq_objs, $seq_obj;

	}#end foreach


	#BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'tblastx',
                      -database => $blast_db,
                      -output   => $blast_output,
                      -e        => .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                while(my $hit = $result->next_hit)
		{

			my $hit_name = $hit->name;
			my $hsp = $hit->next_hsp;

			#my @hit_summary = ($hit->significance, $hit->query_length, $hit-> length, $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'));
			my @hit_summary = ($hit->significance, $hit->query_length, $hit->length,  $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'), $hit->start('hit'), $hit->end('hit'), $hit->raw_score, $hit->description);

                	$blast_results{$result->query_name}{$hit_name} = \@hit_summary;

		}#end if

        }#end while


        return \%blast_results;


}#end tblastx_vs_userDB

#sub blastx_vs_userDB - Takes as input: 1) blast_db		- Database to blast againast
#					 2) $blast_output	- A file to hold blast output
#					 3) r_genes		- A hash where the keys are genes and values are there sequences
#		       - Returns results of the blasting
sub blastx_vs_userDB
{

	my ($blast_db, $blast_output, $r_genes) = @_;

	my %genes = %$r_genes;
	my @seq_objs;
	

	#CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
	foreach my $gene (keys %genes)
	{

		my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},                        
                          		    -display_id => $gene,                        
                          		    -alphabet => "dna" );

		push @seq_objs, $seq_obj;

	}#end foreach


	#BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastx',
                      -database => $blast_db,
                      -output   => $blast_output,
                      -e        => .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                while(my $hit = $result->next_hit)
		{

			my $hit_name = $hit->name;
			my $hsp = $hit->next_hsp;

			my @hit_summary = ($hit->significance, $hit->query_length, $hit->length,  $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'), $hit->start('hit'), $hit->end('hit'), $hit->raw_score, $hit->description);

                	$blast_results{$result->query_name}{$hit_name} = \@hit_summary;

		}#end if

        }#end while


        return \%blast_results;


}#end blastx_vs_userDB


#sub find_long_hits - Takes as input: 1) blast_db	- A database of proteins
#				      2) r_genes	- A hash where the keys are ids and values are coding sequence
#				      3) length_cutoff	- Minimum fraction of length of best hit to be flaged
#		    - Returns a hash where the keys are genes with long hits or short hits and the values are the genes hit
sub find_long_hits
{
	my ($blast_db, $r_genes, $length_cutoff, $identity_cutoff) = @_;

        my %genes = %$r_genes;
        my @seq_objs;
	my %gene_call_errors;


        #CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        foreach my $gene (keys %genes)
        {

                my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},
                                            -display_id => $gene,
                                            -alphabet => "dna" );

                push @seq_objs, $seq_obj;

        }#end foreach


        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastx',
		      -output	=> 'test.long.blast',
                      -database => $blast_db,
                      -e        => .01,
		      -v	=> 250);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %query2hits;
	my %hit2queries;
	my @hit_genes = ();


        while( my $result = $blast_report->next_result )
	{

		my %strain_hits;

		my $q_name = $result->query_name;
		my $q_length = ($result->query_length / 3);

                while(my $hit = $result->next_hit)
                {

			my $h_name = $hit->name;
			my $h_length = $hit->length;

			my $hsp = $hit->next_hsp;
			my $h_frame = ($hsp->query->frame + 1) * $hsp->strand('query');

			my $frac_id = $hsp->frac_identical;


			#SKIP HITS IN SAME GENOME
			my $q_strain = $q_name;
			$q_strain =~ s/^(.*):.*$/$1/;

			my $h_strain = $h_name;
			$h_strain =~ s/^(.*):.*$/$1/;

			if($q_strain eq $h_strain )
			{

				next;

			}#end if


			#CHECK FRAME OF QUERY, IF NOT +1 THEN FLAG AS POTENTIAL BAD CALL
			if($h_frame != 1)
			{

				$gene_call_errors{$q_name} = $h_name;
				$strain_hits{$h_strain} = 1;
				next;

			}#end if


			#SEE IF BEST HIT AGAINST GIVEN STRAIN AND IF SO CHECK LENGTH
			if( !exists($strain_hits{$h_strain}) )
			{

				
				if( ((($q_length / $h_length) < $length_cutoff) || (($h_length / $q_length) < $length_cutoff)) && ($frac_id > $identity_cutoff) )
				{

					if( exists($hit2queries{$h_name}) )
					{
				
						push @{$hit2queries{$h_name}}, $q_name;

					}else{

						@{$hit2queries{$h_name}} = ($q_name);

					}#end if


					if( !exists($query2hits{$q_name}) )
					{

						@{$query2hits{$q_name}} = ($h_name);
						push @hit_genes, $h_name;

					}else{

						push @{$query2hits{$q_name}}, $h_name;
						push @hit_genes, $h_name;

					}#end if

				}#end if

			}#end if


			#NOTE THAT PROTEIN FROM HIT STRAIN HAS BEEN SEEN
			$strain_hits{$h_strain} = 1;


		}#end while

	}#end while



	return (\%query2hits, \%hit2queries, \@hit_genes, \%gene_call_errors);


}#end find_long_hits


#sub find_frsameshifts - Takes as input: 1) blast_db		- Database to blast againast which contains the genome in which you are 
#								  looking for frameshifts
#					 2) r_genes		- A hash where the keys are genes and values are there sequences
#					 3) identity_cutoff	- A percent identity cutoff to consider an HSP
#					 4) gene2candidate	- A hash linking full length genes to the shorter genes in which frame shifts
#								 are being looked for 
#		       - BLASTs genes relative to which a long/short hit was found against the genome where the long/short gene occurs and goes
#			 through HSPs to determine if cause is 1) 5' truncation, 2) 3' truncation, 3) frame-shift and/or 4) contig break
sub find_frameshifts
{

	my ($blast_db, $r_genes, $identity_cutoff, $r_gene2candidate, $dbh) = @_;

	my %genes = %$r_genes;

	my %gene2candidate = %$r_gene2candidate;

	my @seq_objs;
	
	my $r_contig_bordering_genes = &genomes_db_lib::get_genes_bordering_contigs($dbh);
	my %contig_bordering_genes = %$r_contig_bordering_genes;


	#CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
	foreach my $gene (keys %genes)
	{

		my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},                        
                          		    -display_id => $gene,                        
                          		    -alphabet => "dna" );

		push @seq_objs, $seq_obj;

	}#end foreach



	#BLAST AGAINST DATABASE CONTAINING GENOME IN WHICH FRAMESHIFTS ARE BEING LOOKED FOR
        my@params = ( -program  => 'tblastx',
                      -database => $blast_db,
		      -output	=> 'test.frame.out',
                      -e        => .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {
		
		#GET THE COORDINATES AND THE STRAIN OF THE GENE WITH A POTENTIAL FRAMESHIFT
		my $q_name = $result->query_name;


		#FOREACH GENE PROBED BY THIS QUERY CHECK FOR FRAMESHIFTS
		my $r_candidates = $gene2candidate{$q_name};
		

		foreach my $candidate( @$r_candidates)
		{

			#GET COORDINATES OF CANDIDATAE GENE
			my @candidate = ($candidate);
			my $r_gene_coords_hash = &genomes_db_lib::get_gene_coords(\@candidate , $dbh);
			my %gene_coords_hash = %$r_gene_coords_hash;

			my $gene_start;
			my $gene_end;
			my $gene_start_5;
			my $gene_end_5;

			if($gene_coords_hash{$candidate}{'START'} < $gene_coords_hash{$candidate}{'END'})
			{

				$gene_start = $gene_coords_hash{$candidate}{'START'}; 
				$gene_end = $gene_coords_hash{$candidate}{'END'}; 

			}else{

				$gene_start = -1 * $gene_coords_hash{$candidate}{'START'}; 
				$gene_end = -1 * $gene_coords_hash{$candidate}{'END'}; 

			}#end if


			my $strain = $candidate;
			$strain =~ s/^(.*):.*$/$1/;


			#GO THROUGH HSPS AND LOOK FOR EVIDENCE OF FRAME SHIFT OR PREMATURE STOP CODON
                	while(my $hit = $result->next_hit)
			{

				my %hsps;
				my $h_name = $hit->name;


				#CHECK IF HIT IS AGAINST GENOME OF INTEREST
				if( ($h_name !~ /$strain /) && ($h_name !~ /$strain\:/) && ($h_name !~ /$strain$/) && ($h_name !~ /$strain\_/) )
				{

					next;

				}#end if


				while( my $hsp = $hit->next_hsp)
				{

					#GET INFO FOR ALL HSPS WHERE THE HIT IS ON THE CORRECT STRAND/FRAME (+/1), THE HSP PASSES IDENTITY CUTOFF
					#AND THE HSP IS WITHIN 10kb OF THE GENES COORDS
					my $hsp_start;
					my $hsp_end;
					my $frac_id = $hsp->frac_identical;
					my $frame = $hsp->hit->frame * $hsp->strand('hit');

					#ADJUST START AND END DEPENDING ON THE STRAND
					if( $hsp->strand('hit') == 1)
					{
					
						
						$hsp_start = $hsp->start('hit');
						$hsp_end = $hsp->end('hit');

					}else{

						$hsp_start = -1* $hsp->end('hit');
						$hsp_end =-1 *  $hsp->start('hit');

					}#end if


					#DETERMINE IF CURRENT HSP SHOULD BE CONSIDERED
					if( ($hsp->strand('query') == 1) && ($hsp->query->frame == 0) && ($frac_id > $identity_cutoff) && 
					   (  &arith_ops::abs($gene_start - $hsp_start) <  10000) && ( &arith_ops::abs($gene_end - $hsp_end) <  10000 ) )
					{


						$hsps{$hsp_start}{'FRAME'} = $hsp->hit->frame; 
						$hsps{$hsp_start}{'STRAND'} = $hsp->strand('hit'); 
						$hsps{$hsp_start}{'END'} = $hsp_end; 

					}else{

						#print "\t$candidate ($q_name), frac($frac_id), geneS($gene_start_5), geneE($gene_end_5)" . $hsp->strand('query') . ", " . $hsp->query->frame . ", " . $hsp->start('hit') . ", " . $hsp->end('hit') . "\n"; 	

					}#end if				


				}#end while

			
				#DETERMINE WHETHER HSPS INDICATE THAT FRAMESHIFT HAS OCCURED
				my $prev_frame = -1;
				my $prev_end = -1;

				my $first_start = 0;
				my $last_end  = 0;				

				my $frame_shift = 0;
				my $trunc3 = 0;
				my $ext3 = 0;
				my $trunc5 = 0;
				my $ext5 = 0;
				my $contig_border = 0;
				my $bad_gene_call = 0;
			
				#print "$candidate ($q_name): ";	

				foreach	my $start (sort {$a <=> $b} keys %hsps)
				{


					#DETERMINE IF FRAMESHIFT HAS OCCURED
					my $curr_frame = $hsps{$start}{'FRAME'};
					my $curr_strand = $hsps{$start}{'STRAND'};
					#print "($gene_start, $gene_end, " . $start . ", " . $hsps{$start}{'END'} . ", $curr_frame, $curr_strand) , ";				

					if( ($prev_frame != -1) && ($curr_frame != $prev_frame) )
					{

						$frame_shift = 1;

					}#end if


					$prev_frame = $curr_frame;
					$prev_end = $hsps{$start}{'END'};


					#DETERMINE RANGE OF HIT
					if( $first_start == 0 || $start < $first_start)
					{

						$first_start = $start;

					}#end if


					if( ($hsps{$start}{'END'}  > $last_end) || ($last_end == 0) )
					{
	
						$last_end  = $hsps{$start}{'END'};

					}#end if
			

				}#end foreach		

				
				#IF THERE IS NO FRAMESHIFT, DETERMINE IF THERE IS A 5' TRUNCATION/EXTENSION OR A 3' TRUNCATION/EXTENSION
				if(($gene_end < $last_end) && scalar(keys %hsps) > 0)
				{

					$trunc3 = 1;
		
				}#end if

				if(($gene_start > $first_start) && scalar(keys %hsps) > 0){

	
					$trunc5 = 1;

				}#end if

				if(($gene_end > $last_end) && scalar(keys %hsps) > 0){

					$ext3 = 1;	


				}#end if

				if(($gene_start < $first_start) && scalar(keys %hsps) > 0){


					$ext5 = 1;

				}#end if


				#DETERMINE IF GENE BORDERS CONTIG BOUNDARY
				if( exists($contig_bordering_genes{$candidate}) )
				{

					$contig_border = 1;

				}#end if
		

				#STORE RESULTS FOR CURRENT HIT
				my @blast_results = (\%hsps, $frame_shift, $trunc5, $trunc3, $ext5, $ext3,  $contig_border,  $first_start, $last_end, $gene_start, $gene_end);
				$blast_results{$candidate}{$q_name} = \@blast_results;

			}#end while

			#GO THROUGH THE RESULT OF THIS QUERY AGAIN FOR THE NEXT CANDIDATE USING THIS QUERY
			$result->rewind;			


		}#end foreach


        }#end while

        return \%blast_results;


}#end find_frameshifts


#sub get_coords - Takes as input: 1) fasta_file	- A faasta file of nucleotide sequences
#				  2) genome_db	- A database holding a genome
#		- Returns a hash where the keys are the sequence name and START/END and the values 
#		  coordinates of the sequence in the genome
sub get_coords
{
	my ($query_file , $genome_db) = @_;

	my %coords_hash;


        #BLAST PROTEINS AGAINST CDD DATABASE
        my@params = ( -program  => 'blastn',
                      -database => $genome_db,
                      -output   => 'test.blast',
		      -F 	=> 'F',
		      -X	=> 500);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall($query_file);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

		if(my $hit = $result->next_hit)
		{


			my $hsp = $hit->next_hsp;

			#MAKE SURE HSP IS COVERS ENTIRE QUERY
			if( $hsp->hsp_length() == $result->query_length )
			#if( $hsp->hsp_length() > ($result->query_length-10) )
			{

				my $query_name = $result->query_name;

				#CHECK IF HIT IS ON MINUS OR PLUS STRAND
				if($hsp->strand('sbjct') == 1)
				{

					$coords_hash{$query_name}{'START'} = $hsp->start('sbjct');
					$coords_hash{$query_name}{'END'} = $hsp->end('sbjct');
					$coords_hash{$query_name}{'CHR'} = $hit->name;
	
				}else{

					$coords_hash{$query_name}{'START'} = $hsp->end('sbjct');
					$coords_hash{$query_name}{'END'} = $hsp->start('sbjct');
					$coords_hash{$query_name}{'CHR'} = $hit->name;

				}#end if

			}#end if

		}#end if

        }#end while 


	return \%coords_hash;

}#end get_coords


#sub primer_blast - Takes as input: 1) blast_db			- A database holding a genome sequence
#				    2) left_primer_seq		- The sequence of the left primer
#				    3) right_primer_seq		- The sequence of the right primer
#		  - Returns a hash where the key is the amplicon number and the value is a list (left coord, right coord, size)
sub primer_blast
{

	my ($blast_db, $left_primer_seq, $right_primer_seq) = @_;


	#CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        my $left_seq_obj = Bio::Seq->new(-seq => $left_primer_seq,
                      	                 -display_id => "left",
                                         -alphabet => "dna" );


        my $right_seq_obj = Bio::Seq->new(-seq => $right_primer_seq,
                      	                  -display_id => "right",
                                          -alphabet => "dna" );


        my @seq_objs = ($left_seq_obj, $right_seq_obj);


        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastn',
                      -database => $blast_db,
#                     -F                => "F",
                      -e        => .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);

        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

		my $hit_num = 0;

                while(my $hit = $result->next_hit)
                {

                        while(my $hsp = $hit->next_hsp)
			{
				
				$hit_num ++;

                        	my @hit_summary = ($hsp->start('sbjct'), $hsp->end('sbjct'), $hsp->frac_identical, $hsp->length('query'), $hsp->strand('hit'));

                        	$blast_results{$result->query_name}{$hit_num} = \@hit_summary;

			}#end if

                }#end if

        }#end while


	#GO THROUGH EACH PAIR OF PRIMER PRIMER PAIR HITS AND DETERMINE SIZE OF AMPLICON
	my %primer_hits;

	foreach my $left_hit_i (keys %{$blast_results{'left'}})
	{

		my @left_hit = @{$blast_results{'left'}{$left_hit_i}};

		foreach my $right_hit_i(keys %{$blast_results{'right'}})
		{

			my @right_hit = @{$blast_results{'right'}{$right_hit_i}};

			#CHECK THAT: 1) BOTH ANNEAL SUCH THAT PRODUCT WILL BE PRODUCED
			#	     2) PRODUCT SIZE IS < 2000 bp
			#	     3) % ID FOR BOTH IS > 90%
			if(($right_hit[4] == 1) && ($left_hit[4] == -1) && (($left_hit[1] - $right_hit[0]) < 2000) && (($left_hit[1] - $right_hit[0]) > 0) &&
			   ($right_hit[2] > 0.9) && ($left_hit[2] > 0.9))
			{

				$primer_hits{$left_hit[0]} = $left_hit[1] - $right_hit[0] + 1;		

			}elsif(($left_hit[4] == 1) && ($right_hit[4]== -1) && (($right_hit[1] - $left_hit[0]) < 2000) && (($right_hit[1] - $left_hit[0]) > 0) &&
				($right_hit[2] > 0.9) && ($left_hit[2] > 0.9)){

				$primer_hits{$right_hit[0]} = $right_hit[1] - $left_hit[0] + 1;		

			}#end if

			#print "left= @left_hit, right = @right_hit\n";

		}#end foreach


	}#end foreach	
	#print "\n\n";

        return \%primer_hits;

}#end primer_blast


#sub get_top_blast_hit - Takes as input: 1) blast_db            - Database to blast againast
#                                        2) $blast_output       - A file to hold blast output
#                                        3) r_genes             - A hash where the keys are genes and values are there sequences
#                      - Returns results of the blasting
sub get_top_blast_hit
{

        my ($blast_db, $blast_output, $r_genes) = @_;

        my %genes = %$r_genes;
        my @seq_objs;


        #CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        foreach my $gene (keys %genes)
        {

                my $seq_obj = Bio::Seq->new(-seq => $genes{$gene},
                                            -display_id => $gene,
                                            -alphabet => "dna" );

                push @seq_objs, $seq_obj;

        }#end foreach

        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastn',
                      -database => $blast_db,
                      -output   => $blast_output,
#		      -F		=> "F",
#                      -e        => .01,
                      -e        => 1e-10,
		      -b	=> 500);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

                while(my $hit = $result->next_hit)
                {

                        my $hit_name = $hit->name;
                        my $hsp = $hit->next_hsp;
                        
                        $blast_results{$result->query_name}{$hit_name} = $hsp->frac_identical;

			last;

                }#end while


        }#end while

        return \%blast_results;

}#end get_top_blast_hit


#sub get_top_genome_hit - Takes as input: 1) blast_db            - Database to blast againast (genome database)
#                                         2) blast_output       - A file to hold blast output
#                                         3) r_sequences         - A hash where the keys are sequence IDs and the values 
#								   are the sequences 
#					  4) host_genome	- A string describing the host genome, to be used in a regular
#							          expression to filter out self hits
#					  5) num_hits		- The number of hits to the host/other genomes to save
#                        - Returns results of the top hit against the host genome and the top hit against another genome
sub get_top_genome_hit
{

        my ($blast_db, $blast_output, $r_sequences, $host_genome, $num_hits) = @_;

        my %sequences = %$r_sequences;
        my @seq_objs;


        #CREATE LIST OF OBJECTS TO BLAST AGAINST DATABASE
        foreach my $sequence (keys %sequences)
        {

                my $seq_obj = Bio::Seq->new(-seq => $sequences{$sequence},
                                            -display_id => $sequence,
                                            -alphabet => "dna" );

                push @seq_objs, $seq_obj;

        }#end foreach


        #BLAST CONTIGS AGAINST PLASMID DATABASE
        my@params = ( -program  => 'blastn',
                      -database => $blast_db,
                      -output   => $blast_output,
                      -e        => .01);

        my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

        my $blast_report = $factory->blastall(\@seq_objs);


        #PARSE BLAST RESULTS TO DETERMINE BEST ALIGNMENT
        my %blast_results;

        while( my $result = $blast_report->next_result )
        {

		my $host_genome_hit = 0;
		my $other_genome_hit = 0;

                while(my $hit = $result->next_hit)
                {

                        my $hit_name = $hit->name . "(" . $hit->description . ")";


			#MAKE SURE THERE IS AN HSP
			my @hit_summary;

                        if(my $hsp = $hit->next_hsp)
                        {

				@hit_summary = ($hit_name, $hit->significance, $hit->query_length, $hit->length,  $hit->frac_aligned_hit, $hit->frac_aligned_query, $hsp->start('sbjct'), $hsp->end('sbjct'), $hit->start('hit'), $hit->end('hit'));

			}else{

				next;

			}#end if


			#SAVE HIT IF IT IS THE BEST IN THE HOST GENOME OR THE BEST IN ANOTHER GENOME
			if($hit_name =~ /$host_genome/i && $host_genome_hit < $num_hits)
			{

	                        $blast_results{$result->query_name}{"HOST$host_genome_hit"} = \@hit_summary;
				$host_genome_hit ++;
				next;

			}elsif($hit_name !~ /$host_genome/i && $other_genome_hit <$num_hits){

				
	                        $blast_results{$result->query_name}{"OTHER$other_genome_hit"} = \@hit_summary;
				$other_genome_hit ++;
				next;

			}#end if


			#GO TO NEXT RESULT IF HIT TO HOST AND OTHER HAVE BEEN SEEN
			if($other_genome_hit == $num_hits && $host_genome_hit == $num_hits)
			{

				last;

			}#end if


                }#end if

        }#end while


        return \%blast_results;

}#end get_top_genome_hit


#sub blastn_check - Takes as input: 1) blastn_results	- A hash of blastn results
#				    2) in_seqs		- A list of sequence names in which a blast hit should be found
#				    3) out_seqs		- A list of sequence names in which a blast hit should not be found
#				    4) min_frac_aln	- An optional parameter indicating the fraction alignment to be considered
#							  a hit (default = 0.5)
#
#		  - Takes blastn results and confirms that blast hits match presence and absence as inputed.
sub blastn_check
{
	my ($r_blastn_results, $r_in_seqs, $r_out_seqs, $min_frac_aln) = @_;

	my @errors;
	my @fatal_errors;

	my %blastn_results = %$r_blastn_results;
	my @in_seqs = @$r_in_seqs;
	my @out_seqs = @$r_out_seqs;
	my @all_seqs = (@$r_in_seqs, @$r_out_seqs);
	

	my $r_out_seq_blast_hit_hash = &hash_lib::init_hash($r_out_seqs);
        my %out_seq_blast_hit_hash = %$r_out_seq_blast_hit_hash;

        my $r_in_seq_blast_hit_hash = &hash_lib::init_hash($r_in_seqs);
        my %in_seq_blast_hit_hash = %$r_in_seq_blast_hit_hash;

        my $r_all_seq_blast_hit_hash = &hash_lib::init_hash(\@all_seqs);
        my %all_seq_blast_hit_hash = %$r_all_seq_blast_hit_hash;

	
	#CHECK IF FRACTION ALIGNED IS GIVEN, IF NOT THEN SET TO DEFAULT AT 0.5
	if(!$min_frac_aln)
	{

		$min_frac_aln = 0.5;

	}#end if


        #VERIFY PRESENCE AND ABSENCE OF HITS AGAINST CORRECT SEQUENCES
        my @fracs = (0);
        my %blast_hit_hash;
	my %all_seqs_seen;


        foreach my $hit ( keys %blastn_results)
        {


		#NTOE THAT SEQUENCE HAS BEEN SEEN AND PARSE OUT HIT
		$all_seqs_seen{$hit} = 1;

                my ($significance, $query_length, $hit_length,  $frac_aligned_hit, $frac_aligned_query, $sbj_start, $sbj_end, $hit_start, $hit_end, $score, $desc) = @{$blastn_results{$hit}};

                $blast_hit_hash{$hit}{'FRAC_ALN_Q'} = $frac_aligned_query;
                $blast_hit_hash{$hit}{'SIG'} = $significance;

                #IF PRESENT IN MSA, THEN MAKE SURE THERE IS BLAST HIT OF REASONABLE QUALITY
                if($significance < 0.01 && $frac_aligned_query > 0.5 && exists($in_seq_blast_hit_hash{$hit}) )
                {

                        $in_seq_blast_hit_hash{$hit} = 1;
                        $all_seq_blast_hit_hash{$hit} = 1;

                }#end if


                #IF ABSENT IN MSA, AND THERE IS GOOD ALIGNMENT, THEN THROW OUT REGION
                if($significance < 0.01 && $frac_aligned_query > 0.5 && exists($out_seq_blast_hit_hash{$hit}) )
                {


                        push @errors,  "*******Sequence does exist in $hit (E = $significance and fraction aligned = $frac_aligned_query)!*************";

                        $all_seq_blast_hit_hash{$hit} = 1;
                        push @fracs, $frac_aligned_query;

                }#end if

        }#end foreach     


	#REPORT LACK OF HIT IF ANY SEQEUCNE THAT WAS SUPPOSED TO
        foreach my $in_seq (@in_seqs)
        {

                if($in_seq_blast_hit_hash{$in_seq} == 0)
                {

                        push @fatal_errors, print "**********************No significant BLAST hit in $in_seq!******************";

                }#end if

        }#end foreach


        #REPORT WHICH SEQUENECES WERE HIT AND WHICH WERE NOT
        my @hit_seqs;
        my @no_hit_seqs;

        foreach my $seq (keys %all_seqs_seen)
        {

                if(exists($blast_hit_hash{$seq}) && $blast_hit_hash{$seq}{'FRAC_ALN_Q'} > 0.5)
                {

                        push @hit_seqs, "$seq ($blast_hit_hash{$seq}{'FRAC_ALN_Q'})";

                }else{

                        push @no_hit_seqs, "$seq ($blast_hit_hash{$seq}{'FRAC_ALN_Q'})";

                }#end if

        }#end foreach


	my $max_frac = &list_ops::max(\@fracs);

	return (\@hit_seqs, \@no_hit_seqs, \@errors, \@fatal_errors, $max_frac, \%blast_hit_hash);

}#end 


#get_kegg_gene - Takes as input: gene_hash	- A hash linking gene id to protein sequence
#	       - Returns a list of gene names and annotations
#
sub get_kegg_gene
{
	my ($r_gene_hash) = @_;


	#SET UP KEGG VARIABLES
	my $kegg_db = "/Users/snitkines/Documents/BLAST_databases/KEGG/genes_gammaProt";

	my $kegg_code2name = "/Users/snitkines/Documents/BLAST_databases/KEGG/taxonomy_code2name";
	my $r_org_code2full = &hash_lib::dict_to_hash($kegg_code2name);
	my %org_code2full = %$r_org_code2full;


	#BLAST GENES AGAINST KEGG DATABASE
	($r_kegg_blast_results, $r_kegg_blast_results_list) = &blastp_vs_userDB($kegg_db, 'temp.blastOutput_' . rand(5), $r_gene_hash);
	my %kegg_blast_results_list = %$r_kegg_blast_results_list;

	my %orgs;
	my %gene_names;
	my %gene_symbols;

	foreach my $gene_name (keys %$r_gene_hash)
	{


	        #DECLARE SOME VARIABLES
                my @orgs;
                my @descs;
                my $rep_gene_symbol;
                my $rep_gene_symbol_flag = 0;

                my @result_list = @{$kegg_blast_results_list{$gene_name}};


                #GO THROUGH TOP N BLAST HITS, PARSE AND REPORT
                for($i = 0; $i < &arith_ops::min(5,scalar(@result_list)); $i++)
                {

                        my ($hit_name, $significance, $query_length, $hit_length,  $frac_aligned_hit, $frac_aligned_query, $q_start, $q_end, $hit_start, $hit_end, $hit_desc) = split /\t/, $result_list[$i];

                        #GET THE ORG NAME (KEGG 3-LETTER CODE), KEEP TRACK SO KEY CAN BE PRINTED AT END
                        my $org = $hit_name;
                        $org =~ s/^([^:]+):.*$/$1/;


                        #REMOVE KO INFO FROM THE HIT DESCRIPTION
                        $hit_desc =~ s/^(.+?); K\d+.+$/$1/;
                        $hit_desc =~ s/^([^\)]+)\(.*$/$1/;
                        my $gene_symbol = $hit_desc;


                        #PARSE OUT GENE SYMBOL (IF EXISTS) AND GENE DESCRIPTION
                        if($hit_desc =~ /;/)
                        {

                                $hit_desc =~ s/^([^;]+); ([^;]+)$/$2/;
                                $gene_symbol =~ s/^([^;]+); ([^;]+)$/$1/;

                                #GET GENE SYMBOL FROM THE MOST SIGNIFICANT OF THE TOP N GENES
                                if($rep_gene_symbol_flag == 0)
                                {

                                        $rep_gene_symbol = "$gene_symbol/$significance ($hit_desc)";
                                        $rep_gene_symbol_flag = 1;

                                }#end if

                        }else{

                                $gene_symbol = "N\\A";

                                #SET GENE SYMBOL AS MOST SIGNIFICANT DESCRIPTION, IN CASE NONE OF TOP N HAS A SYMBOL
                                if($i == 0)
                                {

                                        $rep_gene_symbol = "N\\A/$significance ($hit_desc)";

                                }#end if

                        }#end if


                        #ONLY SAVE THE TOP 3 FOR OUTPUT
                        if($i < 5)
                        {

                                push @orgs, "$org ($gene_symbol/$significance)";
                                push @descs, "$gene_symbol/$significance ($hit_desc)";

                        }#end if

                }#end for

                $orgs{$gene_name} = \@orgs;
                $genes{$gene_name} = \@descs;
		$gene_symbols{$gene_name} = $rep_gene_symbol;

        }#end foreach


	return (\%gene_symbols, %genes, %orgs);


}#end get_kegg_gene

1;
