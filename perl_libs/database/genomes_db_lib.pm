package genomes_db_lib;

use DBI;
use db_lib;

use Bio::Seq;

use gene_annotation_db_lib;
use arith_ops;
use seq_stats;
use seq_ops;
use seq_io;
use hash_lib;
use mummer;

use Bio::SimpleAlign;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities;
use Bio::AlignIO;
use my_Utilities;
use blast;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::CodonTable;


###################################################### DB CONNECTION #########################################################

#sub genome_DB_connect - Establishes a database connection and returns a handle
sub genome_DB_connect
{
	my ($db_name) = @_;

        my $data_source = "dbi:mysql:$db_name;host=biobase";

        my $user = "snitkines";
        my $password ="ruffian";

        #CONNECT TO DB
        my $dbh = &db_lib::DB_connect($data_source, $user, $password);


        return $dbh;

}#end genome_DB_connect



#################################################### INSERTION ROUTINES ######################################################

#sub insert_newbler_organism - Takes as input: 1) r_newbler_hash - A hash whose keys are N50ContigSize, Q39MinusBases, Q40PlusBases	
#							   	  avgContigSize, largestContigSize, numberOfBases and numberOfContigs.
#		     	      - Makes insertion into organism table
sub insert_newbler_organism
{
	my ($r_newbler_hash, $species, $strain,  $dbh) = @_;

	my %newbler_hash = %$r_newbler_hash;

	#IF NEWBLER HASH EXISTS THEN INSERT ALL METRICS, OTHERWISE JUST INSERT BASICS
	if($r_newbler_hash != 0)
	{

		#PREPARE INSERT STATEMENT
	        my $insert_org = "INSERT INTO organism (Q39_minus, Q40_plus, N50_contig, num_contigs, avg_contig_size, largest_contig_size,
							N50_scaffold, num_scaffolds, avg_scaffold_size, largest_scaffold_size, species, strain )
	                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";


       		#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_org)
        	                 or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($newbler_hash{'Q39MinusBases'},$newbler_hash{'Q40PlusBases'}, $newbler_hash{'N50ContigSize'}, 
					$newbler_hash{'numberOfContigs'}, $newbler_hash{'avgContigSize'}, $newbler_hash{'largestContigSize'},
					$newbler_hash{'N50ScaffoldSize'},$newbler_hash{'numberOfScaffolds'}, $newbler_hash{'avgScaffoldSize'}, 
					$newbler_hash{'largestScaffoldSize'}, $species, $strain )
               	         or die "Can't execute statment: $DBI::errstr";

	}else{

		#PREPARE INSERT STATEMENT
	        my $insert_org = "INSERT INTO organism(species, strain)
	                          VALUES (?, ?)";


       		#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_org)
        	                 or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($species, $strain )
                	         or die "Can't execute statment: $DBI::errstr";

	}#end if


}#end insert_newbler_organism


#sub insert_A5_organism - Takes as input: 1) r_a5_hash - A hash whose keys are column names in A5 assembly stats file
#		     	- Makes insertion into organism table
sub insert_A5_organism
{
	my ($r_a5_hash, $species, $strain,  $dbh) = @_;

	my %a5_hash = %$r_a5_hash;

	#PREPARE INSERT STATEMENT
        my $insert_org = "INSERT INTO organism (Q40_plus, N50_contig, num_contigs,
						num_scaffolds, largest_scaffold_size, species, strain, num_reads )
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?)";


     	#EXECUTE THE STATEMENT
       	my $sth = $dbh -> prepare($insert_org)
       	                 or die "Can't prepare statment: $DBI::errstr";


       	my $rc = $sth -> execute($a5_hash{'bases >= Q40'}, $a5_hash{'N50'}, 
				$a5_hash{'Contigs'}, $a5_hash{'Scaffolds'},  
				$a5_hash{'Longest Scaffold'}, $species, $strain, $a5_hash{'Raw reads'})
              	         or die "Can't execute statment: $DBI::errstr";

}#end insert_a5_organism

#sub insert_mira_organism - Takes as input: 1) r_mira_hash - A hash whose keys are from mira  assembly stats file
#		     	  - Makes insertion into organism table
sub insert_mira_organism
{
	my ($r_mira_hash, $species, $strain,  $dbh) = @_;

	my %mira_hash = %$r_mira_hash;

	#PREPARE INSERT STATEMENT
        my $insert_org = "INSERT INTO organism (N50_contig, num_contigs,
						largest_contig_size, species, strain, num_reads )
                          VALUES (?, ?, ?, ?, ?, ?)";


     	#EXECUTE THE STATEMENT
       	my $sth = $dbh -> prepare($insert_org)
       	                 or die "Can't prepare statment: $DBI::errstr";


       	my $rc = $sth -> execute($mira_hash{'N50 contig size'}, $mira_hash{'Number of contigs'}, 
				$mira_hash{'Largest contig'},  
				$species, $strain, $mira_hash{'Num. reads assembled'} )
              	         or die "Can't execute statment: $DBI::errstr";

}#end insert_mira_organism


#sub insert_contig - Takes as input: 1) r_contig_hash		- A 2D hash whose keys are contig_id and QUALITY_ALL/DEPTH_ALL
#				     2) r_contig_coord_hash	- A hash whose keys are START and END
#				     3) r_contig_desc_hash	- A hash whose key is the contig and value is a description
#		   - Makes insertion into contig table
sub insert_contig
{
	my ($r_contig_hash, $r_contig_coord_hash, $r_contig_desc_hash,  $dbh) = @_;

	my %contig_hash = %$r_contig_hash;
	my %contig_coord_hash = %$r_contig_coord_hash;
	my %contig_desc_hash = %$r_contig_desc_hash;


	foreach my $contig (sort keys %contig_coord_hash)
	{

		#GET THE ID OF THE CHROMOSOME THE CONTIG IS ON
		my $chr_id = &get_chr_id($contig_coord_hash{$contig}{'CHR'}, $dbh);

		#PREPARE INSERT STATEMENT
        	my $insert_contig = "INSERT INTO contig (`contig_name`, `chr_id`, `start`, `end`, `quality`, `depth`, `description`)
        	                     VALUES (?,?, ?, ?, ?, ?, ?)";


        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_contig)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($contig, $chr_id, $contig_coord_hash{$contig}{'START'}, $contig_coord_hash{$contig}{'END'}, 
					$contig_hash{$contig}{'QUALITY_ALL'}, $contig_hash{$contig}{'DEPTH_ALL'}, 
					$contig_desc_hash{$contig})
        	                 or die "Can't execute statment: $DBI::errstr";

	}#end foreach

}#end insert_contig


#sub insert_contig_wQual - Takes as input: 1) r_quality_hash		- A hash where the keys are contig names and the values are space delimited quality scores
#					   2) r_contig_coord_hash	- A hash whose keys are START and END
#					   3) r_contig_desc_hash	- A hash whose key is the contig and value is a description
#		   - Makes insertion into contig table
sub insert_contig_wQual
{
	my ($r_quality_hash, $r_contig_coord_hash, $r_contig_desc_hash,  $dbh) = @_;

	my %contig_coord_hash = %$r_contig_coord_hash;
	my %contig_desc_hash = %$r_contig_desc_hash;
	my %quality_hash = %$r_quality_hash;


	foreach my $contig (sort keys %contig_coord_hash)
	{

		#GET THE ID OF THE CHROMOSOME THE CONTIG IS ON
		my $chr_id = &get_chr_id($contig_coord_hash{$contig}{'CHR'}, $dbh);

		#PREPARE INSERT STATEMENT
        	my $insert_contig = "INSERT INTO contig (`contig_name`, `chr_id`, `start`, `end`, `quality`, `description`)
        	                     VALUES (?, ?, ?, ?, ?, ?)";


        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_contig)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($contig, $chr_id, $contig_coord_hash{$contig}{'START'}, $contig_coord_hash{$contig}{'END'}, 
					$quality_hash{$contig}, $contig_desc_hash{$contig})
        	                 or die "Can't execute statment: $DBI::errstr";

	}#end foreach

}#end insert_contig_wQual

#sub insert_contig_noMetrics - Takes as input: 1) r_contig_coord_hash	- A hash whose keys are START and END
#			      	     	       2) r_contig_desc_hash	- A hash whose key is the contig and value is a description
#		   	     - Makes insertion into contig table
sub insert_contig_noMetrics
{
	my ($r_contig_coord_hash, $r_contig_desc_hash,  $dbh) = @_;

	my %contig_coord_hash = %$r_contig_coord_hash;
	my %contig_desc_hash = %$r_contig_desc_hash;


	foreach my $contig (sort keys %contig_coord_hash)
	{

		print "$contig\n";
		#GET THE ID OF THE CHROMOSOME THE CONTIG IS ON
		my $chr_id = &get_chr_id($contig_coord_hash{$contig}{'CHR'}, $dbh);

		#PREPARE INSERT STATEMENT
        	my $insert_contig = "INSERT INTO contig (contig_name, chr_id, start, end, description)
        	                     VALUES (?, ?, ?, ?, ?)";


        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_contig)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($contig, $chr_id, $contig_coord_hash{$contig}{'START'}, $contig_coord_hash{$contig}{'END'},  $contig_desc_hash{$contig})
        	                 or die "Can't execute statment: $DBI::errstr";

	}#end foreach

}#end insert_contig_noMetrics


#sub insert_chromosome - Takes as input: 1) seq		- The sequence of the chromosome
#					 2) strain	- An strain in the organism table
#					 3) seq_class	- The class of the sequence (e.g. plasmid or chromosome)
#			- Makes insertion into the chromosome table
sub insert_chromosome
{
	my ($seq, $strain, $seq_class, $seq_name, $dbh) = @_;

	my $org_id = "";


	#GET THE ORGANISM ID CORRESPONDING TO THE STRAIN
	my $select_org_id = "SELECT org_id
                             FROM organism
                             WHERE strain = ?";


        my $sth = $dbh -> prepare($select_org_id)
                          or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($strain)
                          or die "Can't execute statment: $DBI::errstr";

        $org_id =  $sth -> fetchrow_array;


	if($org_id == "")
	{

		print "\n\nStrain not in organism table, can't insert into chromosme!\n\n";

		exit;

	}#end if


	#PREPARE INSERT STATEMENT
        my $insert_chrom = "INSERT INTO chromosome (org_id, class, sequence, name)
                            VALUES (?, ?, ?, ?)";


        #EXECUTE THE STATEMENT
        $sth = $dbh -> prepare($insert_chrom)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($org_id, $seq_class, $seq, $seq_name) 
                         or die "Can't execute statment: $DBI::errstr";


}#end insert_chromosome


#sub insert_genes - Takes as input: 1) r_gene_hash	- A hash where the keys are START, END, CLASS
#				    2) chr_id		- An ID from the chromosome table
#		  - Makes insertion into gene table
sub insert_genes
{
	my ($r_gene_hash, $chr_id, $dbh) = @_;

	my %gene_hash = %$r_gene_hash;


	#IF COORDINATES ARE RELATIVE TO CONTIG STARTS, THEN GET CONTIG POSITIONS ON CHROMOSOME
	my %contig2start;
	my %contig2end;

	
        #GET CONTIGS ON CURRENT CHROMSOME
        my $select_contigs = "SELECT contig_name, start, end
                              FROM contig
                              WHERE chr_id = ?";


        my $sth = $dbh -> prepare($select_contigs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($contig_name, $start, $end) =  $sth -> fetchrow_array)
        {

		$contig2start{$contig_name} = $start;
		$contig2end{$contig_name} = $end;


	}#end while


	#GO THROUGH HASH AND MAKE INSERTIONS
	foreach my $gene (sort keys %gene_hash)
	{

		#ADJUST START IF COORDINATES ARE RELATIVE TO CONTIG START
		if( $gene_hash{$gene}{'CONTIG'} !~ /^$/ && exists($contig2start{$gene_hash{$gene}{'CONTIG'}}))
		{

			if($contig2start{$gene_hash{$gene}{'CONTIG'}} < $contig2end{$gene_hash{$gene}{'CONTIG'}})
			{

				$gene_hash{$gene}{'START'} += $contig2start{$gene_hash{$gene}{'CONTIG'}} - 1;
				$gene_hash{$gene}{'END'} += $contig2start{$gene_hash{$gene}{'CONTIG'}} - 1;

			}else{

				$gene_hash{$gene}{'START'} = $contig2start{$gene_hash{$gene}{'CONTIG'}} - $gene_hash{$gene}{'START'} + 1;
				$gene_hash{$gene}{'END'} = $contig2start{$gene_hash{$gene}{'CONTIG'}} - $gene_hash{$gene}{'END'} + 1;

			}#end if


		}elsif($gene_hash{$gene}{'CONTIG'} !~ /^$/){

			print "Contig ($gene_hash{$gene}{'CONTIG'}) not in table!!!!";

		}else{


			print "Contig ($gene_hash{$gene}{'CONTIG'}) does not exist for $gene!!!!";

		}#end if


		#PREPARE INSERT STATEMENT
        	my $insert_gene = "INSERT INTO gene (name, chr_id, start, end, type)
        	                   VALUES (?, ?, ?, ?, ?)";

        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_gene)
        	                 or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($gene, $chr_id, $gene_hash{$gene}{'START'}, $gene_hash{$gene}{'END'}, $gene_hash{$gene}{'TYPE'}) 
                	         or die "Can't execute statment: $DBI::errstr";


	}#end foreach


}#end insert_genes


#sub insert_annotations - Takes as input: 1) annotation_file - A file containing annotations one per line 
#							       (Annotation class | gene name | Hit Name | Annotation | E value)
#					  2) strain	     - The name of the strain for which annotations are being inserted
sub insert_annotations
{
	my ($annotation_file, $strain, $dbh) = @_;

	print "\n\n($annotation_file, $strain, $dbh)\n\n";
	my $genome_name = $strain . "_genome";


	#GET HASH LINKING GENE NAMES TO IDS
	my $r_gene_name2id = &get_all_gene_ids($dbh);
	my %gene_name2id = %$r_gene_name2id;


	#GO THROUGH THE FILE MAKING INSERTIONS
	open INPUT, $annotation_file;


	foreach my $line (<INPUT>)
	{

		chomp $line;

		my ($source, $gene_name, $hit, $annotation, $e_value, $frac_aligned_q, $frac_aligned_h) = split /\t/, $line;

		#IF GENE NAME DOES NOT EXIST, THEN MAKE GENOME NAME PREFIX
		if(!exists($gene_name2id{$gene_name}))
		{
	
			$gene_name = "$genome_name:$gene_name";
			print "$gene_name\n";
 
		}else{

			print "$gene_name\n";

		}#endif


		#INSERT ANNOTATION
        	my $insert_annot = "INSERT INTO annotation (gid, annotation, source, hit, e_value, frac_aligned_query ,frac_aligned_hit)
        	                   VALUES (?, ?, ?, ?, ?,?,?)";


        	#EXECUTE THE STATEMENT
        	$sth = $dbh -> prepare($insert_annot)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	$rc = $sth -> execute($gene_name2id{$gene_name}, $annotation, $source, $hit, $e_value, $frac_aligned_q, $frac_aligned_h) 
                	         or die "Can't execute statment: $DBI::errstr";


	}#end foreach


	close INPUT;

}#end insert_annotations


#sub insert_ortholog_groups - Takes as input: 1) r_ortho_hash - A hash whose keys are GENES< ORGS, NUM_GENES, NUM_ORGS
#			    - Makes insertion into the ortholog_group table
sub insert_ortholog_groups
{
	my ($r_ortho_hash, $dbh) = @_;

	my %ortho_hash = %$r_ortho_hash;


	#INSERT EACH ORTHOLOG GROUP
	foreach my $ortho_id (keys %ortho_hash)
	{

		#INSERT ANNOTATION
        	my $insert_og = "INSERT INTO ortholog_group (gene_names, org_ids, num_orgs, num_genes)
        	                 VALUES (?, ?, ?, ?)";


        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_og)
        	                 or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($ortho_hash{$ortho_id}{'GENES'}, $ortho_hash{$ortho_id}{'ORGS'}, 
				         $ortho_hash{$ortho_id}{'NUM_ORGS'}, $ortho_hash{$ortho_id}{'NUM_GENES'}) 
                	         or die "Can't execute statment: $DBI::errstr";


	}#end foreach

}#end insert_ortholog_groups


#sub insert_og_msa - Takes as input: 1) msa	- A multiple sequence alignment of an orthologous group
#				     2) og_id	- An orthologous group id
#                  - Makes insertion into the ortholog_group table
sub insert_og_msa
{
        my ($msa, $og_id, $dbh) = @_;

	my $update_og = "UPDATE ortholog_group
                         SET msa = ? 
                         WHERE og_id = ?";

        #EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($update_og)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa, $og_id)
                         or die "Can't execute statment: $DBI::errstr";

	

}#end insert_og_msa


#sub insert_pseudo_gene - Takes as input: 1) pseudo_gene_file	- A 3 column flle (gid og_id mutation)
#			- Makes insertion into the pseudo_gene table
sub insert_pseudo_gene
{
	my ($pseudo_gene_file, $dbh) = @_;


	open PSEUDO, $pseudo_gene_file;

	foreach my $line (<PSEUDO>)
	{

		chomp $line;
		my ($gid, $og_id, $mut_type, $position) = split /\t/, $line;

		my $insert_pg = "INSERT INTO pseudogene (gid, og_id, mut_type, position)
                                 VALUES (?, ?, ?, ?)";


                #EXECUTE THE STATEMENT
                my $sth = $dbh -> prepare($insert_pg)
                                 or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($gid, $og_id, $mut_type, $position)
                                 or die "Can't execute statment: $DBI::errstr";


	}#end foreach


	close PSEUDO;

}#end insert_pseudo_gene


#sub insert_nucmer - Takes as input: 1) coords_file	- An output file from show-coords
##				     2) $chr
sub insert_nucmer
{
	my ($coords_file, $dbh) = @_;


	#PARSE COORDS FILE
	my $r_coords_hash = &mummer::parse_coords_file($coords_file);
	my %coords_hash = %$r_coords_hash;


	#GET CHROMSOME IDS
	my $r_chr_ids = &get_chr_id_all($dbh);
	my %chr_ids = %$r_chr_ids;

	#PARSE FILE AND MAKE INSERTIONS INTO DATABASE
	foreach my $aln_id (sort {$a<=>$b} keys %coords_hash)
	{

		#INSERT ALIGNMENT
        	my $insert_os = "INSERT INTO pairwise_alignment (chr_id1, start1, end1, chr_id2, start2, end2)
                       		 VALUES (?, ?, ?, ?, ?, ?)";


        	#EXECUTE THE STATEMENT
        	my $sth = $dbh -> prepare($insert_os)
                       		 or die "Can't prepare statment: $DBI::errstr";


        	my $rc = $sth -> execute($chr_ids{$coords_hash{$aln_id}{'SEQ1'}}, $coords_hash{$aln_id}{'S1'}, $coords_hash{$aln_id}{'E1'}, $chr_ids{$coords_hash{$aln_id}{'SEQ2'}}, $coords_hash{$aln_id}{'S2'}, $coords_hash{$aln_id}{'E2'})
                        		or die "Can't execute statment: $DBI::errstr";



	}#end foreach


}#end insert_nucmer


#sub insert_mauve - Takes as input: 1) mauve_file 	- The .backbone file from rogressive mauve 3.xx
#				    2) r_chr_ids	- A reference to a list of ids from the chromosome table, in the same
#							  order as they appear in  the columns of the mauve_file	 
#		  - Makes insertion into the pairwise_alignment DB
sub insert_mauve
{
	my ($mauve_file, $r_chr_ids, $dbh) = @_;

	my @chr_ids = @$r_chr_ids;


	#PARSE FILE AND MAKE INSERTIONS INTO DATABASE
	open INPUT, $mauve_file;

	<INPUT>;

	foreach my $line (<INPUT>)
	{

		chomp $line;

		my (@coords) = split /\t/, $line;


		while ( @coords > 2)
		{

			my $start1 = shift @coords;
			my $end1 = shift @coords;
			my $chr1 = $chr_ids[ (@chr_ids - (@coords/2) - 1) ];

			for($i = 0;$i < @coords;$i+=2)
			{

				my $start2 = $coords[$i];
				my $end2 = $coords[$i + 1];
				my $chr2 = $chr_ids[ (@chr_ids - (@coords/2) + ($i/2) ) ];


				#INSERT ANNOTATION
        			my $insert_os = "INSERT INTO pairwise_alignment (chr_id1, start1, end1, chr_id2, start2, end2)
        	                		 VALUES (?, ?, ?, ?, ?, ?)";


        			#EXECUTE THE STATEMENT
        			my $sth = $dbh -> prepare($insert_os)
        	                		 or die "Can't prepare statment: $DBI::errstr";


        			my $rc = $sth -> execute($chr1, $start1, $end1, $chr2, $start2, $end2)
                	         		or die "Can't execute statment: $DBI::errstr";


			}#end for


		}#end while


	}#end foreach


}#end insert_mauve



#sub insert_mauve_MSA - Takes as input: 1) mauve_backbone_file	- A backbone file produced by Mauve
#					2) chr_names		- A list of chromosome names corresponding to the columns in the alignment
#					3) aln_desc		- A description of the MSA (e.g. method, sequences, etc.)
#		      - Inserts Mauve MSA into the multiple alignment table
sub insert_mauve_MSA
{
	my ($mauve_backbone_file, $r_chr_names, $aln_desc, $dbh) = @_;


	#CONCATENATE CHROMSOME NAMES
	my $chr_names = join "|", @$r_chr_names; 


	#INSERT MSA INFO AND GET ID

	#PREPARE INSERT STATEMENT
        my $insert_msa_info = "INSERT INTO multiple_alignment_info (chrs, description)
                               VALUES (?, ?)";

        #EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($insert_msa_info)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_names, $aln_desc)
                         or die "Can't execute statment: $DBI::errstr";

	#GET THE ID OF THE INSERTED MSA
	my $msa_id = $sth->{'mysql_insertid'};


	#GO THROUGH BACKBONE COLUMN AND INSERT ALIGNMENT
	open MAUVE, $mauve_backbone_file;

	#GET RID OF HEADER
	<MAUVE>;

	#GO THROUGH EACH LINE
	foreach my $line (<MAUVE>)
	{

		#SPLIT UP LINE
		chomp $line;
		my @line = split /\t/, $line;


		#GO THROUGH GETTING STARTS AND ENDS
		my @starts = ();
		my @ends = ();

		for(my $c = 0;$c < @line; $c=$c+2)
		{

			push @starts, $line[$c];
			push @ends, $line[$c+1];

		}#end for

		my $starts = join "|", @starts;
		my $ends = join "|", @ends;


		#INSERT INTO MSA TABLE
		#PREPARE INSERT STATEMENT
        	my $insert_msa = "INSERT INTO multiple_alignment (msa_id, starts, ends)
                	          VALUES (?, ?, ?)";

	        #EXECUTE THE STATEMENT
	        my $sth = $dbh -> prepare($insert_msa)
	                         or die "Can't prepare statment: $DBI::errstr";
	
	        my $rc = $sth -> execute($msa_id, $starts, $ends)
	                         or die "Can't execute statment: $DBI::errstr";

	}#end foreach


	return $msa_id;

}#end insert_mauve_MSA


#sub insert_msa	- Takes as input: 1) chr_names	- A list of all chromosomes in the alignment
#				  2) aln_desc	- A description of the alignment
#				  3) msa_hash	- A 2D hash where the keys are block_id and CHRS/STARTS/ENDS and
#				  		  the values are the corresponding values for the current alignment block
#		- Inserts the alignment into the multiple_alignment table
sub insert_msa
{
	my ($r_chr_names, $aln_desc, $r_msa_hash, $dbh) = @_;

	my @chr_names = @$r_chr_names;
	my %msa_hash = %$r_msa_hash;


	#IF NO CHROMOSMES, THEN ABORT AND REPORT ERROR
	if(@chr_names == 0)
	{

		print "\n\nNo alignment provided to insert_msa!!!!\n\n";

		return;

	}#end if


	#CONCATENATE CHROMSOME NAMES
	my $chr_names = join "|", @chr_names; 


	#INSERT MSA INFO AND GET ID
	#PREPARE INSERT STATEMENT
        my $insert_msa_info = "INSERT INTO multiple_alignment_info (chrs, description)
                               VALUES (?, ?)";

        #EXECUTE THE STATEMENT
        my $sth = $dbh -> prepare($insert_msa_info)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_names, $aln_desc)
                         or die "Can't execute statment: $DBI::errstr";

	#GET THE ID OF THE INSERTED MSA
	my $msa_id = $sth->{'mysql_insertid'};


	#GO THROUGH EACH ENTRY IN THE MSA HASH AND INSERT INTO DATABASE
	foreach my $block_id (sort {$a <=> $b} keys %msa_hash)
	{

		#INITIALIZE HASH WITH ALL CHROMOSOME IDs, AND PUT IN START AND END POSIITONS
		my %starts;
		@starts{@chr_names} = (0) x @chr_names;
		@starts{split(/\|/, $msa_hash{$block_id}{'CHRS'})} = split(/\|/, $msa_hash{$block_id}{'STARTS'});
		my $starts = join "\|", @starts{@chr_names};

		my @ends;
		@ends{@chr_names} = (0) x @chr_names;
		@ends{split(/\|/, $msa_hash{$block_id}{'CHRS'})} = split(/\|/, $msa_hash{$block_id}{'ENDS'});
		my $ends = join "\|", @ends{@chr_names};


		#INSERT INTO MSA TABLE
		#PREPARE INSERT STATEMENT
        	my $insert_msa = "INSERT INTO multiple_alignment (msa_id, starts, ends)
                	          VALUES (?, ?, ?)";

	        #EXECUTE THE STATEMENT
	        my $sth = $dbh -> prepare($insert_msa)
	                        or die "Can't prepare statment: $DBI::errstr";
	
	        my $rc = $sth -> execute($msa_id, $starts, $ends)
	                         or die "Can't execute statment: $DBI::errstr";

	}#end foreach


	return $msa_id;

}#end insert_MSA

#sub insert_mauve_MSA_snps - Takes as input: 1) mauve_snp_file  - A snp file produced by Mauve
#                                            2) msa_id		- The ID corresponding to the MSA from which the snps were retrieved from
#                          - Inserts Mauve SNPs into the database
sub insert_mauve_MSA_snps
{
	my ($mauve_snp_file, $msa_id, $dbh) = @_;


	#GO THROUGH SNP COLUMN AND INSERT SNPS
	open MAUVE, $mauve_snp_file;

	#GET RID OF HEADER
	<MAUVE>;

	#GO THROUGH EACH LINE
	foreach my $line (<MAUVE>)
	{

		#SPLIT UP LINE
		chomp $line;

		my ($snps, @positions) = split /\t/, $line;
		my $positions = join "|", @positions;
			

		#INSERT INTO MSA TABLE
		#PREPARE INSERT STATEMENT
        	my $insert_msa = "INSERT INTO multiple_alignment_snp (msa_id, snps, positions)
                	          VALUES (?, ?, ?)";

	        #EXECUTE THE STATEMENT
	        my $sth = $dbh -> prepare($insert_msa)
	                         or die "Can't prepare statment: $DBI::errstr";
	
	        my $rc = $sth -> execute($msa_id, $snps, $positions)
	                         or die "Can't execute statment: $DBI::errstr";

	}#end foreach


}#end insert_mauve_MSA_snps


#sub insert_diNucSignature - Takes as input: 1) chr_id	- The chromosome on which di-nuc signature should be calculated for genes 
#			   - Calculates the di-nucleotide signature for each ORF, averaged over all possible 6-ORF clusters of
#			     which it is a member
sub insert_diNucSignature
{
	my ($chr_id, $dbh) = @_;


	#GET ALL GENES ON CHROMOSOME IN ORDER
	my $r_chr_genes = &get_chr_genes($chr_id, $dbh);

	my $r_start2gene = &get_gene_start_5to3($r_chr_genes, $dbh);
	my %start2gene = %$r_start2gene;

	my @sorted_gene_starts = sort {$a <=> $b} keys %start2gene;

	my $r_gene2seq = &get_gene_sequences($r_chr_genes, $dbh);
	my %gene2seq = %$r_gene2seq;

	
	#CALCULATE THE DINUCLEOTIDE FREQUENCIES AND SIGNATURE ACROSS ALL ORFS ON THE CHROMSOOME
	my $r_all_orf_monoNucCounts = &seq_stats::init_monoNucHash();
	my $r_all_orf_diNucCounts = &seq_stats::init_diNucHash();


	foreach my $start (@sorted_gene_starts)
	{

		#CREATE OBJECT FOR CURRENT ORF
		my $seq_obj = Bio::Seq->new(-seq 		=> $gene2seq{$start2gene{$start}},
					    -display_id 	=> $start2gene{$start});

		my $rev_seq_obj = $seq_obj->revcom;

 
		#TABULATE MONO-NUCLEOTIDE FREQUENCIES
		$r_all_orf_monoNucCounts = &seq_stats::get_monoNucCount($r_all_orf_monoNucCounts, $seq_obj);		
		$r_all_orf_monoNucCounts = &seq_stats::get_monoNucCount($r_all_orf_monoNucCounts, $rev_seq_obj);		


		#TABULATE DI-NUCLEOTIDE FREQUENCIES
		$r_all_orf_diNucCounts = &seq_stats::get_diNucCount($r_all_orf_diNucCounts, $seq_obj);		
		$r_all_orf_diNucCounts = &seq_stats::get_diNucCount($r_all_orf_diNucCounts, $rev_seq_obj);		


	}#end foreach


	#CALCULATE GENOMIC SIGNATURE FOR ALL ORFS IN THE GENOME
	my $r_all_orf_diNucSig = &seq_stats::get_diNucSignature($r_all_orf_monoNucCounts, $r_all_orf_diNucCounts);


	#CALCULATE DI-NUCLEOTIDE FREQUENCIES AND SIGNATURE FOR EACH GENE AND THE THREE ORFS TO THE LEFT AND RIGHT
	#NEXT, ASSIGN A VALUE TO EACH GENE WHICH IS THE MEAN DIFFERENCE IN DINUCLEOTIDE SIGNATURE COMPARED TO
	#THE GENOMIC BACKGROUND
	for(my $s=0; $s < (@sorted_gene_starts); $s++)
	{


		#GET MONO AND DI NUCLEOTIDE COUNTS FOR EACH ORF
		my $r_orf_monoNucCounts = &seq_stats::init_monoNucHash();
        	my $r_orf_diNucCounts = &seq_stats::init_diNucHash();


		for(my $i = -3; $i <= 3; $i++)
		{

			
			#CREATE OBJECT FOR CURRENT ORF
			my $seq_obj = Bio::Seq->new(-seq 	=> $gene2seq{$start2gene{$sorted_gene_starts[($s + $i)%(@sorted_gene_starts)]}},
						    -display_id => $start2gene{$sorted_gene_starts[$s + $i]});
		
			my $rev_seq_obj = $seq_obj->revcom;


			#TABULATE MONO-NUCLEOTIDE FREQUENCIES
			$r_orf_monoNucCounts = &seq_stats::get_monoNucCount($r_orf_monoNucCounts, $seq_obj);
			$r_orf_monoNucCounts = &seq_stats::get_monoNucCount($r_orf_monoNucCounts, $rev_seq_obj);

			#TABULATE DI-NUCLEOTIDE FREQUENCIES
			$r_orf_diNucCounts = &seq_stats::get_diNucCount($r_orf_diNucCounts, $seq_obj);
			$r_orf_diNucCounts = &seq_stats::get_diNucCount($r_orf_diNucCounts, $rev_seq_obj);

		}#end for


		#CALCULATE THE DI-NUCLEOTIDE SIGNATURE FOR THE CURRENT ORF-BLOCK
		my $r_orf_diNucSig = &seq_stats::get_diNucSignature($r_orf_monoNucCounts, $r_orf_diNucCounts);


		#CALCULATE THE AVERAGE DIFFERENCE BETWEEN THE DI-NUCLEOTIDE SIGNATURE IN THE CURRENT BLOCK AND THE GENOME BACKGROUND
		my $diNucSigDiff = &arith_ops::rnd(&seq_stats::get_diNucSignatureDiff($r_orf_diNucSig, $r_all_orf_diNucSig),4);		

		#print "$start2gene{$sorted_gene_starts[$s]}\t$diNucSigDiff\n";


		#UPDATE GENE ENTRY WITH DI-NUCLEOTIDE SIGNATURE
		my $update_gene = "UPDATE gene
                          	   SET di_nuc_signature = ? 
                        	   WHERE name = ?";

	        #EXECUTE THE STATEMENT
	        my $sth = $dbh -> prepare($update_gene)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($diNucSigDiff, $start2gene{$sorted_gene_starts[$s]})
                	         or die "Can't execute statment: $DBI::errstr";


	}#end for


}#end insert_diNucSignature


#sub insert_GCcontent - Takes as input: 1) chr_id  - The chromosome on which GC content should be calculated for genes 
#                     - Calculates the GC content for each gene and inserts in the database
sub insert_GCcontent
{
        my ($chr_id, $dbh) = @_;


        #GET ALL GENES ON CHROMOSOME
        my $r_chr_genes = &get_chr_genes($chr_id, $dbh);

        my $r_gene2seq = &get_gene_sequences($r_chr_genes, $dbh);
        my %gene2seq = %$r_gene2seq;


	#CALCULATE GC CONTENT AND INSERT INTO DATABASE
	foreach my $gene (sort keys %gene2seq)
	{

		#CALCULATE BASE COMPOSITION
		my $seq_obj = Bio::Seq->new(-seq        => $gene2seq{$gene},
                                            -display_id => $gene);

		my $r_orf_monoNucCounts = &seq_stats::init_monoNucHash();
		$r_orf_monoNucCounts = &seq_stats::get_monoNucCount($r_orf_monoNucCounts, $seq_obj);


		#CALCULATE GC CONTENT
		my $gene_gc = &arith_ops::rnd(&seq_stats::get_GC($r_orf_monoNucCounts),2);		


		#UPDATE DATABASE ENTRY WITH GC CONTENT
                my $update_gene = "UPDATE gene
                                   SET gc = ? 
                                   WHERE name = ?";

                #EXECUTE THE STATEMENT
                my $sth = $dbh -> prepare($update_gene)
                                 or die "Can't prepare statment: $DBI::errstr";

		#print "$gene\t$gene_gc\n";

                my $rc = $sth -> execute($gene_gc, $gene)
                                or die "Can't execute statment: $DBI::errstr";




	}#end foreach

}#end insert_GCcontent


#sub insert_mobile_element - Takes as input: 1) mobile_hash	- A hash whose first key is an element count and second key is CHR, START, END, CLASS, P_VAL
#			   - Inserts into mobile_element table
sub insert_mobile_element
{
	my ($r_mobile_hash, $dbh) = @_;

        my %mobile_hash = %$r_mobile_hash;


        #INSERT EACH ORTHOLOG GROUP
        foreach my $mobile_count (keys %mobile_hash)
        {

		#GET MAPPING BETWEEN CHROMSOME NAME AND ID
		my $chr_id = &get_chr_id($mobile_hash{$mobile_count}{'CHR'}, $dbh);	


                #INSERT ANNOTATION
                my $insert_mge = "INSERT INTO mobile_element (chr_id, start, end, class, p_value)
                                  VALUES (?, ?, ?, ?, ?)";


                #EXECUTE THE STATEMENT
                my $sth = $dbh -> prepare($insert_mge)
                                 or die "Can't prepare statment: $DBI::errstr";


                my $rc = $sth -> execute($chr_id, $mobile_hash{$mobile_count}{'START'},
                                         $mobile_hash{$mobile_count}{'END'}, $mobile_hash{$mobile_count}{'CLASS'},
				   	 $mobile_hash{$mobile_count}{'P_VAL'})
                                 or die "Can't execute statment: $DBI::errstr";


        }#end foreach

}#end insert_mobile_element


#sub insert_snps - Takes as input: 1) snp_file	- A tab delimited output file produced by the mummer program, show-snps
#						  The parameters used were: show-snps -CTlr -x 10
#		 - Parses file and inserts into snp table
sub insert_snps
{
	my ($snp_file, $dbh) = @_;


	#GET HASH LINKING CHROMOSOME NAME'S TO ID'S
	my $r_chr_name2id = &get_chr_id_all($dbh);
	my %chr_name2id = %$r_chr_name2id;
	

	#OPEN FILE
        open INPUT, $snp_file;


	#GO THROUGH FILE AND MAKE INSERTIONS INTO TABLE
        foreach my $line (<INPUT>)
        {

		#CHECK IF HEADER SECTION IS OVER
		if( $line !~ /^[0-9]+\t/ )
		{
			
			next;

		}#end if


		#PARSE LINE AND INSERT INTO DATABASE
                chomp $line;
                my ($position1, $base1, $base2, $position2, $buffer, $dist, $length1, $length2, $seq1, $seq2, $frm1, $frm2, $chr1, $chr2) = split /\t/, $line;


                #INSERT SNP IF POSITION IS ON CONTIG AND NOT LINKER SEQUENCE
		 #my ($contig1, $contig_start1, $contig_end1, $contig_strand1, $contig_desc1) = &chr2contig_coords($chr_name2id{$chr1}, $position1, $position1, $dbh);
		 #my ($contig2, $contig_start2, $contig_end2, $contig_strand2, $contig_desc2) = &chr2contig_coords($chr_name2id{$chr2}, $position2, $position2, $dbh);

		#if(($contig_start1 != 0) && ($contig_start2 != 0))
		#{

                	my $insert_snp = "INSERT INTO snp (chr_id1, position1, base1, chr_id2, position2, base2)
                        	           VALUES (?, ?, ?, ?, ?, ?)";


                	#EXECUTE THE STATEMENT
                	$sth = $dbh -> prepare($insert_snp)
                	                 or die "Can't prepare statment: $DBI::errstr";

                	$rc = $sth -> execute($chr_name2id{$chr1}, $position1, $base1, $chr_name2id{$chr2}, $position2, $base2)
                        	         or die "Can't execute statment: $DBI::errstr";

		#}#end if

        }#end foreach


        close INPUT;


}#end insert_snps


#sub insert_snps_vcf - Takes as input: 1) vcf_tab_file	- Variant file produced by samtools
#				       2) ref_chr	- Chromosome to which reads were aligned
#				       3) aln_chr	- Chromosome from which reads were aligned
#				       4) method	- Method for read alignment and variant calling
#		     - Parses file and inserts into snp table
sub insert_snps_vcf
{
	my ($vcf_tab_file, $ref_chr, $aln_chr, $method, $dbh) = @_;

	#GET HASH LINKING CHROMOSOME NAME'S TO ID'S
	my $r_chr_name2id = &get_chr_id_all($dbh);
	my %chr_name2id = %$r_chr_name2id;
	

	#OPEN FILE
        open INPUT, $vcf_tab_file;


	#PARSE HEADER
	my $header = <INPUT>;
	chomp $header;
	my @header = split /\t/, $header;


	#GO THROUGH FILE AND MAKE INSERTIONS INTO TABLE
        foreach my $line (<INPUT>)
        {

		#PARSE LINE AND INSERT INTO DATABASE
                chomp $line;
		my @line = split /\t/, $line;


		#MAKE SURE ALN IN HAPLOID FORM
		$line[3] =~ s/^([^\/]+)\/.*$/$1/;	


                #INSERT SNP
                my $insert_snp = "INSERT INTO snp (chr_id1, position1, base1, chr_id2, position2, base2, method)
                                  VALUES (?, ?, ?, ?, ?, ?)";


                #EXECUTE THE STATEMENT
                $sth = $dbh -> prepare($insert_snp)
                                 or die "Can't prepare statment: $DBI::errstr";

                $rc = $sth -> execute($chr_name2id{$ref_chr}, $line[1], $line[2], $chr_name2id{$aln_chr}, '', $line[3], $method)
                                 or die "Can't execute statment: $DBI::errstr";

        }#end foreach


        close INPUT;


}#end insert_snps_vcf


#sub insert_snps_varscan - Takes as input: 1) vs_hash	- A hash containing varScan SNP and indel data (hash{pos} = ref\tvar)
#					   2) ref_chr	- Chromosome to which reads were aligned
#				           3) aln_chr	- Chromosome from which reads were aligned
#				           4) method	- Method for read alignment and variant calling
#		         - Parses file and inserts into snp table
sub insert_snps_varscan
{
	my ($r_vs_hash, $ref_chr, $aln_chr, $method, $dbh) = @_;

	#GET HASH LINKING CHROMOSOME NAME'S TO ID'S
	my $r_chr_name2id = &get_chr_id_all($dbh);
	my %chr_name2id = %$r_chr_name2id;
	

	#GO THROUGH FILE AND MAKE INSERTIONS INTO TABLE
	my %vs_hash = %$r_vs_hash;

        foreach my $pos (keys %vs_hash)
        {

		#PARSE LINE
		my ($ref, $var) = split /\t/, $vs_hash{$pos};


                #INSERT SNP
                my $insert_snp = "INSERT INTO snp (chr_id1, position1, base1, chr_id2, position2, base2, method)
                                  VALUES (?, ?, ?, ?, ?, ?)";


                #EXECUTE THE STATEMENT
                $sth = $dbh -> prepare($insert_snp)
                                 or die "Can't prepare statment: $DBI::errstr";

                $rc = $sth -> execute($chr_name2id{$ref_chr}, $pos, $ref, $chr_name2id{$aln_chr}, '', $var, $method)
                                 or die "Can't execute statment: $DBI::errstr";

        }#end foreach


}#end insert_snps_varscan


#sub insert_repeats - Takes as input: 1) repeat_hash 	- A hash where the keys are alingmnet IDs and the values are descriptors of
#						     	 the repeat
#				     2) chr_name	- The name of the chromosome in which the repeats were found
#				     3) program		- The program used to generate the repeats
#		   - Inserts into repeat table
sub insert_repeats
{
	my ($r_repeat_hash, $chr_name, $program, $dbh) = @_;

	my %repeat_hash = %$r_repeat_hash;
	my $chr_id = &get_chr_id($chr_name, $dbh);
	my %insert_hash;

	#GO THROUGH HASH AND MAKE INSERTIONS INTO TABLE
	foreach my $rep_id (sort {$a<=>$b} keys %repeat_hash)
	{

		#CHECK IF REPEAT ALREADY INSERTED, IF NOT, THEN INSERT
		my $repeat = "$chr_id, $repeat_hash{$rep_id}{'S1'}, $repeat_hash{$rep_id}{'E1'}, $repeat_hash{$rep_id}{'S2'}, $repeat_hash{$rep_id}{'E2'}";

		if(!exists($insert_hash{$repeat}))
		{

                	#INSERT REPEAT
                	my $insert_repeats = "INSERT INTO repeats (chr_id, start1, end1, start2, end2, identity, significance, program)
                                     VALUES (?, ?, ?, ?, ?, ?, ?, ?)";


               		#EXECUTE THE STATEMENT
                	my $sth = $dbh -> prepare($insert_repeats)
                	                 or die "Can't prepare statment: $DBI::errstr";

                	my $rc = $sth -> execute($chr_id, $repeat_hash{$rep_id}{'S1'}, $repeat_hash{$rep_id}{'E1'}, $repeat_hash{$rep_id}{'S2'}, $repeat_hash{$rep_id}{'E2'}, $repeat_hash{$rep_id}{'ID'}, $repeat_hash{$rep_id}{'SIG'}, $program)
                        	         or die "Can't execute statment: $DBI::errstr";
			

			$insert_hash{$repeat} = 1;

		}#end if

	}#end foreach

}#end insert_repeats


#sub insert_tandem_repeat - Takes as input: 1) repeat_hash 	- A hash where the keys are alingmnet IDs and the values are descriptors of
#						     	 	the repeat
#				     	    2) chr_name		- The name of the chromosome in which the repeats were found
#				     	    3) program		- The program used to generate the repeats
#		   	  - Inserts into repeat table
sub insert_tandem_repeat
{
	my ($r_repeat_hash, $chr_name, $program, $dbh) = @_;

	my %repeat_hash = %$r_repeat_hash;
	my $chr_id = &get_chr_id($chr_name, $dbh);


	#GO THROUGH HASH AND MAKE INSERTIONS INTO TABLE
	foreach my $rep_id (sort {$a<=>$b} keys %repeat_hash)
	{

                #INSERT REPEAT
                my $insert_repeat = "INSERT INTO repeat_tandem (chr_id, start, unit_length, copy_number, program)
                                     VALUES (?, ?, ?, ?, ?)";


                #EXECUTE THE STATEMENT
                my $sth = $dbh -> prepare($insert_repeat)
                                 or die "Can't prepare statment: $DBI::errstr";

                my $rc = $sth -> execute($chr_id, $repeat_hash{$rep_id}{'START'}, $repeat_hash{$rep_id}{'UNIT_LEN'}, $repeat_hash{$rep_id}{'NUM_COPIES'}, $program)
                                 or die "Can't execute statment: $DBI::errstr";
			

	}#end foreach

}#end insert_tandem_repeat


#################################################### GENE SELECTION ROUTINES ##################################################

#sub get_unique_genes - Takes as input: 1) strain1	- The strain in which genes should be present
#					2) strain2	- The strain in which genes should not be present
#		      - Returns a list of unique genes
sub get_unique_genes
{
	my ($strain1, $strain2, $dbh) = @_;

	my @unique_genes;
	my %og_genes;


	#SELECT ALL ORTHOLOG GROUPS INVOLVING ONE ORG AND NOT THE OTHER
        my $select_ogs = "SELECT gene_names, org_ids
			  FROM `ortholog_group` 
			  WHERE (org_ids LIKE \'%|$strain1\' OR org_ids LIKE \'%|$strain1|%\') AND 
				org_ids NOT LIKE \'%|$strain2|%\' AND 
				org_ids NOT LIKE \'%|$strain2\'";


        my $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


	#EXTRACT THE GENES FROM THE DESIRED ORGANISM
        while(my ($gene_ids, $org_ids) =  $sth -> fetchrow_array)
        {
	

                my @gene_ids = split /\|/, $gene_ids;
                my @org_ids = split /\|/, $org_ids;


		for(my $i=0; $i < @org_ids;$i++)
		{

			if($org_ids[$i] eq $strain1)
			{

				push @unique_genes, $gene_ids[$i];

			}#end if

		}#end for


        }#end while
	

	#GET ALL ORTHOLOG GROUPS INVOLVING A STRAIN
	$select_ogs = "SELECT gene_names, org_ids
                          FROM `ortholog_group` 
                          WHERE (org_ids LIKE '%|$strain1' OR org_ids LIKE '%|$strain1|%')";


        $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #EXTRACT THE GENES FROM THE DESIRED ORGANISM
        while(my ($gene_ids, $org_ids) =  $sth -> fetchrow_array)
        {


                my @gene_ids = split /\|/, $gene_ids;
                my @org_ids = split /\|/, $org_ids;


                for(my $i=0; $i < @org_ids;$i++)
                {

                        if($org_ids[$i] eq $strain1)
                        {

                                $og_genes{$gene_ids[$i]} = 1;

                        }#end if

                }#end for


        }#end while

	#GET ALL GENES FOR THE GIVEN STRAIN
	my $r_all_genes = &get_org_genes($strain1, $dbh);	


	#IDENTIFY THOSE GENES NOT IN ANY ORTHOLOG GROUP
	foreach my $gene (@$r_all_genes)
	{

		if( !(exists($og_genes{$gene})) )
		{

			push @unique_genes, $gene;

		}#end if

	}#end foreach

	
	#print "There are " .  scalar(@unique_genes) . " unique genes found\n\n";


	#GET SEQUENCES OF GENES
	my $r_unique_genes = &get_gene_sequences(\@unique_genes, $dbh);


	return $r_unique_genes;

}#end get_unique_genes


#sub get_org_genes - Takes as input: 1) strain	- The name of a strain from the organism table
#		   - Returns a list of gene names present in the desired strain
sub get_org_genes
{
	my ($strain, $dbh) = @_;

	my @genes = ();


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_genes = "SELECT DISTINCT g.name
                         FROM gene g, chromosome c, organism o
                         WHERE g.chr_id = c.chr_id AND
			       c.org_id = o.org_id AND
			       g.type = \'protein\' AND
			       o.strain = ?";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($strain)
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($gene) =  $sth -> fetchrow_array)
        {

		push @genes, $gene;

	}#end while


	return \@genes;

}#end get_org_genes


#sub get_chr_genes - Takes as input: 1) chr_id  - An ID from the chromosome table
#                  - Returns a list of gene names present on the desired chromosome
sub get_chr_genes
{
        my ($chr_id, $dbh) = @_;

        my @genes = ();


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_genes = "SELECT DISTINCT g.name
                         FROM gene g, chromosome c
                         WHERE g.chr_id = ? AND
                               g.type = \'protein\'";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($gene) =  $sth -> fetchrow_array)
        {

                push @genes, $gene;

        }#end while


        return \@genes;

}#end get_chr_genes


#sub get_chr_region_genes -Takes as input: chr_id	- An ID from the chromsome table
#					   start	- A start position on the chromosome
#					   end		- An end position on the chromosme
#			  - Returns a list of genes completely in the given region
sub get_chr_region_genes
{
	my ($chr_id, $start, $end, $dbh) = @_;

        my @genes = ();


        #GET ALL GENES PRESENT IN GIVEN REGION
        $select_genes = "SELECT DISTINCT g.name
                         FROM gene g, chromosome c
                         WHERE g.chr_id = ? AND
                               g.start > ? AND
			       g.end < ?";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id, &arith_ops::min($start,$end), &arith_ops::max($start,$end))
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($gene) =  $sth -> fetchrow_array)
        {

                push @genes, $gene;

        }#end while


        return \@genes;

}#end get_chr_region_genes


#sub get_gene_ids - Takes as input: 1) r_gene_list	- A reference to a list of names from the gene table
#		  - Returns a hash linking gene names to gids
sub get_gene_ids
{
	my ($r_gene_list, $dbh) = @_;
	
        my %gene2gid;


        #GET gids CORRESPONDING TO LIST OF GENE NAMES
        my $select_genes = "SELECT gid, name
                             FROM gene 
                             WHERE name IN (";


        $select_genes = &db_lib::add_place_holders($select_genes , scalar(@$r_gene_list));

        $select_genes = $select_genes . ")";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_list)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $gene_name) =  $sth -> fetchrow_array)
        {

                $gene2gid{$gene_name} = $gid;

        }#end while


	return \%gene2gid;

}#end get_gene_ids


#sub get_all_gene_ids - Takes as input: 
#		      - Returns a hash linking all gene names to gids
sub get_all_gene_ids
{
	my ($dbh) = @_;
	
        my %gene2gid;


        #GET gids CORRESPONDING TO LIST OF GENE NAMES
        my $select_genes = "SELECT gid, name
                             FROM gene";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $gene_name) =  $sth -> fetchrow_array)
        {

                $gene2gid{$gene_name} = $gid;

        }#end while


	return \%gene2gid;

}#end get_all_gene_ids


#sub get_gene_annotations - Takes as input: 1) genes	- A reference to a list of gene names in the gene table
#			  - Returns annotations for those genes in a 2D hash {gene}{class}
sub get_gene_annotations
{
	my ($r_gene_list, $dbh) = @_;

	my %gid2gene;
	my %gene_annots;
	
	my $ga_dbh = &gene_annotation_db_lib::gene_annotation_db_connect();
	my $r_cog_info = &gene_annotation_db_lib::get_cog_info($ga_dbh);	
	my %cog_info = %$r_cog_info;


	#IF NO GENES THEN GET OUT	
	if(@$r_gene_list == 0)
	{

		return \%gene_annots;

	}#end if

	
        #GET gids CORRESPONDING TO LIST OF GENE NAMES
        my $select_genes = "SELECT gid, name
                             FROM gene 
                             WHERE name IN (";


        $select_genes = &db_lib::add_place_holders($select_genes , scalar(@$r_gene_list));

        $select_genes = $select_genes . ")";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_list)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $gene_name) =  $sth -> fetchrow_array)
        {

		$gid2gene{$gid} = $gene_name;

        }#end while




	#GET ANNOTATIONS FOR GIDS RECOVERTED
        my $select_annots = "SELECT gid, annotation, source, hit, e_value
			     FROM annotation 
			     WHERE  gid IN (";


	$select_annots = &db_lib::add_place_holders($select_annots , scalar(keys %gid2gene));

        $select_annots = $select_annots . ")";


        $sth = $dbh -> prepare($select_annots)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(keys %gid2gene)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $annotation, $source, $hit, $e_value) =  $sth -> fetchrow_array)
        {

		my $gene = $gid2gene{$gid};

		#IF NR ANNOTATION THEN TRUNCATE
		if($source eq 'NR')
		{
                        $annotation =~ s/^([^\[]+) \[.*$/$1/;

		}#end if


		#CREATE LIST OF COG CODES ASSOCIATED WITH GENE
		if($source eq 'COG')
		{

			my @hits = split /\s?,\s?/, $hit;
			my %codes;

			foreach my $hit (@hits)
			{

				my $codes = $cog_info{$hit}{'FUNC'};;
				my @codes = split //, $codes;

				foreach my $code (@codes)
				{

					$codes{$code} = 1;

				}#end foreach
	
			}#end foreach

			@{$gene_annots{$gene}{'COG CODE'}} = keys %codes;

		}#end if


		#STORE ANNOTATION
		$gene_annots{$gene}{$source} = $annotation;

        }#end while
	

	return \%gene_annots;

}#end get_gene_annotations


#sub get_gene_annotation_list - Takes as input: 1) r_gene_list - A list of genes from the given database
#			      - Returns a list of NCBI/NR annotations for the inputted list of genes
sub get_gene_annotation_list
{
	my ($r_gene_list, $dbh) = @_;

	my %gid2gene;
	my %gene2annots;
	
	my $ga_dbh = &gene_annotation_db_lib::gene_annotation_db_connect();
	my $r_cog_info = &gene_annotation_db_lib::get_cog_info($ga_dbh);	
	my %cog_info = %$r_cog_info;

	
        #GET gids CORRESPONDING TO LIST OF GENE NAMES
        my $select_genes = "SELECT gid, name
                             FROM gene 
                             WHERE name IN (";


        $select_genes = &db_lib::add_place_holders($select_genes , scalar(@$r_gene_list));

        $select_genes = $select_genes . ")";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_list)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $gene_name) =  $sth -> fetchrow_array)
        {

		$gid2gene{$gid} = $gene_name;

        }#end while




	#GET ANNOTATIONS FOR GIDS RECOVERTED
        my $select_annots = "SELECT gid, annotation, source, hit, e_value
			     FROM annotation 
			     WHERE  gid IN (";


	$select_annots = &db_lib::add_place_holders($select_annots , scalar(keys %gid2gene));

        $select_annots = $select_annots . ")";


        $sth = $dbh -> prepare($select_annots)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(keys %gid2gene)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gid, $annotation, $source, $hit, $e_value) =  $sth -> fetchrow_array)
        {

		my $gene = $gid2gene{$gid};

		$gene2annots{$gene}{$source} = $annotation;

        }#end while
	

	#GO THROUGH GENES AND MAKE LIST OF ANNOTATIONS
	my @region_annots;

	foreach my $gene(@$r_gene_list)
	{

		#USE NCBI ANNOTATION IF EXISTS, OTHERWISE USE TRUNCATED NR
		my $annot = "";

		if(exists($gene2annots{$gene}{'NCBI'}))
                {

                        $annot = $gene2annots{$gene}{'NCBI'};

                }elsif(exists($gene2annots{$gene}{'NR'})){

                        $annot = $gene2annots{$gene}{'NR'};
                        $annot =~ s/^([^\[]+)\s\[.*$/$1/;

		}else{

                        $annot = $gene2annots{$gene}{'KEGG'};

                	#REMOVE KO INFO FROM THE HIT DESCRIPTION                                                                                    
                	$annot =~ s/^(.+?); K\d+.+$/$1/;                                                                                         
                	$annot =~ s/^([^\)]+)\(.*$/$1/;                                                                                          
                                                                                                                                            
                                                                                                                                            
                	#PARSE OUT GENE SYMBOL (IF EXISTS) AND GENE DESCRIPTION                                                                     
                	if($annot =~ /;/)                                                                                                        
                	{                                                                                                                           
                                                                                                                                            
                	        $annot =~ s/^([^;]+); ([^;]+)$/$2/;                                                                              
                                                                                                                                            
                	}#end if  


                }#end if

                push @region_annots, $annot;

	}#end foreach


	return \@region_annots;

}#end get_gene_annotation_list


#sub get_gene_COGs - Takes as input: 1) r_genes	- A list of genes
#		   - Returns 1) a 2D hash where the first key is the gene name and the second is a COG 1 letter code
#		     to which the gene belongs AND 2) a 2D hash where the first key is the gene name and the second
#		     is a COG ID to which it belongs.
sub get_gene_COGs
{
	my ($r_genes, $dbh) = @_;

	my %cog_hash;
	my %cog_desc_hash;
	my %cog_code_hash;


	#GET ANNOTATIONS FOR GIDS RECOVERTED
	my $select_annots = "SELECT g.name, a.hit, a.annotation
		FROM annotation a, gene g
		WHERE g.gid = a.gid AND
		a.source = 'COG' AND
		g.name IN (";


				$select_annots = &db_lib::add_place_holders($select_annots , scalar(@$r_genes));

				$select_annots = $select_annots . ")";


				$sth = $dbh -> prepare($select_annots)
				or die "Can't prepare statment: $DBI::errstr";

				$rc = $sth -> execute(@$r_genes)
				or die "Can't execute statment: $DBI::errstr";


	#GO THROUGH RESULTS SOTIRING COGS AND CODES
	my $ga_dbh = &gene_annotation_db_lib::gene_annotation_db_connect();
	my $r_cog_info = &gene_annotation_db_lib::get_cog_info($ga_dbh);	
	my %cog_info = %$r_cog_info;

	$r_cog_code_info = &gene_annotation_db_lib::get_cog_codes($ga_dbh);


        while(my ($gene, $hits, $descs) =  $sth -> fetchrow_array)
        {

		$cog_hash{$gene} = $hits;
		$cog_desc_hash{$gene} = $descs;
		my @hits = split /\s?,\s?/, $hits;

		foreach my $hit(@hits)
		{
			
			my $codes = $cog_info{$hit}{'FUNC'};;
			my @codes = split //, $codes;

			foreach my $code (@codes)
			{

				$cog_code_hash{$gene}{$code} = 1;

			}#end foreach
			
		}#end foreach

        }#end while

	
	return (\%cog_hash, \%cog_code_hash, \%cog_desc_hash, $r_cog_code_info);

}#end get_gene_COGs


#sub get_gene_KEGG - Takes as input: 1) r_genes	- A list of genes
#		   - Returns 1) a 2D hash where the first key is the gene name and the second is 1) KEGG symbol
#		     2)     
sub get_gene_KEGG
{
	my ($r_genes, $dbh) = @_;

	my %KEGG_hash;


	#GET ANNOTATIONS FOR GIDS RECOVERTED
	my $select_annots = "SELECT g.name, a.annotation
		FROM annotation a, gene g
		WHERE g.gid = a.gid AND
		a.source = 'KEGG' AND
		g.name IN (";


	$select_annots = &db_lib::add_place_holders($select_annots , scalar(@$r_genes));

	$select_annots = $select_annots . ")";


	$sth = $dbh -> prepare($select_annots)
			or die "Can't prepare statment: $DBI::errstr";

	$rc = $sth -> execute(@$r_genes)
			or die "Can't execute statment: $DBI::errstr";


	#GO THROUGH RESULTS SOTIRING KEGG HITS AND ANNOTATIONS
        while(my ($gene, $hit_desc) =  $sth -> fetchrow_array)
        {


		#REMOVE KO INFO FROM THE HIT DESCRIPTION
		$hit_desc =~ s/^(.+?); K\d+.+$/$1/;
		$hit_desc =~ s/^([^\)]+)\(.*$/$1/;
		my $gene_symbol = $hit_desc;


		#PARSE OUT GENE SYMBOL (IF EXISTS) AND GENE DESCRIPTION
		if($hit_desc =~ /;/)
		{

			$hit_desc =~ s/^([^;]+); ([^;]+)$/$2/;
			$gene_symbol =~ s/^([^;]+); ([^;]+)$/$1/;

		}else{                                                                                                              

			$gene_symbol = "N\\A";                                                                                      

		}#end if

		$kegg_hash{$gene}{'ANNOT'} = $hit_desc;
		$kegg_hash{$gene}{'SYMBOL'} = $gene_symbol;

        }#end while

	
	return (\%kegg_hash);

}#end get_gene_COGs

#sub get_gene_desc - Takes as input: 1) r_genes - A list of genes
#		   - Retruns a hash where the keys are gene names and 'FRAME SHIFT', 'STOP CODON', 'CONTIG',
#		     and 'LENGTH'
sub get_gene_desc
{
	my ($r_genes, $dbh) = @_;

	my %gene_desc;


        #GET GENE NAMES AND OG GROUPINGS
        my $select_genes = "SELECT g.name, p.mut_type, p.position, g.chr_id
                             FROM pseudogene p, gene g
                             WHERE p.gid = g.gid AND
				   g.name IN (";

	$select_genes = &db_lib::add_place_holders($select_genes , scalar(@$r_genes));

        $select_genes = $select_genes . ")";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $mut_type, $position, $chr_id) =  $sth -> fetchrow_array)
        {

		if($mut_type eq 'FRAME SHIFT')
		{

			#CHECK IF FRAME SHIFT MUTATION IS LIKELY SEQUENCEING ERROR (e.g. poly A or poly T)
			my $seq = "";

			if($position != 0)
			{

				$seq = &extract_sequence($chr_id, abs($position - 1 - 10), abs($position - 1), $dbh);

			}#end if


			if($seq !~ /AAAA|TTTT|GGGG|CCCC|NNNN/i && $seq !~ /[atgc]/ && $seq ne "")
			{

				$gene_desc{$name}{'FRAME SHIFT'} = 1;
				print "$seq\n";

			}else{

				$gene_desc{$name}{'BAD FRAME SHIFT'} = 1;

			}#end if

		}elsif($mut_type eq 'STOP CODON'){

			$gene_desc{$name}{'STOP CODON'} = 1;

		}elsif($mut_type eq 'CONTIG BREAK' ){

			$gene_desc{$name}{'CONTIG BREAK'} = 1;

		}#end if


	}#end while


        #GET GENE NAMES AND OG GROUPINGS
        $select_genes = "SELECT g.name, g.start, g.end
                             FROM gene g
                             WHERE g.name IN (";

        $select_genes = &db_lib::add_place_holders($select_genes , scalar(@$r_genes));

        $select_genes = $select_genes . ")";

        $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $start, $end) =  $sth -> fetchrow_array)
        {

                my $length = &arith_ops::abs(($start-$end));

                $gene_desc{$name}{'LENGTH'} = $length;

        }#end while



	return \%gene_desc;

}#end get_gene_desc


#get_pseudo_genes - Takes as input: 
#		  - Returns a reference to a hash where the keys are gene names and values are potential
#		    orthologous groups that they may fit into
sub get_pseudo_genes
{
	my ($dbh) = @_;

	my %pseudo_hash;


        #GET GENE NAMES AND OG GROUPINGS
        my $select_genes = "SELECT g.name, p.og_id
                             FROM pseudogene p, gene g
                             WHERE p.gid = g.gid";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $og_id) =  $sth -> fetchrow_array)
        {

		$pseudo_hash{$name} = $og_id;

	}

	return \%pseudo_hash;


}#end get_pseudo_genes


#sub get_gene_sequences - Takes as input: 1) r_gene_names	- A list of gene names from the gene table
#			- Returns a hash where the key is the gene name and the value is the sequence
sub get_gene_sequences
{
	my ($r_gene_names, $dbh) = @_;

	my %gene2seq_hash;

	
	#SELECT THE START AND END COORDINATES OF EACH GENE
        my $select_coords = "SELECT name, chr_id, start, end
                             FROM gene
                             WHERE name IN (";


        $select_coords = &db_lib::add_place_holders($select_coords , scalar(@$r_gene_names));

        $select_coords = $select_coords . ")";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_names)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene_name, $chr_id, $start, $end) =  $sth -> fetchrow_array)
        {

		my $seq = &extract_sequence($chr_id, $start, $end, $dbh);

		$gene2seq_hash{$gene_name} = $seq;

        }#end while


	return \%gene2seq_hash;

}#end get_gene_sequence


#sub get_gene_coords- Takes as input: 1) gene - The name of a gene
#		    - Returns the start and end coordinates of the gene
sub get_gene_coords
{
        my ($r_gene_names, $dbh) = @_;

        my %gene2coords_hash;


        #SELECT THE START AND END COORDINATES OF EACH GENE
        my $select_coords = "SELECT name, chr_id, start, end
                             FROM gene
                             WHERE name IN (";


        $select_coords = &db_lib::add_place_holders($select_coords , scalar(@$r_gene_names));

        $select_coords = $select_coords . ")";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_names)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene_name, $chr_id, $start, $end) =  $sth -> fetchrow_array)
        {

                $gene2coords_hash{$gene_name}{'CHR_ID'} = $chr_id;
                $gene2coords_hash{$gene_name}{'START'} = $start;
                $gene2coords_hash{$gene_name}{'END'} = $end;

        }#end while


        return \%gene2coords_hash;

}#end get_gene_coords


#sub get_gene_start_5to3 - Takes as input: 1) r_genes - A list of genes
#			 - Returns a hash where the keys are the positions of the beginning of 
#			   genes when scanning accross the chromosome 5' to 3'
sub get_gene_start_5to3
{
	my ($r_gene_names, $dbh) = @_;


	my %start2gene_hash;


        #SELECT THE START AND END COORDINATES OF EACH GENE
        my $select_coords = "SELECT name, start, end
                             FROM gene
                             WHERE name IN (";


        $select_coords = &db_lib::add_place_holders($select_coords , scalar(@$r_gene_names));

        $select_coords = $select_coords . ")";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_gene_names)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene_name,  $start, $end) =  $sth -> fetchrow_array)
        {

		my $start_coord = &arith_ops::min($start, $end);

                $start2gene_hash{$start_coord} = $gene_name;

        }#end while


        return \%start2gene_hash;


}#end get_gene_start_5to3



#sub get_gene2strain - Takes as input:
#		     - Returns a hash linking gene names to the strain which they are present in
sub get_gene2strain
{
	my ($dbh) = @_;

	my %gene2strain;


       #GET GENE NAMES AND OG GROUPINGS
        my $select_genes = "SELECT g.name, o.strain
                             FROM chromosome c, organism o, gene g
                             WHERE g.chr_id = c.chr_id AND
				   c.org_id = o.org_id";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $strain) =  $sth -> fetchrow_array)
        {

                $gene2strain{$name} = $strain;

        }

	return \%gene2strain;

}#end get_gene2strain


#sub get_proximate_genes - Takes as input: 1) gene_name		- The name of the gene of interest
#					   2) loc_flag		- L, R or B depending on whether genes left right or both should be retrieved
#					   3) num_genes_L	- The number of genes upstream to retrieve
#					   4) num_genes_R	- The number of genes downstream to retireve
#					   5) chr_id		- An id from the choromosome table
#			 - Returns a hash where the keys are start positions and the values are gene names and annotations
sub get_proximate_genes
{
	my ($gene_name, $loc_flag, $num_genes_L, $num_genes_R, $chr_id, $dbh) = @_;

	my %all_gene_hash;
	my %result_gene_hash;
	
	my $my_start;


	#GET ALL START POSITIONS AND GENE NAMES ON THE GIVEN CHROMOSOME
        my $select_genes = "SELECT DISTINCT g.name, g.start, g.end,  a.annotation, a.source
                            FROM gene g, annotation a
                            WHERE g.gid = a.gid AND
				  g.chr_id = ?";



        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene, $start, $end, $annotation, $source) =  $sth -> fetchrow_array)
        {

		$all_gene_hash{$start}{'GENE'} = $gene;
		$all_gene_hash{$start}{$source} = $annotation;
		$all_gene_hash{$start}{'END'} = $end;
		$all_gene_hash{$start}{'SEQ'} = &extract_sequence($chr_id, $start, $end, $dbh);


		#NOTE THE ORIENTATION OF THE GENE
		if($start < $end)
		{

			$all_gene_hash{$start}{'ORIENTATION'} = "+";
	
		}else{

			$all_gene_hash{$start}{'ORIENTATION'} = "-";

		}#end if


		#NOTE THE START OF INPUTTED GENE
		if ($gene eq $gene_name)
		{

			$my_start = $start;

		}#end if


        }#end while


	#GO THROUGH SORTED LIST AND IDENTIFY DESIRED GENES
	my @left_genes = ();
	my @right_genes = ();

	my $found_gene_flag = 0;
	

	if($loc_flag eq "L")
	{

		foreach my $start( sort {$a<=>$b} keys %all_gene_hash)
		{

			if($start == $my_start)
			{

				last;

			}#end if

			#HAVEN't COME ACROSS MY GENE YET
			if( @left_genes == $num_genes_L)
			{
			
				shift @left_genes;
				push @left_genes, $start;

			}else{

				push @left_genes, $start;

			}#end if


		}#end foreach


	}elsif($loc_flag eq "R"){


		foreach my $start( sort {$a<=>$b} keys %all_gene_hash)
		{


			if( ($found_gene_flag > 0) && ($found_gene_flag <= $num_genes_R) )
			{
			
				#SAVE GENE
				push @right_genes, $start;
				$found_gene_flag ++;

			}elsif($found_gene_flag > 0){

				#FOUND ALL GENES
				last;

			}#end if


			if($start == $my_start)
			{

				$found_gene_flag =1;

			}#end if


		}#end foreach


	}elsif($loc_flag eq "B"){

		foreach my $start( sort {$a<=>$b} keys %all_gene_hash)
		{

                        if($start == $my_start)
                        {

                                $found_gene_flag =1;
                        
                        }#end if



			if( ($found_gene_flag > 0) && ($found_gene_flag <= $num_genes_R) )
			{
			
				#SAVE GENE
				push @right_genes, $start;
				$found_gene_flag ++;

			}elsif($found_gene_flag > 0){

				#FOUND ALL GENES
				last;

			}else{

				#HAVEN'T COME ACROSS MY GENE YET
				if( @left_genes == $num_genes_L)
				{
			
					shift @left_genes;
					push @left_genes, $start;

				}else{

					push @left_genes, $start;

				}#end if

			}#end if


		}#end foreach

	}#end if

	
	#GET RESULTS HASH BASED ON MARKED GENES
	$result_gene_hash{$my_start} = $all_gene_hash{$my_start};


	foreach my $start( (@right_genes,@left_genes) )
	{
	
		$result_gene_hash{$start} = $all_gene_hash{$start};

	}#end foreach


	return \%result_gene_hash;


}#end get_proximate_genes



#sub get_proximate_gene_seqs - Takes as input: 1) gene_name		- The name of the gene of interest
#					   2) loc_flag		- L, R or B depending on whether genes left right or both should be retrieved
#					   3) num_genes_L	- The number of genes upstream to retrieve
#					   4) num_genes_R	- The number of genes downstream to retireve
#					   5) chr_id		- An id from the choromosome table
#			      - Returns a hash where the keys are start positions and the values are gene names and annotations
sub get_proximate_gene_seqs
{
	my ($gene_name, $loc_flag, $num_genes_L, $num_genes_R, $chr_id, $dbh) = @_;


	my %result_gene_hash;
	

	#GET THE START OF MY GENE
        my $select_gene = "SELECT g.start
                           FROM gene g
                           WHERE g.name = ? AND
				 g.chr_id = ?";



        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($gene_name, $chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        my $gene_start =  $sth -> fetchrow_array;

	
	
	#GET INFO FOR GENES ON THE LEFT
        my $select_genes = "SELECT DISTINCT g.name, g.start, g.end
                            FROM gene g
                            WHERE g.chr_id = ? AND
				  g.start <= ? 
				  ORDER BY g.start DESC
				  LIMIT ?";


        $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr_id, $gene_start, ($num_genes_L+1))
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene, $start, $end) =  $sth -> fetchrow_array)
        {

		$result_gene_hash{$start}{'GENE'} = $gene;
		$result_gene_hash{$start}{'END'} = $end;
		$result_gene_hash{$start}{'SEQ'} = &extract_sequence($chr_id, $start, $end, $dbh);


		#NOTE THE ORIENTATION OF THE GENE
		if($start < $end)
		{

			$result_gene_hash{$start}{'ORIENTATION'} = "+";
	
		}else{

			$result_gene_hash{$start}{'ORIENTATION'} = "-";

		}#end if

        }#end while


	#GET INFO FOR GENES ON THE RIGHT
        $select_genes = "SELECT DISTINCT g.name, g.start, g.end
                         FROM gene g
                         WHERE g.chr_id = ? AND
			 g.start > ?
			 ORDER BY g.start ASC
			 LIMIT ?";




        $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr_id, $gene_start, $num_genes_R)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene, $start, $end) =  $sth -> fetchrow_array)
        {

		$result_gene_hash{$start}{'GENE'} = $gene;
		$result_gene_hash{$start}{'END'} = $end;
		$result_gene_hash{$start}{'SEQ'} = &extract_sequence($chr_id, $start, $end, $dbh);


		#NOTE THE ORIENTATION OF THE GENE
		if($start < $end)
		{

			$result_gene_hash{$start}{'ORIENTATION'} = "+";
	
		}else{

			$result_gene_hash{$start}{'ORIENTATION'} = "-";

		}#end if

	}#end while


	return \%result_gene_hash;


}#end get_proximate_gene_seqs



#sub get_intervening_genes - Takes as input: 1) gene1	- The name of a gene
#					     2) gene2	- The name of a second gene
#					     2) chr_id	- An id from the choromosome table
#			      - Returns a hash where the keys are start positions and the values are gene names and annotations of genes
#				located in between the inputted genes
sub get_intervening_genes
{
	my ($gene1, $gene2, $chr_id, $dbh) = @_;

	my @starts;
	my %gene_hash;

	#GET THE START OF MY GENE
        my $select_genes = "SELECT g.start
                                   FROM  gene g
                                   WHERE  g.chr_id = ? ";


        #ADD PLACE HOLDERS FOR GENES                                                 
        $select_genes = $select_genes .  "AND g.name IN (";

        $select_genes = &db_lib::add_place_holders($select_genes , 2);

        $select_genes = $select_genes . ")";

        $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr_id, $gene1, $gene2) 
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($start) =  $sth -> fetchrow_array)
        {

		push @starts, $start;

	}#end while


	#GET INFORMATION FOR GENES IN BETWEN INPUTTED GENES
        $select_genes = "SELECT DISTINCT g.name, g.start, g.end,  a.annotation, a.source
                            FROM gene g, annotation a
                            WHERE g.gid = a.gid AND
				  g.chr_id = ? AND
				  g.start >= ? AND
			          g.start <= ?";



        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id, &arith_ops::min($starts[0], $starts[1]), &arith_ops::max($starts[0], $starts[1]))
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($gene, $start, $end, $annotation, $source) =  $sth -> fetchrow_array)
        {

		$gene_hash{$start}{'GENE'} = $gene;
		$gene_hash{$start}{$source} = $annotation;
		$gene_hash{$start}{'END'} = $end;
		$gene_hash{$start}{'SEQ'} = &extract_sequence($chr_id, $start, $end, $dbh);

	}#end while

	return \%gene_hash;

}#end get_intervening_genes


#sub ncbiID2geneName - Takes as input: 1) seq_id 	- An ncbi sequence ID that contains gene coordinates
#				       2) strain_name	- The name of the strain
#		     - Returns the gene name of the gene with the appropriate coordinaates in the given strain
sub ncbiID2geneName
{
	my ($seq_id, $strain_name,  $dbh) = @_;


	#GET THE COORDINATES FROM THE GENE NAME
	my $seq_coords = $seq_id;
	$seq_coords =~ s/^.*:c?(.*)$/$1/;
	
	my @seq_coords = split /-/, $seq_coords;


	#GET THE GENE NAME FROM THE DATABASE
        my $select_gene = "SELECT g.name
                           FROM chromosome c, organism o, gene g
                           WHERE g.chr_id = c.chr_id AND
                                 c.org_id = o.org_id AND
				 o.strain = ? AND
				 g.start = ? AND
				 g.end = ?";

        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($strain_name, $seq_coords[0], $seq_coords[1])
                         or die "Can't execute statment: $DBI::errstr";


        my ($gene_name) =  $sth -> fetchrow_array;


	return $gene_name;

}#end ncbiID2geneName


#sub get_genes_bordering_contigs - Takes as input:
#				 - Returns a list of genes bordering contigs
sub get_genes_bordering_contigs
{
	my ($dbh) = @_;

	my %contig_genes;
	my %contig_dist;


	#GET DISTANCES OF GENES FROM CONTIG STARTS
        my $select_dist = "SELECT c.contig_id, g.name, c.start, c.end, g.start, g.end
			   FROM gene g, contig c 
			   WHERE g.chr_id = c.chr_id";


        my $sth = $dbh -> prepare($select_dist)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



	#DETERMINE THE GENES AT THE BORDER OF EACH CONTIG
        while(my ($contig_id, $gene, $contig_start, $contig_end, $gene_start, $gene_end) =  $sth -> fetchrow_array)
        {

		if(!exists($contig_dist{$contig_id}))
		{
		
			$contig_dist{$contig_id}{'5p'} = 1000000000;
                        $contig_dist{$contig_id}{'5p_gene'} = "NA";
			$contig_dist{$contig_id}{'3p'} = 1000000000;
                        $contig_dist{$contig_id}{'3p_gene'} = "NA";

		}#end if


		if($contig_start < $contig_end)
		{

			if( (($gene_start - $contig_start) < $contig_dist{$contig_id}{'5p'} && (($gene_start - $contig_start) > 0 || ($gene_end - $contig_start > 0))) )
			{

				$contig_dist{$contig_id}{'5p'} = &arith_ops::min(($gene_start - $contig_start), ($gene_end - $contig_start));
				$contig_dist{$contig_id}{'5p_gene'} = $gene;

			}#end if

			if( (($contig_end - $gene_end) < $contig_dist{$contig_id}{'3p'} && (($contig_end - $gene_end) > 0 || ($contig_end - $gene_start > 0))))
			{

				$contig_dist{$contig_id}{'3p'} = &arith_ops::min(($contig_end - $gene_end), ($contig_end - $gene_start));
				$contig_dist{$contig_id}{'3p_gene'} = $gene;

			}#end if

		}else{


			if( (($gene_start - $contig_end) < $contig_dist{$contig_id}{'5p'} && (($gene_start - $contig_end) > 0 || ($gene_end - $contig_end > 0))) )
			{

				$contig_dist{$contig_id}{'5p'} = &arith_ops::min(($gene_start - $contig_end), ($gene_end - $contig_end));
				$contig_dist{$contig_id}{'5p_gene'} = $gene;

			}#end if

			if( (($contig_start - $gene_end) < $contig_dist{$contig_id}{'3p'} && (($contig_start - $gene_end) > 0 || ($contig_start - $gene_start > 0))))
			{

				$contig_dist{$contig_id}{'3p'} = &arith_ops::min(($contig_start - $gene_end), ($contig_start - $gene_start));
				$contig_dist{$contig_id}{'3p_gene'} = $gene;

			}#end if


		}#end if


	}#end while


	#GET GENES THAT BORDER A CONTIG
	foreach my $contig_id (sort {$a<=>$b} keys %contig_dist)
	{

		$contig_genes{$contig_dist{$contig_id}{'3p_gene'}}{'3p_C'} = $contig_id;	
		$contig_genes{$contig_dist{$contig_id}{'3p_gene'}}{'3p_D'} = $contig_dist{$contig_id}{'3p'} ;	

		$contig_genes{$contig_dist{$contig_id}{'5p_gene'}}{'5p_C'} = $contig_id;	
		$contig_genes{$contig_dist{$contig_id}{'5p_gene'}}{'5p_D'} = $contig_dist{$contig_id}{'5p'} ;	

	}#end freach 


	return \%contig_genes;

}#end get_genes_bordering_contigs


#sub is_gene_at_contig_barrier - Takes as input: 1) gene_name	- The name of a gene
#						 2) dist	- Distance from contig barrier
#			       - Returns a 0 if gene is internal to any contig/scaffold and 1 otherwise. Definition
#				 of internal is determined by dist
sub is_gene_at_contig_barrier
{
	my ($gene_name, $dist, $dbh) = @_;


	#GET GENE START AND END
        my $select_gene = "SELECT g.start, g.end
                           FROM gene g
                           WHERE g.name= ?";



        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($gene_name)
                         or die "Can't execute statment: $DBI::errstr";


	my ($gene_start, $gene_end) = $sth -> fetchrow_array;



	#GET ALL CONTIG STARTS AND ENDS ON THE CHROMOSOME WITH THE GENE
        my $select_contig = "SELECT c.start, c.end
                             FROM gene g, contig c
                             WHERE g.chr_id = c.chr_id AND
                                   g.name= ?
			     ORDER BY c.start ASC";



        $sth = $dbh -> prepare($select_contig)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($gene_name)
                         or die "Can't execute statment: $DBI::errstr";


	
	 while(my ($contig_start, $contig_end) =  $sth -> fetchrow_array)
        {

		#CHECK IF GENE FITS LEGALLY IN CURRENT CONTIG

		#CHECK I ON MIMUS OR PLUS STRAND
		if($gene_start < $gene_end)
		{
		
			#CHECK IF CONTIG IS ON PLUS OR MINUS STRAIND
			if($contig_start < $contig_end)
			{

				if( (($gene_start - $contig_start + 1) > $dist) && (($contig_end - $gene_end + 1) > $dist) )
				{
	
					return 0;

				}#end if

			}else{

				if( (($gene_start - $contig_end + 1) > $dist) && (($contig_start - $gene_end + 1) > $dist) )
				{
	
					return 0;

				}#end if

			}#end if


		}else{

			 #CHECK IF CONTIG IS ON PLUS OR MINUS STRAIND
                        if($contig_start < $contig_end)
                        {

				if( (($gene_end - $contig_start + 1) > $dist) && (($contig_end - $gene_start + 1) > $dist) )
                        	{

                        	        return 0;

                        	}#end if

			}else{

				if( (($gene_end - $contig_end + 1) > $dist) && (($contig_start - $gene_start + 1) > $dist) )
                        	{

                        	        return 0;

                        	}#end if

			}#end if

		}#end if

	}#end while


	#GENE IS NOT INTERNAL TO ANY CONTIG
	return 1;

}#end is_gene_at_contig_barrier


#sub is_gene_at_contig_barrier2 - Takes as input: 1) gene_name	- The name of a gene
#						 2) dist	- Distance from contig barrier
#			       - Returns a 0 if gene is internal to any contig and 1 otherwise. Definition
#				 of internal is determined by dist.
sub is_gene_at_contig_barrier2
{
	my ($gene_name, $dist, $dbh) = @_;


	#GET GENE START AND END
        my $select_gene = "SELECT g.start, g.end
                           FROM gene g
                           WHERE g.name= ?";



        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($gene_name)
                         or die "Can't execute statment: $DBI::errstr";


	my ($gene_start, $gene_end) = $sth -> fetchrow_array;



	#GET ALL CONTIG STARTS AND ENDS ON THE CHROMOSOME WITH THE GENE
        my $select_contig = "SELECT c.start, c.end
                             FROM gene g, contig c
                             WHERE g.chr_id = c.chr_id AND
                                   g.name= ?
			     ORDER BY c.start ASC";



        $sth = $dbh -> prepare($select_contig)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($gene_name)
                         or die "Can't execute statment: $DBI::errstr";


	
	 while(my ($contig_start, $contig_end) =  $sth -> fetchrow_array)
        {

		#CHECK IF GENE FITS LEGALLY IN CURRENT CONTIG

		#CHECK I ON MIMUS OR PLUS STRAND
		if($gene_start < $gene_end)
		{
		
			#CHECK IF CONTIG IS ON PLUS OR MINUS STRAIND
			if($contig_start < $contig_end)
			{

				if( (($gene_start < $contig_end && $gene_start > $contig_start) || ($gene_end < $contig_end && $gene_end > $contig_start)) && 
				    ((($gene_start - $contig_start + 1) < $dist) || (($contig_end - $gene_end + 1) < $dist)) )
				{
	
					return 1;

				}#end if

			}else{

				if( (($gene_start > $contig_end && $gene_start < $contig_start) || ($gene_end > $contig_end && $gene_end < $contig_start)) && 
				    ((($gene_start - $contig_end + 1) < $dist) || (($contig_start - $gene_end + 1) < $dist)) )
				{
	
					return 1;

				}#end if

			}#end if


		}else{

			 #CHECK IF CONTIG IS ON PLUS OR MINUS STRAIND
                        if($contig_start < $contig_end)
                        {

				if( (($gene_start < $contig_end && $gene_start > $contig_start) || ($gene_end < $contig_end && $gene_end > $contig_start)) &&
				     ((($gene_end - $contig_start + 1) < $dist) || (($contig_end - $gene_start + 1) < $dist)) )
                        	{

                        	        return 1;

                        	}#end if

			}else{

				if( (($gene_start > $contig_end && $gene_start < $contig_start) || ($gene_end > $contig_end && $gene_end < $contig_start)) &&
				     ((($gene_end - $contig_end + 1) < $dist) || (($contig_start - $gene_start + 1) < $dist)) )
                        	{

                        	        return 1;

                        	}#end if

			}#end if

		}#end if

	}#end while


	#GENE IS NOT INTERNAL TO ANY CONTIG
	return 0;

}#end is_gene_at_contig_barrier



#sub is_dubious_gene - Takes as input: 1) gene 		- THe name of a gene from the database
#				       2) min_amb	- The minimum fraction of ambiguous bases to call gene dubious
#		     - Returns 0 if gene looks OK, and otherwise a description of what is wrong with the gene
sub is_dubious_gene
{
	my ($gene, $min_amb, $dbh) = @_;

	#GET SEQUENCE OF GENE
	my @genes = ($gene);

	my $r_gene_seq_hash = &get_gene_sequences(\@genes, $dbh);
	

	#THROW FLAG IF IT IS TOO SHORT
        my $gene_length = length($r_gene_seq_hash->{$gene});

	if($gene_length < 150)
	{

		return "GENE IS SMALLER THAN 100 BP";

	}#end if


	#THROW FLAG IF IT CONTAINS TOO MANY AMBIGUOUS BASES
        my $gene_seq_obj = Bio::Seq->new(-seq => $r_gene_seq_hash->{$gene});

        my $r_gene_nuc_counts = &seq_stats::init_monoNucHash();
        $r_gene_nuc_counts = &seq_stats::get_monoNucCount($r_gene_nuc_counts, $gene_seq_obj);
        my %gene_nuc_counts = %$r_gene_nuc_counts;

        my $N_count = $gene_length - ($gene_nuc_counts{'A'} + $gene_nuc_counts{'C'} + $gene_nuc_counts{'G'}  + $gene_nuc_counts{'T'});

	if(($N_count / $gene_length) > $min_amb)
	{

		return "GENE IS MORE THAN 50% AMBIGUOUS BASES";

	}#end if


	#PASSED FILTERS
	return "OK";

}#end is_dubious_gene


############################################################ GENE PROPERTY SELECTION ROTUINES #############################################


#sub is_gene_phage - Takes as input: 1) r_genes - A reference to a list of genes in the database
#		   - Returns a hash where the keys are genes that are either in phage regions or have phage in their description
sub is_gene_phage
{
	my ($r_genes, $dbh) = @_;

	my $r_phage_genes = &hash_lib::init_hash($r_genes);
	my %phage_genes = %$r_phage_genes;


	#IF NO GENES GIVEN, THEN RETURN EMPTY LIST
	if(@$r_genes == 0)
	{	

		return $r_phage_genes;

	}#end if


	#DETERMINE IF REGION OVERLAPS WITH PHAGE REGION
	my $select_phage_region = "SELECT g.name
                            	   FROM mobile_element m, gene g
				   WHERE m.class LIKE '%phage%' AND 
      					 m.chr_id = g.chr_id AND
					 (( m.end < g.start  AND g.start < m.start) OR
     		 			 (m.start < g.start AND g.start < m.end)  ) ";


        #ADD PLACE HOLDERS FOR GENES                                                 
        $select_phage_region = $select_phage_region .  "AND g.name IN (";

        $select_phage_region = &db_lib::add_place_holders($select_phage_region , scalar(@$r_genes));

        $select_phage_region = $select_phage_region . ")";


	#EXECUTE QUERY
        my $sth = $dbh -> prepare($select_phage_region)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";


	#PUT GENES OVERLAPPING PHAGE INTO HASH
	while(my $phage_gene = $sth -> fetchrow_array)
	{

		$phage_genes{$phage_gene} = 1;

	}#end while


	#DETERMINE IF GENES HAVE 'PHAGE' IN THEIR NR DESCRIPTION
        my $select_phage_genes = "SELECT g.name
                                  FROM gene g, annotation a
				  WHERE a.gid = g.gid AND
                                        a.source = 'NR' AND
				        a.annotation LIKE '%phage%' ";


        #ADD PLACE HOLDERS FOR GENES                                                 
        $select_phage_genes = $select_phage_genes .  "AND g.name IN (";

        $select_phage_genes = &db_lib::add_place_holders($select_phage_genes , scalar(@$r_genes));

        $select_phage_genes = $select_phage_genes . ")";


        $sth = $dbh -> prepare($select_phage_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";


        while(my $phage_gene = $sth -> fetchrow_array)
        {

		$phage_genes{$phage_gene} = 1;

        }#end while


	return \%phage_genes;

}#end is_gene_phage


#sub get_gene_GC - Takes as input: 1) r_genes	- A list of genes
#		 - Returns the GC content of those genes from the gene table
sub get_gene_GC
{
	my ($r_genes, $dbh) = @_;

	my %gc_hash;


	#GET THE GC OF GIVEN GENES
	my $select_gc =   "SELECT g.name, g.gc
                            FROM gene g
                            WHERE ";


	#ADD PLACE HOLDERS FOR GENES				
	$select_gc = $select_gc . "g.name IN (";

	$select_gc = &db_lib::add_place_holders($select_gc , scalar(@$r_genes));

       	$select_gc = $select_gc . ")";


        $sth = $dbh -> prepare($select_gc)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";

	while(my ($gene, $gc) = $sth -> fetchrow_array)
	{

		$gc_hash{$gene} = $gc;

	}#end while


	return \%gc_hash;

}#end get_gene_GC


#sub get_gene_length - Takes as input: 1) gene	- A gene name
#		      - Returns the length of the gene
sub get_gene_length
{
	my ($gene, $dbh) = @_;


	#GET THE LENGTH OF THE GIVEN GENES
	my $select_length =   "SELECT ABS(g.start - g.end) 
                               FROM gene g
                               WHERE g.name = ?";


        $sth = $dbh -> prepare($select_length)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($gene)
                         or die "Can't execute statment: $DBI::errstr";

	my $length = $sth -> fetchrow_array;

	return ($length+1);

}#end get_gene_length


#sub get_gene_lengths - Takes as input: 1) r_genes	- A list of genes
#		      - Returns the length of the genes
sub get_gene_lengths
{
	my ($r_genes, $dbh) = @_;

	my %length_hash;


	#GET THE LENGTH OF THE GIVEN GENES
	my $select_length =   "SELECT g.name, ABS(g.start - g.end) 
                               FROM gene g
                               WHERE ";


	#ADD PLACE HOLDERS FOR GENES				
	$select_length = $select_length . "g.name IN (";

	$select_length = &db_lib::add_place_holders($select_length , scalar(@$r_genes));

       	$select_length = $select_length . ")";


        $sth = $dbh -> prepare($select_length)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";

	while(my ($gene, $length) = $sth -> fetchrow_array)
	{

		$length_hash{$gene} = $length + 1;

	}#end while


	return \%length_hash;

}#end get_gene_lengths


#sub get_gene_diNuc - Takes as input: 1) r_genes	- A list of genes
#		    - Returns the di-nucleotide signature of those genes from the gene table
sub get_gene_diNuc
{
	my ($r_genes, $dbh) = @_;

	my %di_nuc_hash;


	#GET THE DI-NUC SIG OF GIVEN GENES
	my $select_dns =   "SELECT g.name, g.di_nuc_signature
                            FROM gene g
                            WHERE ";


	#ADD PLACE HOLDERS FOR GENES				
	$select_dns = $select_dns . "g.name IN (";

	$select_dns = &db_lib::add_place_holders($select_dns , scalar(@$r_genes));

       	$select_dns = $select_dns . ")";


        $sth = $dbh -> prepare($select_dns)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";

	while(my ($gene, $di_nuc_sig) = $sth -> fetchrow_array)
	{

		$di_nuc_hash{$gene} = $di_nuc_sig;

	}#end while


	return \%di_nuc_hash;

}#end get_gene_diNuc


#sub is_diNuc_gt1SD - Takes as input: 1) r_genes	- A list of genes on a single chromsome
#				      2) chr_id		- The ID of the chromosome that the genes are on
#		    - Returns a hash where the keys are those genes that have a di-nucleotide signature
#		      greater than 1 SD above the mean for genes on the chromosme
sub is_diNuc_gt1SD
{
	my ($r_genes, $chr_id, $dbh) = @_;

	my %di_nuc_hash;


	#DETERMINE THE MEAN AND STANDARD DEVIATION OF DI-NUCLEOTIDE SIGNATURE ON THE CHROMSOME
	my $select_di_nuc_sd =   "SELECT AVG(di_nuc_signature), STD(di_nuc_signature)
                   	          FROM gene g
                         	  WHERE g.chr_id = ?";



        $sth = $dbh -> prepare($select_di_nuc_sd)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";


        my ($di_nuc_mean, $di_nuc_sd) = $sth -> fetchrow_array;


	#DETERMINE THE NUMBER OF GENES IN REGION WHOSE SIGNATURE IS GREATER THAN 1 SD FROM MEAN
	my $select_gt1sd =   "SELECT g.name
                              FROM gene g
                              WHERE g.di_nuc_signature > ?";


	#IF GENES IN REGIONS THEN ADD PLACE HOLDERS FOR THEM							
	$select_gt1sd = $select_gt1sd . "AND g.name IN (";

	$select_gt1sd = &db_lib::add_place_holders($select_gt1sd , scalar(@$r_genes));

       	$select_gt1sd = $select_gt1sd . ")";


	#EXECTURE QUERY
        $sth = $dbh -> prepare($select_gt1sd)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute( ($di_nuc_mean + $di_nuc_sd) , @$r_genes)
                         or die "Can't execute statment: $DBI::errstr";


	while(my $gene = $sth -> fetchrow_array)
	{

		$di_nuc_hash{$gene} = 1;

	}#end while


	return \%di_nuc_hash

}#end is_diNuc_gt1SD


#sub is_gene_MGE - Takes as input: 1) r_genes	- A list of genes
#		  - Returns a hash where the keys are genes whose annotation contains a keyword for a mge
sub is_gene_MGE
{
	my ($r_genes , $dbh) = @_;

	my $r_mge_hash = &hash_lib::init_hash($r_genes);
	my %mge_hash = %$r_mge_hash;


	#IF NO GENES GIVEN, THEN RETURN EMPTY LIST
        if(@$r_genes == 0)
        {

                return $r_mge_hash;

        }#end if


	#GET GENE ANNOTATIONS
	my $r_gene2annots = &get_gene_annotations($r_genes, $dbh);
	%gene2annots = %$r_gene2annots;


	#CHECK IF CURRENT GENE BELONGS TO GIVEN CLASSES, AND IF SO START COUNT OVER
	foreach my $gene(@$r_genes)
	{

		if( $gene2annots{$gene}{'NR'} =~ /intI|integron|integrase|transposase|recombinase/i || $gene2annots{$gene}{'CDD'} =~ /intI|integron|integrase|transposase|recombinase/i || $gene2annots{$gene}{'NCBI'} =~ /intI|integron|integrase|transposase|recombinase/i )
		{

			$mge_hash{$gene} = 1;

		}#end if

	}#end foreach


	return \%mge_hash;

}#end is_gene_MGE


#sub get_gene_numSNPs - Takes as input: 1) r_genes - A list of genes
#		      - Returns a 2D hash where the first key is a gene name, the second key is a organism name
#			and the value is the number of SNPs in the given gene compared with the given chromosome
sub get_gene_numSNPs
{
	my ($r_genes, $dbh) = @_;

	my $snp_hash;


	#DETERMINE THE NUMBER OF SNPS IN EACH GENE WITH OTHER CHROMSOMES 
	my $select_gene_snps =  "SELECT g.name, o1.strain, o2.strain,  COUNT(*)
                                 FROM gene g, snp s, chromosome c1, chromosome c2, organism o1, organism o2
                                 WHERE ((g.start < s.position1 AND g.end > s.position1 AND g.chr_id = s.chr_id1 AND c1.chr_id = s.chr_id1 AND c2.chr_id = s.chr_id2) OR
				       (g.start > s.position1 AND g.end < s.position1 AND g.chr_id = s.chr_id1 AND c1.chr_id = s.chr_id1 AND c2.chr_id = s.chr_id2) OR
				       (g.start < s.position2 AND g.end > s.position2 AND g.chr_id = s.chr_id2 AND c1.chr_id = s.chr_id2 AND c2.chr_id = s.chr_id1) OR 
				       (g.start > s.position2 AND g.end < s.position2 AND g.chr_id = s.chr_id2 AND c1.chr_id = s.chr_id2 AND c2.chr_id = s.chr_id1)) AND
				        s.base1 <> \'.\' AND
				  	s.base2 <> \'.\' AND
					s.base1 <> \'N\' AND
					s.base2 <> \'N\' AND
					c1.org_id = o1.org_id AND
					c2.org_id = o2.org_id\n";


	#IF GENES IN REGIONS THEN ADD PLACE HOLDERS FOR THEM							
	$select_gene_snps = $select_gene_snps . "AND g.name IN (";

	$select_gene_snps = &db_lib::add_place_holders($select_gene_snps , scalar(@$r_genes));

       	$select_gene_snps = $select_gene_snps . ")";

	$select_gene_snps = $select_gene_snps . "\nGROUP BY g.name, o1.strain, o2.strain";


	#EXECTURE QUERY
        $sth = $dbh -> prepare($select_gene_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";


	while(my ($gene, $org1, $org2, $num_snps) = $sth -> fetchrow_array)
	{

		if(exists($snp_hash{$gene}{$org2}))
		{

			$snp_hash{$gene}{$org2} += $num_snps;

		}else{

			$snp_hash{$gene}{$org2} = $num_snps;
			$snp_hash{$gene}{$org1} = 0;

		}#end if

	}#end while


	return \%snp_hash;

}#end get_gene_numSNPs


#sub get_gene_dN_dS - Takes as input: 1) r_genes	- A list of genes
#		    - Returns a 3D hash where the first key is the gene, the second an org, and the third either S,N,nS,nN. Note
#		      that the statistics are the sum across all genes from a given org in the query genes orthologous group
sub get_gene_dN_dS
{
	my ($r_genes, $dbh) = @_;

	my %dnds_hash;


	#GET THE MSA FOR EACH GENE
	my $select_og_msa =   "SELECT g.name, o.strain, og.gene_names, og.org_ids, og.msa
                               FROM gene g, ortholog_group og, organism o, chromosome c
                               WHERE og.gene_names LIKE concat('%', g.name, '%') AND
				     g.chr_id = c.chr_id AND
				     o.org_id = c.org_id ";


	#ADD PLACE HOLDERS FOR GENES			
	$select_og_msa = $select_og_msa . "AND g.name IN (";

	$select_og_msa = &db_lib::add_place_holders($select_og_msa , scalar(@$r_genes));

       	$select_og_msa = $select_og_msa . ")";


	#EXECTURE QUERY
        my $sth = $dbh -> prepare($select_og_msa)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute(@$r_genes)
                         or die "Can't execute statment: $DBI::errstr";


	while(my ($q_gene, $q_org, $genes, $orgs, $msa) = $sth -> fetchrow_array)
	{
	
		#IF ONLY ONE GENE THEN SKIP
		my @orgs = split /\|/, $orgs;
		my @genes = split /\|/, $genes;

		if(@orgs == 1 || $msa =~ /X/)
		{
		
			next;

		}#end if

	
		#BREAK UP MSA AND CREATE ALIGNMENT OBJECT
		my @msa = split />/, $msa;
		shift @msa;

		#GO THROUGH EACH SEQUENCE AND CREATE SEQUENCE OBJECTS
		my @seq_objs;

		foreach my $entry(@msa)
		{

        		my ($id, @seq) = split /\n/, $entry;

        		my $seq = join "", @seq;

			my $seq_obj = Bio::LocatableSeq->new(-seq => $seq,
                    					     -id  => $id,
                    					     -start => 1,
                    					     -end   => length($seq));

		
			push @seq_objs, $seq_obj;

		}#end foreaach
	

		#MAKE ALIGNMENT OBJECT AND CONVERT TO NUCLEOTIDE
		my $r_dna_seq_hash = &get_gene_sequences(\@genes, $dbh);

		my $aa_aln_obj = Bio::SimpleAlign->new(-seqs => \@seq_objs);

		my $dna_aln_obj = &my_Utilities::aa_to_dna_aln($aa_aln_obj, $r_dna_seq_hash);


		#GET DNA SEQUENCE STATS FOR CURRENT MSA	
		my $stats = Bio::Align::DNAStatistics->new();


		#CALCULATE dN/dS RELATIVE TO EACH ORTHOLOGOUS GENE
		for(my $i = 0;$i <(@genes); $i ++)
		{

			my $gene = $genes[$i];
			my $org = $orgs[$i];


			#CHECK IF GENE IS QUERY, IF SO SKIP
			if($gene eq $q_gene)
			{
				next;

			}#end if


			#CHECK IF ANOTHER GENE FROM SAME ORG IN OG, IF SO SKIP CURRENT MSA
			if($org eq $q_org)
			{

				last;

			}#end if


			#COMPUTE dN/dS FOR CURRENT GENE PAIR
			my $r_dnds_results = $stats-> calc_KaKs_pair($dna_aln_obj, $q_gene, $gene);
			my %dnds_results = %{@$r_dnds_results[0]};


			#GET STATS, BUT IF MULTIPLE GENES FROM ORG IN OG THEN DON'T COUNT
			# S_d - Number of synonymous mutations between the 2 sequences.
    			# N_d - Number of non-synonymous mutations between the 2 sequences.
     			# S - Mean number of synonymous sites in both sequences.
     			# N - mean number of synonymous sites in both sequences.
     			# P_s - proportion of synonymous differences in both sequences given by P_s = S_d/S.
     			# P_n - proportion of non-synonymous differences in both sequences given by P_n = S_n/S.
     			# D_s - estimation of synonymous mutations per synonymous site (by Jukes-Cantor).
     			# D_n - estimation of non-synonymous mutations per non-synonymous site (by Jukes-Cantor).
     			# D_n_var - estimation of variance of D_n .
     			# D_s_var - estimation of variance of S_n.
     			# z_value - calculation of z value.Positive value indicates D_n > D_s, negative value indicates D_s > D_n.
			if(!exists($dnds_hash{$q_gene}{$org}))
			{	

				$dnds_hash{$q_gene}{$org}{'Sd'} = $dnds_results{'S_d'}; 
				$dnds_hash{$q_gene}{$org}{'Nd'} = $dnds_results{'N_d'}; 
				$dnds_hash{$q_gene}{$org}{'NdpSd'} = $dnds_results{'N_d'} + $dnds_results{'S_d'};
				$dnds_hash{$q_gene}{$org}{'dS'} = $dnds_results{'D_s'}; 
				$dnds_hash{$q_gene}{$org}{'dN'} = $dnds_results{'D_n'}; 

				if($dnds_results{'D_s'} > 0)
				{

					$dnds_hash{$q_gene}{$org}{'dNdS'} = $dnds_results{'D_n'} / $dnds_results{'D_s'};
				}else{

					$dnds_hash{$q_gene}{$org}{'dNdS'} = "NaN";

				}#end if

			}else{

				$dnds_hash{$q_gene}{$org}{'Sd'} = "NaN"; 
				$dnds_hash{$q_gene}{$org}{'Nd'} = "NaN"; 
				$dnds_hash{$q_gene}{$org}{'NdpSd'} = "NaN"; 
				$dnds_hash{$q_gene}{$org}{'dS'} = "NaN"; 
				$dnds_hash{$q_gene}{$org}{'dN'} = "NaN"; 

			}#end if

		}#end foreaach

	}#end while


	return \%dnds_hash;

}#end get_gene_dN


#sub get_DNAsim_from_AAaln - Takes as input: 1) genes	- A list of gene names
#			   - Returns the average AA/DNA similarity between all pairs of genes based upon an amino acid alignment
sub get_DNAsim_from_AAaln
{
	my ($r_genes, $dbh) = @_;


	#GET GENE SEQUENCES, TRANSLATE AND CREATE SEQUENCE OBJECTS
        my $r_dna_seq_hash = &get_gene_sequences($r_genes, $dbh);
	my %dna_seq_hash = %$r_dna_seq_hash;

        my @aa_seq_objs;
        my @dna_seq_objs;
	my %dna_seq_hash_mod;
	my @len;

        foreach my $gene (@$r_genes)
        {

		#MAKE UNIQUE ID FOR EACH GENE, IN CASE GENES ARE BEING ALIGNED AGAINST THEMSELVES
		$c ++;
		my $id = $gene . "_$c";
		$id =~ s/\:/_/;
		$dna_seq_hash_mod{$id} = $dna_seq_hash{$gene};

		#TRANSLATE SEQUENCE AND CREATE OBJECT
		my $aa_seq = &seq_ops::translate($dna_seq_hash{$gene});

        	my $aa_seq_obj = Bio::Seq->new(-seq => $aa_seq,
                                               -display_id => $id);

        	my $dna_seq_obj = Bio::Seq->new(-seq => $dna_seq_hash{$gene},
                                                -display_id => $id);

		#print length($aa_seq) . "\n";
		push @len, length($aa_seq);

                push @aa_seq_objs, $aa_seq_obj;
                push @dna_seq_objs, $dna_seq_obj;

        }#end foreaach


        #MAKE ALIGNMENT OBJECT AND CONVERT TO NUCLEOTIDE
        my $factory = new Bio::Tools::Run::Alignment::Clustalw(-quiet => 1);

        #my $dna_aln_obj = $factory->align(\@dna_seq_objs);

        my $aa_aln_obj = $factory->align(\@aa_seq_objs);

	my $aa2dna_aln_obj = &my_Utilities::aa_to_dna_aln($aa_aln_obj, \%dna_seq_hash_mod);

	
	#DETERMINE SIMILARITY AT AA AND NUC LEVEL
	$aa_id = $aa_aln_obj->percentage_identity;
	#$nuc_id = $dna_aln_obj->percentage_identity;
	$nuc_id = $aa2dna_aln_obj->percentage_identity;
	#print "1: " . $dna_aln_obj->percentage_identity . "vs" . $aa2dna_aln_obj->percentage_identity . "vs" . $aa_aln_obj->percentage_identity;

	return ($aa_id, $nuc_id);
	
}#end get_DNAsim_from_AAdna_aln_obj

############################################################# CHROMSOME SELECTION ROUTINES ###################################################

#sub chr_names2ids - Takes as input:
#		   - Returns a hash where the keys are chromsome names and the values are their IDs in the 
#		     chromsome table
sub chr_names2ids
{
	my ($dbh) = @_;

        my %chr_name2id;


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_chr = "SELECT c.chr_id, c.name
                       FROM chromosome c";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($chr_id, $chr_name) =  $sth -> fetchrow_array)
        {

                $chr_name2id{$chr_name} = $chr_id;

        }#end while


        return \%chr_name2id;


}#end chr_names2ids


#sub org_names2ids - Takes as input:
#		   - Returns a hash where the keys are organism names and the values are their IDs in the 
#		     organism table
sub org_names2ids
{
	my ($dbh) = @_;

        my %org_name2id;


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_org = "SELECT o.org_id, o.strain
                       FROM organism o";


        my $sth = $dbh -> prepare($select_org)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($org_id, $org_name) =  $sth -> fetchrow_array)
        {

                $org_name2id{$org_name} = $org_id;

        }#end while


        return \%org_name2id;

}#end org_names2ids


#sub get_org_info - Takes as input: 
#		  - Returns a 2D hash where the first key is the strain name and the second keys are attributes of the assembly
sub get_org_info
{
	my ($dbh) = @_;

	my %org_hash;

        #GET ALL STRAINS AND ASSOCIATED INFO
        $select_org = "SELECT strain, N50_contig, num_contigs, avg_contig_size, largest_contig_size, N50_scaffold, num_scaffolds
                       FROM organism";


        my $sth = $dbh -> prepare($select_org)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
	my @cols = ('strain', 'N50_contig', 'num_contigs', 'avg_contig_size', 'largest_contig_size', 'N50_scaffold', 'num_scaffolds');

        while(my @row  =  $sth -> fetchrow_array)
        {

		for(my $i = 1; $i < (@cols); $i ++)
		{

			$org_hash{$row[0]}{$cols[$i]} = $row[$i];

		}#end for

	}#end while

	return \%org_hash;

}#end get_org_info


#sub get_chr_info - Takes as input:
#		  - Returns a hash where the keys are chromosome ids and the values are the chromsome name, class and gc content
sub get_chr_info
{
	my ($dbh) = @_;


	my %chr_hash;


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_chr = "SELECT c.chr_id, c.name, c.sequence, c.class
                       FROM chromosome c";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($chr_id, $chr_name, $chr_sequence, $chr_class) =  $sth -> fetchrow_array)
        {

		$chr_hash{$chr_id}{'NAME'} = $chr_name;
		$chr_hash{$chr_id}{'CLASS'} = $chr_class;

		
		my $chr_obj = Bio::Seq->new(-seq => $chr_sequence);

        	my $r_chr_nuc_counts = &seq_stats::init_monoNucHash();
        	$r_chr_nuc_counts = &seq_stats::get_monoNucCount($r_chr_nuc_counts, $chr_obj);
		

        	$chr_GC = &seq_stats::get_GC($r_chr_nuc_counts);

		$chr_hash{$chr_id}{'GC'} = $chr_GC;

        }#end while


        return \%chr_hash;

}#end get_chr_info


#sub get_chr_seq - Takes as input: 1) chr_name - The name of the chromsome of interest
#		 - Returns the sequence of the chromosome
sub get_chr_seq
{
	my ($chr_name , $dbh) = @_;


        #GET SEQUENCE OF GIVEN CHROMSOME
        $select_chr = "SELECT c.sequence
                       FROM chromosome c
		       WHERE c.name = ?";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_name)
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        my $chr_seq =  $sth -> fetchrow_array;

	return $chr_seq;

}#end if


#get_chr_desc - Takes as input: 1) chr_id - An ID from the chromosme database
#	      - Returns a 2D hash where the first key is the chromosome name and
#		the second keys are various descriptors of the chromosome 
sub get_chr_desc
{
	my ($chr_id, $dbh) = @_;


	my %chr_hash;


        #GET INFO ABOUT DESIRED CHROMSOME
        $select_chr = "SELECT c.name, c.sequence, c.class, o.strain, o.species
                       FROM chromosome c, organism o
		       WHERE c.chr_id = ? AND
			     c.org_id = o.org_id";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($chr_name, $chr_sequence, $chr_class, $chr_strain,  $chr_species) =  $sth -> fetchrow_array)
        {

		$chr_hash{'NAME'} = $chr_name;
		$chr_hash{'SEQ'} = $chr_sequence;
		$chr_hash{'CLASS'} = $chr_class;
		$chr_hash{'STRAIN'} = $chr_strain;
		$chr_hash{'SPECIES'} = $chr_species;

        }#end while


        return \%chr_hash;


}#end get_chr_desc


#get_chr_desc_all - Takes as input: 1) chr_id - An ID from the chromosme database
#	          - Returns a 2D hash where the first key is the chromosome name and
#		     the second keys are various descriptors of the chromosome 
sub get_chr_desc_all
{
	my ($dbh) = @_;


	my %chr_hash;


        #GET ALL GENES PRESENT IN DESIRED STRAIN
        $select_chr = "SELECT c.chr_id, c.name, c.sequence, c.class, o.strain, o.species
                       FROM chromosome c, organism o
		       WHERE c.org_id = o.org_id";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($chr_id,$chr_name, $chr_sequence, $chr_class, $chr_strain,  $chr_species) =  $sth -> fetchrow_array)
        {

		$chr_hash{$chr_id}{'NAME'} = $chr_name;
		$chr_hash{$chr_id}{'SEQ'} = $chr_sequence;
		$chr_hash{$chr_id}{'CLASS'} = $chr_class;
		$chr_hash{$chr_id}{'STRAIN'} = $chr_strain;
		$chr_hash{$chr_id}{'SPECIES'} = $chr_species;

        }#end while


        return \%chr_hash;

}#end get_chr_desc_all


#sub get_chr_id - Takes as input: 1) chr_name	- The name od a chromsome from the chromsome table
#		- Returns the id corresponding to the name in the table
sub get_chr_id
{
	my ($chr_name, $dbh) = @_;

	#EXTRACT THE SEQUENCE
        my $select_chr = "SELECT chr_id
                          FROM chromosome
                          WHERE name = ?";


        my $sth = $dbh -> prepare($select_chr)
                          or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_name)
                          or die "Can't execute statment: $DBI::errstr";



        my $chr_id  =  $sth -> fetchrow_array;


	return $chr_id;

}#end get_chr_id


#sub get_chr_id_all - Takes as input: 
#		    - Returns a hash linking chromsome names to ids
sub get_chr_id_all
{
	my ($dbh) = @_;

	my %chr2id;

	#EXTRACT THE NAMES AND IDs
        my $select_chr = "SELECT chr_id, name
                          FROM chromosome";


        my $sth = $dbh -> prepare($select_chr)
                          or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                          or die "Can't execute statment: $DBI::errstr";



        while(my ($id, $name ) =  $sth -> fetchrow_array)
	{

		$chr2id{$name} = $id;

	}#end while


	return \%chr2id;

}#end get_chr_id_all


#sub get_chr_orient - Takes as input: 1) chr		- The name of a chromsoome
#				      2) position	- A position in the chromosome
#				      3) msa_id		- An ID of a multiple sequence alignment
#
#		    - Returns a hash where the key is the chromsome name and the value is whether
#		      the chromsome is aligned in forward or reverse complement orientation relative
#		      to the inputted chromsome and the inputted position
sub get_chr_orient
{
	my ($chr, $position, $msa_id, $dbh) = @_;

	my %chr_orient_hash;

	#GET CHROMOSOMES IN MSA
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";



        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs =  $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;


	#GO THROUGH CHROMSOMES AND GET INDEX OF REFERENCE CHROMSOME
	my $my_chr_i = -1000;

	for(my $i=0; $i<(@chrs); $i++)
	{

		if( $chrs[$i] eq $chr )
		{

			$my_chr_i = $i;
			last;

		}#end if

	}#end for


        #GET ALL REGIONS IN THE MSA TABLE WITH THE GIVEN ID
        my $select_regions = "SELECT ma.starts, ma.ends
                              FROM multiple_alignment ma, multiple_alignment_info mai
                              WHERE ma.msa_id = mai.msa_id AND
				    mai.msa_id = ?";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($starts, $ends) =  $sth -> fetchrow_array)
	{

		#CHECK IF ALIGNED REGION ENCOMPASSES POSITION OF INTEREST 
		my @starts = split /\|/, $starts;
		my @ends = split /\|/, $ends;

		if( (abs($starts[$my_chr_i]) > $position) || (abs($ends[$my_chr_i]) < $position) )  
		{

			next;

		}#end if


		for(my $i=0; $i<(@chrs); $i++)
		{

			if(($starts[$i] * $starts[$my_chr_i]) > 0)
			{

				$chr_orient_hash{$chrs[$i]} = 1;

			}else{

				$chr_orient_hash{$chrs[$i]} = -1;

			}#end if

		}#end for

		last;

	}#emd while


	return \%chr_orient_hash;

}#end get_chr_orient


#sub get_chr_close_to_ref - Takes as input: 1) ref_chr	 	- The genome relative to which other genomes should be identified 
#					    2) perc_ID		- The cutoff for percent identity to the reference by which other genomes are selected to align
#
#			  - Returns a list of chromosomes who are within the given identity threshold over shared sequence
sub get_chr_close_to_ref
{
	my ($ref_chr, $perc_ID_cutoff, $length_cutoff, $dbh) = @_;

        #GET CHROMOSOMES IN DATABASE
        my $r_chr_hash = &chr_names2ids($dbh);


        #GET PAIRWISE SNP COUNTS
        my $r_pw_snp_counts = &get_pw_snp_counts($dbh);
        my %pw_snp_counts = %$r_pw_snp_counts;

        my $r_chr2id = &get_chr_id_all($dbh);
        my %chr2id = %$r_chr2id;


	#GET LENGTH OF REFERENCE SEQUENCE
	my $r_ref_chr_info = &get_chr_desc($chr2id{$ref_chr}, $dbh);
	my $ref_length = length($r_ref_chr_info -> {'SEQ'});


        #REMOVE CHROMOSOMES WITHOUT SNPS IN THE DATABASE
        foreach my $chr(keys %$r_chr_hash)
        {

                #GET FRACTION OF CHROMOSOMES ALIGNED TO REFERENCE
                $r_aln_pos = &get_aln_pos_pw($ref_chr, $chr, $dbh);
                %aln_pos = %$r_aln_pos;

                my $perc_id;

                if(keys(%aln_pos) > 0)
                {

                        $perc_id = &arith_ops::max(($pw_snp_counts{$chr2id{$chr}}{$chr2id{$ref_chr}} /  keys(%aln_pos)), ($pw_snp_counts{$chr2id{$ref_chr}}{$chr2id{$chr}}/ keys(%aln_pos)));

                }else{ 

                        $perc_id = 1;

                }#end if

                #THROW OUT CHROMOSOMES WITH GREATER THAN 5% NUCLEOTIDE IDENTITY TO REFERENCE
                print "$chr : aln pos = " . keys(%aln_pos) . ", num snps = " . $pw_snp_counts{$chr2id{$ref_chr}}{$chr2id{$chr}} . ", percent ID = " . $perc_id .  "\n";

                if( (!exists($pw_snp_counts{$chr2id{$chr}}{$chr2id{$ref_chr}}) && !exists($pw_snp_counts{$chr2id{$ref_chr}}{$chr2id{$chr}})) || ($perc_id > $perc_ID_cutoff) || (abs($ref_length - keys(%aln_pos)) > $length_cutoff) )
                {

			print "\tDELETE $chr\n";
                        delete $r_chr_hash->{$chr};

                }#end if

        }#end foreach


        #ADD REFERENCE CHROMOSOME TO BEGINNING
        @chrs = ($ref_chr, sort keys %$r_chr_hash);


	return \@chrs;

}#end get_chr_close_to_ref


#sub get_org_chr - Takes as input: 1) strain	- The name of the strain whose chromosome sequences are desired
#		 - Returns a hash linking chromosome names to sequences
sub get_org_chr
{
	my ($strain, $dbh) = @_;

	my %chr2seq;


        #GET CHROMOSOMES AND THEIR SEQUENCES
        $select_chr = "SELECT c.chr_id, c.sequence, c.class, c.name
                       FROM chromosome c, organism o
                       WHERE c.org_id = o.org_id AND
                             o.strain = ?";


        my $sth = $dbh -> prepare($select_chr)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($strain)
                         or die "Can't execute statment: $DBI::errstr";


        #GET A LIST OF ALL THE GENES
        while(my ($chr_id, $seq, $class, $chr_name) =  $sth -> fetchrow_array)
        {
	
                $chr2seq{$chr_id}{'SEQ'} = $seq;
                $chr2seq{$chr_id}{'NAME'} = $chr_name;

        }#end while
	

	return \%chr2seq;

}#end get_org_chr


#sub get_org_gene_count	- Takes as input: strains	- THe name of a strain in the organism database
#			- Returns the number of protein coding genes associated with the strain in the gene database
sub get_org_gene_count
{
	my ($strain, $dbh) = @_;

        #SELECT THE NUMBER OF GENES FOR THE GIVEN STRAIN
        my $select_genes = "SELECT COUNT(*)
                            FROM organism o, chromosome c, gene g
                            WHERE o.strain = ? AND
			    c.org_id = o.org_id AND
			    g.chr_id = c.chr_id AND
			    g.type = 'protein'";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($strain)
                         or die "Can't execute statment: $DBI::errstr";


	my $gene_count = $sth -> fetchrow_array;


	return $gene_count;	

}#end get_org_gene_count


#sub get_chr_gene_summary - Takes as input: 1) chr_id - An ID from the chromosme table
#			  - Returns a 2D hash where the first keys are the names of genes on a chromosome
#			    and the second keys are attributes of the genes
sub get_chr_gene_summary
{
	my ($chr_id, $dbh) = @_;

	my %gene_hash;


	#SELECT ATTRIBUTES FOR GENES ON GIVEN CHROMOSOME
        my $select_genes = "SELECT name, start, end, type
                            FROM gene
                            WHERE chr_id = ?";

        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while( my ($gene_name, $start, $end, $type) =  $sth -> fetchrow_array)
	{

		#ASSIGN START/END/STRAND DEPENDING ON STRAND WHICH GENE RESIDES ON
		if($start < $end)
		{
			$gene_hash{$gene_name}{'START'} = $start;
			$gene_hash{$gene_name}{'END'} = $end;
			$gene_hash{$gene_name}{'STRAND'} = 1;

		}else{

			$gene_hash{$gene_name}{'START'} = $end;
			$gene_hash{$gene_name}{'END'} = $start;
			$gene_hash{$gene_name}{'STRAND'} = -1;

		}#end if


		#ASSIGN SEQUENCE TYPE TP BE BE GENBANK STANDARD
		if($type eq 'protein')
		{

			$gene_hash{$gene_name}{'TYPE'} = 'CDS';
			

		}else{

			$gene_hash{$gene_name}{'TYPE'} = 'gene';

		}#end if


	}#end while


	return \%gene_hash;


}#end get_chr_gene_summary


#sub chr2contig_coords - Takes as input: 1) chr_id	- The name of a chromosome in the chomosome table
#					 2) start	- The start position of the sequence of interest in chr
#					 3) end		- The end position of the sequence of interest in chr
#		       - Returns the name of a contig which is mapped to the given chromosomal region, and the coordinates in the contig of 
#			 the sequence of interest
sub chr2contig_coords
{
	my ($chr_id, $start, $end, $dbh) = @_;


	#GET THE CONTIG WHICH COMPRISES THE GIVEN CHROMOSOMAL SEQUENCE
	my $select_coords = "SELECT c.contig_name, c.start, c.end, c.description
                             FROM contig c
                             WHERE c.chr_id = ? AND
				   ( (c.start <= ? AND c.start <= ? AND c.end >= ? AND c.end >= ?) OR (c.start >= ? AND c.start >= ? AND c.end <= ? AND c.end <= ?)) AND
				   (contig_name NOT LIKE 'scaffold%')";


	my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id, $start, $end, $start, $end, $start, $end, $start, $end)
                         or die "Can't execute statment: $DBI::errstr";


	my ($contig, $contig_chr_start, $contig_chr_end, $contig_desc) =  $sth -> fetchrow_array;
	#print "($contig, $contig_chr_start, $contig_chr_end, $contig_desc)\n";

	#MAP COORDINATES RELATIVE TO CONTIG
	my ($contig_start, $contig_end, $contig_strand);	

	
	if( $contig =~ /\w+/)
	{

		#CHECK STRAND WHICH CONTIG MAPS TO
		if($contig_chr_start < $contig_chr_end)
		{

			$contig_start = $start - $contig_chr_start + 1;
			$contig_end = $end - $contig_chr_start + 1;
			$contig_strand = 1;

		}else{

			$contig_start = $contig_chr_start - $start + 1;
			$contig_end = $contig_chr_start - $end + 1;
			$contig_strand = -1;

		}#end if

	}else{

		$contig_start = 0;
		$contig_end = 0;
		$contig = "Spans multiple contigs";
		$contig_strand = 0;
		$contig_desc = "NA";

	}#end	


	return ($contig, $contig_start, $contig_end, $contig_strand, $contig_desc);

}#end chr2contig_coords


#sub chr2contigs - Takes as input: 1) chr_id	- The id of a chromosome in the chomosome table
#				   2) start	- The start position of the sequence of interest in chr
#				   3) end	- The end position of the sequence of interest in chr
#		       - Returns the name of a contig which is mapped to the given chromosomal region, and the coordinates in the contig of 
#			 the sequence of interest
sub chr2contigs
{
	my ($chr_id, $start, $end, $dbh) = @_;
	my %contig_coords;

	
	#MAKE SURE THAT START IS LESS THAN END
	if($end < $start)
	{

		my $temp = $start;
		$start = $end;
		$end = $temp;

	}#end if


	#GET THE CONTIG WHICH COMPRISES THE GIVEN CHROMOSOMAL SEQUENCE
	my $select_coords = "SELECT c.contig_name, c.start, c.end, c.description
                             FROM contig c
                             WHERE c.chr_id = ? AND
				   ( (c.start >= ?  AND c.start <= ?) OR  (c.end >= ?  AND c.end <= ?) OR (c.start <= ? AND c.end >= ?) OR (c.start >= ? AND c.end <= ?))";


	my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id, $start, $end, $start, $end, $start, $start, $start, $start)
                         or die "Can't execute statment: $DBI::errstr";



        while( my ($contig, $contig_chr_start, $contig_chr_end, $contig_desc) =  $sth -> fetchrow_array)
	{
	

		#MAP COORDINATES RELATIVE TO CONTIG
		if(1)
		#if( $contig =~ /^contig/)
		{

			#CHECK STRAND WHICH CONTIG MAPS TO
			if( ($contig_chr_start < $contig_chr_end) && ($contig_chr_start <= $start) )
			{

				$contig_coords{$contig}{'START'}= $start - $contig_chr_start + 1;
				$contig_coords{$contig}{'END'} = &arith_ops::min( ($end - $contig_chr_start + 1) , ($contig_chr_end - $contig_chr_start + 1));
				$contig_coords{$contig}{'STRAND'} = 1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}elsif(($contig_chr_start < $contig_chr_end) && ($contig_chr_start > $start) && ($contig_chr_end <= $end)){

				$contig_coords{$contig}{'START'}= 1;
				$contig_coords{$contig}{'END'} = $contig_chr_end - $contig_chr_start + 1;
				$contig_coords{$contig}{'STRAND'} = 1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}elsif(($contig_chr_start < $contig_chr_end) && ($contig_chr_start > $start) && ($contig_chr_end >= $end) ){

				$contig_coords{$contig}{'START'}= 1;
				$contig_coords{$contig}{'END'} = $end - $contig_chr_start + 1;
				$contig_coords{$contig}{'STRAND'} = 1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}elsif( ($contig_chr_start > $contig_chr_end) && ($contig_chr_start >= $end) ){

				$contig_coords{$contig}{'START'}= $contig_chr_start - $end + 1;
				$contig_coords{$contig}{'END'} = &arith_ops::min( ($contig_chr_start - $start + 1) , ($contig_chr_start - $contig_chr_end + 1));
				$contig_coords{$contig}{'STRAND'} = -1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}elsif(($contig_chr_start > $contig_chr_end) && ($contig_chr_end >= $start) && ($contig_chr_start < $end)){

				$contig_coords{$contig}{'START'}= 1;
				$contig_coords{$contig}{'END'} = $contig_chr_start - $contig_chr_end + 1;
				$contig_coords{$contig}{'STRAND'} = -1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}elsif(($contig_chr_start > $contig_chr_end) && ($contig_chr_end <= $start) && ($contig_chr_start < $end) ){

				$contig_coords{$contig}{'START'}= 1;
				$contig_coords{$contig}{'END'} = $contig_chr_start - $start + 1;
				$contig_coords{$contig}{'STRAND'} = -1;
				$contig_coords{$contig}{'DESC'} = $contig_desc;

			}#end if

		}#end	

	}#end while


	return (\%contig_coords);

}#end chr2contigs


#sub chr2scaffold_coords - Takes as input: 1) chr_id	- The name of a chromosome in the chomosome table
#					   2) start	- The start position of the sequence of interest in chr
#					   3) end	- The end position of the sequence of interest in chr
#		       - Returns the name of a scaffold which is mapped to the given chromosomal region, and the coordinates in the scaffold of 
#			 the sequence of interest
sub chr2scaffold_coords
{
	my ($chr_id, $start, $end, $dbh) = @_;


	#GET THE SCAFFOLD WHICH COMPRISES THE GIVEN CHROMOSOMAL SEQUENCE
	my $select_coords = "SELECT c.contig_name, c.start, c.end, c.description
                             FROM contig c, chromosome chr
                             WHERE c.chr_id = ? AND
				   ( (c.start <= ? AND c.end  >= ? AND c.start <= ? AND c.end >= ?) OR (c.end <= ? AND c.start  >= ? AND c.end <= ? AND c.start >= ?)) AND
				   contig_name LIKE 'scaffold%'";


	my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id, $start, $end, $end, $start, $start, $end, $end, $start)
                         or die "Can't execute statment: $DBI::errstr";



        my ($scaffold, $scaffold_chr_start, $scaffold_chr_end, $scaffold_desc) =  $sth -> fetchrow_array;
	

	#MAP COORDINATES RELATIVE TO SCAFFOLD
	my ($scaffold_start, $scaffold_end, $scaffold_strand);	

	#print "$scaffold, $scaffold_chr_start, $scaffold_chr_end, $chr, $start, $end\n";

	if( $scaffold =~ /^scaffold/)
	{

		#CHECK STRAND WHICH SCAFFOLD MAPS TO
		if($scaffold_chr_start < $scaffold_chr_end)
		{

			$scaffold_start = $start - $scaffold_chr_start + 1;
			$scaffold_end = $end - $scaffold_chr_start + 1;
			$scaffold_strand = 1;

		}else{

			$scaffold_start = $scaffold_chr_start - $start + 1;
			$scaffold_end = $scaffold_chr_start - $end + 1;
			$scaffold_strand = -1;

		}#end if

	}else{


		$scaffold_start = 0;
		$scaffold_end = 0;
		$scaffold = "Spans multiple contigs";
		$scaffold_strand = 0;
		$scaffold_desc = "NA";

	}#end	


	return ($scaffold, $scaffold_start, $scaffold_end, $scaffold_strand, $scaffold_desc);

}#end chr2scaffold_coords


#sub contig2chr_coords - Takes as input:
#		       - Returns a hash linking chromosomes:contigs to start, end and strand information
sub contig2chr_coords
{
	my ($dbh) = @_;
	
	my %contig_hash;


	#SELECT FROM CONTIG TABLE
	my $select_coords = "SELECT chr.name, c.start, c.end, c.contig_name
                             FROM contig c, chromosome chr
                             WHERE c.chr_id = chr.chr_id";


	my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while( my ($chr, $start, $end, $contig) =  $sth -> fetchrow_array)
	{

		$contig_hash{"$chr:$contig"}{'START'} = $start;
		$contig_hash{"$chr:$contig"}{'END'} = $end;

	}#end while
	

	return \%contig_hash;

}#end contig2chr_coords


#sub contig2chr - Takes as input: 
#		- Returns a hash linking contig names to chromosome names
sub contig2chr
{
	my ($dbh) = @_;
	
	my %contig_hash;


	#SELECT FROM CONTIG TABLE
	my $select_coords = "SELECT chr.name, c.contig_name
                             FROM contig c, chromosome chr
                             WHERE c.chr_id = chr.chr_id";


	my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while( my ($chr, $contig) =  $sth -> fetchrow_array)
	{

		$contig_hash{$contig} = $chr;

	}#end while
	

	return \%contig_hash;


}#end contig2chr


#sub get_num_contigs - Takes as input:
#		     - Returns a hash linking chromosome name to number of contigs
sub get_num_contigs
{
	my ($dbh) = @_;

	my %num_contigs;


	#SELECT FROM CONTIG TABLE
	my $select_contigs = "SELECT chr.name, COUNT(*)
                              FROM contig c, chromosome chr
                              WHERE c.chr_id = chr.chr_id
		              GROUP BY chr.name";


	my $sth = $dbh -> prepare($select_contigs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while( my ($chr, $num_contigs) =  $sth -> fetchrow_array)
	{

		$num_contigs{$chr} = $num_contigs;

	}#end while

	return \%num_contigs;

}#end get_num_contigs


#################################################################### SEQUENCE SELECTION ROUTINES #############################################

#sub extract_sequence - Takes as input: 1) chr_id	- A chr_id from the database
#					2) start	- The sequence start position
#					3) enb		- The sequence end position
#		      - Returns the desired sequence
sub extract_sequence
{
	my ($chr_id, $start, $end, $dbh) = @_;

	my $seq;


	#DETERMINE IF REV COM NEEDS TO BE TAKES
	if($end >= $start)
	{

        	#EXTRACT THE SEQUENCE
        	my $select_seq = "SELECT SUBSTRING(sequence, ?, ?)
        	                  FROM chromosome
        	                  WHERE chr_id = ?";


        	my $sth = $dbh -> prepare($select_seq)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($start, ($end-$start+1), $chr_id)
        	                 or die "Can't execute statment: $DBI::errstr";



        	$seq  =  $sth -> fetchrow_array;


	}else{

		#EXTRACT THE SEQUENCE
                my $select_seq = "SELECT SUBSTRING(sequence, ?, ?)
                                  FROM chromosome
                                  WHERE chr_id = ?";


                my $sth = $dbh -> prepare($select_seq)
                                 or die "Can't prepare statment: $DBI::errstr";

                my $rc = $sth -> execute($end, ($start-$end+1), $chr_id)
                                 or die "Can't execute statment: $DBI::errstr";



                $seq  =  $sth -> fetchrow_array;

		my $seq_obj = Bio::Seq->new(-seq => $seq);
		$rev_seq_obj = $seq_obj->revcom;

		$seq = $rev_seq_obj->seq;

	}#end if


	return $seq;

}#end extract_sequence


#################################################################### GENOME PROPERTIES SELECTION ROTUINES  #############################################


#sub get_gene_seq_positions - Takes as input: 1) chr - The name of a chromsome in the database
#			    - Returns a hash whose keys are all positions on the chromsome that
#			      are in genes
sub get_gene_seq_positions
{
	my ($chr, $dbh) = @_;

	my %gene_coord_hash;


	#SELECT THE START AND END COORDINATES OF EACH GENE
        my $select_coords = "SELECT start, end
                             FROM gene g, chromosome c
                             WHERE g.chr_id = c.chr_id AND
				   c.name = ?";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
        {

		my $first = &arith_ops::min($start, $end);
		my $last = &arith_ops::max($start, $end);


		for(my $p = $first; $p <= $last; $p++)
		{

			$gene_coord_hash{$p} = 0;

		}#end 

        }#end while


	return \%gene_coord_hash;

}#end get_gene_seq_positions


#sub get_amplicons	- Takes as input: 1) primers - A file of primers
#			- Returns a hash with the results of running ePCR using the input primers
#			  against the database of chromosomes
sub get_amplicons
{
	my ($primers, $dbh) = @_;


	#CREATE A FASTA FILE OF GENOMES IN THE DATABASE
	my $r_chr_desc_hash = &get_chr_desc_all($dbh);
	my %chr_desc_hash = %$r_chr_desc_hash;
	
	my %chr_hash;
	
	foreach my $chr_id (keys %chr_desc_hash)
	{
	
	        $chr_hash{$chr_desc_hash{$chr_id}{'NAME'}} = $chr_desc_hash{$chr_id}{'SEQ'};
	
	}#end foreach
	
	&seq_io::hash2fasta('temp_chr.fa', \%chr_hash);
	
	
	#PERFORM E-PCR
	`e-PCR -w9 -f 1 -m100 -u+ $primers D=100-400 temp_chr.fa V=+ O=temp_pcr N=1 G=1 T=3`;
	
	
	#GET THE NAMES OF THE PRIMERS
	my $r_primer_hash = &hash_lib::dict_to_hash($primers);
	my @primers = sort keys(%$r_primer_hash);
	
	
	#PARSE OUTPUT FILE
	my %pcr_results;
	
	open PCR, "temp_pcr";                                                                                                                                                                                                     
	                                                                                                                                                                                                                          
	foreach my $line(<PCR>)                                                                                                                                                                                                   
	{                                                                                                                                                                                                                         
	                                                                                                                                                                                                                          
	        #SPLIT UP LINE                                                                                                                                                                                                    
	        chomp $line;                                                                                                                                                                                                      
	                                                                                                                                                                                                                          
	        my ($chr, $primer, $strand, $start, $end, $size, @junk) = split /\t/, $line;                                                                                                                                      
	        my @size = split /[-\/]/, $size;                                                                                                                                                                                  
	                                                                                                                                                                                                                          
	                                                                                                                                                                                                                          
	        #MAKE SURE AMPLICON MATCHES EXPECTED SIZE                                                                                                                                                                         
	        if(abs($size[1] - $size[2]) < 10)                                                                                                                                                                                 
	        {                                                                                                                                                                                                                 
	                                                                                                                                                                                                                          
	                $pcr_results{$chr}{$primer} = "$start-$end";                                                                                                                                                                          
	                                                                                                                                                                                                                          
	        }else{                                                                                                                                                                                                            
	                                                                                                                                                                                                                          
	                print "$chr/$primer SIZE DISCREPANCY ($size[1] vs. $size[2])!!!\n\n";                                                                                                                                     
	                $pcr_results{$chr}{$primer} = "$start-$end";                                                                                                                                                                          
	                                                                                                                                                                                                                          
	        }#end if                                                                                                                                                                                                          
	                                                                                                                                                                                                                          
	}#end foreach                                                                                                                                                                                                             
	close PCR;  
	
	#`rm temp_pcr`;
	`rm temp_chr.fa`;
	
	
	return \%pcr_results;

}#end get_amplicons


#sub get_gene_4fold_degenerate_positions - Takes as input: 1) chr - The name of a chromsome in the database
#			                 - Returns a hash whose keys are all positions on the chromsome that
#			                   are 4fold degenerate positions in genes
sub get_gene_4fold_degenerate_positions
{
	my ($chr, $dbh) = @_;

	my %four_fold_pos_hash;


	#GET SEQUENCES OF GENES ON CHROMSOME
	my $r_chr_genes = &get_chr_genes(&get_chr_id($chr, $dbh), $dbh);

        my $r_gene2seq = &get_gene_sequences($r_chr_genes, $dbh);
	my %gene2seq = %$r_gene2seq;

	my $r_gene_coords = &get_gene_coords($r_chr_genes, $dbh);
	my %gene_coords = %$r_gene_coords;


	#GO THROUGH EACH SEQUENCE AND GET 4-FOLD DEGENERATE SITES
	foreach my $gene(keys %gene2seq)
	{

		#GET 4-FOLD POSITIONS IN GENE
		my $r_four_fold_pos_in_gene  = &seq_stats::get_4fold_degenerate_positions($gene2seq{$gene});

		#CONVERT TO CHROMSOMAL COORDINATES
		if($gene_coords{$gene}{'START'} < $gene_coords{$gene}{'END'} )
		{

			foreach my $pos(@$r_four_fold_pos_in_gene)
			{

				$four_fold_pos_hash{($pos + $gene_coords{$gene}{'START'} - 1)} = 0;

			}#end foreach

		}else{

			foreach my $pos(@$r_four_fold_pos_in_gene)
			{

				$four_fold_pos_hash{($gene_coords{$gene}{'START'} - $pos + 1)} = 0;

			}#end foreach

		}#end if


	}#end foreach


	return \%four_fold_pos_hash;


}#end get_gene_4fold_degen_positions



############################################################ ORTHOLOGY SELECTION ROUTINES ##########################################################

#sub get_singleton_genes - Takes as input: 1) r_strains		- Reference to a list of strains for which singleton genes should
#								  be retrieved
#			 - Returns a list of singlteon genes
sub get_singleton_genes
{
	my ($r_strains, $dbh) = @_;

	my @all_genes;
	my %og_genes;
	my %singleton_genes;


	#GET ALL GENES IN GIVEN STRAINS
	foreach my $strain (@$r_strains)
	{

		my $r_strain_genes = &get_org_genes($strain, $dbh);	

		push @all_genes, @$r_strain_genes;

	}#end foreach


	#GET ALL GENES IN ORTHOLOG GROUPS
        my $select_ogs = "SELECT gene_names
                       FROM `ortholog_group`"; 

        my $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";


        #EXTRACT THE GENES FROM THE DESIRED ORGANISM
        while(my ($genes) =  $sth -> fetchrow_array)
        {

		$genes =~ s/^\|//;
                my @genes = split /\|/, $genes;


		foreach my $gene (@genes)
		{

                	$og_genes{$gene} = 1;

		}#end foreach

        }#end while


	#COMPILE LIST OF GENES FROM STRAINS NOT IN ORTHOLOG GROUP
	foreach my $gene (@all_genes)
	{

		if( !exists($og_genes{$gene}) )
		{

			push @singleton_genes, $gene;

		}#end if

	}#end foreach


	return \@singleton_genes;	

}#end get_singleton_genes


#sub get_ortholog_groups - Takes as input:
#			 - Returns a list of all the gene name strings from the ortholog_group table for the 
#			   desired organism
sub get_ortholog_groups
{
	my ($dbh) = @_;

	my %ortholog_groups;


	#SELECT ORTHOLOG GROUPS
        my $select_ogs = "SELECT og_id, gene_names
			  FROM ortholog_group";



        my $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($og_id, $og) =  $sth -> fetchrow_array)
        {

		$ortholog_groups{$og_id} = $og;


        }#end while


	return \%ortholog_groups;

}#end get_ortholog_groups


#sub get_gene_ortholog_groups - Takes as input: 
#               	      - Returns a hash linking gene names to the ortholog groups to which they belong
sub get_gene_ortholog_groups
{
        my ($dbh) = @_;

        my %gene2og;


        #SELECT ORTHOLOG GROUPS
        my $select_ogs = "SELECT og_id, gene_names
                          FROM ortholog_group";



        my $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($og_id, $og) =  $sth -> fetchrow_array)
        {

		$og =~ s/^\|//;
		my @genes = split /\|/, $og;

		foreach my $gene(@genes)
		{

                	$gene2ogs{$gene} = $og_id;

		}#end foreach

        }#end while


        return \%gene2ogs;

}#end get_gene_ortholog_groups


#sub get_ortholog_group_orgs - Takes as input: 
#                            - Returns a string of strain IDs present in each orthologous group
#			       keyed by the table id
sub get_ortholog_group_orgs
{
        my ($dbh) = @_;

        my %ortholog_groups;


        #SELECT ORTHOLOG GROUPS
        my $select_ogs = "SELECT og_id, org_ids
                          FROM ortholog_group";



        my $sth = $dbh -> prepare($select_ogs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($og_id, $org_ids) =  $sth -> fetchrow_array)
        {

                $ortholog_groups{$og_id} = $org_ids;


        }#end while


        return \%ortholog_groups;

}#end get_ortholog_group_orgs


############################################################# GENOMIC REGION SELECTION ROUTINES ###################################################

#sub get_msa_id	- Takes as input: 1) chrs - 	A reference to a list of chromosomes

#		- Returns an msa_id from the multiple_alignment_info table containing ALL chromosomes in the list
sub get_msa_id
{
	my ($r_chrs, $dbh) = @_;
	my @chrs = @$r_chrs;


	#GET ALIGNMENTS FROM DATABASE
        my $select_msa = "SELECT msa_id, chrs
                          FROM multiple_alignment_info
			  ORDER BY msa_id ASC";


        my $sth = $dbh -> prepare($select_msa)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";



	#GO THROUGH ALIGNMENTS AND SELCT THE ONE CONTAINING ALL DESIRED CHROMOSOMES, AND AS FEW OTHERS AS POSSIBLE
	my $select_aln = -1;
	my $select_chr_count = 100000;

        ALN: while(my ($msa_id, $chrs) =  $sth -> fetchrow_array)
	{

		my @aln_chrs = split /\|/, $chrs;

		#VERIFY THAT ALL CHROMOSOMES ARE IN THE LIST
		foreach my $chr (@chrs)
		{

			if(($chrs !~ /\|$chr/) && ($chrs !~ /$chr\|/) && ($chrs !~ /\|$chr\|/))
			{
					
				next ALN;

			}#end if

		}#end foreach
	

		#IF FEWER CHROMOSOMES IN THIS ALIGNMENT, THEN GO WITH IT
		if(@aln_chrs <= $select_chr_count)
		{

			$select_aln = $msa_id;
			$select_chr_count = @aln_chrs;

		}#end if

	}#end while


	return $select_aln;

}#end get_msa_id


#sub get_unique_regions - Takes as input: 1) chr1	- An ID from the chromosome table
#					  2) chr2	- An ID from the chromsoome table
#					  3) min_size	- The minimum size to consider a unique region
#
#			- Returns a hash where the keys are start and end coordinates of regions unique to
#			  chr1 that are greater than min_size
sub get_unique_regions
{
	my ($chr1, $chr2, $min_size, $dbh) = @_;

	my %unique_regions;
	my $region_counter = 1;


        #GET ALL REGIONS UNIQUE TO CHR1, WHERE CHR1 IS THE FIRST ID IN THE TABLE
        my $select_regions = "SELECT o.start1, o.end1
                              FROM pairwise_alignment o
                              WHERE o.chr_id1 = ? AND
			 	    o.chr_id2 = ? AND
				    o.start2 = 0 AND
				    o.end2 = 0 AND
				    (o.end1 - o.start1) > ?
                              ORDER BY o.start1 ASC";



        my $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2, $min_size)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		$unique_regions{$region_counter}{'START'} = $start;
		$unique_regions{$region_counter}{'END'} = $end;

		$region_counter ++;

	}#end while


        #GET ALL REGIONS UNIQUE TO CHR1, WHERE CHR1 IS THE SECOND ID IN THE TABLE
        $select_regions = "SELECT o.start2, o.end2
                              FROM pairwise_alignment o
                              WHERE o.chr_id2 = ? AND
			 	    o.chr_id1 = ? AND
				    o.start1 = 0 AND
				    o.end1 = 0 AND
                                    (o.end2 - o.start2) > ?
                              ORDER BY o.start2 ASC";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr1, $chr2, $min_size)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		$unique_regions{$region_counter}{'START'} = $start;
		$unique_regions{$region_counter}{'END'} = $end;

		$region_counter ++;

	}#end while


	return \%unique_regions;

}#end get_unique_regions

#sub get_unique_bases - Takes as input: 1) chr1		- The ID of a chromsome
#					2) chr2		- The ID of a chromosome
#					3) min_size	- The minimum size to consider a unique region
#
#			- Returns a hash where the keys are bases unique to chromosome 1
sub get_unique_bases
{
	my ($chr1, $chr2, $min_size, $dbh) = @_;

	my %aln_bases;
	my %unique_bases;

	#GET THE LENGTH OF THE CHROMSOME
        my $r_chr_desc_hash = &get_chr_desc($chr1, $dbh);
        my $genome_length = length($r_chr_desc_hash->{'SEQ'});


        #GET ALL REGIONS UNIQUE TO CHR1, WHERE CHR1 IS THE FIRST ID IN THE TABLE
        my $select_regions = "SELECT o.start1, o.end1
                              FROM pairwise_alignment o
                              WHERE o.chr_id1 = ? AND
			 	    o.chr_id2 = ? AND
				    (o.end1 - o.start1) > ?
                              ORDER BY o.start1 ASC";



        my $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2, $min_size)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		if($start < $end)
		{

			@aln_bases{$start..$end} = 1;

		}else{

			@aln_bases{$end..$start} = 1;

		}#end if

	}#end while


        #GET ALL REGIONS UNIQUE TO CHR1, WHERE CHR1 IS THE SECOND ID IN THE TABLE
        $select_regions = "SELECT o.start2, o.end2
                              FROM pairwise_alignment o
                              WHERE o.chr_id2 = ? AND
			 	    o.chr_id1 = ? AND
                                    (o.end2 - o.start2) > ?
                              ORDER BY o.start2 ASC";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr1, $chr2, $min_size)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		if($start < $end)
		{

			@aln_bases{$start..$end} = 1;

		}else{

			@aln_bases{$end..$start} = 1;

		}#end if

	}#end while

	#DETERMINE UNIQUE BASES
	my @genome = 1..$genome_length;

	my $r_unique_bases = &hash_lib::init_hash(\@genome);
	%unique_bases = %$r_unique_bases;

	delete @unique_bases{keys(%aln_bases)};


	return \%unique_bases;

}#end get_unqiue_bases


#sub get_diff_genomic_regions - Takes as input: 1) chr_list1	- A list of chromsome names
#						2) chr_list2	- A second list of chrosomsome names
#						3) min_size	- The minimum sized region to consider
#
#			      - Returns the set of regions from the multiple_alignment table that are
#				present in chr_list1 and not in chr_list2 that are at least min_size
sub get_diff_genomic_regions
{
	my ($r_chr_list1, $r_chr_list2, $min_size, $msa_id, $dbh) = @_;

	my %diff_regions;
	my $region_counter = 1;

	my $r_chr_in_hash = &hash_lib::init_hash($r_chr_list1);
	my %chr_in_hash = %$r_chr_in_hash;
	my @chr_in_inds;

	my $r_chr_out_hash = &hash_lib::init_hash($r_chr_list2);
	my %chr_out_hash = %$r_chr_out_hash;
	my @chr_out_inds;


	#GET CHROMOSOMES IN MSA
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";



        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs =  $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;


	#GO THROUGH CHROMSOMES AND NOTE INDICIES OF IN AND OUT CHROMSOMES
	for(my $i=0; $i<(@chrs); $i++)
	{

		if( exists($chr_in_hash{$chrs[$i]}) )
		{

			push @chr_in_inds, $i;

		}elsif( exists($chr_out_hash{$chrs[$i]}) ){

			push @chr_out_inds, $i;

		}#end if

	}#end for


        #GET ALL REGIONS IN THE MSA TABLE WITH THE GIVEN ID
        my $select_regions = "SELECT ma.starts, ma.ends
                              FROM multiple_alignment ma, multiple_alignment_info mai
                              WHERE ma.msa_id = mai.msa_id AND
				    mai.msa_id = ?";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($starts, $ends) =  $sth -> fetchrow_array)
	{

		my @starts = split /\|/, $starts;
		my @ends = split /\|/, $ends;

		#CHECK IF REGION PRESENT IN CORRECT CHROMSOMES AND GREATER THAN MINIMUM SIZE
		if((join("|", @starts[@chr_in_inds]) =~ /^(\|?\-?[1-9]\d*\|?)+$/) && (join("|", @starts[@chr_out_inds]) =~ /^(\|?0\|?)+$/) && (abs($ends[$chr_in_inds[0]] - $starts[$chr_in_inds[0]]) > $min_size))
		{


			#GET INFO ABOUT IN AND OUT CHROMSOMES
			my @out_inds = ();
			my @in_inds = ();

			for(my $i=0; $i < (@starts); $i++)
			{

				if($starts[$i] == 0)
				{

					push @out_inds, $i;

				}else{

					push @in_inds, $i;

					#CHECK STRAND, AND ADJUST START AND END FROM MAUVE FORMAT, TO MY OWN
					if($starts[$i] < 0)
					{

						my $temp = $starts[$i];
						$starts[$i] = abs($ends[$i]);
						$ends[$i] = abs($temp);

					}#end if

				}#end if

			}#end for	


			#STORE INFO
			$diff_regions{$region_counter}{'OUT_CHRS'} = join " , ", @chrs[@out_inds];

			$diff_regions{$region_counter}{'IN_CHRS'} = join " , ", @chrs[@in_inds];
			$diff_regions{$region_counter}{'IN_STARTS'} = @starts[@in_inds];
			$diff_regions{$region_counter}{'IN_ENDS'} = @ends[@in_inds];

			$diff_regions{$region_counter}{'CHR'} = $chrs[$chr_in_inds[0]];
			$diff_regions{$region_counter}{'START'} = $starts[$chr_in_inds[0]];
			$diff_regions{$region_counter}{'END'} = $ends[$chr_in_inds[0]];

			$region_counter ++;

		}#end if

	}#end while
	
	return \%diff_regions;

}#end get_diff_genomic_regions


#sub get_group_diff_genomic_regions - Takes as input: 1) chr_list1	- A list of chromsome names
#						      2) min_size	- The minimum sized region to consider
#
#			      - Returns the set of regions from the multiple_alignment table that are
#				not present in all of the inputed chromosomes
sub get_group_diff_genomic_regions
{
	my ($r_chr_list, $min_size, $msa_id, $dbh) = @_;

	my @chr_list = @$r_chr_list;

	my %diff_regions;
	my %region_seqs;
	my %chr_starts;
	my $region_counter = 1;

	my $r_chr_hash = &hash_lib::init_hash($r_chr_list);
	my %chr_hash = %$r_chr_hash;

	my @chr_inds;
	my %chr2ind;


	#GET CHROMOSOMES IN MSA
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";



        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs =  $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;


	#GO THROUGH CHROMSOMES AND NOTE INDICIES INPUTTED CHROMSOMES
	for(my $i=0; $i<(@chrs); $i++)
	{

		if( exists($chr_hash{$chrs[$i]}) )
		{

			push @chr_inds, $i;
			$chr2ind{$chrs[$i]} = $i;

		}#end if

	}#end for


	#CHECK IF ALL CHROMOSOMES ARE IN MSA, IF NOT THEN THROW ERROR
	if((@chr_inds) != (@chr_list))
	{

		print "\n\nNot all inputted chromosomes are in the MSA!!!\n\n";

		exit;

	}#end if


        #GET ALL REGIONS IN THE MSA TABLE WITH THE GIVEN ID
        my $select_regions = "SELECT ma.starts, ma.ends
                              FROM multiple_alignment ma, multiple_alignment_info mai
                              WHERE ma.msa_id = mai.msa_id AND
				    mai.msa_id = ?";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";



        INDEL:while(my ($starts, $ends) =  $sth -> fetchrow_array)
	{

		my @starts = split /\|/, $starts;
		my @ends = split /\|/, $ends;


		#GET REGIONS THAT ARE 1) NOT PRESENT IN ALL CHROMOSOMES AND  2) NOT ABSENT IN ALL CHROMOSOMES
		if((join("|", @starts[@chr_inds]) !~ /^(\|?\-?[1-9]\d*\|?)+$/) && (join("|", @starts[@chr_inds]) !~ /^(\|?0\|?)+$/))
		{

			#GET INFO ABOUT IN AND OUT CHROMSOMES
			my @out_inds = ();
			my @in_inds = ();
			my $pa_vector = "";

			my $max_indel = 0;
			my $max_indel_i;

			#STORE INDS OF IN AN OUT CHROMSOMES, AND MAKE SURE CURRENT INDEL IS LARGE ENOUGH TO CONSIDER			
			foreach my $ind (@chr2ind{@chr_list})
			{

				if( $starts[$ind] == 0) 
				{

					push @out_inds, $ind;
					$pa_vector .= 0;

				}else{


					#NOTE IND OF THE CHROMOSOME WITH THE MAXIMUM INDEL TO USE AS REFERENCE
					if(abs($ends[$ind] - $starts[$ind]) > $max_indel)
					{

						$max_indel = abs($ends[$ind] - $starts[$ind]);
						$max_indel_i = $ind;
	
					}#end if

					#CHECK STRAND OF REGION, AND ADJUST FROM MAUVE COORDS TO MY OWN
					if($starts[$ind] < 0)
					{

						my $temp = $starts[$ind];
						$starts[$ind] = abs($ends[$ind]);
						$ends[$ind] = abs($temp);

					}#end if


					push @in_inds, $ind;
					$pa_vector .= 1;

				}#end if

			}#end for	


			#CHECK THAT INDEL IS BIGGER THAN MIN SIZE
			if($max_indel < $min_size)
			{

				next INDEL;

			}#end if


			#STORE INFO
			$diff_regions{$region_counter}{'OUT_CHRS'} = join " , ", @chrs[@out_inds];

			$diff_regions{$region_counter}{'IN_CHRS'} = join " , ", @chrs[@in_inds];
			@{$diff_regions{$region_counter}{'IN_STARTS'}} = @starts[@in_inds];
			@{$diff_regions{$region_counter}{'IN_ENDS'}} = @ends[@in_inds];

			@{$diff_regions{$region_counter}{'ALL_STARTS'}} = @starts[@chr2ind{@chr_list}];
			@{$diff_regions{$region_counter}{'ALL_ENDS'}} = @ends[@chr2ind{@chr_list}];

			$diff_regions{$region_counter}{'CHR'} = $chrs[$max_indel_i];
			$diff_regions{$region_counter}{'START'} = $starts[$max_indel_i];
			$diff_regions{$region_counter}{'END'} = $ends[$max_indel_i];
			
			$diff_regions{$region_counter}{'PA_VECTOR'} = $pa_vector;


        		#GET THE SEQUENCE OF THE REGION AND STORE IF LESS THAN 50% AMBIGUOS
			my $region_name = "$chrs[$max_indel_i]\_$starts[$max_indel_i]\_$ends[$max_indel_i]";

        		my $region_seq = &genomes_db_lib::extract_sequence(&get_chr_id($chrs[$max_indel_i], $dbh), $starts[$max_indel_i], $ends[$max_indel_i], $dbh);
			
			if(&seq_stats::get_N($region_seq) < .5)
        		{

				$region_seqs{$region_name} = $region_seq;
				$chr_starts{$starts[$max_indel_i]} = $region_counter;

			}#end if


			$region_counter ++;

		}#end if

	}#end while

	return (\%diff_regions, \%region_seqs, \%chr_starts);

}#end get_group_diff_genomic_regions


#sub get_recombined_regions - Takes as input: 1) chr1		- The name of a chromsome
#			   		      2) chr2		- The name of a second chromsome
#					      3) min_max	- The minimum value during a walk for which a recombinant region is declared
#               			      4) max_drop	- The maximum amount by which the walk can drop before the region ends
#					      5) msa_id		- The multiple sequence alignment to use to get snps/indels
#				              6) min_dist	- The minimum allowed distance between SNPs before filtering
#
#			    - Returns the set of regions with a SNP density higher than some threshold. THe parameters min_max and max_drop
#			      control the minimum enrichment for a region to be flagged and the minimum decrease from the maximum to end
#			      a given region
sub get_recombined_regions
{
	my ($chr1, $chr2, $min_max, $max_drop, $msa_id, $min_dist, $dbh) = @_;

	my @chrs = ($chr1, $chr2);

	print "$chr1 vs. $chr2\n";

	#GET SNPS BETWEEN PAIRS OF CHROMSOMES
	my $r_snp_hash = &get_all_snps($chr1, $chr2, $dbh);
	my $r_filt_snps = &filter_snps_ref($chr1, $r_snp_hash, \@chrs, $min_dist, $dbh);
	my %snp_hash = %$r_filt_snps;


	#CREATE VECTOR WHERE A SNP POSITION IS A POSITIVE VALUE, NON-SNP NEGATIVE AND INDEL ZERO
	my @my_snps = keys %snp_hash;

	my $r_chr_desc_hash = &get_chr_desc(&get_chr_id($chr1, $dbh), $dbh);
	my $genome_length = length($r_chr_desc_hash->{'SEQ'});

	my @walk_vector = ();
	
	my $snp_inc = 1 / scalar(@my_snps);
	my $non_snp_dec = -1 / ($genome_length - scalar(@my_snps));

	#ASSIGN VALUE DEPENDING ON WHETHER POSITION HAS SNP OR NOT
	@walk_vector[0 .. ($genome_length-1)] = (($non_snp_dec) x $genome_length);
	@walk_vector[@my_snps] = (($snp_inc) x scalar(@my_snps));


	#ZERO OUT POSITIONS WHICH DO NOT ALIGN
	my @chr1 = ($chr1);
	my @chr2 = ($chr2);
	my $r_diff_regions = &get_diff_genomic_regions(\@chr1, \@chr2, 100, $msa_id, $dbh);
	my %diff_regions = %$r_diff_regions;

	foreach my $region(keys %diff_regions)
	{

		#CHECK THE STRAND OF THE INDEL, AND ADJUST STATEMENT ACCORDINGLY
		if($diff_regions{$region}{'START'} < $diff_regions{$region}{'END'})
		{

		        @walk_vector[$diff_regions{$region}{'START'} .. $diff_regions{$region}{'END'}] = ((0) x ($diff_regions{$region}{'END'} - $diff_regions{$region}{'START'} + 1));

		}else{

		        @walk_vector[$diff_regions{$region}{'END'} .. $diff_regions{$region}{'START'}] = ((0) x ($diff_regions{$region}{'START'} - $diff_regions{$region}{'END'} + 1));

		}#end if

	}#end foreach


	#IDENTIFY PUTITIVE RECOMBINED REGIONS BASED ON HIGHER THAN EXPECTED SNP DENSITY
	my @running_walk_vector;
	$running_walk_vector[0] = $walk_vector[0];

	my $start = 1;
	my $curr_max = 1;
	my $curr_max_pos = 0;
	my $in_rec = 0;
	my %rec_regions;
	my $pos = 1;

	while($pos < scalar(@walk_vector))
	{

		#DETERMINE RUNNING SUM AT CURRENT POSITION 
	        $running_walk_vector[$pos] = &arith_ops::max(($running_walk_vector[$pos-1] + $walk_vector[$pos]), 0);
		#print "Pos($pos), start($start), tally($running_walk_vector[$pos]), curr_max($curr_max), curr_max_pos($curr_max_pos), in_rec($in_rec)\n";


	        #IF NEWLY ABOVE THRESHOLD, THEN SET VARIABLE
	        if(($in_rec == 0) && ($running_walk_vector[$pos] > $min_max))
	        {

	                $in_rec = 1;
	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;

	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF MAX SHOULD BE INCREMENTED
	        if(($in_rec == 1) && ($running_walk_vector[$pos] > $curr_max))
	        {

	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;
	
	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF HAVE GONE FAR ENOUGH BELOW MAX TO END
	        if(($in_rec == 1) && ($running_walk_vector[$pos] < ($curr_max-$max_drop)))
	        {
	
	                $rec_regions{$start}{'END'} = $curr_max_pos;
	                $rec_regions{$start}{'MAX'} = &arith_ops::rnd($curr_max,3);

	                $start = $curr_max_pos + 1;
	                $pos = $curr_max_pos + 1;
	                $in_rec = 0;
	                $curr_max = 0;
	                $curr_max_pos = 0;
	                $running_walk_vector[$pos] = &arith_ops::max($walk_vector[$pos], 0);

	        }#end if

	        #IF AT ZERO, THEN MOVE START AHEAD
	        if($running_walk_vector[$pos] == 0)
	        {
	
        	        $start = $pos;

        	}#end if


	        $pos ++;

	}#end while

	#TALLY LAST REGION IF ENDED IN THE MIDST OF IT
	if($in_rec == 1)
	{

	        $rec_regions{$start}{'END'} = $curr_max_pos;
	        $rec_regions{$start}{'MAX'} = $curr_max;

	}#end if


	return \%rec_regions;


}#end get_recombined_regions


#sub get_recombined_regions_snpInput - Takes as input: 1) chr1		- The name of a chromsome
#			   			       2) chr2		- The name of a second chromsome
#						       3) min_max	- The minimum value during a walk for which a recombinant region is declared
#               				       4) max_drop	- The maximum amount by which the walk can drop before the region ends
#						       5) msa_id	- The multiple sequence alignment to use to get snps/indels
#				        	       6) r_snps	- A list of positions in chr1
#
#			    	      - Returns the set of regions with a SNP density higher than some threshold. THe parameters min_max and max_drop
#			      		control the minimum enrichment for a region to be flagged and the minimum decrease from the maximum to end
#			      		a given region
sub get_recombined_regions_snpInput
{
	my ($chr1, $chr2, $min_max_temp, $max_drop, $msa_id, $r_snps, $dbh) = @_;

	my @chrs = ($chr1, $chr2);
	my @my_snps = @$r_snps;

	print "$chr1 vs. $chr2\n";

	#IF NO SNPS, THEN ABORT
	if(@my_snps == 0)
	{

		return 0;	

	}#end if

	#CREATE VECTOR WHERE A SNP POSITION IS A POSITIVE VALUE, NON-SNP NEGATIVE AND INDEL ZERO
	my $r_chr_desc_hash = &get_chr_desc(&get_chr_id($chr1, $dbh), $dbh);
	my $genome_length = length($r_chr_desc_hash->{'SEQ'});

	my @walk_vector = ();
	
	my $snp_inc = 1 / scalar(@my_snps);
	my $non_snp_dec = $min_max_temp * (-1 / ($genome_length - scalar(@my_snps)));

	###########
	$min_max = $snp_inc * &arith_ops::min($min_max_temp, (0.99*scalar(@my_snps)));
	$max_drop = $snp_inc * $max_drop;
	#print "$min_max, $max_drop, $snp_inc, $non_snp_dec, " . scalar(@indel_pos)  . "\n";
	#$min_max = $snp_inc * $min_max_temp;
	#$max_drop = $snp_inc * $max_drop;
	#print "$min_max, $max_drop\n";
	###########

	#ASSIGN VALUE DEPENDING ON WHETHER POSITION HAS SNP OR NOT
	@walk_vector[0 .. ($genome_length-1)] = (($non_snp_dec) x $genome_length);

	#ZERO OUT POSITIONS WHICH DO NOT ALIGN
	my $r_unique_bases = &get_unique_bases(&get_chr_id($chr1, $dbh), &get_chr_id($chr2, $dbh), 100, $dbh);
	my %unique_bases = %$r_unique_bases;

        @walk_vector[keys(%unique_bases)] = ((0) x scalar(keys(%unique_bases)));

	#ADJUST SNP POSITIONS
	@walk_vector[@my_snps] = (($snp_inc) x scalar(@my_snps));

	#IDENTIFY PUTITIVE RECOMBINED REGIONS BASED ON HIGHER THAN EXPECTED SNP DENSITY
	my @running_walk_vector;
	$running_walk_vector[0] = $walk_vector[0];

	my $start = 1;
	my $curr_max = 1;
	my $curr_max_pos = 0;
	my $in_rec = 0;
	my %rec_regions;
	my $pos = 1;

	while($pos < scalar(@walk_vector))
	{

		#DETERMINE RUNNING SUM AT CURRENT POSITION 
	        $running_walk_vector[$pos] = &arith_ops::max(($running_walk_vector[$pos-1] + $walk_vector[$pos]), 0);

##########
		if($chr2 eq "Efae_AC78_121410_genome")
		{

##			print "Pos($pos), start($start), tally($running_walk_vector[$pos]), curr_max($curr_max), curr_max_pos($curr_max_pos), in_rec($in_rec)\n";

		}
###########
	        #IF NEWLY ABOVE THRESHOLD, THEN SET VARIABLE
	        if(($in_rec == 0) && ($running_walk_vector[$pos] > $min_max))
	        {

	                $in_rec = 1;
	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;

	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF MAX SHOULD BE INCREMENTED
	        if(($in_rec == 1) && ($running_walk_vector[$pos] > $curr_max))
	        {

	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;
	
	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF HAVE GONE FAR ENOUGH BELOW MAX TO END
	        if(($in_rec == 1) && ($running_walk_vector[$pos] < ($curr_max-$max_drop)))
	        {
	
	                $rec_regions{$start}{'END'} = $curr_max_pos;
	                $rec_regions{$start}{'MAX'} = &arith_ops::rnd($curr_max,3);

	                $start = $curr_max_pos + 1;
	                $pos = $curr_max_pos + 1;
	                $in_rec = 0;
	                $curr_max = 0;
	                $curr_max_pos = 0;
	                $running_walk_vector[$pos] = &arith_ops::max($walk_vector[$pos], 0);

	        }#end if

	        #IF AT ZERO, THEN MOVE START AHEAD
	        if($running_walk_vector[$pos] == 0)
	        {
	
        	        $start = $pos;

        	}#end if


	        $pos ++;

	}#end while

	#TALLY LAST REGION IF ENDED IN THE MIDST OF IT
	if($in_rec == 1)
	{

	        $rec_regions{$start}{'END'} = $curr_max_pos;
	        $rec_regions{$start}{'MAX'} = $curr_max;

	}#end if


	return \%rec_regions;


}#end get_recombined_regions_snpInput
 

#sub get_recombined_regions_snpInput_indelInput - Takes as input: 1) ref_chr		- THe name of the chrosome relative to which SNP posiitons are given
##							   	  2) r_snps        	- A list of positions in the reference chromsome
#			   			       	   	  3) r_indel_pos   	- A list of positions relative to the reference chromsome which are 
#			   			       	   	  			  not present in all aligned chromsomes	
#						           	  4) min_max	    	- The minimum value during a walk for which a recombinant region is 
#						           	   			  declared
#               				       	   	  5) max_drop		- The maximum amount by which the walk can drop before the region ends
#
#			    	      		- Returns the set of regions with a SNP density higher than some threshold. THe parameters min_max and 
#			    	      		  max_drop control the minimum enrichment for a region to be flagged and the minimum decrease from the 
#			    	      		  maximum to end a given region
#
sub get_recombined_regions_snpInput_indelInput
{
	my ($ref_chr, $r_snps, $r_indel_pos, $min_max_temp, $max_drop, $dbh) = @_;

	my @chrs = ($ref_chr);
	my @my_snps = @$r_snps;
	my @indel_pos = @$r_indel_pos;


	#CREATE VECTOR WHERE A SNP POSITION IS A POSITIVE VALUE, NON-SNP NEGATIVE AND INDEL ZERO
	my $r_chr_desc_hash = &get_chr_desc(&get_chr_id($ref_chr, $dbh), $dbh);
	my $genome_length = length($r_chr_desc_hash->{'SEQ'});

	my @walk_vector = ();
	
	my $snp_inc = 1 / scalar(@my_snps);
	my $non_snp_dec = -1 / ($genome_length - scalar(@my_snps));

	###########
	$min_max = $snp_inc * &arith_ops::min($min_max_temp, (0.99*scalar(@my_snps)));
	$max_drop = $snp_inc * $max_drop;
	#print "$min_max, $max_drop, $snp_inc, $non_snp_dec, " . scalar(@indel_pos)  . "\n";
	###########

	#ASSIGN VALUE DEPENDING ON WHETHER POSITION HAS SNP OR NOT
	@walk_vector[0 .. ($genome_length-1)] = (($non_snp_dec) x $genome_length);

	#ZERO OUT POSITIONS WHICH DO NOT ALIGN
	@walk_vector[@indel_pos] = ((0) x scalar(@indel_pos));

	#ADJUST SNP POSITIONS
	@walk_vector[@my_snps] = (($snp_inc) x scalar(@my_snps));


	#IDENTIFY PUTITIVE RECOMBINED REGIONS BASED ON HIGHER THAN EXPECTED SNP DENSITY
	my @running_walk_vector;
	$running_walk_vector[0] = $walk_vector[0];

	my $start = 1;
	my $curr_max = 1;
	my $curr_max_pos = 0;
	my $in_rec = 0;
	my %rec_regions;
	my $pos = 1;

	while($pos < scalar(@walk_vector))
	{

		#DETERMINE RUNNING SUM AT CURRENT POSITION 
	        $running_walk_vector[$pos] = &arith_ops::max(($running_walk_vector[$pos-1] + $walk_vector[$pos]), 0);

	        #IF NEWLY ABOVE THRESHOLD, THEN SET VARIABLE
	        if(($in_rec == 0) && ($running_walk_vector[$pos] > $min_max))
	        {

	                $in_rec = 1;
	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;

	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF MAX SHOULD BE INCREMENTED
	        if(($in_rec == 1) && ($running_walk_vector[$pos] > $curr_max))
	        {

	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;
	
	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF HAVE GONE FAR ENOUGH BELOW MAX TO END
	        if(($in_rec == 1) && ($running_walk_vector[$pos] < ($curr_max-$max_drop)))
	        {
	
	                $rec_regions{$start}{'END'} = $curr_max_pos;
	                $rec_regions{$start}{'MAX'} = &arith_ops::rnd($curr_max,3);

	                $start = $curr_max_pos + 1;
	                $pos = $curr_max_pos + 1;
	                $in_rec = 0;
	                $curr_max = 0;
	                $curr_max_pos = 0;
	                $running_walk_vector[$pos] = &arith_ops::max($walk_vector[$pos], 0);

	        }#end if

	        #IF AT ZERO, THEN MOVE START AHEAD
	        if($running_walk_vector[$pos] == 0)
	        {
	
        	        $start = $pos;

        	}#end if


	        $pos ++;

	}#end while

	#TALLY LAST REGION IF ENDED IN THE MIDST OF IT
	if($in_rec == 1)
	{

	        $rec_regions{$start}{'END'} = $curr_max_pos;
	        $rec_regions{$start}{'MAX'} = $curr_max;

	}#end if


	return \%rec_regions;

}#end get_recombined_regions_snpInput_indelInput


#sub get_recombined_regions_msa - Takes as input: 1) chr1	- The name of a chromsome
#			   		          2) chr2	- The name of a second chromsome
#					          3) min_max	- The minimum value during a walk for which a recombinant region is declared
#               			          4) max_drop	- The maximum amount by which the walk can drop before the region ends
#					          5) msa_id	- The multiple sequence alignment to use to get snps/indels
#					          6) min_dist	- The minimum distance between two SNPs, for them to be filtered out
#								  because of assumed sequencing/assembly/alingment error
#					          7) rep_chrs	- The chromsomes in which repetative sequences are used to filter SNPs
#
#			        - Returns the set of regions with a SNP density higher than some threshold. THe parameters min_max and max_drop
#				  control the minimum enrichment for a region to be flagged and the minimum decrease from the maximum to end
#				  a given region
sub get_recombined_regions_msa
{
	my ($chr1, $chr2, $min_max, $max_drop, $msa_id, $min_dist, $r_rep_chrs,  $dbh) = @_;

	my @chrs = ($chr1, $chr2);

	print "$chr1 vs. $chr2\n";

	#GET SNPS BETWEEN PAIR OF CHROMSOMES
	my ($r_snp_hash, $r_chr_names, $r_snp_msa_hash, $r_snp2pos) = &genomes_db_lib::get_all_snps_msa($msa_id, \@chrs, $dbh);
        $r_filt_snp_hash = &genomes_db_lib::filter_snps_msa($r_snp2pos, \@chrs, $r_rep_chrs, $msa_id, $min_dist, $dbh);

	#my %snp_hash = %$r_snp_msa_hash;
	my %snp_hash = %$r_filt_snp_hash;


	#CREATE VECTOR WHERE A SNP POSITION IS A POSITIVE VALUE, NON-SNP NEGATIVE AND INDEL ZERO
	my @my_snps = keys %snp_hash;

	my $r_chr_desc_hash = &get_chr_desc(&get_chr_id($chr1, $dbh), $dbh);
	my $genome_length = length($r_chr_desc_hash->{'SEQ'});

	my @walk_vector = ();
	
	my $snp_inc = 1 / scalar(@my_snps);
	my $non_snp_dec = -1 / ($genome_length - scalar(@my_snps));
	print "inc = $snp_inc, dec = $non_snp_dec\n";

	#ASSIGN VALUE DEPENDING ON WHETHER POSITION HAS SNP OR NOT
	@walk_vector[0 .. ($genome_length-1)] = (($non_snp_dec) x $genome_length);
	@walk_vector[@my_snps] = (($snp_inc) x scalar(@my_snps));


	#ZERO OUT POSITIONS WHICH DO NOT ALIGN
	my @chr1 = ($chr1);
	my @chr2 = ($chr2);
	my $r_diff_regions = &get_diff_genomic_regions(\@chr1, \@chr2, 100, $msa_id, $dbh);
	my %diff_regions = %$r_diff_regions;

	foreach my $region(keys %diff_regions)
	{

		#CHECK THE STRAND OF THE INDEL, AND ADJUST STATEMENT ACCORDINGLY
		if($diff_regions{$region}{'START'} < $diff_regions{$region}{'END'})
		{

	        	@walk_vector[$diff_regions{$region}{'START'} .. $diff_regions{$region}{'END'}] = ((0) x ($diff_regions{$region}{'END'} - $diff_regions{$region}{'START'} + 1));

		}else{

	        	@walk_vector[$diff_regions{$region}{'END'} .. $diff_regions{$region}{'START'}] = ((0) x ($diff_regions{$region}{'START'} - $diff_regions{$region}{'END'} + 1));

		}#end if

	}#end foreach


	#IDENTIFY PUTITIVE RECOMBINED REGIONS BASED ON HIGHER THAN EXPECTED SNP DENSITY
	my @running_walk_vector;
	$running_walk_vector[0] = $walk_vector[0];

	my $start = 1;
	my $curr_max = 1;
	my $curr_max_pos = 0;
	my $in_rec = 0;
	my %rec_regions;
	my $pos = 1;

	while($pos < scalar(@walk_vector))
	{

		#DETERMINE RUNNING SUM AT CURRENT POSITION 
	        $running_walk_vector[$pos] = &arith_ops::max(($running_walk_vector[$pos-1] + $walk_vector[$pos]), 0);
		#print "Pos($pos), start($start), tally($running_walk_vector[$pos]), curr_max($curr_max), curr_max_pos($curr_max_pos), in_rec($in_rec)\n";


	        #IF NEWLY ABOVE THRESHOLD, THEN SET VARIABLE
	        if(($in_rec == 0) && ($running_walk_vector[$pos] > $min_max))
	        {

	                $in_rec = 1;
	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;

	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF MAX SHOULD BE INCREMENTED
	        if(($in_rec == 1) && ($running_walk_vector[$pos] > $curr_max))
	        {

	                $curr_max = $running_walk_vector[$pos];
	                $curr_max_pos = $pos;
	
	        }#end if

	        #IF IN RECOMBINANT REGION, THEN CHECK IF HAVE GONE FAR ENOUGH BELOW MAX TO END
	        if(($in_rec == 1) && ($running_walk_vector[$pos] < ($curr_max-$max_drop)))
	        {
	
	                $rec_regions{$start}{'END'} = $curr_max_pos;
	                $rec_regions{$start}{'MAX'} = &arith_ops::rnd($curr_max,3);

	                $start = $curr_max_pos + 1;
	                $pos = $curr_max_pos + 1;
	                $in_rec = 0;
	                $curr_max = 0;
	                $curr_max_pos = 0;
	                $running_walk_vector[$pos] = &arith_ops::max($walk_vector[$pos], 0);

	        }#end if

	        #IF AT ZERO, THEN MOVE START AHEAD
	        if($running_walk_vector[$pos] == 0)
	        {
	
        	        $start = $pos;

        	}#end if


	        $pos ++;

	}#end while

	#TALLY LAST REGION IF ENDED IN THE MIDST OF IT
	if($in_rec == 1)
	{

	        $rec_regions{$start}{'END'} = $curr_max_pos;
	        $rec_regions{$start}{'MAX'} = $curr_max;

	}#end if


	return \%rec_regions;


}#end get_recombined_regions_msa


#sub get_recombined_regions_rel_ref_phyloDisc - Takes as input: 1) chrs       - The list of chromosomes, with the first being the reference
#                                          	      		2) min_max    - The minimum value during a walk for which a recombinant region is declared
#                                                     		3) max_drop   - The maximum amount by which the walk can drop before the region ends
#                                                     		4) msa_id     - The multiple sequence alignment to use to get snps/indels
#						      		5) snp_hash   - A 2D hash, where the first key is the chromosome in which SNPs are being
#								      		found in, and the second key is the position in the reference chromsome
#								      		***REFERENCE MUST BE THE FIRST CHROMOSOME IN THE chrs array
#
#                                   	      - Returns the set of regions with a SNP density higher than some threshold, relative to a reference. THe 
#                                   	        parameters min_max and max_drop control the minimum enrichment for a region to be flagged and the minimum 
#                                   	        decrease from the maximum to end a given region. 
sub get_recombined_regions_rel_ref_phyloDisc
{
	my ($r_chrs, $min_max, $max_drop, $msa_id, $r_phylo_hash, $dbh) = @_;


	#SEPARATE SNPS BY PHYLOGENETIC PATTERN
	my @chrs = @$r_chrs;
	my %phylo_hash = %$r_phylo_hash;
	my %snp_hash;
	my %pcodes;

	foreach my $pos (keys %phylo_hash)
	{

		if(exists($pcodes{$phylo_hash{$pos}}))
		{
			#PCODE SEEN, ADD SNP TO LIST
			push @{$snp_hash{$phylo_hash{$pos}}}, $pos;

		}else{

			#PCODE NOT SEEN, ADD SNP TO LIST AND PCODE TO LIST
			@{$snp_hash{$phylo_hash{$pos}}} = ($pos);

			$pcodes{$phylo_hash{$pos}} = 1;

		}#end if

	}#end foreach
	

	#GET GENOMIC SEQUENCE PRESENT IN REFERENCE, BUT NOT IN ALL OTHER CHROMOSOMES IN ALIGNMENT
	my @indel_pos;
	my %indel_pos;

	for(my $c = 1; $c < @chrs; $c++)
	{

		#ZERO OUT POSITIONS WHICH DO NOT ALIGN
		my $r_unique_bases = &get_unique_bases(&get_chr_id($chrs[0], $dbh), &get_chr_id($chrs[$c], $dbh), 100, $dbh);
		my %unique_bases = %$r_unique_bases;

        	@indel_pos{keys %unique_bases} = ((1) x scalar(keys(%unique_bases)));

	}#end foreach

	@indel_pos = sort {$a<=>$b} keys %indel_pos;


	#GET BASES IN FIRST CHROMSOMES THAT ARE PART OF RECOMBINE REGION RELATIVE TO ANY OTHER CHROMSOME
	my @rec_bases = ();

	open DEBUG, ">>temp_snp_pos_sort_3";
	print DEBUG "\n-----------------------------------------------------\n";
	print DEBUG "PHYLO DISC ROUTINE\n";
	print DEBUG "MIN_MAX = $min_max, MAX_DROP = $max_drop\n";

	for my $pcode (keys %pcodes)
	{

		print "$pcode\n";
		print DEBUG "---\n";
		print DEBUG "$pcode\n";

		my @snps_unsorted = @{$snp_hash{$pcode}};
		my @snps = sort{$a <=> $b} @snps_unsorted;
		my @temp_rec_bases = ();


		my $r_rec_regions;
		my %rec_regions;

		if(@snps > 1)
		{
	
		        $r_rec_regions = &get_recombined_regions_snpInput_indelInput($chrs[0], \@snps, \@indel_pos, $min_max, $max_drop, $dbh) ;
		        %rec_regions = %$r_rec_regions;

		}#end if

	        foreach my $start(sort {$a <=> $b} keys %rec_regions)
	        {

	                print DEBUG "$pcode: $start-$rec_regions{$start}{'END'}\n";
	                push @rec_bases, $start..$rec_regions{$start}{'END'};
	                push @temp_rec_bases, $start..$rec_regions{$start}{'END'};

	        }#end foreach

		print DEBUG "REC BASES = " . @temp_rec_bases . "\n\n";
		print DEBUG join("\n", @snps) . "\n";

	}#end for

	close DEBUG;


	#PARTITION INTO DISCRETE REGIONS
	my $r_rec_bases_hash = &hash_lib::init_hash(\@rec_bases);
	my %rec_bases_hash = %$r_rec_bases_hash;

	my %rec_region_hash;

	my $start_pos = -1;
	my $prev_pos = -1;

	#print "Total DNA in recombination regions is: " . keys(%rec_bases_hash) . "\n";

	foreach my $pos (sort {$a<=>$b} keys %rec_bases_hash)
	{

	        #SEE IF AT THE START OF A NEW REGION, AND IF SO SAVE PREVIOUS ONE
	        if($pos > ($prev_pos+1) && $start_pos != -1)
	        {

	                #print "$start_pos - $prev_pos\n";
	                $rec_region_hash{$start_pos}{'END'} = $prev_pos;
	                $start_pos = $pos;

	        }elsif($start_pos == -1){

	                $start_pos = $pos;

        	}#end if


	        #KEEP TRACK OF POSITION
	        $prev_pos = $pos

	}#end foreach

	#SAVE LAST REGION
	$rec_region_hash{$start_pos}{'END'} = $prev_pos;


	return (\%rec_region_hash, \%rec_bases_hash);

}#end get_recombined_regions_rel_ref_phyloDisc



#sub get_recombined_regions_rel_ref_phyloDisc_iter - Takes as input: 1) chrs       - The list of chromosomes, with the first being the reference
#                                          	      		     2) min_max    - The minimum value during a walk for which a recombinant region is declared
#                                                     		     3) max_drop   - The maximum amount by which the walk can drop before the region ends
#                                                     		     4) msa_id     - The multiple sequence alignment to use to get snps/indels
#						      		     5) snp_hash   - A 2D hash, where the first key is the chromosome in which SNPs are being
#								      		found in, and the second key is the position in the reference chromsome
#								      		***REFERENCE MUST BE THE FIRST CHROMOSOME IN THE chrs array
#
#                                   	      - Returns the set of regions with a SNP density higher than some threshold, relative to a reference. THe 
#                                   	        parameters min_max and max_drop control the minimum enrichment for a region to be flagged and the minimum 
#                                   	        decrease from the maximum to end a given region. 
sub get_recombined_regions_rel_ref_phyloDisc_iter
{
	my ($r_chrs, $min_max, $max_drop, $msa_id, $r_phylo_hash, $dbh) = @_;


	#SEPARATE SNPS BY PHYLOGENETIC PATTERN
	my @chrs = @$r_chrs;
	my %phylo_hash = %$r_phylo_hash;
	my %snp_hash;
	my %pcodes;

	foreach my $pos (keys %phylo_hash)
	{

		if(exists($pcodes{$phylo_hash{$pos}}))
		{
			#PCODE SEEN, ADD SNP TO LIST
			push @{$snp_hash{$phylo_hash{$pos}}}, $pos;

		}else{

			#PCODE NOT SEEN, ADD SNP TO LIST AND PCODE TO LIST
			@{$snp_hash{$phylo_hash{$pos}}} = ($pos);

			$pcodes{$phylo_hash{$pos}} = 1;

		}#end if

	}#end foreach
	

	#GET GENOMIC SEQUENCE PRESENT IN REFERENCE, BUT NOT IN ALL OTHER CHROMOSOMES IN ALIGNMENT
	my @indel_pos;
	my %indel_pos;

	for(my $c = 1; $c < @chrs; $c++)
	{

		#ZERO OUT POSITIONS WHICH DO NOT ALIGN
		my $r_unique_bases = &get_unique_bases(&get_chr_id($chrs[0], $dbh), &get_chr_id($chrs[$c], $dbh), 100, $dbh);
		my %unique_bases = %$r_unique_bases;

        	@indel_pos{keys %unique_bases} = ((1) x scalar(keys(%unique_bases)));

	}#end foreach

	@indel_pos = sort {$a<=>$b} keys %indel_pos;


	#GET BASES IN FIRST CHROMSOMES THAT ARE PART OF RECOMBINE REGION RELATIVE TO ANY OTHER CHROMSOME
	my @rec_bases = ();

	#open DEBUG, ">>temp_snp_pos_sort_3";
	#print DEBUG "\n-----------------------------------------------------\n";
	#print DEBUG "PHYLO DISC ROUTINE\n";
	#print DEBUG "MIN_MAX = $min_max, MAX_DROP = $max_drop\n";

	for my $pcode (keys %pcodes)
	{

		#print DEBUG "---\n";
		#print DEBUG "$pcode\n";

		my @snps_unsorted = @{$snp_hash{$pcode}};
		my @snps = sort{$a <=> $b} @snps_unsorted;
		my @all_snps = @snps;
		my @temp_rec_bases = ();


		#ITERATIVELY FIND REC REGIONS, REMOVING THOSE FOUND IN PREVIOUS ITERATION
		my $num_rec_regions = 1;
		my $iter = 1;

		while ($num_rec_regions > 0)
		{

			#FIND REC REGIONS AND SAVE POSITIONS
			my $r_rec_regions;
			my %rec_regions;

			if(@snps > 1)
			{
	
			        $r_rec_regions = &get_recombined_regions_snpInput_indelInput($chrs[0], \@snps, \@indel_pos, $min_max, $max_drop, $dbh) ;
			        %rec_regions = %$r_rec_regions;

			}#end if

	        	foreach my $start(sort {$a <=> $b} keys %rec_regions)
	        	{

	        	        #print DEBUG "$pcode ($iter): $start-$rec_regions{$start}{'END'}\n";
	                	push @rec_bases, $start..$rec_regions{$start}{'END'};
	                	push @temp_rec_bases, $start..$rec_regions{$start}{'END'};

	        	}#end foreach

			#REMOVE SNPS IN PREVIOUSLY IDENTIFIED REC REGIONS
			my $lc = List::Compare->new(\@snps, \@temp_rec_bases);
			@snps = $lc->get_unique;

			#UPDATE NUMBER OF REC REGIONS FOUND IN THIS ITERATION
			$num_rec_regions = scalar(keys %rec_regions);
			$iter ++;

		}#end while

		#print DEBUG "REC BASES = " . @temp_rec_bases . "\n\n";
		#print DEBUG join("\n", @all_snps) . "\n";

	}#end for

	#close DEBUG;


	#PARTITION INTO DISCRETE REGIONS
	my $r_rec_bases_hash = &hash_lib::init_hash(\@rec_bases);
	my %rec_bases_hash = %$r_rec_bases_hash;

	my %rec_region_hash;

	my $start_pos = -1;
	my $prev_pos = -1;

	#print "Total DNA in recombination regions is: " . keys(%rec_bases_hash) . "\n";

	foreach my $pos (sort {$a<=>$b} keys %rec_bases_hash)
	{

	        #SEE IF AT THE START OF A NEW REGION, AND IF SO SAVE PREVIOUS ONE
	        if($pos > ($prev_pos+1) && $start_pos != -1)
	        {

	                #print "$start_pos - $prev_pos\n";
	                $rec_region_hash{$start_pos}{'END'} = $prev_pos;
	                $start_pos = $pos;

	        }elsif($start_pos == -1){

	                $start_pos = $pos;

        	}#end if


	        #KEEP TRACK OF POSITION
	        $prev_pos = $pos

	}#end foreach

	#SAVE LAST REGION
	$rec_region_hash{$start_pos}{'END'} = $prev_pos;


	return (\%rec_region_hash, \%rec_bases_hash);

}#end get_recombined_regions_rel_ref_phyloDisc_iter


#sub get_recombined_regions_rel_ref - Takes as input: 1) chrs       - The list of chromosomes, with the first being the reference
#                                          	      2) min_max    - The minimum value during a walk for which a recombinant region is declared
#                                                     3) max_drop   - The maximum amount by which the walk can drop before the region ends
#                                                     4) msa_id     - The multiple sequence alignment to use to get snps/indels
#						      5) snp_hash   - A 2D hash, where the first key is the chromosome in which SNPs are being
#								      found in, and the second key is the position in the reference chromsome
#								      ***REFERENCE MUST BE THE FIRST CHROMOSOME IN THE chrs array
#
#                                   - Returns the set of regions with a SNP density higher than some threshold, relative to a reference. THe parameters 
#                                     min_max and max_drop control the minimum enrichment for a region to be flagged and the minimum decrease from the 
#                                     maximum to end a given region. 
sub get_recombined_regions_rel_ref
{
	my ($r_chrs, $min_max, $max_drop, $msa_id, $r_snp_hash, $dbh) = @_;


	#GET RECOMBINED REGIONS AMONG THE CHROMSOMES
	my @rec_bases = ();

	#GET BASES IN FIRST CHROMSOMES THAT ARE PART OF RECOMBINE REGION RELATIVE TO ANY OTHER CHROMSOME
	my @chrs = @$r_chrs;
	my %snp_hash = %$r_snp_hash;

	open DEBUG, ">>temp_snp_pos_sort_2";
	print DEBUG "\n-----------------------------------------------------\n";
	print DEBUG "MIN_MAX = $min_max, MAX_DROP = $max_drop\n";

	for(my $c = 1; $c < @chrs; $c++)
	{

		print DEBUG "$chrs[0] vs. $chrs[$c]\n";
		print "$chrs[0] vs. $chrs[$c]\n";

		my @snps = sort{$a <=> $b} keys(%{$snp_hash{$chrs[$c]}});


	        my $r_rec_regions = &get_recombined_regions_snpInput($chrs[0], $chrs[$c], $min_max, $max_drop, $msa_id, \@snps, $dbh);
	        my %rec_regions = %$r_rec_regions;
		my @temp_rec_bases = ();


	        foreach my $start(sort {$a <=> $b} keys %rec_regions)
	        {

	                print DEBUG "$chrs[$c]: $start-$rec_regions{$start}{'END'}\n";
	                push @rec_bases, $start..$rec_regions{$start}{'END'};
	                push @temp_rec_bases, $start..$rec_regions{$start}{'END'};

	        }#end foreach

		print DEBUG "REC BASES = " . @temp_rec_bases . "\n\n";
		print DEBUG join("\n", @snps) . "\n";

	}#end for
	close DEBUG;


	#PARTITION INTO DISCRETE REGIONS
	my $r_rec_bases_hash = &hash_lib::init_hash(\@rec_bases);
	my %rec_bases_hash = %$r_rec_bases_hash;

	my %rec_region_hash;

	my $start_pos = -1;
	my $prev_pos = -1;

	#print "Total DNA in recombination regions is: " . keys(%rec_bases_hash) . "\n";

	foreach my $pos (sort {$a<=>$b} keys %rec_bases_hash)
	{

	        #SEE IF AT THE START OF A NEW REGION, AND IF SO SAVE PREVIOUS ONE
	        if($pos > ($prev_pos+1) && $start_pos != -1)
	        {

	                #print "$start_pos - $prev_pos\n";
	                $rec_region_hash{$start_pos}{'END'} = $prev_pos;
	                $start_pos = $pos;

	        }elsif($start_pos == -1){

	                $start_pos = $pos;

        	}#end if


	        #KEEP TRACK OF POSITION
	        $prev_pos = $pos

	}#end foreach

	#SAVE LAST REGION
	$rec_region_hash{$start_pos}{'END'} = $prev_pos;


	return (\%rec_region_hash, \%rec_bases_hash);

}#end get_recombined_regions_rel_ref


#sub get_recombined_regions_msa_rel_ref - Takes as input: 1) chrs       - The list of chromosomes, with the first being the reference
#                                          		  3) min_max    - The minimum value during a walk for which a recombinant region is declared
#                                                 	  4) max_drop   - The maximum amount by which the walk can drop before the region ends
#                                                 	  5) msa_id     - The multiple sequence alignment to use to get snps/indels
#                                                 	  6) min_dist   - The minimum distance between two SNPs, for them to be filtered out
#                                                     	                  because of assumed sequencing/assembly/alingment error
#                                                 	  7) rep_chrs   - The chromsomes in which repetative sequences are used to filter SNPs
#
#                               	- Returns the set of regions with a SNP density higher than some threshold, relative to a reference. THe parameters min_max and max_drop
#                                 	control the minimum enrichment for a region to be flagged and the minimum decrease from the maximum to end
#                                 	a given region. 
sub get_recombined_regions_msa_rel_ref
{
	my ($r_chrs, $min_max, $max_drop, $msa_id, $min_dist, $r_rep_chrs, $dbh) = @_;


	#GET RECOMBINED REGIONS AMONG THE CHROMSOMES
	my @rec_bases = ();

	#GET BASES IN FIRST CHROMSOMES THAT ARE PART OF RECOMBINE REGION RELATIVE TO ANY OTHER CHROMSOME
	my @chrs = @$r_chrs;
	my @rep_chrs = @$r_rep_chrs;

	for(my $c = 1; $c < @chrs; $c++)
	{

	        #my $r_rec_regions = &genomes_db_lib::get_recombined_regions($chr_names[0], $chr_names[$c], .01, .005, $msa_id, $dbh);
	        my $r_rec_regions = &get_recombined_regions_msa($chrs[0], $chrs[$c], $min_max, $max_drop, $msa_id, $min_dist, \@rep_chrs, $dbh);
	        my %rec_regions = %$r_rec_regions;


	        foreach my $start(sort {$a <=> $b} keys %rec_regions)
	        {

	                print "$start-$rec_regions{$start}{'END'}\n";
	                push @rec_bases, $start..$rec_regions{$start}{'END'};

	        }#end foreach

	}#end for


	#PARTITION INTO DISCRETE REGIONS
	my $r_rec_bases_hash = &hash_lib::init_hash(\@rec_bases);
	my %rec_bases_hash = %$r_rec_bases_hash;

	my %rec_region_hash;

	my $start_pos = -1;
	my $prev_pos = -1;

	#print "Total DNA in recombination regions is: " . keys(%rec_bases_hash) . "\n";

	foreach my $pos (sort {$a<=>$b} keys %rec_bases_hash)
	{

	        #SEE IF AT THE START OF A NEW REGION, AND IF SO SAVE PREVIOUS ONE
	        if($pos > ($prev_pos+1) && $start_pos != -1)
	        {

	                #print "$start_pos - $prev_pos\n";
	                $rec_region_hash{$start_pos}{'END'} = $prev_pos;
	                $start_pos = $pos;

	        }elsif($start_pos == -1){

	                $start_pos = $pos;

        	}#end if


	        #KEEP TRACK OF POSITION
	        $prev_pos = $pos

	}#end foreach

	#SAVE LAST REGION
	$rec_region_hash{$start_pos}{'END'} = $prev_pos;


	return (\%rec_region_hash, \%rec_bases_hash);

}#end get_recombined_regions_msa_rel_ref

#sub get_genomic_region_summary - Take as input: 1) chr_id	- An ID from the chromsome table
#						 2) start	- A start location
#						 3) end		- An end location
#				- Returns a hash containing summary info about the genomic region
sub get_genomic_region_summary
{
	my ($chr_id, $start, $end, $dbh) = @_;

	my %region_summary;
	my %gene2annots;

	my $num_integrase = 0;
	my $num_transposase = 0;
	my $num_intI = 0;
	my $num_recombinase = 0;

	my %integrase_dist;
	my %transposase_dist;
	my %intI_dist;
	my %recombinase_dist;

	my %cog_summary;

	my $perc_phage_genes = 0;
	my $region_GC = 0;
	my $chr_GC = 0;
	my $region_di_nuc_gt1sd_perc = 0;


	#DETERMINE WHICH CONTIG THE REGION IS PART OF
	my $r_contig_info = &chr2contigs($chr_id, $start, $end, $dbh);
	my %contig_info = %$r_contig_info;

	my @contigs = sort keys %contig_info;
	my @contig_lengths = ();
	my @contig_descs = ();


	foreach my $contig (@contigs)
	{

		push @contig_lengths, abs($contig_info{$contig}{'START'} - $contig_info{$contig}{'END'});
		push @contig_descs, $contig_info{$contig}{'DESC'};

	}#end foreach	


	my $contigs = join ", ", @contigs;
	my $contig_lengths = join ", ", @contig_lengths;
	my $contig_descs = join ", ", @contig_descs;


	#GET GENES AND ANNOTATIONS FROM GENOMIC REGION
	my $r_region_genes = &get_chr_region_genes($chr_id, $start, $end, $dbh);
	my $r_chr_genes = &get_chr_genes($chr_id, $dbh);

	my $r_gene2annots = &get_gene_annotations($r_chr_genes, $dbh);
	%gene2annots = %$r_gene2annots;

	
	#ONLY PROCEED IF THERE ARE GENES IN THE REGION
	if( scalar(@$r_region_genes) == 0)
	{

		print "NO GENES IN REGION!!\n\n";

		return \%region_summary;

	}#end if


	
	#GET NR ANNOTATIONS FOR GENES IN REGION
	my $r_region_annots = &get_gene_annotation_list($r_region_genes, $dbh);
	my @region_annots = @$r_region_annots;


	#DETERMINE COG FUNCTIONAL SUMMARY FOR GENE LIST
	my $ga_dbh = &gene_annotation_db_lib::gene_annotation_db_connect();

	my $r_cog_code_info = &gene_annotation_db_lib::get_cog_codes($ga_dbh);
	my %cog_code_info = %$r_cog_code_info;

	my @cog_codes = keys %cog_code_info;

	my $r_cog_summary = &hash_lib::init_hash(\@cog_codes);
	%cog_summary = %$r_cog_summary;

	foreach my $gene (@$r_region_genes)
	{

		if(exists($gene2annots{$gene}{'COG CODE'}) && $gene2annots{$gene}{'COG CODE'} ne "")
		{

			foreach my $code ( @{$gene2annots{$gene}{'COG CODE'}} )
			{

				$cog_summary{$code} ++;

			}#end foreach

		}#end if	

	}#end foreach


	#COMPUTE THE PERCENTAGE OF PAHGE GENES IN REGION (OVERLAP PHAGE FINDER REGION OR HAVE PHAGE IN ANNOTATION)
	my $r_phage_genes = &is_gene_phage($r_region_genes, $dbh);
	my %phage_genes = %$r_phage_genes;

	my $num_phage_genes = &hash_lib::hash_sum(\%phage_genes);

        $perc_phage_genes = &arith_ops::rnd($num_phage_genes / scalar(@$r_region_genes), 2);


        #COMPUTE THE PERCENTAGE OF GENES WHOSE DI-NUC SIGNATURE IS GREATER THAN 1 SD ABOVE MEAN IN GENOME
	my $r_di_nuc_hash = &is_diNuc_gt1SD($r_region_genes, $chr_id, $dbh);
       	$region_di_nuc_gt1sd_perc = &arith_ops::rnd(keys(%$r_di_nuc_hash) / scalar(@$r_region_genes), 2);


	#DETERMINE REGION GC CONTENT
	my $r_chr_info = &get_chr_info($dbh);
	my %chr_info = %$r_chr_info;

	my $region_seq = &extract_sequence($chr_id, $start, $end, $dbh);
	my $region_seq_obj = Bio::Seq->new(-seq => $region_seq);

	my $r_region_nuc_counts = &seq_stats::init_monoNucHash();
	$r_region_nuc_counts = &seq_stats::get_monoNucCount($r_region_nuc_counts, $region_seq_obj);
	my %region_nuc_counts = %$r_region_nuc_counts;

	my $region_length = abs($end - $start);
	my $N_count = $region_length - ($region_nuc_counts{'A'} + $region_nuc_counts{'C'} + $region_nuc_counts{'G'}  + $region_nuc_counts{'T'});
	my $LQ_count = ($region_seq =~ tr/actg/ACTG/);

	if( ($N_count/$region_length) < 0.50 )
	{
	
		$region_GC = &arith_ops::rnd(&seq_stats::get_GC($r_region_nuc_counts), 2);

	}else{
	
		$region_GC = "NA";

		print "SCAFFOLDING REGION WITH GREATER THAN 50% N'S!!\n\n";
		return \%region_summary;

	}#end if


	#DETERMINE DISTANCE TO NEAREST TRANSPOSASE/INTEGRASE
	my $r_start2gene = &get_gene_start_5to3($r_chr_genes, $dbh);	
	my %start2gene = %$r_start2gene;

	my $integrase_3p_flag = 0;
	my $intI_3p_flag = 0;
	my $transposase_3p_flag = 0;
	my $recombinase_3p_flag = 0;
	my $integrase_5p_flag = 0;
	my $intI_5p_flag = 0;
	my $transposase_5p_flag = 0;
	my $recombinase_5p_flag = 0;

	my $gene_counter = 0;
	my $integrase_dist_from_start = 0;
	my $integrase_dist_from_end = 0;
	my $intI_dist_from_start = 0;
	my $intI_dist_from_end = 0;
	my $transposase_dist_from_start = 0;
	my $transposase_dist_from_end = 0;
	my $recombinase_dist_from_start = 0;
	my $recombinase_dist_from_end = 0;

	$intI_dist{'5_count'} = 0;
	$intI_dist{'3_count'} = 0;
	$transposase_dist{'5_count'} = 0;
	$transposase_dist{'3_count'} = 0;
	$integrase_dist{'5_count'} = 0;
	$integrase_dist{'3_count'} = 0;
	$recombinase_dist{'5_count'} = 0;
	$recombinase_dist{'3_count'} = 0;


	foreach my $gene_start (sort {$a <=> $b} keys %start2gene)
	{

		#KEEP TRACK OF GENE NUMBER SO LOCATION OF FIRST GENE TYPE CAN BE NOTED FOR PURPOSES OF EDGE EFFECT
		$gene_counter ++;

		#KEEP TRACK OF GENES 5' OF REGION
		if($gene_start < $start)
		{

			#CHECK IF CURRENT GENE BELONGS TO GIVEN CLASSES, AND IF SO START COUNT OVER
			if( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /intI|integron/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /intI|integron/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /intI|integron/i)
			{

				$intI_dist{'5_count'} = 0;
				$intI_5p_flag = 1;
				
				if($intI_dist_from_start == 0)
				{
				
					$intI_dist_from_start = $gene_counter;

				}#end if


			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /transposase/i){

				$transposase_dist{'5_count'} = 0;
				$transposase_5p_flag = 1;

				if($transposase_dist_from_start == 0)
                                {
                              
				          $transposase_dist_from_start = $gene_counter;

                                }#end if

			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /integrase/i){

				$integrase_dist{'5_count'} = 0;
				$integrase_5p_flag = 1;

				if($integrase_dist_from_start == 0)
				{
				
					$integrase_dist_from_start = $gene_counter;

				}#end if


			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /recombinase/i){

                                $recombinase_dist{'5_count'} = 0;
                                $recombinase_5p_flag = 1;

                                if($recombinase_dist_from_start == 0)
                                {

                                        $recombinase_dist_from_start = $gene_counter;

                                }#end if


			}#end if

			
			#INCREMENT DISTANCE
			$integrase_dist{'5_count'} ++;
			$transposase_dist{'5_count'} ++;
			$intI_dist{'5_count'} ++;
			$recombinase_dist{'5_count'} ++;

		}#end if

		
		#KEEP TRACK OF GENES IN REGION
		if($gene_start > $start && $gene_start < $end)
                {

                        #CHECK IF CURRENT GENE BELONGS TO GIVEN CLASSES, AND IF SO START COUNT OVER
			#print "GENE ANNOT = " . $gene2annots{$start2gene{$gene_start}}{'NR'} . "\n";
                        if( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /intI|integron/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /intI|integron/i ||  $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /intI|integron/i )
                        {

                                $num_intI ++;

                        }elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /transposase/i ){
				
                                $num_transposase ++;

                        }elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /integrase/i){

                                $num_integrase ++;

			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /recombinase/i){

                                $num_recombinase ++;

                        }#end if


                }#end if



		#KEEP TRACK OF CLOSEST GENE 3' OF THE REGION
		if($gene_start > $end)
		{

			#INCREMENT GENE DISTANCE IF GENE OF GIVEN CLASS NOT SEEN YET
			if($intI_3p_flag == 0)
			{
	
				$intI_dist{'3_count'} ++;
	
			}#end if

			if($transposase_3p_flag == 0)
			{
	
				$transposase_dist{'3_count'} ++;
	
			}#end if

			if($integrase_3p_flag == 0)
			{
	
				$integrase_dist{'3_count'} ++;
	
			}#end if

			if($recombinase_3p_flag == 0)
			{
	
				$recombinase_dist{'3_count'} ++;
	
			}#end if


			#CHECK IF CURRENT GENE IS OF GIVEN CLASS
			if( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /intI|integron/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /intI|integron/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /intI|integron/i)
			{

				$intI_3p_flag = 1;
                                $intI_dist_from_end = 0;

			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /transposase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /transposase/i){

				$transposase_3p_flag = 1;
                                $transposase_dist_from_end = 0;

			}elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /integrase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /integrase/i){

				$integrase_3p_flag = 1;
                                $integrase_dist_from_end = 0;

                        }elsif( $gene2annots{$start2gene{$gene_start}}{'NR'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'COG'} =~ /recombinase/i || $gene2annots{$start2gene{$gene_start}}{'NCBI'} =~ /recombinase/i){

                                $recombinase_3p_flag = 1;
                                $recombinase_dist_from_end = 0;

			}#end if

			
			$intI_dist_from_end ++;
			$transposase_dist_from_end ++;
			$integrase_dist_from_end ++;
			$recombinase_dist_from_end ++;

		}#end if


	}#end foreach


	#TAKE CARE OF POSSIBLE EDGE EFFECT IF NO GENE 3' OR 5' OF REGION IN LINEAR FORM
	if($intI_5p_flag == 0)
	{

		$intI_dist{'5_count'} = $intI_dist{'5_count'} + $intI_dist_from_end;

	}#end if

	if($intI_3p_flag == 0)
	{

		$intI_dist{'3_count'} = $intI_dist{'3_count'} + $intI_dist_from_start;

	}#end if

	if($transposase_5p_flag  == 0)
	{

		$transposase_dist{'5_count'} = $transposase_dist{'5_count'} + $transposase_dist_from_end;

	}#end if

	if($transposase_3p_flag == 0)
	{

		$transposase_dist{'3_count'} = $transposase_dist{'3_count'} + $transposase_dist_from_start;

	}#end if

	if($integrase_5p_flag  == 0)
	{

		$integrase_dist{'5_count'} = $integrase_dist{'5_count'} + $integrase_dist_from_end;

	}#end if

	if($integrase_3p_flag == 0)
	{

		$integrase_dist{'3_count'} = $integrase_dist{'3_count'} + $integrase_dist_from_start;

	}#end if

	if($recombinase_5p_flag  == 0)
	{

		$recombinase_dist{'5_count'} = $recombinase_dist{'5_count'} + $recombinase_dist_from_end;

	}#end if

	if($integrase_3p_flag == 0)
	{

		$recombinase_dist{'3_count'} = $recombinase_dist{'3_count'} + $recombinase_dist_from_start;

	}#end if

	#PRINT OVERALL STATISTICS AND GET SUMMARY OF REGIONS
	my $region_genes = join " , ", @$r_region_genes;
	my $region_annots = join " , ", @region_annots;

	print "Number of genes in region: " . scalar(@$r_region_genes) . "\n";
	print "Genes in region: $region_genes\n\n";
	print "Annotations in region: $region_annots\n\n";
	print "Transposases (5' gene dist, 3' gene dist, interior count): $transposase_dist{'5_count'}, $transposase_dist{'3_count'}, $num_transposase\n";
	print "Integrasases (5' gene dist, 3' gene dist, interior count): $integrase_dist{'5_count'}, $integrase_dist{'3_count'}, $num_integrase\n";
	print "IntI's (5' gene dist, 3' gene dist, interior count): $intI_dist{'5_count'}, $intI_dist{'3_count'}, $num_intI\n";
	print "Recombinases (5' gene dist, 3' gene dist, interior count): $recombinase_dist{'5_count'}, $recombinase_dist{'3_count'}, $num_recombinase\n";
	print "Frequency of genes in phage regions: $perc_phage_genes\n";	
	print "Frequency of genes with di-nucleotide signature > 1 SD: $region_di_nuc_gt1sd_perc\n";
	print "Region,chromosome GC: $region_GC, $chr_info{$chr_id}{'GC'}\n";
	print "Nucleotide frequency: A($region_nuc_counts{'A'}), T($region_nuc_counts{'T'}), G($region_nuc_counts{'G'}), C($region_nuc_counts{'C'}), N($N_count)\n";
	print "Contig info: $contigs (lengths = $contig_lengths, classes = $contig_descs) \n\n";

	my @TE_5 = ($transposase_dist{'5_count'}, $integrase_dist{'5_count'}, $intI_dist{'5_count'});
	my @TE_3 = ($transposase_dist{'3_count'}, $integrase_dist{'3_count'}, $intI_dist{'3_count'});


	$region_summary{'Number of Genes'} = scalar(@$r_region_genes);
	$region_summary{'TE: Number of Phage Genes'} = $num_phage_genes;
	$region_summary{'TE: Fraction of Genes With di-nuc sig > 1 SD'} = $region_di_nuc_gt1sd_perc;
	$region_summary{'TE: Number IS/Int/Tr in Region'} = $num_transposase + $num_integrase + $num_intI;
	$region_summary{'TE: Closest 5p IS/Int/Tr'} = &list_ops::min(\@TE_5);
	$region_summary{'TE: Closest 3p IS/Int/Tr'} = &list_ops::min(\@TE_3);
	$region_summary{'Base Content: GC'} = $region_GC;
	$region_summary{'Base Content: Ambiguous'} = &arith_ops::rnd($N_count/$region_length,2);
	$region_summary{'Base Content: Low Quality'} = &arith_ops::rnd($LQ_count/$region_length,2);
	

	#PRINT COG FUCNTIONS
	print "COG Functions\n";
	
	foreach my $class (sort keys %cog_summary)
	{

		my $key = " COG_$cog_code_info{$class}{'DESC'} ($class)";

		$region_summary{$key} = $cog_summary{$class};

		print "$cog_code_info{$class}{'DESC'} ($class): $cog_summary{$class}\n";

	}#end foreach

	return \%region_summary;

}#end get_genomic_region_summary


#sub get_multiple_genomic_region_summary - Takes as input: region_hash  - A reference to a hash where the first key is the start position
#									  of the region and among the second keys is 'END'
#							   chr		- The name of a chromsome from the database
#							   summary_file	- A file to write one line summaries of each region
#
#					 - Produces detailed description of each region to standard out and a one line summary is written
#					   to summary_file
sub get_multiple_genomic_region_summary
{
	my ($r_region_hash, $chr, $summary_file, $dbh) = @_;


	open SUMMARY, ">$summary_file";


	#OPTIONALLY BLAST REGIONS AGAINST GENOME DATABASE
	#GET SEQUENCES OF RECOMBINED REGIONS
	#my %rec_seq_hash;
	#
	#foreach my $start (sort {$a <=> $b} keys %rec_region_hash)
	#{
	#
	#        my $region_name =  "$chrs[0]($start-$rec_region_hash{$start}{'END'})";
	#        $rec_seq_hash{$region_name} = &genomes_db_lib::extract_sequence(&genomes_db_lib::get_chr_id($chrs[0], $dbh), $start, $rec_region_hash{$start}{'END'}, $dbh);
	#
	#}#end foreach
	#
	#BLAST
	#my $blast_output = "temp.blast";
	#
	#my $r_bac_blast_results = &blast::get_top_genome_hit($bac_genome_db, $blast_output, \%rec_seq_hash, "acinetobacter");
	#my %bac_blast_results = %$r_bac_blast_results;
	#
	#`rm $blast_output`;



	#DECLARE SOME VARIABLES
	my %region_hash = %$r_region_hash;

	my $region_counter = 0;
	my $total_rec_seq = 0;
	my $first_region = 0;
	my %chr_rec_summary;
	my $r_chr_rec_summary;

	my @mean_stats = ('TE: Closest 5p IS/Int/Tr', 'TE: Closest 3p IS/Int/Tr', 'Base Content: GC', 'Base Content: Ambiguous', 'TE: Fraction of Genes With di-nuc sig > 1 SD');
	my $r_mean_stats_hash = &hash_lib::init_hash(\@mean_stats);
	my %mean_stats_hash = %$r_mean_stats_hash;


	#GO THROUGH EACH REGION AND 1) PRINT DETAILED SUMMARY TO STANDARD OUT AND 2) ONE LINE SUMMARY TO summary_file
	foreach my $start (sort {$a <=> $b} keys %region_hash)
	{


	        #PRINT OUT SUMMMARY OF REGION
	        #PRINT OUT HEADER
	        my $region_length = $region_hash{$start}{'END'} - $start;
		my $region_name;

		if(!exists($region_hash{$start}{'NAME'}))
		{

		        $region_name =  "$chr($start-$region_hash{$start}{'END'})";

		}else{

		        $region_name = $region_hash{$start}{'NAME'};

		}#end if


	        print "-----------------------------------------------------------------------------------------------------------\n\n";
	        print "Region $region_counter: $region_name\n\n";

	        #GET SUMMARY OF REGION
	        my $r_region_summary = &genomes_db_lib::get_genomic_region_summary(&genomes_db_lib::get_chr_id($chr, $dbh), $start, $region_hash{$start}{'END'}, $dbh);
	        my %region_summary = %$r_region_summary;


	        #ADD INFO TO CUMULATIVE STATS AND PRINT ONE LINE SUMMARY
	        if($first_region == 0 && scalar(keys(%region_summary)) != 0)
	        {

	                print SUMMARY "Region Name\tRegion Length\tStart\tEnd";

	                foreach my $desc (sort keys %region_summary)
	                {

	                        print SUMMARY "\t$desc";

	                }#end foreach

	               print SUMMARY "\n";

	                $first_region = 1;

	        }#end if

	        if(scalar(keys(%region_summary)) != 0)
	        {

	                $r_chr_rec_summary = &hash_lib::hash_add(\%chr_rec_summary, \%region_summary);
	                %chr_rec_summary = %$r_chr_rec_summary;

	                print SUMMARY "$region_name\t$region_length\t$start\t$region_hash{$start}{'END'}";

	                foreach my $desc (sort keys %region_summary)
	                {

	                        print SUMMARY "\t$region_summary{$desc}";

	                }#end foreach
	
	                print SUMMARY "\n";

	                $total_length += abs($region_length);
	                $region_counter ++;

	        }#end if


		#IF BLASTIONG DONE THEN REPORT RESULTS HERE
	        #GO THROUGH BEST HITS, RETRIEVING 1) THE BEST 3 Acinetobacter HITS AND 2) THE BEST 3 NON-ACINETOBACTER HIT
	        #print "\n";
	        #
	        #for(my $i = 0; $i < 3; $i++)
	        #{
	        #
	        #       my ($host_name, $host_significance, $host_query_length, $host_hit_length,  $host_frac_aligned_hit, $host_frac_aligned_query, $host_sbj_start, $host_sbj_end, $host_hit_start, $host_hit_end) = @{$bac_blast_results{$region_name}{"HOST$i"}};
	        #
	        #       my ($other_name, $other_significance, $other_query_length, $other_hit_length,  $other_frac_aligned_hit, $other_frac_aligned_query, $other_sbj_start, $other_sbj_end, $other_hit_start, $other_hit_end) = @{$bac_blast_results{$region_name}{"OTHER$i"}};
	        #
	        #       print "\n";
	        #       print "Host hit $i: $host_name, $host_significance, $host_frac_aligned_query\n";
	        #       print "Other hit $i: $other_name, $other_significance, $other_frac_aligned_query\n";
	        #
	        #}#end for


	        print "\n-----------------------------------------------------------------------------------------------------------\n";


	}#end foreach

	#PRINT OUT FINAL LINE TO SUMMARY FILE
	print SUMMARY "TOTAL\t" . $total_length/$region_counter . "\tNA\tNA";

	foreach my $desc (sort keys %chr_rec_summary)
	{

	        if(!exists($mean_stats_hash{$desc}))
	        {

	                print SUMMARY "\t$chr_rec_summary{$desc}";

	        }else{

	                print SUMMARY "\t" . &arith_ops::rnd($chr_rec_summary{$desc}/$region_counter,2);

	        }#end if

	}#end foreach
	

	close SUMMARY;
	
}#end get_multiple_genomic_region_summary


#sub get_msa_shared_seq - Takes as input: 1) msa_id	- An ID from the multiple_alignment_info table
#					  2) min_size	- The minimum size of an aligned region to be included
#					  3) phage	- Indicates whether phage regions should be included (1) or excluded (0)
#			- Returns a 2D hash where the keys are chromsoome names and the values are the length of
#			  sequence alignmed between them in the msa
sub get_msa_shared_seq
{
	my ($msa_id, $min_size, $phage, $dbh) = @_;

	my %shared_seq_hash;
	my %shared_gene_hash;

	my $r_chr2id = &get_chr_id_all($dbh);
        my %chr2id = %$r_chr2id;


	#GET CHROMSOMES PART OF THE ALIGNMENT
	#SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $chrs = $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;


	#GET ALL ALIGNED REGIONS
        my $select_regions = "SELECT ma.starts, ma.ends
                              FROM multiple_alignment ma, multiple_alignment_info mai
                              WHERE ma.msa_id = mai.msa_id AND
				    mai.msa_id = ?";



        $sth = $dbh -> prepare($select_regions)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";



	#GO THROUGH EACH ALIGNED REGION
        REGION:while(my ($starts, $ends) =  $sth -> fetchrow_array)
	{

		my $num_genes = 0;
		my %gene_hash;

		#SPLIT UP STARTS AND ENDS
		my @starts = split /\|/, $starts;
		my @ends = split /\|/, $ends;

		
		#MAKE SURE THAT EACH REGION PASSES SIZE REQUIREMENTS AND IF DESIRED, IS NOT ANNOTATED AS PHAGE
		for(my $c = 0;$c < (@chrs); $c++)
		{

			#ONLY CHECK FOR CHROMOSOMES THAT ARE PART OF ALIGNMENT
			if($starts[$c] == 0)
			{

				next;

			}#end if


			#CHECK IF SIZE IS LESS THAN MINIMUM
			my $size = abs($ends[$c] - $starts[$c]);

			if($size < $min_size)
			{

				#print "REGION SIZE: $size\n";
				next REGION;

			}#end if
	
	
			#CHECK IF SEQUENCE IS AMBIGUOUS
			my $region_seq = &extract_sequence($chr2id{$chrs[$c]}, $starts[$c], $ends[$c], $dbh);
        		my $region_seq_obj = Bio::Seq->new(-seq => $region_seq);

        		my $r_region_nuc_counts = &seq_stats::init_monoNucHash();
        		$r_region_nuc_counts = &seq_stats::get_monoNucCount($r_region_nuc_counts, $region_seq_obj);
        		my %region_nuc_counts = %$r_region_nuc_counts;

        		my $N_count = $size - ($region_nuc_counts{'A'} + $region_nuc_counts{'C'} + $region_nuc_counts{'G'}  + $region_nuc_counts{'T'});

			if(($N_count / $size) > 0.5 )
			{

				#print "AMBIG BASES: $chrs[$c], $N_count, $size, $starts\n";
				next REGION;

			}#end if


			#CHECK IF THERE ARE ANY GENES IN THE REGION
			#CHECK THE STRAND AND ADJUST STATEMENT ACCORDINGLY
			my $r_region_genes;

			if($starts[$c] < $ends[$c])
			{

				#$r_region_genes = &get_chr_region_genes($chr2id{$chrs[$c]}, $starts[$c], $ends[$c], $dbh);

			}else{

				#$r_region_genes = &get_chr_region_genes($chr2id{$chrs[$c]}, $ends[$c], $starts[$c], $dbh);

			}#end if

			$gene_hash{$chrs[$c]} = $r_region_genes;
			$num_genes = scalar(@$r_region_genes);

			if( $num_genes == 0 )
			{

				#print "NUMBER OF GENES in $chrs[$c] ($size):" . @$r_region_genes . ", $starts\n";
				#next REGION;

			}#end if


			#CHECK IF REGION IS ANNOTATED AS PHAGE, IF GREATER THAN 50% OF GENES IN REGION ARE, THEN SKIP
			if($phage)
			{

				my $r_phage_genes = &is_gene_phage($r_region_genes, $dbh);
			        my %phage_genes = %$r_phage_genes;

        			my $num_phage_genes = &hash_lib::hash_sum(\%phage_genes);

				if(($num_phage_genes / scalar(@$r_region_genes)) > 0.4)
				{

					next REGION;

				}#end if

			}#end if

		}#end for


		#FOREACH PAIR OF CHROMSOMES TALLLY SHARED SEQUENCE
		for(my $c1 = 0;$c1 < (@chrs); $c1++)
		{

			#UPDATE TALLY WITHIN CURRENT GENOME
			if($starts[$c1] != 0)
			{

				$shared_seq_hash{$chrs[$c1]}{$chrs[$c1]} += abs($starts[$c1] - $ends[$c1]);

				push @{$shared_gene_hash{$chrs[$c1]}{$chrs[$c1]}}, @{$gene_hash{$chrs[$c1]}};

			}#end if


			#KEEP TRACK OF ALIGNED SEQUENCE FOR PAIRS
			for(my $c2 = 0; $c2 < $c1; $c2++)
			{

				if( ($starts[$c1] != 0) && ($starts[$c2] != 0) )
				{

					$shared_seq_hash{$chrs[$c1]}{$chrs[$c2]} += abs($starts[$c1] - $ends[$c1]);
					$shared_seq_hash{$chrs[$c2]}{$chrs[$c1]} += abs($starts[$c2] - $ends[$c2]);

					push @{$shared_gene_hash{$chrs[$c1]}{$chrs[$c2]}}, @{$gene_hash{$chrs[$c1]}};
					push @{$shared_gene_hash{$chrs[$c2]}{$chrs[$c1]}}, @{$gene_hash{$chrs[$c2]}};

				}#end if

			}#end for

		}#end for		

	}#end while


	foreach my $chr1(@chrs)
	{
		print $chr1;

		foreach my $chr2(@chrs)
		{

			print "\t$shared_seq_hash{$chr1}{$chr2}";

		}#end foreach

		print "\n";
	}#end foreach

	return (\%shared_seq_hash, \%shared_gene_hash);

}#get_msa_shared_seq

############################################################################### SNP SELECTION ROUTINES #############################################################


#sub get_snps - Takes as input: 1) chr1	- A name from the chromsome table
#				2) chr2	- A second name from the chromosome table
#	      - Returns a hash whose first key is the position in chr1 of the SNP and second keys are:
# 		1) position in chromosome 2, 2) base in chromsome 1 and 3) base in chromosome 2
sub get_snps
{
	my ($chr1, $chr2, $dbh) = @_;

	my %snp_info;


	#GET SNPs FROM THE snp TABLE
	my $select_snps = "SELECT base1, position1, base2, position2, SUBSTRING(c1.sequence, position1 - 5, 11), SUBSTRING(c2.sequence, position2 - 5, 11), SUBSTRING(c1.sequence, position1 - 20, 41), SUBSTRING(c2.sequence, position2 - 20, 41)
                           FROM snp s, chromosome c1, chromosome c2, contig c
                           WHERE c1.name = ? AND 
			         c2.name = ? AND
				 s.chr_id1 = c1.chr_id AND
				 s.chr_id2 = c2.chr_id AND
			  	 c.chr_id = c1.chr_id AND
       				 ((c.start < position1 AND c.end> position1) OR (c.start > position1 AND c.end < position1)) AND
      				 (c.contig_name NOT LIKE 'scaffold%')
			    ORDER BY position1 ASC";

        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2)
                         or die "Can't execute statment: $DBI::errstr";


	#CHECK IF QUERY YIELDED RESULTS, IF NOT, INVERT QUERY
	if($sth->rows == 0)
	{

		$select_snps = "SELECT base2, position2, base1, position1, SUBSTRING(c2.sequence, position1 - 5, 11), SUBSTRING(c1.sequence, position2 - 5, 11), SUBSTRING(c2.sequence, position2 - 20, 41), SUBSTRING(c1.sequence, position1 - 20, 41)
	                           FROM snp s, chromosome c1, chromosome c2, contig c
	                           WHERE c1.name = ? AND 
				         c2.name = ? AND
					 s.chr_id1 = c1.chr_id AND
					 s.chr_id2 = c2.chr_id AND
				  	 c.chr_id = c1.chr_id AND
       					 ((c.start < position2 AND c.end> position2) OR (c.start > position2 AND c.end < position2)) AND
      					 (c.contig_name NOT LIKE 'scaffold%')
				    ORDER BY position1 ASC";

        	$sth = $dbh -> prepare($select_snps)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	$rc = $sth -> execute($chr2, $chr1)
                	         or die "Can't execute statment: $DBI::errstr";


	}#end if


	my $indel_pos = 0;
	my $prev_pos = 0;
	my $indel_seq = "";
	my $in_seq = 0;
	my $indel = 0;

        while(my ($base1, $position1, $base2, $position2, $sub_seq1, $sub_seq2, $ext_sub_seq1, $ext_sub_seq2) =  $sth -> fetchrow_array)
        {

		#IGNORE SNP IF 1) AFTER HOMOPOLYMER, 2) LOW QUALITY BASE OR 3) AMBIGOUS NUCLEOTIDE
		#if( $sub_seq1 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i || $sub_seq2 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i)
		if( ($sub_seq1 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i && ($base1 eq '.' || $base2 eq '.')) || 
		    ($sub_seq2 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i && ($base1 eq '.' || $base2 eq '.')) || 
		    $sub_seq1 =~ /^.{3}.*[actgN].*.{3}$/ || 
		    $sub_seq2 =~ /^.{3}.*[actgN].*.{3}$/)
		{

			next;

		}#end if


		#CHECK IF SNP PART OF INDEL
		if( $indel > 0 && (($in_seq == 1 && $position1 == ($prev_pos + 1) && $base2 eq '.' ) || ($in_seq == 2 && $position2 == ($prev_pos + 1) && $base1 eq '.') ) )
		{


			#INCREMENT SIZE OF INDEL, NOTE THE CURRENT POSITION IN THE INSERTED SEQUENCE AND MOVE TO NEXT SNP
			$indel ++;

			if($in_seq == 1)
			{

				$prev_pos = $position1;
				$indel_seq = $indel_seq . $base1;

			}else{

				$prev_pos = $position2;
				$indel_seq = $indel_seq . $base2;

			}#end if

			next;
	
		}elsif( $indel > 0 &&  ($base1 eq '.' || $base2 eq '.') ){

			#THIS IS A NEW INDEL, SO SAVE THE PREVIOUS ONE
			$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
			$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;


			#SET UP VARIABLES FOR NEW INDEL
			$indel = 1;

			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if
	


		}elsif( $indel == 0 && ($base1 eq '.' || $base2 eq '.') ){

			#SET UP VARIABLES DESCRIBING INDEL
			$indel = 1;


			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if
	

		}elsif( $indel > 0){

			#SAVE INFO ON PREVIOUS INDEL
			$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
			$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;

			$indel = 0;
			
		}#end if


		#RECORD INFO ON SNP
		$snp_info{$position1}{'BASE1'} = $base1;
		$snp_info{$position1}{'BASE2'} = $base2;
		$snp_info{$position1}{'POS2'} = $position2;
		$snp_info{$position1}{'SEQ1'} = $sub_seq1;
		$snp_info{$position1}{'SEQ2'} = $sub_seq2;
		$snp_info{$position1}{'EXT_SEQ1'} = $ext_sub_seq1;
		$snp_info{$position1}{'EXT_SEQ2'} = $ext_sub_seq2;


	}#end while


	#STORE THE LAST INDEL IF ONE EXISTS
	if( $indel > 0 )
	{

		$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
		$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;

	}#end if

	return \%snp_info;

}#end get_snps


#sub get_pw_snp_counts - Takes as input 
#		       - Returns pairwise snp counts for each pair of chromosomes in the database
sub get_pw_snp_counts
{
	my ($dbh) = @_;


	#GET PAIRWISE SNP COUNTS FROM TABLE
	my $select_counts = "SELECT chr_id1, chr_id2, COUNT(*) 
			     FROM snp 
			     WHERE base1 <> '.' AND base2 <> '.'
			     GROUP BY chr_id1, chr_id2";

        my $sth = $dbh -> prepare($select_counts)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute()
                         or die "Can't execute statment: $DBI::errstr";

	#GET RESULTS
	my %snp_counts;

	while(my ($chr1, $chr2, $num_snps) = $sth -> fetchrow_array)
	{

		#print "$chr1, $chr2, $num_snps\n";
		$snp_counts{$chr1}{$chr2} = $num_snps;

	}#end while


	return \%snp_counts;

}#end get_pw_snp_counts


#sub pw_snp_count_chr_reduce - Takes as input: 1) chrs 	   - An initial set of chromsomes
#					       2) min_snps - The minimum number of SNPs to already selected chromosome
#							     to be included in list
#			     - Returns a reduced list of chromosomes based on the criteria that no two chromosmes have less
#			      than a desired number of snps
sub pw_snp_count_chr_reduce
{
	my ($r_chrs, $min_snps, $dbh) = @_;


	#GET CHROMOSOME IDS
	my $r_chr_ids = &get_chr_id_all($dbh);
	my %chr_ids = %$r_chr_ids;
	#my $r_chr_desc = &get_chr_desc_all($dbh);
	#my %chr_desc = %$r_chr_desc;


	#GET PAIRWISE SNP COUNTS
	my $r_pw_snp_counts = &get_pw_snp_counts($dbh);
	my %pw_snp_counts = %$r_pw_snp_counts;


	#GO THROUGH CHROMOSOMES, DISCARDING THOSE WITH LESS THAN THRESHOLD SNPS TO ALREADY SELECTED CHROMSOME
	my @chr_red = ();

	CHR:foreach $chr (@$r_chrs)
	{

		foreach $sel_chr (@chr_red)
		{

			if((defined($pw_snp_counts{$chr_ids{$chr}}{$chr_ids{$sel_chr}}) && $pw_snp_counts{$chr_ids{$chr}}{$chr_ids{$sel_chr}} < $min_snps) || 
			   (defined($pw_snp_counts{$chr_ids{$sel_chr}}{$chr_ids{$chr}}) && $pw_snp_counts{$chr_ids{$sel_chr}}{$chr_ids{$chr}} < $min_snps) )
			{			

				next CHR;

			}#end if

		}#end foreach
		
		push @chr_red, $chr;
		print "$chr\n";

	}#end for

	return \@chr_red;

}#end pw_snp_count_chr_reduce


#sub pw_snp_count_cluster    - Takes as input: 1) pw_snp_counts	- A 2D hash where the keys are chromosome names and the values are SNP differents
#						             	between the chromosomes
#					       2) min_snps 	- The minimum number of SNPs to already selected chromosome
#							     	to be included in list
#			     - 
#
sub pw_snp_count_cluster
{
	my ($r_pw_snp_counts, $min_snps) = @_;

	my %pw_snp_counts = %$r_pw_snp_counts;
	my @chrs = keys %pw_snp_counts;


	#GO THROUGH CHROMOSOMES, MAKE INITIAL ASSIGNMENTS TO CLUSTERS
	my %clusters;
	my %chr2clusters;
	my $cluster_num = 0;

	CHR:foreach $chr (@chrs)
	{

		#CHECK IF CURRENT CHROMOSOME IS CLOSE TO ANY CHROMOSOME ALREADY ASSIGNED TO A CLUSTER
		foreach my $clust_chr (keys %chr2clusters)
		{

			if($pw_snp_counts{$chr}{$clust_chr} <= $min_snps)
			{

				$chr2clusters{$chr} = $chr2clusters{$clust_chr};
				push @{$clusters{$chr2clusters{$clust_chr}}}, $chr;
				next CHR;

			}#end if
			

		}#end foreach		

		#CHROSOME NOT CLOSE TO OTHERS SEEN, SO ASSIGN TO OWN CLUSTER
		$cluster_num ++;

		$chr2clusters{$chr} = $cluster_num;
		@{$clusters{$cluster_num}} = ($chr);

	}#end foreach


	#MERGE CLUSTERS
	#GO THROUGH EACH PAIR OF CHROMOSOMES AND MERGE CLUSTERS IF THEY ARE LESS THAN MINIUMUM DISTANCE APART
	foreach my $chr1 (@chrs)
	{

		foreach my $chr2 (@chrs)
		{

			#IF TWO CHROMOSOMES ARE CLOSER THAN MINIMUM DISTANCE AND IN DIFFERENT CLUSTERS THEN MERGE
			if(($pw_snp_counts{$chr1}{$chr2} <= $min_snps) && ($chr2clusters{$chr1} != $chr2clusters{$chr2}))
			{

				$cluster1 = $chr2clusters{$chr1};
				$cluster2 = $chr2clusters{$chr2};

				@{$clusters{$cluster1}} = (@{$clusters{$cluster1}}, @{$clusters{$cluster2}});

				@chr2clusters{@{$clusters{$cluster2}}} = $cluster1;

				delete $clusters{$cluster2};

			}#end if

		}#end foreach

	}#end foreach
		

	return \%clusters;

}#end pw_snp_count_cluster


#sub get_all_snps - Takes as input: 1) chr1	- A name from the chromsome table
#	  			    2) chr2	- A second name from the chromosome table
#	          - Returns a hash whose first key is the position in chr1 of the SNP and second keys are:
# 	  	    1) position in chromosome 2, 2) base in chromsome 1 and 3) base in chromosome 2
sub get_all_snps
{
	my ($chr1, $chr2, $dbh) = @_;

	my %snp_info;


	#GET SNPs FROM THE snp TABLE
	my $select_snps = "SELECT base1, position1, base2, position2, SUBSTRING(c1.sequence, position1 - 5, 11), SUBSTRING(c2.sequence, position2 - 5, 11), SUBSTRING(c1.sequence, position1 - 20, 41), SUBSTRING(c2.sequence, position2 - 20, 41)
                             FROM snp s, chromosome c1, chromosome c2
                             WHERE c1.name = ? AND 
				   c2.name = ? AND
				   s.chr_id1 = c1.chr_id AND
				   s.chr_id2 = c2.chr_id
			    ORDER BY position1 ASC";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2)
                         or die "Can't execute statment: $DBI::errstr";


        #CHECK IF QUERY YIELDED RESULTS, IF NOT, INVERT QUERY
        if($sth->rows == 0)
        {

                $select_snps = "SELECT base2, position2, base1, position1, SUBSTRING(c2.sequence, position2 - 5, 11), SUBSTRING(c1.sequence, position1 - 5, 11), SUBSTRING(c2.sequence, position2 - 20, 41), SUBSTRING(c1.sequence, position1 - 20, 41)
                                   FROM snp s, chromosome c1, chromosome c2, contig c
                                   WHERE c1.name = ? AND 
                                         c2.name = ? AND
                                         s.chr_id1 = c1.chr_id AND
                                         s.chr_id2 = c2.chr_id AND
                                         c.chr_id = c1.chr_id AND
                                         ((c.start < position2 AND c.end> position2) OR (c.start > position2 AND c.end < position2)) AND
                                         (c.contig_name NOT LIKE 'scaffold%')
                                    ORDER BY position2 ASC";

                $sth = $dbh -> prepare($select_snps)
                                 or die "Can't prepare statment: $DBI::errstr";

                $rc = $sth -> execute($chr2, $chr1)
                                 or die "Can't execute statment: $DBI::errstr";


        }#end if


	my $indel_pos = 0;
	my $prev_pos = 0;
	my $indel_seq = "";
	my $in_seq = 0;
	my $indel = 0;

        while(my ($base1, $position1, $base2, $position2, $sub_seq1, $sub_seq2, $ext_sub_seq1, $ext_sub_seq2) =  $sth -> fetchrow_array)
        {


		#CHECK IF SNP PART OF INDEL
		if( $indel > 0 && (($in_seq == 1 && $position1 == ($prev_pos + 1) && $base2 eq '.' ) || ($in_seq == 2 && $position2 == ($prev_pos + 1) && $base1 eq '.') ) )
		{

			#INCREMENT SIZE OF INDEL, NOTE THE CURRENT POSITION IN THE INSERTED SEQUENCE AND MOVE TO NEXT SNP
			$indel ++;

			if($in_seq == 1)
			{

				$prev_pos = $position1;
				$indel_seq = $indel_seq . $base1;

			}else{

				$prev_pos = $position2;
				$indel_seq = $indel_seq . $base2;

			}#end if

			next;
	
		}elsif( $indel > 0 &&  ($base1 eq '.' || $base2 eq '.') ){

			#THIS IS A NEW INDEL, SO SAVE THE PREVIOUS ONE
			$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
			$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;


			#SET UP VARIABLES FOR NEW INDEL
			$indel = 1;

			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if
	


		}elsif( $indel == 0 && ($base1 eq '.' || $base2 eq '.') ){

			#SET UP VARIABLES DESCRIBING INDEL
			$indel = 1;


			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if
	

		}elsif( $indel > 0){

			#SAVE INFO ON PREVIOUS INDEL
			$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
			$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;

			$indel = 0;
			
		}#end if


		#RECORD INFO ON SNP
		$snp_info{$position1}{'BASE1'} = $base1;
		$snp_info{$position1}{'BASE2'} = $base2;
		$snp_info{$position1}{'POS2'} = $position2;
		$snp_info{$position1}{'SEQ1'} = $sub_seq1;
		$snp_info{$position1}{'SEQ2'} = $sub_seq2;
		$snp_info{$position1}{'EXT_SEQ1'} = $ext_sub_seq1;
		$snp_info{$position1}{'EXT_SEQ2'} = $ext_sub_seq2;


	}#end while


	#STORE THE LAST INDEL IF ONE EXISTS
	if( $indel > 0 )
	{

		$snp_info{$indel_key}{'INDEL_LENGTH'} = $indel;
		$snp_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;

	}#end if

	return \%snp_info;

}#end get_all_snps


#sub check_pos_across_chrs - Takes as input: 1) chr	- A chromsome in the database
#					     2) pos	- A list of positions in the chromosome
#			   - Returns a hash where the first key are positions in the reference chromosome and the second key are chromosomes 
#			     which have alignments covering that position
sub check_pos_across_chrs
{
	my ($chr, $r_pos, $dbh) = @_;

	my %pos_pres;

	my $r_chrid2name = &get_chr_desc_all($dbh);
	my %chrid2name = %$r_chrid2name;
	my $ref_chr_id = &get_chr_id($chr, $dbh);

	foreach my $pos (@$r_pos)
	{

		#GET ALIGNMENTS FROM DATABASE COVERING POSITIONS OF INTEREST
		my $select_aln = "SELECT p.chr_id1, p.chr_id2
				  FROM pairwise_alignment p 
				  WHERE (chr_id1 = ? AND
              				((start1 <= ? AND end1 >= ?) OR (start1 >= ? AND end1 <= ?))) 
              				OR
              				(chr_id2 = ? AND
              				((start2 <= ? AND end2 >= ?) OR (start2 >= ? AND end2 <= ?)))";
		#my $select_aln = "SELECT c1.name, c2.name
		#		  FROM pairwise_alignment p, chromosome c1, chromosome c2 
		#		  WHERE (c1.chr_id = p.chr_id1 AND c2.chr_id = p.chr_id2 AND
               	#			 c1.name = ? AND
              	#			((start1 <= ? AND end1 >= ?) OR (start1 >= ? AND end1 <= ?))) 
              	#			OR
              	#			(c1.chr_id = p.chr_id2 AND c2.chr_id = p.chr_id1 AND
               	#			c1.name = ? AND
              	#			((start2 <= ? AND end2 >= ?) OR (start2 >= ? AND end2 <= ?)))";

        	my $sth = $dbh -> prepare($select_aln)
                         or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($ref_chr_id, $pos, $pos, $pos, $pos, $ref_chr_id, $pos, $pos, $pos, $pos)
                         or die "Can't execute statment: $DBI::errstr";
	
		#KEEP TRACK OF CHROMOSOMES COVERING DIFFERENT POSITIONS
        	while(my ($ref_chr, $aln_chr) =  $sth -> fetchrow_array)
        	{

			if($chrid2name{$ref_chr}{'NAME'} eq $chr)
			{

				$pos_pres{$pos}{$chrid2name{$aln_chr}{'NAME'}} = 1;

			}else{

				$pos_pres{$pos}{$chrid2name{$ref_chr}{'NAME'}} = 1;

			}#end if

		}#end while

	}#end foreach

	return \%pos_pres;

}#end check_pos_across_chrs


#sub get_aln_pos - Takes as input: 1) ref_chr	- A chromsome in the database
#				   2) aln_chrs	- Chromosomes to check for alignment against reference chromosme
#			   - Returns a hash where the first key are positions in the reference chromosome and the second key are chromosomes 
#			     which have alignments covering that position
sub get_aln_pos
{
	my ($ref_chr, $r_aln_chrs, $dbh) = @_;

	my %pos_pres;

	my $r_chrid2name = &get_chr_desc_all($dbh);
	my %chrid2name = %$r_chrid2name;
	my $ref_chr_id = &get_chr_id($ref_chr, $dbh);

	my $r_chr_name2id = &get_chr_id_all($dbh);
        my %chr_name2id = %$r_chr_name2id;

	foreach my $chr (@$r_aln_chrs)
	{

		#GET ALIGNMENTS FROM DATABASE COVERING POSITIONS OF INTEREST
		my $select_aln = "SELECT p.chr_id1, p.start1, p.end1, p.chr_id2, p.start2, p.end2
				  FROM pairwise_alignment p 
				  WHERE (chr_id1 = ? AND chr_id2 = ?) OR
					(chr_id1 = ? AND chr_id2 = ?)";

        	my $sth = $dbh -> prepare($select_aln)
                         or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($ref_chr_id, $chr_name2id{$chr}, $chr_name2id{$chr}, $ref_chr_id)
                         or die "Can't execute statment: $DBI::errstr";
	
		#KEEP TRACK OF CHROMOSOMES COVERING DIFFERENT POSITIONS
        	while(my ($chr_id1, $start1, $end1, $chr_id2, $start2, $end2) =  $sth -> fetchrow_array)
        	{

			if($chrid2name{$chr_id1}{'NAME'} eq $ref_chr)
			{
		
				for(my $i = &arith_ops::min($start1, $end1); $i <= &arith_ops::max($start1, $end1); $i++)
				{

					$pos_pres{$i}{$chrid2name{$chr_id2}{'NAME'}} = 1;

				}#end for
				
			}else{

				for(my $i = &arith_ops::min($start2, $end2); $i <= &arith_ops::max($start2, $end2); $i++)
				{

					$pos_pres{$i}{$chrid2name{$chr_id1}{'NAME'}} = 1;

				}#end for

			}#end if

		}#end while

	}#end foreach


	return \%pos_pres;

}#end get_aln_pos


#sub get_aln_pos_pw - Takes as input: 1) ref_chr	- A chromsome in the database
#				      2) aln_chr	- Chromosome to check for alignment against reference chromosme
#			   - Returns a hash where the key are positions in the reference chromosome
sub get_aln_pos_pw
{
	my ($ref_chr, $aln_chr, $dbh) = @_;

	my %pos_pres;

	my $r_chrid2name = &get_chr_desc_all($dbh);
	my %chrid2name = %$r_chrid2name;
	my $ref_chr_id = &get_chr_id($ref_chr, $dbh);
	my $aln_chr_id = &get_chr_id($aln_chr, $dbh);



	#GET ALIGNMENTS FROM DATABASE COVERING POSITIONS OF INTEREST
	my $select_aln = "SELECT p.chr_id1, p.start1, p.end1, p.chr_id2, p.start2, p.end2
			  FROM pairwise_alignment p 
			  WHERE (chr_id1 = ? AND chr_id2 = ?) OR
				(chr_id1 = ? AND chr_id2 = ?)";

       	my $sth = $dbh -> prepare($select_aln)
                  or die "Can't prepare statment: $DBI::errstr";

       	my $rc = $sth -> execute($ref_chr_id, $aln_chr_id, $aln_chr_id, $ref_chr_id)
                   or die "Can't execute statment: $DBI::errstr";
	
	#KEEP TRACK OF CHROMOSOMES COVERING DIFFERENT POSITIONS
       	while(my ($chr_id1, $start1, $end1, $chr_id2, $start2, $end2) =  $sth -> fetchrow_array)
       	{

		if($chrid2name{$chr_id1}{'NAME'} eq $ref_chr)
		{
		

			if($start1 < $end1)
			{

				@pos_pres{$start1..$end1} = 1;

			}else{

				@pos_pres{$end1..$start1} = 1;

			}#end if
				
		}else{

			if($start2 < $end2)
			{

				@pos_pres{$start2..$end2} = 1;

			}else{

				@pos_pres{$end2..$start2} = 1;

			}#end if

		}#end if

	}#end while



	return \%pos_pres;

}#end get_aln_pos_pw


#sub get_non_snp_pos - Takes as input: 1) ref_chr	- A chromsome in the database
#				   	2) aln_chrs	- Chromosomes to get non-snp positions
#			   - Returns a hash where the key is positions in the reference chromosome that : 1) have alignments
#			     to all chromosomes at that position and 2) do not have SNPs at any position  
sub get_non_snp_pos
{
	my ($ref_chr, $r_aln_chrs, $dbh) = @_;


	my $r_chrid2name = &get_chr_desc_all($dbh);
	my %chrid2name = %$r_chrid2name;
	my $ref_chr_id = &get_chr_id($ref_chr, $dbh);

	my $r_aln_chrs_hash = &hash_lib::init_hash($r_aln_chrs);
	my %aln_chrs = %$r_aln_chrs_hash;

	my $r_chr_name2id = &get_chr_id_all($dbh);
        my %chr_name2id = %$r_chr_name2id;

	my %snp_pos;


	#GET POSITIONS IN DIFFERENT CHROMOSOMES ALIGNED TO REFERENCE
	my $r_pos_pres = &get_aln_pos($ref_chr, $r_aln_chrs, $dbh);
	my %pos_pres = %$r_pos_pres;


	#GET ALL SNPS RELATIVE TO REFERENCE
	foreach my $chr (@$r_aln_chrs)
	{

		my $select_aln = "SELECT chr_id1, position1, chr_id2, position2
				  FROM snp
				  WHERE (chr_id1 = ?  AND chr_id2 = ?) OR (chr_id2 = ? AND chr_id1 = ?)";

	       	my $sth = $dbh -> prepare($select_aln)
        	                or die "Can't prepare statment: $DBI::errstr";
	
       		my $rc = $sth -> execute($ref_chr_id, $chr_name2id{$chr}, $ref_chr_id, $chr_name2id{$chr})
        	                or die "Can't execute statment: $DBI::errstr";
	
	
       		while(my ($chr_id1, $pos1, $chr_id2, $pos2) =  $sth -> fetchrow_array)
		{

			if($chr_id1 == $ref_chr_id)
			{
	
				$snp_pos{$pos1} = 1;

			}else{

				$snp_pos{$pos2} = 1;

			}#end if

		}#end while


	}#end foreach


	#GO THROUGH ALIGNED POSITIONS AND SAVE THOSE PRESENT ACROSS ALL CHRS AND LACKING A SNP
	foreach my $pos (keys %pos_pres)
	{

		my @pres_chrs = keys(%{$pos_pres{$pos}});
		
		if(!exists($snp_pos{$pos}) & (@pres_chrs  == @$r_aln_chrs) )
		{

			$non_snp_pos{$pos} = 1;

		}#end if


	}#end foreach


	return \%non_snp_pos;

}#end get_non_snp_pos



#sub pw_aln_in_db - Takes as input: 1) chrs - A list of chromosmes, with the first being the reference
#				    
#		  - Returns chromsoomes with pairwise alignments to the reference in the database
sub pw_aln_in_db
{
	my ($r_chrs, $dbh) = @_;

	my @old_chrs = @$r_chrs;
	my @chrs = ($old_chrs[0]);

	my $ref_chr_id = &get_chr_id($chrs[0], $dbh);


	#KEEP CHROMSOMES WITH ALIGNMENTS TO REFERENCE
	for(my $c = 1;$c < @old_chrs; $c++)
	{

		my $chr = $old_chrs[$c];

		#GET ALIGNMENTS FROM DATABASE COVERING POSITIONS OF INTEREST
		my $select_aln = "SELECT p.chr_id1
				  FROM pairwise_alignment p 
				  WHERE (chr_id1 = ? AND chr_id2 = ?) OR
				   (chr_id1 = ? AND chr_id2 = ?)";

        	my $sth = $dbh -> prepare($select_aln)
                         or die "Can't prepare statment: $DBI::errstr";

        	my $rc = $sth -> execute($ref_chr_id, &get_chr_id($chr, $dbh), &get_chr_id($chr, $dbh), $ref_chr_id)
                         or die "Can't execute statment: $DBI::errstr";
	
		if ($sth->rows > 0)
		{

			push @chrs, $chr;

		}else{

			print "\n\n$chr has no alignment in the database!!\n\n";

		}#end if

	}#end foreach	

	return \@chrs;

}#end pw_aln_in_db


#sub filter_snps_ref - Takes as input: 1) ref_chr       - A reference chromsome relative to which SNPs are detected
#				       2) snps		- A hash where the fist key is a chromsome, the second key is a position in
#							  the reference chromsome and the value is the position on the key chromsome
#							  that has a SNP with the given position on the reference
#				       3) rep_chrs	- A list of chromsomes which repative elements are found and filtered
#				       4) min_dist	- The minimum distance allowed between two SNPs before they are filtered
sub filter_snps_ref
{
	my ($ref_chr, $r_snps, $r_rep_chrs, $min_dist, $dbh) = @_;

	my %snps = %$r_snps;
	my @rep_chrs = @$r_rep_chrs;
	my @my_chrs = ($ref_chr, keys %snps);


	my %filt_snps;


	#GET LISTS OF SNPS THAT SHOULD BE FILTERED OUT
	my %snp_filter_function;
	my %snp_filter_contig;

	my @snp_filter_cluster;
	my @snp_filter_context;
	my %snp_filter_cluster;
	my %snp_filter_context;

	my @snp_filter_large_indels;
	my %snp_filter_large_indels;

	my %snp_filter_repeat;
	my %snp_filter_tandem_repeat;
	
	my %snp_filter_depth;
	my %snp_filter_quality;


	#DETERMINE POSITIONS IN EACH CHROMSOME THAT ARE AT END OF CONTIG OR HAVE FUNCTION OF REPATIVE ELEMENT
	foreach my $chr (@my_chrs)
	{

		print "1\n";
		$snp_filter_function{$chr} = &filter_snps_function($chr, $dbh);
		print "2a\n";
		$snp_filter_contig{$chr} = &filter_snps_contig($chr, $dbh);
		print "2b\n";
		#$snp_filter_quality{$chr} = &filter_snps_quality($chr, $dbh);
		print "2c\n";
		$snp_filter_depth{$chr} = &filter_snps_depth($chr, $dbh);
		

	}#end foreach

	#DETERMINE SNP POSITIONS IN REFERENCE CHROMSOME THAT ARE 1) ALIGNED TO DUBIOUS SEQUENCE IN ANY CHROMSOME OR
	#2) CLUSTER TOGETHER IN THE REFERENCE CHROMSOME, OR 3) ARE NOT PRESENT IN ANOTHER CHROMOSOME
	foreach my $chr (keys %snps)
	{

		print "3\n";
		my ($r_chr1_pos, $r_chr2_pos) = &filter_snps_cluster($ref_chr, $chr,  $min_dist, $dbh);
		@snp_filter_cluster =  keys(%$r_chr1_pos);
		$snp_filter_cluster{$chr} = &hash_lib::init_hash(\@snp_filter_cluster);

		print "4\n";
		($r_chr1_pos, $r_chr2_pos) = &filter_snps_context($ref_chr, $chr,  $dbh);
		@snp_filter_context = keys(%$r_chr1_pos);
		$snp_filter_context{$chr} = &hash_lib::init_hash(\@snp_filter_context);


	}#end foreach

	#DETERMINE POSITIONS OF REPEATS IN CHROMSOMES THAT HAVE BEEN DESIGNATED TO SEARCH IN
	foreach my $chr (@rep_chrs)
	{

		print "5\n";
		$snp_filter_repeat{$chr} = &filter_snps_repeats($chr,  $dbh);
		print "6\n";
		$snp_filter_tandem_repeat{$chr} = &filter_snps_tandem_repeats($chr,  $dbh);

	}#end foreach


	#SCREEN INPUTTED LIST OF SNPS
	my %filt_snp_pos;

	for(my $c=1; $c < @my_chrs; $c++)
	{

		my %filt_snp_count;
		my $filt_snp_flag = 0;

		foreach my $pos(keys %{$snps{$my_chrs[$c]}})
		{


			#FILTER SNPS BASED ON THERE LOCAL CONTEXT
			if(exists($snp_filter_context{$my_chrs[$c]}{$pos}))
			{
			
				$filt_snp_flag = 1;
				$filt_snp_pos{'CONTEXT'}{$pos} = 1;
				$filt_snp_count{'CONTEXT'} ++;
	
			}#end if

			#FILTER SNPS BASED ON LOCAL DENSITY
			if(exists($snp_filter_cluster{$my_chrs[$c]}{$pos}))
			{
			
				$filt_snp_flag = 1;
				$filt_snp_pos{'CLUSTER'}{$pos} = 1;
				$filt_snp_count{'CLUSTER'} ++;

			}#end if

			#FILTER SNPS BASED ON THEIR BEING IN PHAGE OR REPATIVE ELEMENTS
			if(exists($snp_filter_function{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}})  || exists($snp_filter_function{$my_chrs[0]}{$pos}) )
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'FUNCTION'}{$pos} = 1;
				$filt_snp_count{'FUNCTION'} ++;

			}#end foreach

			#FILTER SNPS BASED ON BEING NEAR THE BEGINNING OR END OF A CONTIG
			if(exists($snp_filter_contig{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_contig{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'CONTIG'}{$pos} = 1;
				$filt_snp_count{'CONTIG'} ++;

			}#end foreach
		
			#FILTER SNPS BASED ON QUALITY
			if(exists($snp_filter_quality{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_quality{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'QUAL'}{$pos} = 1;
				$filt_snp_count{'QUAL'} ++;

			}#end foreach
		
			#FILTER SNPS BASED ON DEPTH
			if(exists($snp_filter_depth{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_depth{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'DEPTH'}{$pos} = 1;
				$filt_snp_count{'DEPTH'} ++;

			}#end foreach
		
			#FILTER SNPS THAT ARE IN REPETATIVE REGIONS
			if(exists($snp_filter_repeat{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_tandem_repeat{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_repeat{$my_chrs[0]}{$pos}) || exists($snp_filter_tandem_repeat{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'REPEAT'}{$pos} = 1;
				$filt_snp_count{'REPEAT'} ++;

			}#end foreach


			#CHECK IF FILTERS PASSED, IF SO PUT CURRENT SNP INTO HASH
			if($filt_snp_flag == 0)
			{

				$filt_snps{$my_chrs[$c]}{$pos} = $snps{$my_chrs[$c]}{$pos};

			}else{

				$filt_snp_flag = 0;

			}#end if

		}#end foreach

		print "\n\n$my_chrs[$c] SNPS filtered out: $filt_snp_count{'CONTEXT'}(context), $filt_snp_count{'CLUSTER'}(cluster), $filt_snp_count{'CONTIG'}(contig), $filt_snp_count{'FUNCTION'}(function), $filt_snp_count{'REPEAT'}(repeat), $filt_snp_count{'DEPTH'}(depth), $filt_snp_count{'QUAL'}(quality)\n SNPs remaining: " . scalar(keys(%{$filt_snps{$my_chrs[$c]}}))  . "\n\n";

	}#end foreach

	
	return \%filt_snps;

}#end filter_snps_ref


#sub filter_snps_ref_w_annot - Takes as input: 1) ref_chr       - A reference chromsome relative to which SNPs are detected
#					       2) snps		- A hash where the fist key is a chromsome, the second key is a position in
#								  the reference chromsome and the value is the position on the key chromsome
#								  that has a SNP with the given position on the reference
#					       3) rep_chrs	- A list of chromsomes which repative elements are found and filtered
#				       	       4) min_dist	- The minimum distance allowed between two SNPs before they are filtered
sub filter_snps_ref_w_annot
{
	my ($ref_chr, $r_snps, $r_rep_chrs, $min_dist, $dbh) = @_;

	my %snps = %$r_snps;
	my @rep_chrs = @$r_rep_chrs;
	my @my_chrs = ($ref_chr, keys %snps);


	my %filt_snps;


	#GET LISTS OF SNPS THAT SHOULD BE FILTERED OUT
	my %snp_filter_function;
	my %snp_filter_contig;

	my @snp_filter_cluster;
	my @snp_filter_context;
	my %snp_filter_cluster;
	my %snp_filter_context;

	my @snp_filter_large_indels;
	my %snp_filter_large_indels;

	my %snp_filter_repeat;
	my %snp_filter_tandem_repeat;
	
	my %snp_filter_depth;
	my %snp_filter_quality;


	#DETERMINE POSITIONS IN EACH CHROMSOME THAT ARE AT END OF CONTIG OR HAVE FUNCTION OF REPATIVE ELEMENT
	foreach my $chr (@my_chrs)
	{

		print "1\n";
		$snp_filter_function{$chr} = &filter_snps_function($chr, $dbh);
		print "2a\n";
		$snp_filter_contig{$chr} = &filter_snps_contig($chr, $dbh);
		print "2b\n";
		#$snp_filter_quality{$chr} = &filter_snps_quality($chr, $dbh);
		print "2c\n";
		$snp_filter_depth{$chr} = &filter_snps_depth($chr, $dbh);
		

	}#end foreach

	#DETERMINE SNP POSITIONS IN REFERENCE CHROMSOME THAT ARE 1) ALIGNED TO DUBIOUS SEQUENCE IN ANY CHROMSOME OR
	#2) CLUSTER TOGETHER IN THE REFERENCE CHROMSOME, OR 3) ARE NOT PRESENT IN ANOTHER CHROMOSOME
	foreach my $chr (keys %snps)
	{

		print "3\n";
		my ($r_chr1_pos, $r_chr2_pos) = &filter_snps_cluster($ref_chr, $chr,  $min_dist, $dbh);
		@snp_filter_cluster =  keys(%$r_chr1_pos);
		$snp_filter_cluster{$chr} = &hash_lib::init_hash(\@snp_filter_cluster);

		print "4\n";
		($r_chr1_pos, $r_chr2_pos) = &filter_snps_context($ref_chr, $chr,  $dbh);
		@snp_filter_context = keys(%$r_chr1_pos);
		$snp_filter_context{$chr} = &hash_lib::init_hash(\@snp_filter_context);


	}#end foreach

	#DETERMINE POSITIONS OF REPEATS IN CHROMSOMES THAT HAVE BEEN DESIGNATED TO SEARCH IN
	foreach my $chr (@rep_chrs)
	{

		print "5\n";
		$snp_filter_repeat{$chr} = &filter_snps_repeats($chr,  $dbh);
		print "6\n";
		$snp_filter_tandem_repeat{$chr} = &filter_snps_tandem_repeats($chr,  $dbh);

	}#end foreach


	#SCREEN INPUTTED LIST OF SNPS
	my %filt_snp_pos;
	my %filt_snps_annot;

	for(my $c=1; $c < @my_chrs; $c++)
	{

		my %filt_snp_count;
		my $filt_snp_flag = 0;

		foreach my $pos(keys %{$snps{$my_chrs[$c]}})
		{


			my @filt_cat;

			#FILTER SNPS BASED ON THERE LOCAL CONTEXT
			if(exists($snp_filter_context{$my_chrs[$c]}{$pos}))
			{
			
				$filt_snp_flag = 1;
				$filt_snp_pos{'CONTEXT'}{$pos} = 1;
				$filt_snp_count{'CONTEXT'} ++;
				$filt_annot{$pos}{'CONTEXT'} = 1;

			}#end if

			#FILTER SNPS BASED ON LOCAL DENSITY
			if(exists($snp_filter_cluster{$my_chrs[$c]}{$pos}))
			{
			
				$filt_snp_flag = 1;
				$filt_snp_pos{'CLUSTER'}{$pos} = 1;
				$filt_snp_count{'CLUSTER'} ++;
				$filt_annot{$pos}{'CLUSTER'} = 1;

			}#end if

			#FILTER SNPS BASED ON THEIR BEING IN PHAGE OR REPATIVE ELEMENTS
			if(exists($snp_filter_function{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}})  || exists($snp_filter_function{$my_chrs[0]}{$pos}) )
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'FUNCTION'}{$pos} = 1;
				$filt_snp_count{'FUNCTION'} ++;
				$filt_annot{$pos}{'FUNCTION'} = 1;

			}#end foreach

			#FILTER SNPS BASED ON BEING NEAR THE BEGINNING OR END OF A CONTIG
			if(exists($snp_filter_contig{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_contig{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'CONTIG'}{$pos} = 1;
				$filt_snp_count{'CONTIG'} ++;
				$filt_annot{$pos}{'CONTIG'} = 1;

			}#end foreach
		
			#FILTER SNPS BASED ON QUALITY
			if(exists($snp_filter_quality{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_quality{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'QUAL'}{$pos} = 1;
				$filt_snp_count{'QUAL'} ++;
				$filt_annot{$pos}{'QUAL'} = 1;

			}#end foreach
		
			#FILTER SNPS BASED ON DEPTH
			if(exists($snp_filter_depth{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_depth{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'DEPTH'}{$pos} = 1;
				$filt_snp_count{'DEPTH'} ++;
				$filt_annot{$pos}{'DEPTH'} = 1;

			}#end foreach
		
			#FILTER SNPS THAT ARE IN REPETATIVE REGIONS
			if(exists($snp_filter_repeat{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_tandem_repeat{$my_chrs[$c]}{$snps{$my_chrs[$c]}{$pos}}) || exists($snp_filter_repeat{$my_chrs[0]}{$pos}) || exists($snp_filter_tandem_repeat{$my_chrs[0]}{$pos}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'REPEAT'}{$pos} = 1;
				$filt_snp_count{'REPEAT'} ++;
				$filt_annot{$pos}{'REPEAT'} = 1;

			}#end foreach


			#CHECK IF FILTERS PASSED, IF SO PUT CURRENT SNP INTO HASH
			if($filt_snp_flag == 0)
			{

				$filt_snps{$my_chrs[$c]}{$pos} = $snps{$my_chrs[$c]}{$pos};

			}else{

				$filt_snp_flag = 0;


			}#end if

		}#end foreach

		print "\n\n$my_chrs[$c] SNPS filtered out: $filt_snp_count{'CONTEXT'}(context), $filt_snp_count{'CLUSTER'}(cluster), $filt_snp_count{'CONTIG'}(contig), $filt_snp_count{'FUNCTION'}(function), $filt_snp_count{'REPEAT'}(repeat), $filt_snp_count{'DEPTH'}(depth), $filt_snp_count{'QUAL'}(quality)\n SNPs remaining: " . scalar(keys(%{$filt_snps{$my_chrs[$c]}}))  . "\n\n";

	}#end foreach


	return (\%filt_snps, \%filt_annot);

}#end filter_snps_ref_w_annot



#sub filter_snps_msa - Takes as input: 1) snps		- A hash where the keys are the positions of snps in the first of my_chrs, and 
#							  the values are pipe-delimited lists of positions in all chromsomes in the same 
#							  order as my_chrs
#				       2) my_chrs	- A list of chromsomes relative to which SNPs exist
#				       3) rep_chrs	- A list of chromsomes which repative elements are found and filtered
#				       4) msa_id	- A multiple alignment id, with the alignment being used to look for
#							  clusters of SNPs
#				       5) min_dist	- The minimum distance allowed between two SNPs before they are filtered
sub filter_snps_msa
{
	my ($r_snps, $r_my_chrs, $r_rep_chrs, $msa_id, $min_dist, $dbh) = @_;

	my %snps = %$r_snps;
	my @my_chrs = @$r_my_chrs;
	my @rep_chrs = @$r_rep_chrs;


	my %filt_snps;

	#GET INDICIES OF rep_chrs IN my_chrs
	my %rep2my_chr;

	REP:foreach my $rep_chr(@rep_chrs)
	{

		for(my $i = 0; $i< @my_chrs; $i++)
		{

			if($rep_chr eq $my_chrs[$i])
			{

				$rep2my_chr{$rep_chr} = $i;
				next REP;

			}#end if

		}#end foreach

	}#end foreach

	my @final_rep_chrs = keys %rep2my_chr;


	#GET LISTS OF SNPS THAT SHOULD BE FILTERED OUT
	my %snp_filter_function;
	my %snp_filter_contig;
	my %snp_filter_cluster;
	my %snp_filter_context;
	my %snp_filter_repeat;
	my %snp_filter_tandem_repeat;


	#DETERMINE POSITIONS IN EACH CHROMSOME THAT ARE AT END OF CONTIG OR HAVE FUNCTION OF REPATIVE ELEMENT
	foreach my $chr (@my_chrs)
	{

		print "1\n";
		$snp_filter_function{$chr} = &filter_snps_function($chr, $dbh);
		print "2a\n";
		$snp_filter_contig{$chr} = &filter_snps_contig($chr, $dbh);
		print "2b\n";
		#$snp_filter_quality{$chr} = &filter_snps_quality($chr, $dbh);
		print "2c\n";
		$snp_filter_depth{$chr} = &filter_snps_depth($chr, $dbh);

	}#end foreach

	#DETERMINE SNP POSITIONS IN REFERENCE CHROMSOME THAT ARE 1) ALIGNED TO DUBIOUS SEQUENCE IN ANY CHROMSOME OR
	#2) CLUSTER TOGETHER IN THE REFERENCE CHROMSOME
	print "3\n";
	my $r_snp_filter_cluster = &filter_snps_cluster_msa(\@my_chrs, $msa_id, $min_dist, $dbh);
	%snp_filter_cluster = %$r_snp_filter_cluster;
	
	print "4\n";
	my $r_snp_filter_context = &filter_snps_context_msa(\@my_chrs, $msa_id,  $dbh);
	%snp_filter_context = %$r_snp_filter_context;

	#DETERMINE POSITIONS OF REPEATS IN CHROMSOMES THAT HAVE BEEN DESIGNATED TO SEARCH IN
	foreach my $chr (@rep_chrs)
	{

		print "5\n";
		$snp_filter_repeat{$chr} = &filter_snps_repeats($chr,  $dbh);
		print "6\n";
		$snp_filter_tandem_repeat{$chr} = &filter_snps_tandem_repeats($chr,  $dbh);

	}#end foreach


	#SCREEN INPUTTED LIST OF SNPS
	my %filt_snp_pos;
	my %filt_snp_count;
	

	my $filt_snp_flag = 0;

	foreach my $pos(keys %snps)
	{

		#GET POSITIONS OF SNPS ON DIFFERENT CHROMSOMES
		my @my_chr_pos = split /\|/, $snps{$pos};


		#FILTER SNPS BASED ON THERE LOCAL CONTEXT
		if(exists($snp_filter_context{$pos}))
		{

			$filt_snp_flag = 1;
			$filt_snp_pos{'CONTEXT'}{$pos} = 1;
			$filt_snp_count{'CONTEXT'} ++;

		}#end if

		#FILTER SNPS BASED ON DENSITY
		if(exists($snp_filter_cluster{$pos}) )
		{

			$filt_snp_flag = 1;
			$filt_snp_pos{'CLUSTER'}{$pos} = 1;
			$filt_snp_count{'CLUSTER'} ++;

		}#end if

		#FILTER SNPS BASED ON THEIR BEING IN PHAGE OR REPATIVE ELEMENTS, BEING AT THE END OF A CONTIG, QUALITY OR DEPTH
		for(my $i=0; $i < @my_chrs; $i++)
		{

			if(exists($snp_filter_quality{$my_chrs[$i]}{$my_chr_pos[$i]}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'QUALITY'}{$pos} = 1;
				$filt_snp_count{'QUALITY'} ++;

			}#end if

			if(exists($snp_filter_depth{$my_chrs[$i]}{$my_chr_pos[$i]}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'DEPTH'}{$pos} = 1;
				$filt_snp_count{'DEPTH'} ++;

			}#end if

			if(exists($snp_filter_contig{$my_chrs[$i]}{$my_chr_pos[$i]}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'CONTIG'}{$pos} = 1;
				$filt_snp_count{'CONTIG'} ++;

			}#end if

			if(exists($snp_filter_function{$my_chrs[$i]}{$my_chr_pos[$i]}) )
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'FUNCTION'}{$pos} = 1;
				$filt_snp_count{'FUNCTION'} ++;

			}#end if

		}#end foreach
		
		#FILTER SNPS THAT ARE IN REPETATIVE REGIONS
		foreach my $chr(@final_rep_chrs)
		{

			if(exists($snp_filter_repeat{$chr}{$my_chr_pos[$rep2my_chr{$chr}]}) || exists($snp_filter_tandem_repeat{$chr}{$my_chr_pos[$rep2my_chr{$chr}]}))
			{

				$filt_snp_flag = 1;
				$filt_snp_pos{'REPEAT'}{$pos} = 1;
				$filt_snp_count{'REPEAT'} ++;

			}#end foreach

		}#end foreach

		#CHECK IF FILTERS PASSED, IF SO PUT CURRENT SNP INTO HASH
		if($filt_snp_flag == 0)
		{

			$filt_snps{$pos} = $snps{$pos};

		}else{

			$filt_snp_flag = 0;

		}#end if

	}#end foreach

	#PRINT NUMBER OF SNPS FILTERED BY DIFFERENT CRITERIA
	print "\n\nSNPS filtered out: $filt_snp_count{'CONTEXT'}(context), $filt_snp_count{'CLUSTER'}(cluster), $filt_snp_count{'CONTIG'}(contig), $filt_snp_count{'FUNCTION'}(function), $filt_snp_count{'REPEAT'}(repeat), $filt_snp_count{'QUALITY'}(QUALITY), $filt_snp_count{'DEPTH'}(depth)\n SNPs remaining: " . scalar(keys(%filt_snps))  . "\n\n";

#######################################
#	#PRINT OUT POSITIONS OF FILTERED SNPS
#	print "MSA FILTERED SNP POSITIONS\n";
#	foreach my $class (sort keys %filt_snp_pos)
#	{
#
#		my @pos = keys %{$filt_snp_pos{$class}};
#		print "$class\n" . @pos . "\n";
#
#	}#end foreach
#####################################

	return \%filt_snps;

}#end filter_snps_msa


#sub filter_snps_large_indels - Takes as input: 1) chr1	- A chromosome from the database
#						2) chr2	- A second chromsome from the database
#			      - Returns a hash with positions in each crhomsome that do not align with the other
sub filter_snps_large_indels
{
	($chr1, $chr2, $dbh) = @_;

	my $r_chr1_positions_hash;
	my $r_chr2_positions_hash;

	my @chr1_positions;
	my @chr2_positions;


	#SELECT REPEAT SEQUENCES FROM THE DATABASE
	my $select_repeat = "SELECT c1.name, p.start1, p.end1, c2.name, p.start2, p.end2
			     FROM pairwise_alignment p, chromosome c1, chromosome c2
			     WHERE p.chr_id1 = c1.chr_id AND
				   p.chr_id2 = c2.chr_id AND
				   (p.start1 = 0 OR p.start2 = 0) AND
				   ((c1.name = ? AND c2.name = ?) OR (c1.name = ? AND c2.name = ?))";

        my $sth = $dbh -> prepare($select_repeat)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2, $chr1, $chr2)
                         or die "Can't execute statment: $DBI::errstr";


	#MAKE LIST CONTAINING ALL POSITIONS THAT ARE PART OF A REPEAT
        while(my ($name1, $start1, $end1, $name2, $start2, $end2) =  $sth -> fetchrow_array)
	{

		if( ($name1 eq $chr1) && ($start2 == 0) )
		{

			push @chr1_positions, ($start1..$end1);

		}elsif( ($name1 eq $chr1) && ($start1 == 0) ){

			push @chr2_positions, ($start2..$end2);

		}elsif( ($name2 eq $chr1) && ($start2 == 0) ){

			push @chr2_positions, ($start1..$end1);

		}elsif( ($name2 eq $chr1) && ($start1 == 0) ){

			push @chr1_positions, ($start2..$end2);

		}#end if


	}#end while


	#PUT POSITIONS INTO A HASH
	$r_chr1_positions_hash = &hash_lib::init_hash(\@chr1_positions);
	$r_chr2_positions_hash = &hash_lib::init_hash(\@chr2_positions);


	return ($r_chr1_positions_hash, $r_chr2_positions_hash);

}#end filter_snps_large_indels


#sub filter_snps_function - Takes as input: chr	- THe name of a chromsome from the chrosmome table
#		 	   - Returns a reference to a hash where keys are positions in the genome that
#		   	     are to be filtered out of SNP analysis based on functional characteristics
#			     of genes surrounding the SNP. 1) being in a phage region OR
#		             2) being in insertion sequence
#		   
sub filter_snps_function
{
	my ($chr, $dbh) = @_;

	my $r_chr_positions_hash;
	my @chr_positions;


	#SELECT GENOMIC COORDINATES OF PHAGE AND TRANSPOSABLE ELEMENT GENES
	my $select_coords = "SELECT DISTINCT g.start, g.end
                             FROM gene g, annotation a, chromosome c
                             WHERE c.name = ? AND
				   g.chr_id = c.chr_id AND
				   g.gid = a.gid AND 
				   (a.annotation LIKE '%integrase%' OR
				   a.annotation LIKE '%integron%' OR
				   a.annotation LIKE '%transposase%' OR
				   a.annotation LIKE '%insertion sequence%' OR
				   a.annotation LIKE '%phage%' OR
				   a.annotation LIKE '%IntI%') AND
				   (source = 'COG' OR source = 'NR' OR source = 'NCBI' OR source = 'PC' OR SOURCE = 'KEGG')";


        my $sth = $dbh -> prepare($select_coords)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		if($start < $end)
		{

			push @chr_positions, $start..$end;

		}else{

			push @chr_positions, $end..$start;

		}#end if


	}#end while


	#SELECT GENOMIC REGIONS THAT HAVE BEEN DESIGNATED AS PHAGE
	my $select_phage = "SELECT m.start, m.end
			    FROM mobile_element m, chromosome c
			    WHERE c.name = ? AND
				  c.chr_id = m.chr_id";

        $sth = $dbh -> prepare($select_phage)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		if($start < $end)
		{

			push @chr_positions, $start..$end;

		}else{

			push @chr_positions, $end..$start;

		}#end if


	}#end while


	#CREATE HAS WITH POSITIONS TO BE FILTERED
	$r_chr_positions_hash = &hash_lib::init_hash(\@chr_positions);


	return $r_chr_positions_hash;

}#end filter_snps_function


#sub filter_snps_quality - Takes as input: chr - The name of a chromosome form the chromosome table
#			 - Returns a reference to a hash where the keys are positions in the genome that
#			   are to be filtered out of the SNP analysis based on the SNP being a low quality base
sub filter_snps_quality
{
	my ($chr, $dbh) = @_;

	#GET QUALITY SCORES
	my $r_quality_scores = &get_quality_scores($chr, $dbh);
	my @quality_scores = @$r_quality_scores;


	#CREATE HASH OF BASES TO FILTER BASED ON QUALITY SCORES
	my %low_quality_bases;

	for(my $i = 0; $i < @quality_scores; $i++)
	{

		if($quality_scores[$i] < 50 && $quality_scores[$i] > 0)
		{

			$low_quality_bases{$i} = 1;

		}#end if

	}#end for


	return \%low_quality_bases;

}#end filter_snps_quality


#sub filter_snps_depth	 - Takes as input: chr - The name of a chromosome form the chromosome table
#			 - Returns a reference to a hash where the keys are positions in the genome that
#			   are to be filtered out of the SNP analysis based on the SNP being low depth
sub filter_snps_depth
{
	my ($chr, $dbh) = @_;

	#GET QUALITY SCORES
	my $r_depth = &get_depth($chr, $dbh);
	my @depth = @$r_depth;


	#CREATE HASH OF BASES TO FILTER BASED ON QUALITY SCORES
	my %low_depth_bases;

	for(my $i = 0; $i < @depth; $i++)
	{

		if($depth[$i] < 3 && $depth[$i] > 0)
		{

			$low_depth_bases{$i} = 1;

		}#end if

	}#end for


	return \%low_depth_bases;

}#end filter_snps_depth


#sub filter_snps_contig - Takes as input: chr	- THe name of a chromsome from the chrosmome table
#		 	- Returns a reference to a hash where keys are positions in the genome that
#		   	  are to be filtered out of SNP analysis based on the SNP being internal to 
#			  a contig
#		   
sub filter_snps_contig
{
	my ($chr, $dbh) = @_;

	my $r_chr_positions_hash;
	my @chr_positions;
	my @contig_coords;


	#GET THE COORDINATES OF ALL CONTIGS ASSOCIATED WITH THE GENOME	
	my $select_contig = "SELECT co.start, co.end
			     FROM contig co, chromosome c
			     WHERE c.name = ? AND
				   co.chr_id = c.chr_id 
			     ORDER BY co.start ASC";
	
        my $sth = $dbh -> prepare($select_contig)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($start, $end) =  $sth -> fetchrow_array)
	{

		if($start < $end)
		{

			push @contig_coords, ($start, $end);

		}else{

			push @contig_coords, ($end, $start);

		}#end if
	
	}#end while


	#GO THROUGH CONTIG COORDINATES AND NOT REGIONS 20 BP BEFORE END OF CONTIG TO 20 BP INTO START OF NEXT CONTIG
	for(my $i = 1; $i < (@contig_coords-1); $i+=2)
	{

		push @chr_positions, ($contig_coords[$i] - 19)..($contig_coords[$i+1] + 19);

	}#end for

	#CREATE HAS WITH POSITIONS TO BE FILTERED
	$r_chr_positions_hash = &hash_lib::init_hash(\@chr_positions);


	return $r_chr_positions_hash;

}#end filter_snps_contig


#sub filter_snps_cluster - Takes as input: chr1		- The name of a chromsome from the chrosmome table
#					   chr1         - The name of a second chromsome from the chrosmome table
#					   min_dist	- The minimum distance to another SNP by which a
#					  		  SNP is called dubious
#		 	 - Returns a reference to two hashes where keys are positions in each genomes that
#		   	   are to be filtered out of SNP analysis based on the local density of SNPs.
#			   A high density of SNPs is taken to indicate an alignment or assembly error.
#			   
sub filter_snps_cluster
{
	my ($chr1, $chr2, $min_dist,  $dbh) = @_;

	my %chr1_snps;

	my %chr1_positions_hash;
	my %chr2_positions_hash;


	#GET POSITIONS OF ALL SNPS BETWEEN TWO CHROMSOMES
	my $select_snps = "SELECT c1.name, s.position1, c2.name, s.position2
			   FROM snp s, chromosome c1, chromosome c2
			   WHERE s.chr_id1 = c1.chr_id AND 
				 s.chr_id2 = c2.chr_id AND
				 ((c1.name = ? AND c2.name = ?) OR (c1.name = ? AND c2.name = ?))";	



	#GO THROUGH LIST OF SNPS AND ADD TO FILTERED LIST IF THEY ARE LESS THAN MINIMUM DISTANCE APART
        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2, $chr2, $chr1)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name1, $pos1, $name2, $pos2) =  $sth -> fetchrow_array)
	{

		if($name1 eq $chr1)
		{

			$chr1_snps{$pos1} = $pos2;

		}else{

			$chr1_snps{$pos2} = $pos1;

		}#end if

	}#end while


	#GO THROUGH SNPS AND NOTE ANY THAT ARE WITHIN MINIMUM DISTANCE OF EACH OTHER
	my $prev_pos1 = -1;
	my $prev_pos2 = -1;

	foreach my $pos (sort {$a <=> $b} keys %chr1_snps)
	{

		#MUST GET TO THIRD SNP BEFORE CAN BEGIN	
		if($prev_pos2 == -1)
		{

			$prev_pos2 = $prev_pos1;
			$prev_pos1 = $pos;

			next;

		}#end if


		if( (($prev_pos1- $prev_pos2) <= $min_dist) || (($pos- $prev_pos1) <= $min_dist) )
		{

			$chr1_positions_hash{$prev_pos1} = 1;
			$chr2_positions_hash{$chr1_snps{$prev_pos1}} = 1;

		}#end if

		#INCREMENT VARIABLES
		$prev_pos2 = $prev_pos1;
		$prev_pos1 = $pos;

	}#end foreach


	return (\%chr1_positions_hash, \%chr2_positions_hash);

}#end filter_snps_cluster


#sub filter_snps_cluster_msa - Takes as input: chrs		- A reference to a list of chromsomes
#					       msa_id		- An ID from the multiple sequence alignment table
#					       min_dist		- The minimum distance to another SNP by which a
#					  		  SNP is called dubious
#		 	     - Returns a reference to a hash where the values are positions in the first chromsome in
#			       chrs contianing SNPs that should be filtered out because they are too close to other SNPs,
#			       suggesting a misassembly or poor alignment
#			   
sub filter_snps_cluster_msa
{
	my ($r_my_chrs, $msa_id, $min_dist,  $dbh) = @_;

	my %snp_msa_hash;
	my %chr_positions_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

	my %my_chr2inds;
	my @my_chr_inds;
	my @my_chrs = @$r_my_chrs;


	#SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $chrs = $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;

	#GET INDEX OF DESIRED CHROMSOMES IN LIST
	if(@$r_my_chrs == 0)
	{

		@my_chr_inds = 0..(@chrs-1);
		@my_chrs = @chrs;

	}else{

		for(my $c=0; $c < (@chrs); $c++)
		{

			if(exists($my_chrs_hash{$chrs[$c]}))
			{

				$my_chr2inds{$chrs[$c]} = $c;

			}#end if
	
		}#end foreach

		@my_chr_inds = @my_chr2inds{@my_chrs};

	}#end if


        #MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
        my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
        my %chrs_hash = %$r_chrs_hash;

        foreach my $chr (@my_chrs)
        {

                if( !exists($chrs_hash{$chr}) )
                {

                        print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
                        exit;

                }#end if

        }#end foreach


	#SELECT ALL SNPS
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	#KEEP SNPS THAT ARE NOT IDENTICAL AMONG INPUTTED CHROMSOMES AND DO NOT CONTAIN AMBIGUOS SEQUENCE
        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN REFERENCE GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i || ($snps[$my_chr_inds[0]] =~ /-/))
		{
	
			next;
				
		}#end if


		#PUT CURRENT SNP INTO HASH
		$snp_msa_hash{$positions[$my_chr_inds[0]]} = join("", @snps[@my_chr_inds]);

	}#end while



	#GO THROUGH SNPS AND NOTE ANY THAT ARE WITHIN MINIMUM DISTANCE OF EACH OTHER
	my $prev_pos1 = -1;
	my $prev_pos2 = -1;

	foreach my $pos (sort {$a <=> $b} keys %snp_msa_hash)
	{

		#MUST GET TO THIRD SNP BEFORE CAN BEGIN	
		if($prev_pos2 == -1)
		{

			$prev_pos2 = $prev_pos1;
			$prev_pos1 = $pos;

			next;

		}#end if


		#IF BASE IS WITHIN MIN DISTANCE OF PREVIOUS OR FOLLOWING BASE THEN NOTE IT
		if( (($prev_pos1- $prev_pos2) <= $min_dist) || (($pos- $prev_pos1) <= $min_dist) )
		{

			$chr_positions_hash{$prev_pos1} = 1;

		}#end if


		#INCREMENT VARIABLES
		$prev_pos2 = $prev_pos1;
		$prev_pos1 = $pos;

	}#end foreach


	return (\%chr_positions_hash);

}#end filter_snps_cluster_msa


#sub filter_snps_repeats - Takes as input: chr		- THe name of a chromsome from the chrosmome table
#		 	 - Returns a reference to a hash where keys are positions in the genome that
#		   	   are to be filtered out of SNP analysis based on being in repeat sequence  
sub filter_snps_repeats
{
	my ($chr,  $dbh) = @_;

	my @chr_positions;
	my $r_chr_positions_hash;


	#SELECT REPEAT SEQUENCES FROM THE DATABASE
	my $select_repeat = "SELECT r.start1, r.end1, r.start2, r.end2
			     FROM repeats r, chromosome c
			     WHERE r.chr_id = c.chr_id AND
				   c.name = ?";

        my $sth = $dbh -> prepare($select_repeat)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";


	#MAKE LIST CONTAINING ALL POSITIONS THAT ARE PART OF A REPEAT
        while(my ($start1, $end1, $start2, $end2) =  $sth -> fetchrow_array)
	{

		#PUT BOTH COPIES OF REPEAT IN LIST
		if($start1 < $end1)
		{
		
			push @chr_positions, ($start1..$end1);

		}else{

			push @chr_positions, ($end1..$start1);

		}#end if


		if($start2 < $end2)
		{
		
			push @chr_positions, ($start2..$end2);

		}else{

			push @chr_positions, ($end2..$start2);

		}#end if

	}#end while


	#PUT POSITIONS INTO A HASH
	$r_chr_positions_hash = &hash_lib::init_hash(\@chr_positions);


	return $r_chr_positions_hash;

}#end filter_snps_repeats


#sub filter_snps_tandem_repeats - Takes as input: chr		- THe name of a chromsome from the chrosmome table
#		 	 	- Returns a reference to a hash where keys are positions in the genome that
#		   		   are to be filtered out of SNP analysis based on being in tandem-repeat sequence.
#			   
sub filter_snps_tandem_repeats
{
	my ($chr,  $dbh) = @_;

	my @chr_positions;
	my $r_chr_positions_hash;


	#SELECT REPEAT SEQUENCES FROM THE DATABASE
	my $select_repeat = "SELECT r.start, r.unit_length, r.copy_number
			     FROM repeat_tandem r, chromosome c
			     WHERE r.chr_id = c.chr_id AND
				   c.name = ?";

        my $sth = $dbh -> prepare($select_repeat)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";


	#MAKE LIST CONTAINING ALL POSITIONS THAT ARE PART OF A REPEAT
        while(my ($start, $unit_length, $copy_num) =  $sth -> fetchrow_array)
	{

		push @chr_positions, ($start..($start+($unit_length*$copy_num)));

	}#end while


	#PUT POSITIONS INTO A HASH
	$r_chr_positions_hash = &hash_lib::init_hash(\@chr_positions);


	return $r_chr_positions_hash;

}#end filter_snps_tandem_repeats


#sub filter_snps_context - Takes as input: chr1	- THe name of a chromsome from the chrosmome table
#					   chr2	- THe name of a second chromsome from the chrosmome table
#		 	 - Returns a reference to two hashes where keys are positions in the genome that
#		   	   are to be filtered out of SNP analysis based on the sequence context around
#			   the SNPs in the two genomes
#		   
sub filter_snps_context
{
	my ($chr1, $chr2, $dbh) = @_;


	#GET QUALITY SCORES FOR CHROMSOMES
	my %quality;
	my @chrs = ($chr1, $chr2);

	foreach my $chr (@chrs)
	{

		$quality{$chr} = &get_quality_scores($chr, $dbh);		

	}#end foreach


	#GET SNPS FROM THE SNO TABLE AND SURROUNDING SEQUENCE
	my $select_seq = "SELECT c1.name, s.position1, SUBSTRING(c1.sequence, (s.position1 - 5), 11), c2.name, s.position2, SUBSTRING(c2.sequence, (s.position2 - 5), 11), s.base1, s.base2
                          FROM snp s, chromosome c1, chromosome c2
                          WHERE s.chr_id1 = c1.chr_id AND 
                                s.chr_id2 = c2.chr_id AND
                                ((c1.name = ? AND c2.name = ?) OR (c1.name = ? AND c2.name = ?))";


        my $sth = $dbh -> prepare($select_seq)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2, $chr2, $chr1)
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($name1, $pos1, $seq1, $name2, $pos2, $seq2, $base1, $base2) =  $sth -> fetchrow_array)
	{

		#GET QUALITY VECTORS	
		my @q1 = @{$quality{$name1}}[($pos1-5)..($pos1+5)];
		my @q2 = @{$quality{$name2}}[($pos2-5)..($pos2+5)];


		#IF QUALITY SCORES PROVIDED, THEN DETERMINE THE NUMBER OF BASES IN RANGE BELOW THRESHOLD
		my $q1_c = grep($_ < 50, @q1);
		my $q2_c = grep($_ < 50, @q2);
		my $q1_l = @{$quality{$name1}};
		my $q2_l = @{$quality{$name2}};

		my $q1 = &arith_ops::min($q1_c, $q1_l); 
		my $q2 = &arith_ops::min($q2_c, $q2_l); 
	
	
		#DETERMINE IF POSITIONS SHOULD BE FILTERED
		if($chr1 eq $name1)
		{

			#STORE SNP IF 1) THERE ARE MORE THAN 2 LOW QUALITY/AMBIGUOS BASES IN SURROUNDING 10 BP OR 2) THE SURROUNDING 10 BP HAS A 
			#HOMOPOLYMER RUN OF LENGTH 5 OR 3) EITHER BASE IS LOW QUALITY OR AMBIGUOUS
			#if(($seq1 =~ /^.{5}[actgN].{5}$/) || ($seq2 =~ /^.{5}[actgN].{5}$/) )
			#if( (($seq1 =~ tr/actgN/actgN/) > 2) || (($seq2 =~ tr/actgN/actgN/) > 2)  || ($seq1 =~ /AAAAA|TTTTT/i) || ($seq2 =~ /AAAAA|TTTTT/i) || ($seq1 =~ /^.{5}[N].{5}$/) || ($seq2 =~ /^.{5}[N].{5}$/) )
			#if( $q1 > 2 || $q2 > 2 || (($seq1 =~ tr/actgN/actgN/) > 2) || (($seq2 =~ tr/actgN/actgN/) > 2)  || ($seq1 =~ /AAAAA|TTTTT/i) || ($seq2 =~ /AAAAA|TTTTT/i) || ($seq1 !~ /^.{5}[ACGT].{5}$/) || ($seq2 !~ /^.{5}[ACGT].{5}$/) )
			if( $q1 > 2 || $q2 > 2 || ($seq1 =~ /AAAAA|TTTTT/i && $base2 == '.') || ($seq2 =~ /AAAAA|TTTTT/i && $base1 == '.') || ($seq1 !~ /^.{5}[ACGT].{5}$/i) || ($seq2 !~ /^.{5}[ACGT].{5}$/i) )
			{
				#print "$chr1, $chr2, $base1, $base2, $pos1, $pos2, $q1, $q2, $seq1, $seq2, @q2\n";
				$chr1_positions{$pos1} = 1;
				$chr2_positions{$pos2} = 1;

			}#end if

		}else{

			#STORE SNP IF 1) THERE ARE MORE THAN 2 LOW QUALITY/AMBIGUOS BASES IN SURROUNDING 10 BP OR 2) THE SURROUNDING 10 BP HAS A 
			#HOMOPOLYMER RUN OF LENGTH 5 OR 3) EITHER BASE IS LOW QUALITY OR AMBIGUOUS
			#if( ($seq1 =~ /^.{5}[actgN].{5}$/) || ($seq2 =~ /^.{5}[actgN].{5}$/) )
			#if( (($seq1 =~ tr/actgN/actgN/) > 2) || (($seq2 =~ tr/actgN/actgN/) > 2)  || ($seq1 =~ /AAAAA|TTTTT/i) || ($seq2 =~ /AAAAA|TTTTT/i) || ($seq1 =~ /^.{5}[N].{5}$/) || ($seq2 =~ /^.{5}[N].{5}$/) )
			#if( $q1 > 2 || $q2 > 2 || (($seq1 =~ tr/actgN/actgN/) > 2) || (($seq2 =~ tr/actgN/actgN/) > 2)  || ($seq1 =~ /AAAAA|TTTTT/i) || ($seq2 =~ /AAAAA|TTTTT/i) || ($seq1 !~ /^.{5}[ACGT].{5}$/) || ($seq2 !~ /^.{5}[ACGT].{5}$/) )
			if( $q1 > 2 || $q2 > 2 || ($seq1 =~ /AAAAA|TTTTT/i && $base2 == '.') || ($seq2 =~ /AAAAA|TTTTT/i && $base1 == '.') || ($seq1 !~ /^.{5}[ACGT].{5}$/i) || ($seq2 !~ /^.{5}[ACGT].{5}$/i) )
			{

				#print "$chr1, $chr2, $pos1, $pos2, $q1, $q2, $seq1, $seq2, @q1\n";
				$chr1_positions{$pos2} = 1;
				$chr2_positions{$pos1} = 1;

			}#end if

		}#end if

	}#end while

	return (\%chr1_positions, \%chr2_positions);

}#end filter_snps_context


#sub filter_snps_context_msa - Takes as input: chrs	- A reference to a list of chromsomes
#					       msa_id	- An id from the mutlitple_alignment table 
#		 	     - Returns a reference to a hash where the values are positions in the first chromsome in
#			       chrs contianing SNPs that should be filtered out because the surrounding sequence in any 
#			       of the chromosomes in dubious
#
sub filter_snps_context_msa
{
	my ($r_my_chrs, $msa_id, $dbh) = @_;

	my %chr_positions_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

	my %my_chr2inds;
	my @my_chr_inds;
	my @my_chrs = @$r_my_chrs;


	#GET QUALITY SCORES FOR CHROMSOMES
	my %quality;

	foreach my $chr (@my_chrs)
	{

		$quality{$chr} = &get_quality_scores($chr, $dbh);		

	}#end foreach


	#SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $chrs = $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;


	#GET INDEX OF DESIRED CHROMSOMES IN LIST
	if(@$r_my_chrs == 0)
	{

		@my_chr_inds = 0..(@chrs-1);
		@my_chrs = @chrs;

	}else{

		for(my $c=0; $c < (@chrs); $c++)
		{

			if(exists($my_chrs_hash{$chrs[$c]}))
			{

				$my_chr2inds{$chrs[$c]} = $c;

			}#end if
	
		}#end foreach

		@my_chr_inds = @my_chr2inds{@my_chrs};

	}#end if


	#GET THE SEQUENCE OF ALL CHROMOSOMES, SO THAT THE SEQUENCE SURROUNDING SNPS CAN BE CHECKED
	my $r_chr_desc = &get_chr_desc_all($dbh);
        my %chr_desc_hash = %$r_chr_desc;

	my $r_chr2id = &get_chr_id_all($dbh);
	my %chr2id = %$r_chr2id;


	#SELECT ALL SNPS
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	#KEEP SNPS THAT ARE NOT IDENTICAL AMONG INPUTTED CHROMSOMES AND DO NOT CONTAIN AMBIGUOS SEQUENCE
        SNP:while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN ANY GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/i || join("", @snps[@my_chr_inds]) =~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i)
		{
	
			next;
				
		}#end if


		#CHECK IF SEQUENCE SURROUNDING ANY SNP IS DUBIOUS, IF SO PUT CURRENT SNP INTO HASH
		foreach my $chr (@my_chrs)
		{	

			#GET THE SEQUENCE SURROUNDING THE CURRENT SNP
			my $subseq = substr($chr_desc_hash{$chr2id{$chr}}{'SEQ'}, ($positions[$my_chr2inds{$chr}] - 6), 11);
			#my $subseq = "ATGCATGCATG";

			#GET QUALITY SCORES
			my @q;
			my $q;

			print "$chr:" . @{$quality{$chr}} . "\n";

			if(@{$quality{$chr}} > 0)
			{

				@q = @{$quality{$chr}}[($positions[$my_chr2inds{$chr}]-5)..($positions[$my_chr2inds{$chr}]+5)];
				$q = grep($_ < 50, @q); 
			}else{

				$q = 0;

			}#end if		
	

			#if(($subseq =~ /^.{5}[actgN].{5}$/))
			if( $q > 2 || (($subseq =~ tr/actgN/actgN/) > 2) || ($subseq =~ /AAAAA|TTTTT/i) || ($subseq !~ /^.{5}[ACGT].{5}$/))
			#if( (($subseq =~ tr/actgN/actgN/) > 2) || ($subseq =~ /AAAAA|TTTTT/i) || ($subseq =~ /^.{5}[N].{5}$/))
			{	

				print "$chr. $q\n";
				$chr_positions{$positions[$my_chr_inds[0]]} = 1;

				next SNP;

			}#end if

		}#end foreach

	}#end while


	return (\%chr_positions);

}#end filter_snps_context_msa


#sub validate_snps - Takes as input: chr_list 	- A list of chromsomes
#				     pos_list	- A list of positions in the chromosomes that are aligned
#				     pos_bases	- A string of bases that occur at the given positions in the given chromsomes
#
#	           - Validates the SNPs occur as predicted by performing a multiple sequece alingnment of the surrounding regions
sub validate_snps
{
	my ($r_chr_list, $r_pos_list, $pos_bases, $dbh) = @_;

	my @chr_list = @$r_chr_list;
	my @pos_list = @$r_pos_list;

	$pos_bases =~ tr/[a-z]/[A-Z]/;


	#CREATE CLUSTAL ALIGNMENT OBJECT
	my $factory = new Bio::Tools::Run::Alignment::Clustalw(-quiet => 1);


	#GET THE SEQUENCES SURROUNDING EACH POSITION IN ITS GIVEN CHROMSOME, AND CREATE SEQUENCE OBJECT
	my $r_chr2id = &get_chr_id_all($dbh);
	my %chr2id = %$r_chr2id;

	my @seq_objs;

	#ADD REFERENCE TO OBJECT LIST
	my $seq = &extract_sequence($chr2id{$chr_list[0]}, (abs($pos_list[0]) - 25), (abs($pos_list[0])+25), $dbh);	
	print "$chr_list[0]($pos_list[0]): $seq\n";

	my $seq_obj = Bio::Seq->new(-seq => &extract_sequence($chr2id{$chr_list[0]}, (abs($pos_list[0]) - 25), (abs($pos_list[0])+25), $dbh),
                        	    -id  => $chr_list[0]);

	if( $seq_obj->alphabet ne 'dna' ) 
	{
		
		print "BAD (NOT DNA): $seq\n";

		return 0;

	}#end if

	push @seq_objs, $seq_obj;


	#ADD REST OF SEQUENCES, CHECKING FOR CORRECT STRAND
	for(my $c = 1; $c <@chr_list; $c ++)
	{
		#GET SEQUENCE AND CREATE OBJECT
		my $seq = &extract_sequence($chr2id{$chr_list[$c]}, (abs($pos_list[$c]) - 25), (abs($pos_list[$c])+25), $dbh);	

		my $seq_obj = Bio::Seq->new(-seq => &extract_sequence($chr2id{$chr_list[$c]}, (abs($pos_list[$c]) - 25), (abs($pos_list[$c])+25), $dbh),
                              		    -id  => $chr_list[$c]);


		#VERIFY THAT SEQUENCE ONLY COMPOSED OF NUCLEOTIDES
		my $seq_obj_rc;

		if( $seq_obj->alphabet eq 'dna' ) 
		{
		
			$seq_obj_rc =  $seq_obj->revcom;

		}else{

			print "BAD (NOT DNA): $seq\n";

			return 0;

		}#end if


		#ALIGN TO REFERENCE AND DETERMINE WHETHER TO REVERSE COMPLEMENT FOR MSA
        	my $aln_for = $factory->align([$seq_obj, $seq_objs[0]]);
        	my $aln_rev = $factory->align([$seq_obj_rc, $seq_objs[0]]);
		
		if(($aln_for->percentage_identity) > ($aln_rev->percentage_identity))
        	{
			push @seq_objs, $seq_obj;
			print "$chr_list[$c]:($pos_list[$c]) " .  $seq_obj->seq . "\n";

        	}else{

			push @seq_objs, $seq_obj_rc;
			print "$chr_list[$c]:($pos_list[$c]) " .  $seq_obj_rc->seq . "\n";

        	}#end if

	}#end for
	

	#ALIGN THE SEQUENCES
        my $aln = $factory->align(\@seq_objs);


	#GET THE BASE AT THE MIDDLE POSITION IN EACH SEQUENCE
	my $snp_pos = $aln->column_from_residue_number( $chr_list[0], 26);
	my $snp_slice = $aln->slice($snp_pos, $snp_pos);

	my $cons_aln = $aln->consensus_string(100);

	my @snps;
	my @seqs;

	for(my $c = 0;$c < @chr_list; $c ++)
	{

		#GET THE BASE AT THE APPROPRIATE PLACE IN THE ALIGNMENT
		my $snp = $snp_slice->get_seq_by_id($chr_list[$c]);


		#MAKE SURE NOT INDEL AT POSITION
		if($snp)
		{

			push @snps, $snp->seq;

		}else{

			push @snps, "-";

		}#end if


		#ALO GET FULL SEQUENCE TO REPORT
		push @seqs, ($aln->get_seq_by_id($chr_list[$c]))->seq;

	}#end foreach
	
	
	#VERIFY THAT BASES AT THE GIVEN POSITION MATCH EXPECTATION
	my $snps = join("", @snps);
	print "SNPS = $snps\n";
	my $snps_rc = reverse &seq_ops::rev_com($snps);

	if( (($snps eq $pos_bases) || ($snps_rc eq $pos_bases)) && (($cons_aln =~ tr/?//) < 3))
	{

		#ALIGNMENT AS EXPECTED
		print "GOOD: $snps, $pos_bases, $cons_aln\n\n";
		return 1;

	}else{

		#ALIGNMENT DIFFERENT FROM EXPECTATION
		print "BAD: $snps, $pos_bases, $cons_aln\n\n";
		return 0;

	}#end if


}#end validate_snps


#sub validate_indels - Takes as input: chr_list 	- A list of chromsomes, where the first chromsome has the sequence
#		 		       pos_list		- A list of positions in the chromosomes that are aligned
#				       pos_indel	- A binary list indicating whether the indel sequence should be present
#				       indel_seq	- The sequence of the indel
#
#	           - Validates that indels occur as predicted by performing a multiple sequece alingnment of the surrounding regions
sub validate_indels
{
	my ($r_chr_list, $r_pos_list, $pos_indel, $indel_seq, $dbh) = @_;

	my @chr_list = @$r_chr_list;
	my @pos_list = @$r_pos_list;

	$indel_seq =~ tr/[a-z]/[A-Z]/;
	
	my $indel_length = length($indel_seq);


	#GET THE SEQUENCES SURROUNDING EACH POSITION IN ITS GIVEN CHROMSOME, AND CREATE SEQUENCE OBJECT
	my $r_chr2id = &get_chr_id_all($dbh);
	my %chr2id = %$r_chr2id;

	my @seq_objs;

	for(my $c = 0; $c <@chr_list; $c ++)
	{
		my $seq = &extract_sequence($chr2id{$chr_list[$c]}, (abs($pos_list[$c]) - 25), (abs($pos_list[$c])+25+$indel_length), $dbh);	
		print "$chr_list[$c]: $seq\n";

		my $seq_obj = Bio::Seq->new(-seq => &extract_sequence($chr2id{$chr_list[$c]}, (abs($pos_list[$c]) - 25), (abs($pos_list[$c])+25+$indel_length), $dbh),
                              		    -id  => $chr_list[$c]);

		push @seq_objs, $seq_obj;

	}#end for
	

	#ALIGN THE SEQUENCES
	my $factory = new Bio::Tools::Run::Alignment::Clustalw(-quiet => 1);

        my $aln = $factory->align(\@seq_objs);


	#GET THE CONSENSOUS ALIGNMENT
	my $cons_aln = $aln->consensus_string(100);


	#VERIFY THAT CONSENSOUS FITS EXPECTATION
	if( $cons_aln =~ /^[^?]*\??[^?]+\?{$indel_length}[^?]*\??[^?]+\?{$indel_length}$/ )
	{

		#ALIGNMENT AS EXPECTED
		print "GOOD: $cons_aln\n\n";
		return 1;

	}else{

		#ALIGNMENT DIFFERENT FROM EXPECTATION
		print "BAD: $cons_aln\n\n";
		return 0;

	}#end if


}#end validate_indels

#sub get_indels - Takes as input: 1) chr1	- A name from the chromsome table
#				  2) chr2	- A second name from the chromosome table
#	        - Returns a hash whose first key is the position in chr1 of the indel and second keys are:
# 		 1) position in chromosome 2, 2) chromsome sequence is in and 3) the sequence
sub get_indels
{
	my ($chr1, $chr2, $dbh) = @_;

	my %indel_info;

	#GET SNPs FROM THE snp TABLE
	my $select_snps = "SELECT base1, position1, base2, position2, SUBSTRING(c1.sequence, position1 - 5, 11), SUBSTRING(c2.sequence, position2 - 5, 11), SUBSTRING(c1.sequence, position1 - 20, 41), SUBSTRING(c2.sequence, position2 - 20, 41)
                             FROM snp s, chromosome c1, chromosome c2, contig c
                             WHERE c1.name = ? AND 
				   c2.name = ? AND
				   s.chr_id1 = c1.chr_id AND
				   s.chr_id2 = c2.chr_id AND
			  	   c.chr_id = c1.chr_id AND
       				   ((c.start < position1 AND c.end> position1) OR (c.start > position1 AND c.end < position1)) AND
      				   (c.contig_name LIKE 'contig%' OR c.contig_name LIKE '%chr%' OR c.contig_name LIKE '%_plasmid%%') AND
				   (base1 = '.' OR base2 = '.')
			    ORDER BY position1, position2 ASC";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr1, $chr2)
                         or die "Can't execute statment: $DBI::errstr";


	#MAKE SURE THAT CHROMSOMES ARE NOT IN REVERSED ORDER IN DATABASE, IF SO RE-EXECUTE ACCORDINGLY
	if ($sth->rows == 0) 
	{

		$select_snps = "SELECT base2, position2, base1, position1, SUBSTRING(c2.sequence, position2 - 5, 11), SUBSTRING(c1.sequence, position1 - 5, 11), SUBSTRING(c2.sequence, position2 - 20, 41), SUBSTRING(c1.sequence, position1 - 20, 41)
                	             FROM snp s, chromosome c1, chromosome c2, contig c
                	             WHERE c1.name = ? AND 
					   c2.name = ? AND
					   s.chr_id1 = c1.chr_id AND
					   s.chr_id2 = c2.chr_id AND
				  	   c.chr_id = c1.chr_id AND
       					   ((c.start < position1 AND c.end> position1) OR (c.start > position1 AND c.end < position1)) AND
      					   (c.contig_name LIKE 'contig%' OR c.contig_name LIKE '%chr%' OR c.contig_name LIKE '%_plasmid%%') AND
					   (base1 = '.' OR base2 = '.')
				    ORDER BY position2, position1 ASC";


        	$sth = $dbh -> prepare($select_snps)
        	                 or die "Can't prepare statment: $DBI::errstr";

        	$rc = $sth -> execute($chr2, $chr1)
        	                 or die "Can't execute statment: $DBI::errstr";


	}#end if


	my $indel_pos = 0;
	my $prev_pos = 0;
	my $indel_seq = "";
	my $in_seq = 0;
	my $indel = 0;

        while(my ($base1, $position1, $base2, $position2, $sub_seq1, $sub_seq2, $ext_sub_seq1, $ext_sub_seq2) =  $sth -> fetchrow_array)
        {

		#IGNORE SNP IF 1) AFTER HOMOPOLYMER, 2) LOW QUALITY BASE OR 3) AMBIGOUS NUCLEOTIDE
		#if( $sub_seq1 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i || $sub_seq2 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i)
		if( $sub_seq1 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i || 
		    $sub_seq2 =~ /AAAAA|TTTTT|CCCCC|GGGGG|NNNNN/i || 
		    $sub_seq1 =~ /^.{3}.*[actgN].*.{3}$/ || 
		    $sub_seq2 =~ /^.{3}.*[actgN].*.{3}$/)
		{

			#IF IN THE MIDDLE OF INDEL, THEN END IT
			if($indel > 0)
			{

				$indel_info{$indel_key}{'INDEL_LENGTH'} = $indel;
				$indel_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;
				$indel = 0;

				if((($in_seq == 1 && $position1 == ($prev_pos + 1) && $base2 eq '.' ) || ($in_seq == 2 && $position2 == ($prev_pos + 1) && $base1 eq '.') ))
				{

					print "END INDEL BECAUSE OF BAD SEQUENCE ($position1)!\n";

				}#end if

			}#end if

			next;

		}#end if


		#CHECK IF SNP PART OF CURRENT INDEL
		if( $indel > 0 && (($in_seq == 1 && $position1 == ($prev_pos + 1) && $base2 eq '.' ) || ($in_seq == 2 && $position2 == ($prev_pos + 1) && $base1 eq '.') ) )
		{

			#INCREMENT SIZE OF INDEL, NOTE THE CURRENT POSITION IN THE INSERTED SEQUENCE AND MOVE TO NEXT SNP
			$indel ++;

			if($in_seq == 1)
			{

				$prev_pos = $position1;
				$indel_seq = $indel_seq . $base1;

			}else{

				$prev_pos = $position2;
				$indel_seq = $indel_seq . $base2;

			}#end if

			next;
	
		}elsif( $indel > 0 ){

			#THIS IS A NEW INDEL, SO SAVE THE PREVIOUS ONE
			$indel_info{$indel_key}{'INDEL_LENGTH'} = $indel;
			$indel_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;


			#SET UP VARIABLES FOR NEW INDEL
			$indel = 1;

			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if

		}elsif( $indel == 0 ){

			#NEW INDEL, SET UP VARIABLES DESCRIBING INDEL
			$indel = 1;


			#DETERMINE WHICH SEQUENCE HAS EXTRA BASES
			if($base1 eq '.')
			{

				$in_seq = 2;
				$indel_pos = $position2;
				$prev_pos = $position2;
				$indel_seq = $base2;
				$indel_key = $position1;

			}else{

				$in_seq = 1;
				$indel_pos = $position1;
				$prev_pos = $position1;
				$indel_seq = $base1;
				$indel_key = $position1;

			}#end if
	
		}#end if


		#RECORD INFO ON SNP
		$indel_info{$position1}{'BASE1'} = $base1;
		$indel_info{$position1}{'BASE2'} = $base2;
		$indel_info{$position1}{'POS2'} = $position2;
		$indel_info{$position1}{'SEQ1'} = $sub_seq1;
		$indel_info{$position1}{'SEQ2'} = $sub_seq2;
		$indel_info{$position1}{'EXT_SEQ1'} = $ext_sub_seq1;
		$indel_info{$position1}{'EXT_SEQ2'} = $ext_sub_seq2;


	}#end while


	#STORE THE LAST INDEL IF ONE EXISTS
	if( $indel > 0 )
	{

		$indel_info{$indel_key}{'INDEL_LENGTH'} = $indel;
		$indel_info{$indel_key}{'INDEL_SEQ'} = $indel_seq;

	}#end if

	return \%indel_info;

}#end get_indels


#sub get_all_chr_pair_snps - Takes as input:
#			   - Returns a 2D hash where the keys are chromsome names and the values are references to lists of positions
#			     where SNPs occur between the two chromsomes	 
sub get_all_chr_pair_snps
{
	my ($dbh) = @_;

	my %snp_hash;


	#GET ALL CHRS
	my $r_chr_hash = &get_chr_desc_all($dbh);
	my %chr_hash = %$r_chr_hash;
	
	my @chr_ids = keys %chr_hash;


	#FOREACH PAIR OF CHROMOSOMES, GET PAIRS OF SNPs
	for(my $c1 = 0; $c1 < @chr_ids; $c1++)
	{


		for(my $c2 = $c1+1; $c2 < @chr_ids; $c2++)
		{

			#GET CURRENT CHROMOSOMES
			print "$c1, $c2\n";
			my $chr1 = $chr_ids[$c1];
			my $chr2 = $chr_ids[$c2];


			#SELECT SNPs FROM THE DATABASE
			my $select_snps = "SELECT s.position1, s.chr_id1, s.position2, s.chr_id2
					   FROM snp s
					   WHERE ((s.chr_id1 = ? AND s.chr_id2 = ?) OR (s.chr_id1 = ? AND s.chr_id2 = ?)) AND
      						 s.base1 <> '.' AND
      						 s.base2 <> '.'";

			my $sth = $dbh -> prepare($select_snps)
       	                	  or die "Can't prepare statment: $DBI::errstr";

			my $rc = $sth -> execute($chr1, $chr2, $chr2, $chr1)
        		          or die "Can't execute statment: $DBI::errstr";


			#PUT SNPS IN HASH
        		while(my ($snp_pos1, $snp_chr1, $snp_pos2, $snp_chr2) =  $sth -> fetchrow_array)
        		{

				if(exists($snp_hash{$chr_hash{$snp_chr1}{'NAME'}}{$chr_hash{$snp_chr2}{'NAME'}}))
				{

					push @{$snp_hash{$chr_hash{$snp_chr1}{'NAME'}}{$chr_hash{$snp_chr2}{'NAME'}}}, $snp_pos1;
					push @{$snp_hash{$chr_hash{$snp_chr2}{'NAME'}}{$chr_hash{$snp_chr1}{'NAME'}}}, $snp_pos2;

				}else{

					@{$snp_hash{$chr_hash{$snp_chr1}{'NAME'}}{$chr_hash{$snp_chr2}{'NAME'}}} = ($snp_pos1);
					@{$snp_hash{$chr_hash{$snp_chr2}{'NAME'}}{$chr_hash{$snp_chr1}{'NAME'}}} = ($snp_pos2);

				}#end if

			}#end while

		}#end for

	}#end for


	return \%snp_hash;

}#end get_all_chr_pair_snps


#sub get_all_chr_pair_indels - Takes as input: min_length - THe minimum size of an indel to consider
#			     - Returns a 2D hash where the keys are chromsome names and the values are references to lists of positions
#			       where indels occur between the two chromsomes (e.g. start1, end1, start2, end2, ... startN, endN)	 
sub get_all_chr_pair_indels
{
	my ($min_length, $dbh) = @_;

	my $genome_db = "blast_db/Abau_pe_libs_baseCall_v2.3/AbauALL/Abau_ALL";

	my %indel_hash;
	my %indel_region_hash;
	my %indel_seqs;


	#GET ALL CHRS
	my $r_chr_hash = &get_chr_desc_all($dbh);
	my %chr_hash = %$r_chr_hash;
	
	my @chr_ids = keys %chr_hash;


	#FOREACH PAIR OF CHROMOSOMES, GET PAIRS OF INDELS
	for(my $c1 = 0; $c1 < @chr_ids; $c1++)
	{


		for(my $c2 = $c1+1; $c2 < @chr_ids; $c2++)
		{

			#GET CURRENT CHROMOSOMES
			print "$c1, $c2\n";
			my $chr1 = $chr_ids[$c1];
			my $chr2 = $chr_ids[$c2];

		
			#INITITALIZE LISTS OF COORDINATES
			$indel_hash{$chr_hash{$chr1}{'NAME'}}{$chr_hash{$chr2}{'NAME'}} = ();
			$indel_hash{$chr_hash{$chr2}{'NAME'}}{$chr_hash{$chr1}{'NAME'}} = ();


			#SELECT INDELS FROM THE DATABASE
			my $select_indels = "SELECT pa.chr_id1, pa.start1, pa.end1, pa.chr_id2, pa.start2, pa.end2
					     FROM pairwise_alignment pa
					     WHERE ((pa.chr_id1 = ? AND pa.chr_id2 = ?) OR (pa.chr_id1 = ? AND pa.chr_id2 = ?)) AND
      						   (pa.start1 = 0 OR pa.start2 = 0)";

			my $sth = $dbh -> prepare($select_indels)
       	                	  or die "Can't prepare statment: $DBI::errstr";

			my $rc = $sth -> execute($chr1, $chr2, $chr2, $chr1)
        		          or die "Can't execute statment: $DBI::errstr";


			#PUT INDELS IN HASH
        		while(my ($indel_chr1, $indel_start1, $indel_end1, $indel_chr2, $indel_start2, $indel_end2) =  $sth -> fetchrow_array)
        		{

				#CHECK IF INDEL IS LEGIT BY BLASTING AND VERIFYING THAT SEQUENCE IS NOT PRIMARILY AMBIGUOUS
				my $seq;
				my $gc;
				my $region_name;

	
				#GET INDEL SEQUENCE 
				if($indel_start1 == 0)
				{
					$seq = &extract_sequence($indel_chr2, $indel_start2, $indel_end2, $dbh);

					$region_name = $chr_hash{$indel_chr2}{'NAME'} . ":$indel_start2-$indel_end2";

				}else{

					$seq = &extract_sequence($indel_chr1, $indel_start1, $indel_end1, $dbh);

					$region_name = $chr_hash{$indel_chr1}{'NAME'} . ":$indel_start1-$indel_end1";

				}#end if


				#MAKE SURE SEQUENCE IS LONG ENOUGH
				if(length($seq) < $min_length)
				{

					next;

				}#end if


				#MAKE SURE SEQUENCE IS LESS THAN 50% AMBIGUOS
        			my $region_seq_obj = Bio::Seq->new(-seq => $seq);

        			my $r_region_nuc_counts = &seq_stats::init_monoNucHash();
        			$r_region_nuc_counts = &seq_stats::get_monoNucCount($r_region_nuc_counts, $region_seq_obj);
        			my %region_nuc_counts = %$r_region_nuc_counts;

        			my $region_length = length($seq);
        			my $N_count = $region_length - ($region_nuc_counts{'A'} + $region_nuc_counts{'C'} + $region_nuc_counts{'G'}  + $region_nuc_counts{'T'});

        			if( ($N_count/$region_length) > 0.50 )
				{

					print "C1 = $chr_hash{$indel_chr1}{'NAME'}, C2 = $chr_hash{$indel_chr2}{'NAME'},L = $region_length, N = $N_count\n";
					next;

				}#end if


				#IF LEGIT, ADD INDEL TO LIST
				if($indel_start1 == 0)
				{

					#EXTRA SEQUENCE IN SEQUENCE 2
					@{$indel_region_hash{$region_name}{$chr_hash{$indel_chr1}{'NAME'}}} = ($indel_start2, $indel_end2);
					$indel_seqs{$region_name} = $seq;

				}else{

					#EXTRA SEQUENCE IN SEQUENCE 1
					@{$indel_region_hash{$region_name}{$chr_hash{$indel_chr2}{'NAME'}}} = ($indel_start1, $indel_end1);
					$indel_seqs{$region_name} = $seq;

				}#end if

			}#end while

		}#end for

	}#end for


	#BLAST CHR SEQUENCES TO MAKE SURE THEY ARE NOT IN OTHER CHROMSOME
	my $r_blast_results = &blast::blastn_vs_userDB($genome_db, "temp.blast", \%indel_seqs);
	my %blast_results = %$r_blast_results;				


	#IF EXISTS IN CHROMSOME THAT IS SUPPOSED TO LACK, THEN SKIP
	foreach my $region(keys %indel_region_hash)
	{

		my $in_chr = $region;
		$in_chr =~ s/^([^:]+):.*$/$1/;

		foreach my $out_chr(keys %{$indel_region_hash{$region}})
		{

			if(exists($blast_results{$region}{$out_chr}))
			{
				#BLAST HIT EXISTS, IF NOT SIGNIFICANT ENOUGH, THEN STILL KEEP
				my ($significance, $query_length, $hit_length,  $frac_aligned_hit, $frac_aligned_query, $sbj_start, $sbj_end, $hit_start, $hit_end, $score, $desc) = @{$blast_results{$region}{$out_chr}};

				if($significance < 0.01 && $frac_aligned_query > 0.5)
				{

					print "C1 = $in_chr, C2 = $out_chr, S = $significance, FA = $frac_aligned_query\n";
					next;

				}else{

					push @{$indel_hash{$in_chr}{$out_chr}}, @{$indel_region_hash{$region_name}{$out_chr}};

				}#end if
			
			}else{

				#NO BLAST HIT, KEEP INDEL
				push @{$indel_hash{$in_chr}{$out_chr}}, @{$indel_region_hash{$region_name}{$out_chr}};

			}#end if

		}#end foreach

	}#end foreach



	return \%indel_hash;

}#end get_all_chr_pair_indels


#sub get_all_snps_relative_to_reference - Takes as input: 1) ref_chr 	- THe name of a chromosome in the database to be used as
#									  a reference for SNP detection
#							  2) chrs	- A reference to a list of chromsomes for which SNPs with 
#									  the reference chromosome are to be found
#					- Returns a 2D hash where the first key is the position in the reference genome
#					  the second key is the name of the second chromsome and the value is the base
#					  in the second chrosome. Also returns a hash linking positions in the reference
#					  genome to the base at that position
sub get_all_snps_relative_to_reference
{
	($ref_chr, $r_chrs, $dbh) = @_;

	my @chrs = @$r_chrs;

	my %snp_hash;
	my %ref_base_hash;
	my %indel_hash;
	my %position_hash;


	#SELECT ALL SUBSTITUTION SNPs INVOLVING THE REFERENCE GENOME
	my $select_snps = "SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
      				 s.chr_id2 = c2.chr_id AND
      				 c1.name = ? AND 
				 c2.name IN (";
	
        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";

	$select_snps = $select_snps . 		   
			  "UNION
			   SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
                                 s.chr_id2 = c2.chr_id AND
                                 c2.name = ? AND
				 c1.name IN (";

        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ref_chr, @chrs, $ref_chr, @chrs)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($chr1,$base1, $position1, $chr2, $base2, $position2) =  $sth -> fetchrow_array)
	{
		
		#CHECK WHICH CHROMOSOME IS THE REFERENCE
		if($chr1 eq $ref_chr)
		{
			
			#IF SNP IS INDEL THEN DON'T SAVE AND NOTE THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
			{
				
				$snp_hash{$chr2}{$position1} = $base2;
				$ref_base_hash{$position1} = $base1;
				$position_hash{$chr2}{$position1} = $position2;

			}else{

				$indel_hash{$position1} = 1;

			}#end if

		}else{

			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
                        {

				$snp_hash{$chr1}{$position2} = $base1;
				$ref_base_hash{$position2} = $base2;
				$position_hash{$chr1}{$position2} = $position1;

			}else{

                                $indel_hash{$position2} = 1;

                        }#end if

		}#end if
	
	}#end while


	#GET INDELS RELATIVE TO REFERENCE AND ADD TO INDEL POSITIONS, AND AT THE SAME TIME MAKE SURE THAT SNPS WERE FOUND FOR
	#ALL INPUTTED CHROMSOSOMES
	my @indel_pos = keys %indel_hash;

	foreach my $chr (@chrs)
	{

		#MAKE SURE THAT SNPS WERE FOUND FOR CURRENT CHROMSOME, AND IF NOT THROW ERROR
		if(!exists($snp_hash{$chr}))
		{

			print "\n\nNo pairwise SNPs found between $ref_chr and $chr!!!!\n\n";

			exit;

		}#end if
		

		#GET INDELS AND ADD TO LIST
		my ($r_chr1_pos_hash, $r_chr2_pos_hash) = &filter_snps_large_indels($ref_chr, $chr, $dbh);

		push @indel_pos, keys(%$r_chr1_pos_hash);

	}#end foreach	

	$r_indel_pos = &hash_lib::init_hash(\@indel_pos);
	%indel_hash = %$r_indel_pos;


	return (\%snp_hash, \%ref_base_hash, \%indel_hash, \%position_hash);

}#end get_all_snps_relative_to_reference


#sub get_all_snps_msa  - Takes as input: 1) msa_id	- The id of the MSA to use for snp analysis
#					 2) chrs	- A list of chrosomes in which SNPs should be found
#
#		        - Returns a 2D hash where the first key is the name of a chromsome and
#		          the second key is the position of the SNP in the first chrosmome  and the value 
#		          is the base  in the chrosome.
sub get_all_snps_msa
{
	($msa_id, $r_my_chrs, $dbh) = @_;

	my %snp_hash;
	my %snp_msa_hash;
	my %snp2pos_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

	my %my_chr2inds;
	my @my_chr_inds;
	my @my_chrs = @$r_my_chrs;
	

	#SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $chrs = $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;

	#GET INDEX OF DESIRED CHROMSOMES IN LIST
	if(@my_chrs == 0)
	{

		@my_chr_inds = 0..(@chrs-1);
		@my_chrs = @chrs;

	}else{

		for(my $c=0; $c < (@chrs); $c++)
		{

			if(exists($my_chrs_hash{$chrs[$c]}))
			{

				$my_chr2inds{$chrs[$c]} = $c;

			}#end if
	
		}#end foreach

		@my_chr_inds = @my_chr2inds{@my_chrs};

	}#end if

	
	#MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
	my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
	my %chrs_hash = %$r_chrs_hash;

	foreach my $chr (@my_chrs)
	{

		if( !exists($chrs_hash{$chr}) )
		{

			print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
			exit;

		}#end if
	
	}#end foreach
	

	#SELECT ALL SNPS THAT DO NOT HAVE INDELS AND SAVE IN HASH
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN ANY GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/i || join("", @snps[@my_chr_inds]) =~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i)
		{
	
			next;
				
		}#end if


		#PUT CURRENT SNP INTO HASH
		$snp_msa_hash{$positions[$my_chr_inds[0]]} = join "", @snps[@my_chr_inds];
		$snp2pos_hash{$positions[$my_chr_inds[0]]} = join "\|", @positions[@my_chr_inds];

		foreach my $chr_ind (@my_chr_inds)
		{

                        $snp_hash{$chrs[$chr_ind]}{$positions[$my_chr_inds[0]]} = $snps[$chr_ind];

		}#end for	
	

	}#end while


	return (\%snp_hash, \@my_chrs, \%snp_msa_hash, \%snp2pos_hash);

}#end get_all_snps_msa


#sub get_all_indels_msa  - Takes as input: 1) msa_id	- The id of the MSA to use for snp analysis
#					   2) chrs	- A list of chrosomes in which SNPs should be found
#
#		        - Returns several hashes relating indels to chromsomes and positions
sub get_all_indels_msa
{
	($msa_id, $r_my_chrs, $dbh) = @_;

	my %snp_hash;
	my %snp_msa_hash;
	my %snp2pos_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

	my %my_chr2inds;
	my @my_chr_inds;
	my @my_chrs = @$r_my_chrs;
	

	#SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
	my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	my $chrs = $sth -> fetchrow_array;
	my @chrs = split /\|/, $chrs;

	#GET INDEX OF DESIRED CHROMSOMES IN LIST
	if(@my_chrs == 0)
	{

		@my_chr_inds = 0..(@chrs-1);
		@my_chrs = @chrs;

	}else{

		for(my $c=0; $c < (@chrs); $c++)
		{

			if(exists($my_chrs_hash{$chrs[$c]}))
			{

				$my_chr2inds{$chrs[$c]} = $c;

			}#end if
	
		}#end foreach

		@my_chr_inds = @my_chr2inds{@my_chrs};

	}#end if

	
	#MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
	my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
	my %chrs_hash = %$r_chrs_hash;

	foreach my $chr (@my_chrs)
	{

		if( !exists($chrs_hash{$chr}) )
		{

			print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
			exit;

		}#end if
	
	}#end foreach
	

	#SELECT ALL SNPS THAT DO NOT HAVE INDELS AND SAVE IN HASH
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR THERE IS NO INDEL
		if(join("", @snps[@my_chr_inds]) =~ /N/ || join("", @snps[@my_chr_inds]) !~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(-+)$/i)
		{
	
			next;
				
		}#end if


		#FIND FIRST CHROMOSOME WITH SEQUENCE, AND USE IT TO STORE INDEL
		foreach my $chr_ind (@my_chr_inds)
		{
	
			if($positions[$chr_ind] != 0)
			{
				$indel_pos_hash{$chrs[$chr_ind]}{$positions[$chr_ind]} = join "\|", @positions[@my_chr_inds];
				$indel_base_hash{$chrs[$chr_ind]}{$positions[$chr_ind]} = join "\|", @snps[@my_chr_inds];

				last;

			}#end if

		}#end if

	}#end while


	#ASSEMBLE INDEL POSITIONS INTO SMALL INDELS
	my %small_indel_info;
	my %small_indel_pos;

	my $indel_pos = 0;
	my $prev_pos = 0;
	my %indel_seq;
	my $indel_chr;
	my $indel = 0;

	foreach my $chr (@my_chrs)
	{

		foreach my $pos (sort {$a<=>$b} keys %{$indel_pos_hash{$chr}})
		{

			my @positions = split /\|/, $indel_pos_hash{$chr}{$pos};
			my @bases = split /\|/, $indel_base_hash{$chr}{$pos};


			#CHECK IF 1)CONTINUING INDEL, 2) NEW INDEL OR 3) FIRST INDEL
			if( $indel > 0 &&  $pos == ($prev_pos + 1) )
			{

				#CONTINUING INDEL

				#INCREMENT SIZE OF INDEL, NOTE THE CURRENT POSITION IN THE INSERTED SEQUENCE AND MOVE TO NEXT POSITION
				$indel ++;

				$prev_pos = $pos;

				for(my $i = 0; $i < @my_chrs; $i ++)
				{

					$indel_seq{$my_chrs[$i]} = $indel_seq{$my_chrs[$i]} . $bases[$i];

				}#end foreach

	
			}elsif( $indel > 0 ){

				#NEW INDEL

				#SO SAVE THE PREVIOUS ONE
				my %saved_indel = %indel_seq;
				$small_indel_info{$indel_chr}{$indel_pos} = \%saved_indel;
				$small_indel_pos{"$indel_chr:$indel_pos"} = $indel_pos_hash{$indel_chr}{$indel_pos};
		
				#my $indel_seq = join "\n", @indel_seq{@my_chrs};
				#my $indel_seq_2 = join "\n", values(%{$small_indel_info{$indel_chr}{$indel_pos}});
				#print "$chr, $indel_pos\n$indel_seq_2\n\n";


				#SET UP VARIABLES FOR NEW INDEL
				$indel = 1;
				$indel_pos = $pos;
				$prev_pos = $pos;
				$indel_chr = $chr;

				for(my $i = 0; $i < @my_chrs; $i ++)
				{

					$indel_seq{$my_chrs[$i]} = $bases[$i];

				}#end foreach

			}else{

				#FIRST INDEL

				#INCREMENT SIZE OF INDEL, NOTE POSITION AND INITIALIZE SEQUENCE OF INDEL
				$indel = 1;
				$indel_pos = $pos;
				$prev_pos = $pos;
				$indel_chr = $chr;

				for(my $i = 0; $i < @my_chrs; $i ++)
				{

					$indel_seq{$my_chrs[$i]} = $bases[$i];

				}#end foreach

			}#end if

		}#end foreach

	}#end foreach


	#STORE THE LAST INDEL
	$small_indel_info{$indel_chr}{$indel_pos} = \%indel_seq;
	$small_indel_pos{"$indel_chr:$indel_pos"} = $indel_pos_hash{$indel_chr}{$indel_pos};


	return (\%small_indel_info, \%small_indel_pos);

}#end get_all_indels_msa


#sub get_coding_snps_relative_to_reference - Takes as input: 1) ref_chr    - THe name of a chromosome in the database to be used as
#                                                                            a reference for SNP detection
#                                                            2) chrs       - A reference to a list of chromsomes for which SNPs with 
#                                                                            the reference chromosome are to be found
#					   - Returns a 2D hash where the first key is the position in the reference genome
#					     the second key is the name of the second chromsome and the value is the base
#					     in the second chrosome. Also returns a hash linking positions in the reference
#					     genome to the base at that position
sub get_coding_snps_relative_to_reference
{
	($ref_chr, $r_chrs, $dbh) = @_;

	my @chrs = @$r_chrs;

	my %snp_hash;
	my %ref_base_hash;
	my %indel_hash;
	my %position_hash;


	#SELECT ALL SUBSTITUTION SNPs INVOLVING THE REFERENCE GENOME
	my $select_snps = "SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
      				 s.chr_id2 = c2.chr_id AND
      				 c1.name = ? AND 
				 c2.name IN (";
	
        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";

	$select_snps = $select_snps . 		   
			  "UNION
			   SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
                                 s.chr_id2 = c2.chr_id AND
                                 c2.name = ? AND
				 c1.name IN (";

        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ref_chr, @chrs, $ref_chr, @chrs)
                         or die "Can't execute statment: $DBI::errstr";


	#GET ALL THE CODING SNPS
	my $r_coding_seq_pos = &get_gene_seq_positions($chr, $dbh);
	my %coding_seq_pos = %$r_coding_seq_pos;

        while(my ($chr1,$base1, $position1, $chr2, $base2, $position2) =  $sth -> fetchrow_array)
	{


		#DETERMINE THE EFFECT OF THE SEQUENCE CHANGE AND FILTER FOR CODING SNPS
		if($chr1 eq $ref_chr)
		{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(!exists($coding_seq_pos{$position1}))
			{
	
				next;

			}#end if


			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
			{

				$snp_hash{$chr2}{$position1} = $base2;
				$ref_base_hash{$position1} = $base1;

			}else{

				$indel_hash{$position1} = 1;

			}#end if

		}else{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(!exists($coding_seq_pos{$position2}))
			{
	
				next;

			}#end if

			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
                        {

				$snp_hash{$chr1}{$position2} = $base1;
				$ref_base_hash{$position2} = $base2;

			}else{

                                $indel_hash{$position2} = 1;

                        }#end if

		}#end if
	
	}#end while


	#GET INDELS RELATIVE TO REFERENCE AND ADD TO INDEL POSITIONS, AND AT THE SAME TIME MAKE SURE THAT SNPS WERE FOUND FOR
	# ALL INPUTTED CHROMSOSOMES
	my @indel_pos = keys %indel_hash;

	foreach my $chr (@chrs)
	{

		#MAKE SURE THAT SNPS WERE FOUND FOR CURRENT CHROMSOME, AND IF NOT THROW ERROR
		if(!exists($snp_hash{$chr}))
		{

			print "\n\nNo pairwise SNPs found between $ref_chr and $chr!!!!\n\n";

			exit;

		}#end if
		

		#GET INDELS AND ADD TO LIST
		my ($r_chr1_pos_hash, $r_chr2_pos_hash) = &filter_snps_large_indels($ref_chr, $chr, $dbh);

		push @indel_pos, keys(%$r_chr1_pos_hash);

	}#end foreach	

	$r_indel_pos = &hash_lib::init_hash(\@indel_pos);
	%indel_hash = %$r_indel_pos;


	return (\%snp_hash, \%ref_base_hash, \%indel_hash);

}#end get_coding_snps_relative_to_reference


#sub get_coding_snps_msa  - Takes as input: 1) msa_id	- The id of the MSA to use for snp analysis
#					    2) chrs	- A list of chrosomes in which SNPs should be found
#
#	 	          - Returns a 2D hash where the first key is the name of a chromsome and
#		            the second key is the position of the SNP in the first chrosmome and the value 
#		            is the base  in the chrosome. Only SNPs in coding regions are returned
sub get_coding_snps_msa
{
	($msa_id, $r_my_chrs, $dbh) = @_;

	my %snp_hash;
	my %snp_msa_hash;
	my %snp2pos_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

        my %my_chr2inds;
        my @my_chr_inds;
        my @my_chrs = @$r_my_chrs;


        #SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
        my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs = $sth -> fetchrow_array;
        my @chrs = split /\|/, $chrs;

        #GET INDEX OF DESIRED CHROMSOMES IN LIST
        if(@$r_my_chrs == 0)
        {

                @my_chr_inds = 0..(@chrs-1);
                @my_chrs = @chrs;

        }else{

                for(my $c=0; $c < (@chrs); $c++)
                {

                        if(exists($my_chrs_hash{$chrs[$c]}))
                        {

                                $my_chr2inds{$chrs[$c]} = $c;

                        }#end if
        
                }#end foreach

                @my_chr_inds = @my_chr2inds{@my_chrs};

        }#end if

	
        #MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
        my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
        my %chrs_hash = %$r_chrs_hash;

        foreach my $chr (@my_chrs)
        {

                if( !exists($chrs_hash{$chr}) )
                {

                        print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
                        exit;

                }#end if

        }#end foreach


	#SELECT ALL SNPS THAT DO NOT HAVE INDELS AND SAVE IN HASH
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	#KEEP TRACK OF ONLY THOSE SNPS IN CODING REGIONS
	my $r_coding_seq_pos = &get_gene_seq_positions($chrs[$my_chr_inds[0]], $dbh);
	my %coding_seq_pos = %$r_coding_seq_pos;

        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN ANY GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/ || join("", @snps[@my_chr_inds]) =~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i || !exists($coding_seq_pos{$positions[$my_chr_inds[0]]}))
		{
	
			next;
				
		}#end if


		#PUT CURRENT SNP INTO HASH
                $snp_msa_hash{$positions[$my_chr2inds{$my_chrs[0]}]} = join "", @snps[@my_chr_inds];
		$snp2pos_hash{$positions[$my_chr_inds[0]]} = join "\|", @positions[@my_chr_inds];

                foreach my $chr_ind (@my_chr_inds)
                {

                        $snp_hash{$chrs[$chr_ind]}{$positions[$my_chr_inds[0]]} = $snps[$chr_ind];

                }#end for 
		
	}#end while


	return (\%snp_hash, \@my_chrs, \%snp_msa_hash, \%snp2pos_hash);

}#end get_coding_snps_msa


#sub get_4fold_degen_snps_relative_to_reference - Takes as input: 1) ref_chr    - THe name of a chromosome in the database to be used as
#                                                                            	a reference for SNP detection
#                                                            	  2) chrs       - A reference to a list of chromsomes for which SNPs with 
#                                                                            	the reference chromosome are to be found
#		         			- Returns a 2D hash where the first key is the position in the reference genome
#					          the second key is the name of the second chromsome and the value is the base
#					          in the second chrosome. Also returns a hash linking positions in the reference
#					          genome to the base at that position
sub get_4fold_degen_snps_relative_to_reference
{
	($ref_chr, $r_chrs, $dbh) = @_;

	my @chrs = @$r_chrs;

	my %snp_hash;
	my %ref_base_hash;
	my %indel_hash;
	my %position_hash;


	#SELECT ALL SUBSTITUTION SNPs INVOLVING THE REFERENCE GENOME
	my $select_snps = "SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
      				 s.chr_id2 = c2.chr_id AND
      				 c1.name = ? AND 
				 c2.name IN (";
	
        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";

	$select_snps = $select_snps . 		   
			  "UNION
			   SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
                                 s.chr_id2 = c2.chr_id AND
                                 c2.name = ? AND
				 c1.name IN (";

        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ref_chr, @chrs, $ref_chr, @chrs)
                         or die "Can't execute statment: $DBI::errstr";


	#GET ALL THE CODING SNPS
	my $r_four_fold_seq_pos = &get_gene_4fold_degenerate_positions($chr, $dbh);
	my %four_fold_seq_pos = %$r_four_fold_seq_pos;

        while(my ($chr1,$base1, $position1, $chr2, $base2, $position2) =  $sth -> fetchrow_array)
	{


		#DETERMINE THE EFFECT OF THE SEQUENCE CHANGE AND FILTER FOR CODING SNPS
		if($chr1 eq $ref_chr)
		{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(!exists($four_fold_seq_pos{$position1}))
			{
	
				next;

			}#end if


			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
			{

				$snp_hash{$chr2}{$position1} = $base2;
				$ref_base_hash{$position1} = $base1;

			}else{

				$indel_hash{$position1} = 1;

			}#end if

		}else{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(!exists($four_fold_seq_pos{$position2}))
			{
	
				next;

			}#end if

			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
                        {

				$snp_hash{$chr1}{$position2} = $base1;
				$ref_base_hash{$position2} = $base2;

			}else{

                                $indel_hash{$position2} = 1;

                        }#end if

		}#end if
	
	}#end while


	#GET INDELS RELATIVE TO REFERENCE AND ADD TO INDEL POSITIONS, AND AT THE SAME TIME MAKE SURE THAT SNPS WERE FOUND FOR
	# ALL INPUTTED CHROMSOSOMES
	my @indel_pos = keys %indel_hash;

	foreach my $chr (@chrs)
	{

		#MAKE SURE THAT SNPS WERE FOUND FOR CURRENT CHROMSOME, AND IF NOT THROW ERROR
		if(!exists($snp_hash{$chr}))
		{

			print "\n\nNo pairwise SNPs found between $ref_chr and $chr!!!!\n\n";

			exit;

		}#end if
		

		#GET INDELS AND ADD TO LIST
		my ($r_chr1_pos_hash, $r_chr2_pos_hash) = &filter_snps_large_indels($ref_chr, $chr, $dbh);

		push @indel_pos, keys(%$r_chr1_pos_hash);

	}#end foreach	

	$r_indel_pos = &hash_lib::init_hash(\@indel_pos);
	%indel_hash = %$r_indel_pos;


	return (\%snp_hash, \%ref_base_hash, \%indel_hash);


}#end get_4fold_degen_snps_relative_to_reference


#sub get_4fold_degen_snps_msa  - Takes as input: 1) msa_id	- The id of the MSA to use for snp analysis
#					    	 2) chrs	- A list of chrosomes in which SNPs should be found
#
#	 	          	- Returns a 2D hash where the first key is the name of a chromsome and
#		            	  the second key is the position of the SNP in the first chrosmome and the value 
#		            	  is the base  in the chrosome. Only SNPs at 4-fold degenerate sites are reported.
sub get_4fold_degen_snps_msa
{
	($msa_id, $r_my_chrs, $dbh) = @_;

	my %snp_hash;
	my %snp_msa_hash;
	my %snp2pos_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

        my %my_chr2inds;
        my @my_chr_inds;
        my @my_chrs = @$r_my_chrs;


        #SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
        my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs = $sth -> fetchrow_array;
        my @chrs = split /\|/, $chrs;

        #GET INDEX OF DESIRED CHROMSOMES IN LIST
        if(@$r_my_chrs == 0)
        {

                @my_chr_inds = 0..(@chrs-1);
                @my_chrs = @chrs;

        }else{

                for(my $c=0; $c < (@chrs); $c++)
                {

                        if(exists($my_chrs_hash{$chrs[$c]}))
                        {

                                $my_chr2inds{$chrs[$c]} = $c;

                        }#end if
        
                }#end foreach

                @my_chr_inds = @my_chr2inds{@my_chrs};

        }#end if


        #MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
        my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
        my %chrs_hash = %$r_chrs_hash;

        foreach my $chr (@my_chrs)
        {

                if( !exists($chrs_hash{$chr}) )
                {

                        print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
                        exit;

                }#end if

        }#end foreach


	#SELECT ALL SNPS THAT DO NOT HAVE INDELS AND SAVE IN HASH
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	#KEEP TRACK OF ONLY THOSE SNPS IN CODING REGIONS
	my $r_four_fold_seq_pos = &get_gene_4fold_degenerate_positions($chrs[$my_chr_inds[0]], $dbh);
        my %four_fold_seq_pos = %$r_four_fold_seq_pos;


        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN ANY GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/ || join("", @snps[@my_chr_inds]) =~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i || !exists($four_fold_seq_pos{$positions[$my_chr_inds[0]]}))
		{
	
			next;
				
		}#end if

	
		#PUT CURRENT SNP INTO HASH
                $snp_msa_hash{$positions[$my_chr2inds{$my_chrs[0]}]} = join "", @snps[@my_chr_inds];
		$snp2pos_hash{$positions[$my_chr_inds[0]]} = join "\|", @positions[@my_chr_inds];

                foreach my $chr_ind (@my_chr_inds)
                {

                        $snp_hash{$chrs[$chr_ind]}{$positions[$my_chr_inds[0]]} = $snps[$chr_ind];

                }#end for 
		
	}#end while


	return (\%snp_hash, \@my_chrs, \%snp_msa_hash, \%snp2pos_hash);

}#end get_4fold_degen_snps_msa


#sub get_noncoding_snps_relative_to_reference - Takes as input: 1) ref_chr    - THe name of a chromosome in the database to be used as
#                                                                      	   	a reference for SNP detection
#                                                         	2) chrs       - A reference to a list of chromsomes for which SNPs with 
#                                                                         	the reference chromosome are to be found
#		     			       - Returns a 2D hash where the first key is the position in the reference genome
#					        the second key is the name of the second chromsome and the value is the base
#					        in the second chrosome. Also returns a hash linking positions in the reference
#					        genome to the base at that position
sub get_noncoding_snps_relative_to_reference
{
	($ref_chr, $r_chrs, $dbh) = @_;

	my @chrs = @$r_chrs;

	my %snp_hash;
	my %ref_base_hash;
	my %indel_hash;
	my %position_hash;


	#SELECT ALL SUBSTITUTION SNPs INVOLVING THE REFERENCE GENOME
	my $select_snps = "SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
      				 s.chr_id2 = c2.chr_id AND
      				 c1.name = ? AND 
				 c2.name IN (";
	
        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";

	$select_snps = $select_snps . 		   
			  "UNION
			   SELECT c1.name,  base1, position1, c2.name, base2, position2
                           FROM snp s, chromosome c1, chromosome c2
                           WHERE s.chr_id1 = c1.chr_id AND 
                                 s.chr_id2 = c2.chr_id AND
                                 c2.name = ? AND
				 c1.name IN (";

        $select_snps = &db_lib::add_place_holders($select_snps , scalar(@chrs));

        $select_snps = $select_snps . ")";


        my $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($ref_chr, @chrs, $ref_chr, @chrs)
                         or die "Can't execute statment: $DBI::errstr";


	#GET ALL THE NON-CODING SNPS
	my $r_coding_seq_pos = &get_gene_seq_positions($chr, $dbh);
	my %coding_seq_pos = %$r_coding_seq_pos;

        while(my ($chr1,$base1, $position1, $chr2, $base2, $position2) =  $sth -> fetchrow_array)
	{


		#DETERMINE THE EFFECT OF THE SEQUENCE CHANGE AND FILTER FOR CODING SNPS
		if($chr1 eq $ref_chr)
		{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(exists($coding_seq_pos{$position1}))
			{
	
				next;

			}#end if


			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
			{

				$snp_hash{$chr2}{$position1} = $base2;
				$ref_base_hash{$position1} = $base1;

			}else{

				$indel_hash{$position1} = 1;

			}#end if

		}else{

			#IF NOT CODING, THEN GO TO NEXT SNP	
			if(exists($coding_seq_pos{$position2}))
			{
	
				next;

			}#end if

			#IF SNP IS INDEL THEN DON'T SAVE AND NOT THAT THIS POSITION HAS AN INDEL IN SOME ChROMSOME
			if($base1 ne "." && $base2 ne ".")
                        {

				$snp_hash{$chr1}{$position2} = $base1;
				$ref_base_hash{$position2} = $base2;

			}else{

                                $indel_hash{$position2} = 1;

                        }#end if

		}#end if
	
	}#end while


	#GET INDELS RELATIVE TO REFERENCE AND ADD TO INDEL POSITIONS, AND AT THE SAME TIME MAKE SURE THAT SNPS WERE FOUND FOR
	# ALL INPUTTED CHROMSOSOMES
	my @indel_pos = keys %indel_hash;

	foreach my $chr (@chrs)
	{

		#MAKE SURE THAT SNPS WERE FOUND FOR CURRENT CHROMSOME, AND IF NOT THROW ERROR
		if(!exists($snp_hash{$chr}))
		{

			print "\n\nNo pairwise SNPs found between $ref_chr and $chr!!!!\n\n";

			exit;

		}#end if
		

		#GET INDELS AND ADD TO LIST
		my ($r_chr1_pos_hash, $r_chr2_pos_hash) = &filter_snps_large_indels($ref_chr, $chr, $dbh);

		push @indel_pos, keys(%$r_chr1_pos_hash);

	}#end foreach	

	$r_indel_pos = &hash_lib::init_hash(\@indel_pos);
	%indel_hash = %$r_indel_pos;


	return (\%snp_hash, \%ref_base_hash, \%indel_hash);

}#end get_noncoding_snps_relative_to_reference


#sub get_noncoding_snps_msa  - Takes as input: 1) msa_id	- The id of the MSA to use for snp analysis
#					       2) chrs	- A list of chrosomes in which SNPs should be found
#
#	 	             - Returns a 2D hash where the first key is the name of a chromsome and
#		               the second key is the position of the SNP in the first chrosmome and the value 
#		               is the base  in the chrosome. Only SNPs in non-coding regions are returned
sub get_noncoding_snps_msa
{
	($msa_id, $r_my_chrs, $dbh) = @_;

	my %snp_hash;
	my %snp_msa_hash;
	my %snp2pos_hash;

	my $r_my_chrs_hash = &hash_lib::init_hash($r_my_chrs);
	my %my_chrs_hash = %$r_my_chrs_hash;

        my %my_chr2inds;
        my @my_chr_inds;
        my @my_chrs = @$r_my_chrs;


        #SELECT THE CHROMSOMES THAT SNPS ARE FOUND IN
        my $select_chrs = "SELECT mai.chrs
                           FROM multiple_alignment_info mai
                           WHERE mai.msa_id = ?";


        my $sth = $dbh -> prepare($select_chrs)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


        my $chrs = $sth -> fetchrow_array;
        my @chrs = split /\|/, $chrs;

        #GET INDEX OF DESIRED CHROMSOMES IN LIST
        if(@$r_my_chrs == 0)
        {

                @my_chr_inds = 0..(@chrs-1);
                @my_chrs = @chrs;

        }else{

                for(my $c=0; $c < (@chrs); $c++)
                {

                        if(exists($my_chrs_hash{$chrs[$c]}))
                        {

                                $my_chr2inds{$chrs[$c]} = $c;

                        }#end if
        
                }#end foreach

                @my_chr_inds = @my_chr2inds{@my_chrs};

        }#end if


        #MAKE SURE THAT ALL CHROMSOMES ARE IN THE ALIGNMENT
        my $r_chrs_hash = &hash_lib::init_hash(\@chrs);
        my %chrs_hash = %$r_chrs_hash;

        foreach my $chr (@my_chrs)
        {

                if( !exists($chrs_hash{$chr}) )
                {

                        print "\n\nChromsome $chr is not part of the selected MSA!!!\n\n";
                        exit;

                }#end if

        }#end foreach


	#SELECT ALL SNPS THAT DO NOT HAVE INDELS AND SAVE IN HASH
	my $select_snps = "SELECT mas.snps, mas.positions
                           FROM multiple_alignment_snp mas
                           WHERE mas.msa_id = ?";


        $sth = $dbh -> prepare($select_snps)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($msa_id)
                         or die "Can't execute statment: $DBI::errstr";


	#KEEP TRACK OF ONLY THOSE SNPS IN CODING REGIONS
	my $r_coding_seq_pos = &get_gene_seq_positions($chrs[$my_chr_inds[0]], $dbh);
	my %coding_seq_pos = %$r_coding_seq_pos;

        while(my ($snps, $positions) = $sth -> fetchrow_array)
	{

		#GET THE CURRENT SNPS AND POSITIONS	
		my @snps = split //, $snps;
		my @positions = split /\|/, $positions;


		#SKIP IF SNP IS AMBIGUOUS OR INDEL IN ANY GENOME
		if(join("", @snps[@my_chr_inds]) =~ /N/ || join("", @snps[@my_chr_inds]) =~ /-/ || join("", @snps[@my_chr_inds]) =~ /^(A+|T+|C+|G+)$/i || exists($coding_seq_pos{$positions[$my_chr_inds[0]]}))
		{
	
			next;
				
		}#end if

	
		#PUT CURRENT SNP INTO HASH
                $snp_msa_hash{$positions[$my_chr2inds{$my_chrs[0]}]} = join "", @snps[@my_chr_inds];
		$snp2pos_hash{$positions[$my_chr_inds[0]]} = join "\|", @positions[@my_chr_inds];

                foreach my $chr_ind (@my_chr_inds)
                {

                        $snp_hash{$chrs[$chr_ind]}{$positions[$my_chr_inds[0]]} = $snps[$chr_ind];

                }#end for 
		
	}#end while


	return (\%snp_hash, \@my_chrs, \%snp_msa_hash, \%snp2pos_hash);

}#end get_noncoding_snps_msa


#sub get_snp_class - Takes as input: 1) position	- The position of the snp
#				     2) chr		- The chromsome the snp occurs on
#				     3) base_change	- One of the 4 nucleotides (relative to the forward strand of chr)
#		   - Returns a number indiciating the impact of the base change on the given chromsome at the
#		     given position (1 = non-coding, 2 = synonomous, 3 = non-synonomous)
sub get_snp_class
{
	my ($position, $chr, $base_change, $dbh) = @_;

	my $snp_class =0;
	my $aa_change = "NA";
	my $codon_change = "NA";
	my $gene_pos_range = "NA";


	#GET THE GENE THAT THE SNP OCCURS IN
	my $select_gene = "SELECT g.name, g.start, g.end, g.chr_id
                           FROM gene g, chromosome c
                           WHERE ((g.start < ? AND g.end > ?) OR
				 (g.start > ? AND g.end < ?)) AND
				 c.chr_id = g.chr_id AND
 				 c.name = ?";


        my $sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($position, $position, $position, $position, $chr)
                         or die "Can't execute statment: $DBI::errstr";


	if(($gene_name, $start, $end, $chr_id) = $sth -> fetchrow_array)
	{
		#SNP IS IN A GENE, DETERMINE IF IT IS NON-SYNOMOUS

		#GET THE SEQEUENCE OF THE GENE
		my @gene = ($gene_name);

		my $gene_seq = &extract_sequence($chr_id, $start, $end, $dbh);
		my @gene_seq = split //, $gene_seq;


		#DETERMINE THE CODON CHANGE
		my @codon;
		my @new_codon;
		my $gene_pos;
		my $base_change_gene = $base_change;

		#DETERMINE THE POSITION IN THE GENE OF THE CHANGE
		if($start < $end)
		{

			$gene_pos = $position - $start + 1;
			$gene_pos_range = $gene_pos . " of " .  ($end-$start+1);

		}else{

			$gene_pos = $start - $position + 1;
			$gene_pos_range = $gene_pos . " of " . ($start-$end+1);

			#CHANGE BASE TO IT'S COMPLEMENT IF GENE ON - STRAND
			$base_change_gene = &seq_ops::rev_com($base_change);

		}#end if

		#GET THE OLD AND NEW CODON
		if($gene_pos%3 == 1)
		{

			@codon = @gene_seq[($gene_pos - 1)..($gene_pos + 2 - 1)];

			@new_codon = @codon;
			$new_codon[0] = $base_change_gene;

		}elsif($gene_pos%3 == 2){

			@codon = @gene_seq[($gene_pos - 1 - 1)..($gene_pos + 1 - 1)];

			@new_codon = @codon;
			$new_codon[1] = $base_change_gene;

		}else{

			@codon = @gene_seq[($gene_pos - 2 - 1)..($gene_pos - 1)];

			@new_codon = @codon;
			$new_codon[2] = $base_change_gene;

		}#end if

		
		#TRANSLATE CODONS AND COMPARE
		$trans_codon = &seq_ops::translate((join "", @codon));
		$trans_new_codon = &seq_ops::translate((join "", @new_codon));

		$codon_change = (join "", @codon) . " -> " . (join "", @new_codon);
		$aa_change = "$trans_codon -> $trans_new_codon";


		if($trans_codon eq $trans_new_codon)
		{

			$snp_class = 2;

		}else{

			$snp_class = 3;

		}#end if

	}else{
		#SNP IS NOT CODING
		$snp_class = 1;

	}#end if


	return ($snp_class, $codon_change, $aa_change, $gene_pos_range);

}#end get_snp_class


#sub get_base_info - Takes as input: 1) chr		- A name from the chromosome table
#				     2) position	- A posiiton in chr
#		   - Returns a hash summarizing the position in the chromosome 
sub get_base_info
{
	my ($chr, $position, $dbh) = @_;

	my %base_summary;

	my $genes_5p = 0;	
	my $genes_3p = 0;	


	#GET THE SEQUENCE SURROUNDING THE POSITION
	my $select_seq = "SELECT SUBSTRING(c.sequence, (? - 5), 11)
			  FROM chromosome c
			  WHERE c.name = ?";	

	my $sth = $dbh -> prepare($select_seq)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($position, $chr)
                         or die "Can't execute statment: $DBI::errstr";


        $base_summary{'SEQ'} =  $sth -> fetchrow_array;


	#DETERMINE IF THE POSITION IS INSIDE A GENE, AND IF SO GET THE NAME OF THE GENE
        my $select_gene = "SELECT g.name, g.start, g.end, a.annotation, a.source
                           FROM gene g, chromosome c, annotation a
                           WHERE g.chr_id = c.chr_id AND 
                                 c.name = ? AND
                                 ((g.start < ? AND g.end > ?) OR
                                 (g.start > ? AND g.end < ?)) AND
				 a.gid = g.gid AND
				 (a.source = \'KEGG\' OR a.source = \'NCBI\')
			    LIMIT 1";

	$sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr, $position, $position, $position, $position)
                         or die "Can't execute statment: $DBI::errstr";




        if( my ($gene, $start, $end, $annotation, $source) =  $sth -> fetchrow_array)
	{

		#TRIM ANNOTATION
		if($source =~ /NCBI/)
		{

			$annotation =~ s/^([^\|]+)\|.*$/$1/;
			$annotation =~ s/^(.+\]).+$/$1/;

		}else{

			$annotation =~ s/^(.+?); K\d+.+$/$1/;
                	$annotation =~ s/^([^\)]+)\(.*$/$1/;                                                                                          
                                                                                                                                            
                                                                                                                                            
                        if($annotation=~ /;/)                                                                                                        
                	{
                        
				$annotation =~ s/^([^;]+); ([^;]+)$/$2/;

			}#end id

		}#end if


		$base_summary{'IN_GENE'} = $gene; 
		$base_summary{'IN_ANNOT'} = $annotation;
		$base_summary{'IN_START'} = $start;

		if($start < $end)
		{

			$base_summary{'IN_POS'} = $position - $start + 1;
			$base_summary{'IN_LENGTH'} = $end - $start;
		
		}else{


			$base_summary{'IN_POS'} = $start - $position + 1;
			$base_summary{'IN_LENGTH'} = $start - $end;

		}#end if

		#$base_summary{'IN_ANNOT'} =~ s/^([^\|]+)\|.*$/$1/;
		$base_summary{'IN_ANNOT'} =~ s/^([^\[]+)\[.*$/$1/;


		#GET COG INFORMATION FOR GENE
		my @gene = ($gene);

		my ($r_cog_hash, $r_cog_code_hash, $r_cog_desc_hash, $r_cog_code_info) = &get_gene_COGs(\@gene, $dbh);

		$base_summary{'COG_CODES'} = join "|" , keys(%{$r_cog_code_hash->{$gene}});
		$base_summary{'COGS'} = join "|" , keys(%{$r_cog_hash->{$gene}});
		$base_summary{'COG_DESCS'} = join "|" , keys(%{$r_cog_desc_hash->{$gene}});


	}else{

		$base_summary{'IN_GENE'} = "";
		$base_summary{'IN_POS'} = 0;
		$base_summary{'IN_LENGTH'} = 0;
		$base_summary{'IN_START'} = 0;
		$base_summary{'IN_ANNOT'} = "";

	}#end if

	#DETERMINE THE TWO GENES UPSTREAM AND DOWNSTREAM OF THE POSITION
	my $select_genes = "SELECT g.name, ABS(g.start - ?) AS DIST, (g.start - ?), a.annotation
                            FROM gene g, chromosome c, annotation a
                            WHERE g.chr_id = c.chr_id AND 
                                  c.name = ? AND
				  a.gid = g.gid AND
				  (a.source = \'KEGG\' OR a.source = \'NCBI\')
			    ORDER BY DIST";

	$sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($position, $position, $chr)
                         or die "Can't execute statment: $DBI::errstr";



	#GO THROUGH GENES IN ORDER OF PROXIMITY TO GENE OF INTEREST, UNTIL FINDING TWO CLOSEST UPSTREAM AND DOWNSTREAM GENES
        while(my ($gene, $abs_dist, $dist, $annotation) =  $sth -> fetchrow_array)
	{

		#CHECK IF GENE CONTAINS SNP, IF SO SKIP
		if($gene eq $base_summary{'IN_GENE'})
		{

			next;

		}#end if


		#TRIM ANNOTATION
		$annotation =~ s/^([^\|]+)\|.*$/$1/;
		$annotation =~ s/^(.+\]).+$/$1/;


		if( ($dist < 0) && ($genes_5p < 2) )
		{

			if($genes_5p == 0)
			{

				$base_summary{'5P_GENE1'} = $gene; 
				$base_summary{'5P_ANNOT1'} = $annotation; 
				$genes_5p ++;

			}else{

				$base_summary{'5P_GENE2'} = $gene; 
				$base_summary{'5P_ANNOT2'} = $annotation; 
				$genes_5p ++;

			}#end if


		}elsif( ($dist > 0) && ($genes_3p < 2) ){

			if($genes_3p == 0)
			{

				$base_summary{'3P_GENE1'} = $gene; 
				$base_summary{'3P_ANNOT1'} = $annotation; 
				$genes_3p ++;

			}else{


				$base_summary{'3P_GENE2'} = $gene; 
				$base_summary{'3P_ANNOT2'} = $annotation; 
				$genes_3p ++;

			}#end if

		}elsif( ($genes_5p >= 2) && ($genes_3p >= 2)){


			last;

		}#end if


	}#end while


	return \%base_summary;

}#end get_base_info


#sub get_base_anot - Takes as input: 1) chr		- A name from the chromosome table
#				     2) position	- A posiiton in chr
#		   - Returns a hash summarizing the position in the chromosome 
sub get_base_annot
{
	my ($chr, $position, $dbh) = @_;

	my %base_summary;


	#DETERMINE IF THE POSITION IS INSIDE A GENE, AND IF SO GET THE NAME OF THE GENE
        my $select_gene = "SELECT g.name, a.annotation, a.source
                           FROM gene g, chromosome c, annotation a
                           WHERE g.chr_id = c.chr_id AND 
                                 c.name = ? AND
                                 ((g.start < ? AND g.end > ?) OR
                                 (g.start > ? AND g.end < ?)) AND
				 a.gid = g.gid AND
				 (a.source = \'KEGG\' OR a.source = \'NCBI\')
			    LIMIT 1";

	$sth = $dbh -> prepare($select_gene)
                         or die "Can't prepare statment: $DBI::errstr";

        $rc = $sth -> execute($chr, $position, $position, $position, $position)
                         or die "Can't execute statment: $DBI::errstr";



        if( my ($gene, $annotation, $source) =  $sth -> fetchrow_array)
	{

		#TRIM ANNOTATION
		if($source =~ /NCBI/)
		{

			$annotation =~ s/^([^\|]+)\|.*$/$1/;
			$annotation =~ s/^(.+\]).+$/$1/;

		}else{

			$annotation =~ s/^(.+?); K\d+.+$/$1/;
                	$annotation =~ s/^([^\)]+)\(.*$/$1/;                                                                                          
                                                                                                                                            
                                                                                                                                            
                        if($annotation=~ /;/)                                                                                                        
                	{
                        
				$annotation =~ s/^([^;]+); ([^;]+)$/$2/;

			}#end id

		}#end if


		$base_summary{'IN_GENE'} = $gene; 

		$base_summary{'IN_ANNOT'} = $annotation;
		#$base_summary{'IN_ANNOT'} =~ s/^([^\|]+)\|.*$/$1/;
		$base_summary{'IN_ANNOT'} =~ s/^([^\[]+)\[.*$/$1/;


		#GET COG INFORMATION FOR GENE
		my @gene = ($gene);

		my ($r_cog_hash, $r_cog_code_hash, $r_cog_desc_hash, $r_cog_code_info) = &get_gene_COGs(\@gene, $dbh);

		$base_summary{'COG_CODES'} = join "|" , keys(%{$r_cog_code_hash->{$gene}});
		$base_summary{'COGS'} = join "|" , keys(%{$r_cog_hash->{$gene}});
		$base_summary{'COG_DESCS'} = join "|" , keys(%{$r_cog_desc_hash->{$gene}});


	}else{

		$base_summary{'IN_GENE'} = "";
		$base_summary{'IN_POS'} = 0;
		$base_summary{'IN_LENGTH'} = 0;
		$base_summary{'IN_START'} = 0;
		$base_summary{'IN_ANNOT'} = "";

	}#end if

	return \%base_summary;

}#end get_base_annot



#sub get_base_depth - Takes as input: chr_name 			- The name of the chromosome the base resides on
#				       pos			- The position of the chromosome the base resides on
#				       contig_depth_hash	- A 2D hash where the first key is the contig name for hte
#							          given chromosome and the second key is the position, with
#								  the value being the sequencing depth on that contig at that position	 	
#		    - Returns the depth of sequencing at the given position
sub get_base_depth
{
	my ($chr_name, $pos, $r_contig_depth_hash, $dbh) = @_;

	my $depth = -1;
	my %contig_depth_hash = %$r_contig_depth_hash;


	#CONVERT THE CHROMSOME POSITION TO CONTIG SPACE
	my ($contig, $contig_start, $contig_end, $contig_strand, $contig_desc) = &chr2contig_coords(&get_chr_id($chr_name, $dbh), $pos, $pos, $dbh);


	#GET THE DEPTH
	$depth = $contig_depth_hash{$contig}{$contig_start};


	return $depth; 

}#end get_base_depth


#sub get_quality_scores - Takes as input: 1) chr - The name of a chromosome
#			- Returns a list where the index is chromosome position and the value is base quality
sub get_quality_scores
{
	my ($chr, $dbh) = @_;

	my @quality_scores;


	#GET CONTIG POSITIONS AND QUALITY SCORES
        my $select_qual = "SELECT c.start, c.end, c.quality, c.contig_name
                           FROM chromosome chr, contig c
                           WHERE c.chr_id = chr.chr_id AND
			   	 chr.name = ?";


        my $sth = $dbh -> prepare($select_qual)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";




        while(my ($start, $end, $quality) =  $sth -> fetchrow_array)
        {

		#IF NO QUALITY SCORES FOR CONTIG, THEN ABORT
		if((split / /, $quality) <= 1)
		{

			return 0;

		}#end if


		if($start < $end)
		{

			@quality_scores[$start..$end] = split / /, $quality;

		}else{

			@quality_scores[$end..$start] = split / /, $quality;

		}#end if


	}#end while

	return \@quality_scores;

}#end get_quality_scores


#sub get_depth - Takes as input: 1) chr - The name of a chromosome
#	       - Returns a list where the index is chromosome position and the value is base depth
sub get_depth
{
	my ($chr, $dbh) = @_;

	my @depth;


	#GET CONTIG POSITIONS AND DEPTH
        my $select_depth = "SELECT c.start, c.end, c.depth, c.contig_name
                            FROM chromosome chr, contig c
                            WHERE c.chr_id = chr.chr_id AND
			   	  chr.name = ?";


        my $sth = $dbh -> prepare($select_depth)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr)
                         or die "Can't execute statment: $DBI::errstr";




        while(my ($start, $end, $depth) =  $sth -> fetchrow_array)
        {

		#IF NO DEPTH FOR CONTIG, THEN ABORT
		if((split / /, $depth) <= 1)
		{

			return 0;

		}#end if


		if($start < $end)
		{

			@depth[$start..$end] = split / /, $depth;

		}else{

			@depth[$end..$start] = split / /, $depth;

		}#end if


	}#end while

	return \@depth;

}#end get_depth

##################################################### GFF SELECT QUERIES ###########################################

#sub get_gene_gff - Takes as input: 1) chr_id	- An id from the chromosome table
#		  - Returns a hash with keys corresponding to gff columns
#		    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
sub get_gene_gff
{
	my ($chr_id, $dbh) = @_;
	
	my %gff_hash;


	#GET GENE INFO AND COORDINATES
        my $select_genes = "SELECT name, start, end, type
                             FROM gene 
                             WHERE chr_id = ?";


        my $sth = $dbh -> prepare($select_genes)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $start, $end, $type) =  $sth -> fetchrow_array)
        {

		#SOURCE AND FEATURE
		if($type eq 'protein')
		{

	                $gff_hash{$name}{'source'} = 'Glimmer3';
			$gff_hash{$name}{'feature'} = 'CDS';

		}elsif($type eq 'tRNA'){

	                $gff_hash{$name}{'source'} = 'tRNAscanSE';
			$gff_hash{$name}{'feature'} = 'tRNA';

		}#end if


		#SCORE
		$gff_hash{$name}{'score'} = "NA";


		#START, END AND STRAND
		if($start < $end)
		{

			$gff_hash{$name}{'strand'} = "+";
			$gff_hash{$name}{'start'} = $start;
			$gff_hash{$name}{'end'} = $end;

		}else{

			$gff_hash{$name}{'strand'} = "-";
			$gff_hash{$name}{'end'} = $start;
			$gff_hash{$name}{'start'} = $end;

		}#end if


		#FRAME
		$gff_hash{$name}{'frame'} = "";


        }#end while


	#GET THE ANNOTATIONS FOR EACH GENE TO PUT IN THE ATTRIVUTES COLUMN
	my @genes = keys %gff_hash;

	my $r_annotations = &get_gene_annotations(\@genes, $dbh);
	my %annotations = %$r_annotations;
	
	foreach my $gene (@genes)
	{

		if($gff_hash{$gene}{'feature'} eq 'CDS')
		{

			#SHORTEN THE NR ANNOTATION STRING
			my $nr_annot;

			if(exists($annotations{$gene}{'NCBI'}))
			{

				$nr_annot = $annotations{$gene}{'NCBI'};

			}else{

				$nr_annot = $annotations{$gene}{'NR'};
				$nr_annot =~ s/^([^\[]+)\[.*$/$1/;

			}#end if


			#MAKE GFF ENTRY FOR EACH PROTEIN, COLORING DEPENDING ON FUNCTION
			$gff_hash{$gene}{'attributes'} = "label=$gene;cog=$annotations{$gene}{'COG'};nr=$nr_annot";

			if( exists($annotations{$gene}{'ARDB'}) )
			{

				$gff_hash{$gene}{'attributes'} = $gff_hash{$gene}{'attributes'} . ";ardb=$annotations{$gene}{'ARDB'};color=150 0 0";

			}elsif($annotations{$gene}{'NR'} =~ /IntI|integron/i || $annotations{$gene}{'CDD'} =~ /IntI|integron/i || $annotations{$gene}{'NCBI'} =~ /IntI|integron/i){

				$gff_hash{$gene}{'attributes'} = $gff_hash{$gene}{'attributes'} . ";color=0 200 100";

			}elsif($annotations{$gene}{'NR'} =~ /transposase|integrase|recombinase/i || $annotations{$gene}{'CDD'} =~ /transposase|integrase|recombinase/i || $annotations{$gene}{'NCBI'} =~ /transposase|integrase|recombinase/i ){

				$gff_hash{$gene}{'attributes'} = $gff_hash{$gene}{'attributes'} . ";color=0 250 150";

			}elsif($annotations{$gene}{'NR'} =~ /phage/i || $annotations{$gene}{'NCBI'} =~ /phage/i ){

				$gff_hash{$gene}{'attributes'} = $gff_hash{$gene}{'attributes'} . ";color=100 100 100";

			}else{

				$gff_hash{$gene}{'attributes'} = $gff_hash{$gene}{'attributes'} . ";color=0 0 150";

			}#end if
	

		}elsif($gff_hash{$gene}{'feature'} eq 'tRNA'){
		
			$gff_hash{$gene}{'attributes'} = "label=$gene;function=$annotations{$gene}{'tRNAscanSE'};color=0 150 0";

		}#end if


	}#end if


	return \%gff_hash;

}#end get_gene_gff


#sub get_contig_gff - Takes as input: 1) chr_id	- An id from the chromosome table
#	 	    - Returns a hash with keys corresponding to gff columns
#		     <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
sub get_contig_gff
{
	my ($chr_id, $dbh) = @_;
	
	my %gff_hash;


	#GET GENE INFO AND COORDINATES
        my $select_contig = "SELECT contig_name, start, end, depth
                             FROM contig
                             WHERE chr_id = ? AND
				   contig_name NOT LIKE \'scaffold%\'";


        my $sth = $dbh -> prepare($select_contig)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $start, $end, $depth) =  $sth -> fetchrow_array)
        {

		#print "$name, $start, $end, $depth\n";
		#SOURCE
	        $gff_hash{$name}{'source'} = 'Newbler';


		#FEATURE
		$gff_hash{$name}{'feature'} = 'contig';


		#SCORE
		$gff_hash{$name}{'score'} = "NA";


		#STRAND, START AND END
		if($start < $end)
		{

			$gff_hash{$name}{'strand'} = "+";
			$gff_hash{$name}{'start'} = $start;
			$gff_hash{$name}{'end'} = $end;

		}else{

			$gff_hash{$name}{'strand'} = "-";
			$gff_hash{$name}{'start'} = $end;
			$gff_hash{$name}{'end'} = $start;

		}#end if


		#FRAME
		$gff_hash{$name}{'frame'} = "";


		#ATTRIBUTES
		$gff_hash{$name}{'attributes'} = "label=$name;depth=$depth";;

	}#end while


	return \%gff_hash;

}#end get_contig_gff


#sub get_scaffold_gff - Takes as input: 1) chr_id	- An id from the chromosome table
#	 	       - Returns a hash with keys corresponding to gff columns
#		        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
sub get_scaffold_gff
{
	my ($chr_id, $dbh) = @_;
	
	my %gff_hash;


	#GET GENE INFO AND COORDINATES
        my $select_contig = "SELECT contig_name, start, end, depth
                             FROM contig
                             WHERE chr_id = ? AND
				   contig_name LIKE \'%scaffold%\'";


        my $sth = $dbh -> prepare($select_contig)
                         or die "Can't prepare statment: $DBI::errstr";

        my $rc = $sth -> execute($chr_id)
                         or die "Can't execute statment: $DBI::errstr";



        while(my ($name, $start, $end, $depth) =  $sth -> fetchrow_array)
        {

		#SOURCE
	        $gff_hash{$name}{'source'} = 'Newbler';


		#FEATURE
		$gff_hash{$name}{'feature'} = 'scaffold';


		#SCORE
		$gff_hash{$name}{'score'} = "NA";


		#STRAND, START AND END
		if($start < $end)
		{

			$gff_hash{$name}{'strand'} = "+";
			$gff_hash{$name}{'start'} = $start;
			$gff_hash{$name}{'end'} = $end;

		}else{

			$gff_hash{$name}{'strand'} = "-";
			$gff_hash{$name}{'start'} = $end;
			$gff_hash{$name}{'end'} = $start;

		}#end if


		#FRAME
		$gff_hash{$name}{'frame'} = "";


		#ATTRIBUTES
		$gff_hash{$name}{'attributes'} = "label=$name;depth=$depth";;

	}#end while


	return \%gff_hash;

}#end get_scaffold_gff


############################################################ RECOMBINEERING ROUTINES ######################################################

#sub get_mut_aa_change - Takes as input: 1) chr		- The name of a chromosome
#					 2) gene	- The name of the gene to make the mutation in
#					 3) aa_pos	- THe amino acid position to change
#					 4) new_aa	- The amino acid to change to
#
#		       - Returns the chromosome position(s) and base changes resulting in the desired amino acid change
sub get_mut_aa_change
{
	my ($chr, $gene, $aa_pos, $new_aa, $dbh) = @_;
	
	my @genes = ($gene);

	#GET GENE SEQUENCE AND CODON OF AMINO AICD OF INTEREST
	my $r_gene_seq = &get_gene_sequences(\@genes, $dbh);
	my $gene_seq = $r_gene_seq->{$gene};

	my $old_codon = substr($gene_seq, (($aa_pos-1) * 3), 3);

	#DETERMINE CODON COORDINATES
	#GET GENOMIC COORDINATES OF GENE
        my ($r_gene2coords_hash) = &get_gene_coords(\@genes, $dbh);
        my %gene2coords_hash = %$r_gene2coords_hash;

	#GET POSITION OF CODON
	my @codon_pos;

	if($gene2coords_hash{$gene}{'START'} < $gene2coords_hash{$gene}{'END'})
	{

		$codon_pos[0] = (($aa_pos-1) * 3) + $gene2coords_hash{$gene}{'START'};
		$codon_pos[1] = $codon_pos[0] + 1;
		$codon_pos[2] = $codon_pos[0] + 2;

	}else{

		$codon_pos[0] = $gene2coords_hash{$gene}{'START'} -  (($aa_pos-1) * 3);
		$codon_pos[1] = $codon_pos[0] - 1;
		$codon_pos[2] = $codon_pos[0] - 2;

	}#end if


	#GET CODONS FOR NEW AMINO ACIDS
	my $codon_table_obj = Bio::Tools::CodonTable->new();

	my @new_codons = $codon_table_obj->revtranslate($new_aa);


	#SELECT NEW CODON WHICH HAS MINIMAL DISTANCE FROM CURRENT CODON, AND NOTE MUTATIONS
	my @min_codon;
	my @min_codon_diffs = (2,2,2);


	foreach my $new_codon (@new_codons)
	{

		#DETERMINE THE NUMBER OF SEQUENCE DIFFERENCES BETWEEN THE TWO CODONS
		my $r_codon_diffs = &seq_ops::get_num_seq_diffs($old_codon, $new_codon);

		#KEEP CODON IF MINIMALLY DISTANT FROM ORIGINAL
		if(&list_ops::sum($r_codon_diffs) < &list_ops::sum(\@min_codon_diffs))
		{

			@min_codon = split //, $new_codon;
			@min_codon_diffs = @$r_codon_diffs;

		}#end if
	
	}#end foreach


	#DETERMINE CHANGES BETWEEN OLD AND NEW CODON
	my %base_changes;

	for(my $p = 0; $p <= 2; $p++)
        {

		if($min_codon_diffs[$p])
		{

			$base_changes{$p}{$codon_pos[$p]} = $min_codon[$p];

		}#end if

	}#end for


	return (\%base_changes);

}#end get_mut_aa_change


#sub get_mut_silent_changes - Takes as input: 1) chr	- The name of a chromosome
#					      2) gene	- The name of a gene on that chromosome
#					      3) aa_pos	- The position of the amino acid of interest on the gene
#		
#			    - Returns silent mutations in the desired gene in order of distance from inputted position
sub get_mut_silent_changes
{
	my ($chr, $gene, $aa_pos, $dbh) = @_;

	my @genes = ($gene);

	#GET COORDINATES IN GENE SURROUNDING AMINO ACID OF INTEREST TO LOOK FOR SILENT MUTATIONS

	#GO THROUGH GENE SEQUENCE AND GET ALL SILENT CHANGES
	#GET GENE SEQUENCE
	my $r_gene_seq = &get_gene_sequences(\@genes, $dbh);
	my $gene_seq = $r_gene_seq->{$gene};


	#GET GENOMIC COORDINATES OF GENE
        my ($r_gene2coords_hash) = &get_gene_coords(\@genes, $dbh);
        my %gene2coords_hash = %$r_gene2coords_hash;

	#print "\n$gene2coords_hash{$gene}{'START'}\n$gene_seq\n\n";
	#GET COORDINATES OF SUBSEQUENCE
	my $codon_start = (($aa_pos-1) * 3);

	my $start =  &arith_ops::max(0, $codon_start - 51);
	my $end = &arith_ops::min(length($gene_seq), $codon_start + 51 - 1);

	#GET STRAND GENE IS ON
	my $strand;

	if($gene2coords_hash{$gene}{'START'} < $gene2coords_hash{$gene}{'END'})
	{

		$strand = 1;

	}else{

		$strand = -1;

	}#end if

	#GO THROUGH SUBSEQUENCE NOTING SILENT MUTATIONS AND THEIR DISTANCE FROM THE MUTATION
	my $codon_table_obj = Bio::Tools::CodonTable->new();
	my %base_changes;

	for($gene_p = $start; $gene_p < $end; $gene_p+=3)
	{
	
		#SKIP IF CURRENT CODON IS INPUTTED ONE
		if($gene_p == $codon_start)
		{

			next;

		}#end if

		#GET CURRENT CODON AND AMINO ACID IT ENCODES
		my $old_codon = substr($gene_seq, $gene_p, 3);
		my $aa = $codon_table_obj->translate($old_codon);;

		#DETERMINE ALL SYNONYMOUS CODONS
        	my @new_codons = $codon_table_obj->revtranslate($aa);
		
		#DETERMINE THE SYNONYMOUS CODON WITH THE MAXIMAL NUMBER OF DIFFERENCES
		my @max_codon_diffs = (0, 0, 0);
		my @max_codon;

		foreach my $new_codon (@new_codons)
		{

			#DETERMINE THE NUMBER OF SEQUENCE DIFFERENCES BETWEEN THE TWO CODONS
			$r_codon_diffs = &seq_ops::get_num_seq_diffs($old_codon, $new_codon);

			#KEEP CODON IF MINIMALLY DISTANT FROM ORIGINAL
			if(&list_ops::sum($r_codon_diffs) > &list_ops::sum(\@max_codon_diffs))
			{

				@max_codon = split //, $new_codon;
				@max_codon_diffs = @$r_codon_diffs;
	
			}#end if
	
		}#end foreach


		#NOTE SYNONYMOUS DIFFERENCES FOR CHOSEN CODON
		for(my $codon_p = 0; $codon_p <= 2; $codon_p++)
        	{

			if($max_codon_diffs[$codon_p])
			{
				
				#DETERMINE POSITION ON CHROMOSOME OF DIFFERNCE
				my $chr_pos;

				if($strand == 1)
				{

					$chr_pos = $gene2coords_hash{$gene}{'START'} + $gene_p + $codon_p;

				}else{

					$chr_pos = $gene2coords_hash{$gene}{'START'} - $gene_p - $codon_p;

				}#end if

				my $dist_from_aa = &arith_ops::min(abs($gene_p + $codon_p - $codon_start), abs($gene_p + $codon_p - $codon_start + 2));
				$base_changes{$dist_from_aa}{$chr_pos} = $max_codon[$codon_p];

			}#end if

		}#end for

	}#end for	
	

	return (\%base_changes);

}#end get_mut_silent_changes


#sub get_gene_strand - Takes as input: 1) chr	- The chromosome on which the gene resides
#				       2) gene	- The name of the gene		
#
#		     - Returns the stand on which the gene resides
sub get_gene_strand
{
	my ($chr, $gene, $dbh) = @_;

	#GET POSITIONS OF ORIGIN OF REPLICATION AND TERMINUS
	my $ori = 0;
	my $ter = length(&get_chr_seq($chr, $dbh)) / 2;


	#GET POSITION OF GENE ON CHROMSOME
	my @genes = ($gene);

	my ($r_gene2coords_hash) = &get_gene_coords(\@genes, $dbh);
	my %gene2coords_hash = %$r_gene2coords_hash;


	#DETERMINE WHICH SIDE OF THE ORIGIN THE GENE IS ON (ASSUMES: 1) SEQUENCE IS 5'->3' AND 2) GENE DOES NOT CROSS ORIGIN)
	if($ori < $ter && $gene2coords_hash{$gene}{'START'} > $ori && $gene2coords_hash{$gene}{'START'} < $ter && $gene2coords_hash{$gene}{'START'} < $gene2coords_hash{$gene}{'END'})
	{

		return 1;

	}elsif($ori > $ter && ($gene2coords_hash{$gene}{'START'} > $ori || $gene2coords_hash{$gene}{'START'} < $ter) && $gene2coords_hash{$gene}{'START'} < $gene2coords_hash{$gene}{'END'}){

		return 1;

	}else{

		return -1;

	}#end if

}#end get_gene_strand

#sub design_ss_mut_oligo - Takes as input: 1) chr			- The chromosome for which the oligo is designed
#					   2) gene			- The gene to be mutated
#					   3) oligo_length		- The length of the designed oligo
#					   4) num_prox_mut		- The number of nonsyonymous mutations to make in addition
#									  to the desired coding mutation
#					   5) r_coding_changes		- A 2D hash of coding mutations (second key is chr position)
#					   6) r_non_coding_changes	- A 2D hash of synonymous mutations (first key is rank, second is chr position)
#					   7) oligo_strand		- Indicates whether oligo should be complimented to target desired strand (Y/-1, N/1)
#
#			 - Returns 1) the designed oligo, 2) the original gene sequence, 2) the updated gene sequence
sub design_ss_mut_oligo
{
	my ($chr, $gene, $oligo_length, $num_prox_mut, $r_coding_changes, $r_non_coding_changes, $oligo_strand, $dbh) = @_;


	#GET GENOMIC COORDINATES OF GENE
	my @genes = ($gene);

        my ($r_gene2coords_hash) = &get_gene_coords(\@genes, $dbh);
        my %gene2coords_hash = %$r_gene2coords_hash;

	#EXTRACT SEQUENCE OF GENE OF INTEREST AND FLANKING 100 BP ON EITHER SIDE
	my $start_coord;
	my $end_coord;

	if($gene2coords_hash{$gene}{'START'} < $gene2coords_hash{$gene}{'END'})
	{

		$start_coord = $gene2coords_hash{$gene}{'START'}- 100;
		$end_coord = $gene2coords_hash{$gene}{'END'} + 100;

	}else{

		$start_coord = $gene2coords_hash{$gene}{'START'} + 100;
		$end_coord = $gene2coords_hash{$gene}{'END'} - 100;

	}#end if

	my $orig_sequence = &extract_sequence(&get_chr_id($chr, $dbh), $start_coord, $end_coord, $dbh);
	my $alt_sequence = $orig_sequence;


	#MAKE CODING CHANGES TO OLIGO
	my %coding_changes = %$r_coding_changes;
	my $oligo_mid;

	foreach my $ind (keys %coding_changes)
	{

		foreach my $pos (keys %{$coding_changes{$ind}})
		{

			substr($alt_sequence, abs($pos - $start_coord), 1, $coding_changes{$ind}{$pos});
			$oligo_mid = $pos;

		}#end foreach

	}#end foreach


	#MAKE DESIRED NUMBER OF NONSYNONYMOUS CHANGES TO THE OLIGO
	my %non_coding_changes = %$r_non_coding_changes;

	NC:foreach my $ind (sort {$a<=>$b} keys %non_coding_changes)
	{

		foreach my $pos (keys %{$non_coding_changes{$ind}})
		{

			substr($alt_sequence, abs($pos - $start_coord), 1, $non_coding_changes{$ind}{$pos});

			$num_prox_mut --;

			#IF DESIRED NUMBER OF CHANGES MADE, THEN GET OUT
			if($num_prox_mut == 0)
			{

				last NC;

			}#end if

		}#end foreach

	}#end foreach


	#EXTRACT OLIGO, AND REVERSE COMPLIMENT THEM IF NECCESARY
	my $oligo_seq;

	if($oligo_strand == 1)
	{

		$oligo_seq = substr($alt_sequence, (abs($oligo_mid - $start_coord) - $oligo_length/2), $oligo_length);
	
	}else{

		$oligo_seq = &seq_ops::rev_com(substr($alt_sequence, (abs($oligo_mid - $start_coord) - $oligo_length/2), $oligo_length));

	}#end if


	#EXTRACT MODIFIED GENE SEQUENCE
	my $alt_gene_seq = substr($alt_sequence, abs($gene2coords_hash{$gene}{'START'} - $start_coord), 
					  abs($gene2coords_hash{$gene}{'END'} - $gene2coords_hash{$gene}{'START'}) + 1);

	#EXTRACT ORIGINAL GENE SEQUENCE
	my $orig_gene_seq = substr($orig_sequence, abs($gene2coords_hash{$gene}{'START'} - $start_coord), 
					  abs($gene2coords_hash{$gene}{'END'} - $gene2coords_hash{$gene}{'START'}) + 1);


	return($oligo_seq, $alt_gene_seq, $orig_gene_seq);

}#end design_ss_mut_oligo


1;
