#!/usr/bin/perl -w
#Performs GO enrichment analysis on lists of genes 



#DEFINE PARAMTERS USED TO IDENTIFY DATA IN DATBASE AND ENRICHMENT PARAMTERS

my $exp_date = '9_6';
#my @percentiles = ('one' , 'five' , 'ten' , 'fifteen' , 'twenty' , 'twentyfive' , 'thirty');
#my @perc_sizes = (49 , 243 , 486 , 728 , 971 , 1214, 1456);
#my @perc_sizes = (36 , 178 , 357 , 535 , 714 , 892, 1071);


#POPUKATION FILE FOR ALL GENES
my $all_population_file = '/home/esnitkin/projects/single_KO_analysis/model_comparison/exp.orfs';


#DESIGNATES FILE OF GENE ANNOTATIONS
my $gene_ass_file = '/home/esnitkin/packages/GeneMerge1.2/AssociationFiles/S_cerevisiae.PATH';
#my $gene_ass_file = '/home/esnitkin/packages/GeneMerge1.2/AssociationFiles/E_coli.PATH';
#my $gene_ass_file = '/home/esnitkin/packages/GeneMerge1.2/AssociationFiles/S_cerevisiae.hBP';
#my $gene_ass_file = '/home/esnitkin/packages/GeneMerge1.2/AssociationFiles/S_cerevisiae.hMF';

#DESIGNATES FILE WITH ANNOTATION DESCRIPTIONS
my $gene_desc_file = '/home/esnitkin/packages/GeneMerge1.2/DescriptionFiles/S_cerevisiae.PATH.use';
#my $gene_desc_file = '/home/esnitkin/packages/GeneMerge1.2/DescriptionFiles/E_coli.PATH.use';
#my $gene_desc_file = '/home/esnitkin/packages/GeneMerge1.2/DescriptionFiles/GO.BP.use';
#my $gene_desc_file = '/home/esnitkin/packages/GeneMerge1.2/DescriptionFiles/GO.MF.use';


#DIRECTORIES WHERE GENE MERGE OUTPUT AND PARSED GENE MERGE OUTPUT ARE STORED
my $gene_merge_out_dir = "/home/esnitkin/projects/single_KO_analysis/model_comparison/enrichment_analysis/iND750_exp0p0_model0p4/iND750_exp0p0model0p4_KEGG\_gene_merge_out/";

my $parsed_gene_merge_out_dir = "/home/esnitkin/projects/single_KO_analysis/model_comparison/enrichment_analysis/iND750_exp0p0_model0p4/iND750_exp0p0_model0p4_KEGG\_parsed_gene_merge_out/";


#THE COLUMNS IN THE GENE MERGE OUTPUT TO BE SAVED
my @gene_merge_cols = (2,3,4,5,6);
#my @gene_merge_cols = (4,6);


#THE COLUMN CONTAINING THE RAW E-SCORE IN THE GENE MERGE OUTPUT
my $rawE_col = 4;


#THE DIRECTORY WHERE THE INPUT FILES RESIDE IF THE BATCH METHOD USED
my $input_dir = "/home/esnitkin/projects/single_KO_analysis/model_comparison/enrichment_analysis/iND750_exp0p0_model0p4/";


#MAKE DIRECTORIES FOR GENE MERGE OUTPUT AND PARSED OUTPUT
`mkdir $gene_merge_out_dir`;
`mkdir $parsed_gene_merge_out_dir`;



#PERFORM GENE MERGE ANALYSIS ON FILE IN SELECTED DIRECTORY
&gene_merge_batch($input_dir , $gene_merge_out_dir , $parsed_gene_merge_out_dir , 
		  $gene_ass_file, $gene_desc_file, $all_population_file, \@gene_merge_cols, $rawE_col);



#############################################################SUBROUTINES########################################################

#sub gene_merge_batch - Takes as input a directory
#		      - Performs gene merge analysis on the files in the directory
sub gene_merge_batch
{
	my ($input_dir , $gene_merge_out_dir , $parsed_gene_merge_out_dir , 
	    $gene_ass_file, $gene_desc_file, $population_file, $r_gene_merge_cols, $rawE_col) = @_;
	my @gene_merge_cols = @$r_gene_merge_cols;
	
	#GET THE NAMES OF THE FILE IN THE DIRECTORY TO PERFORM ANALYSIS ON 
	@input_files = `ls $input_dir`;


	#PERFORM GENE MERGE ANALYSIS ON THE INPUT FILES
	foreach my $input_file (@input_files)
	{

		#ADD DIRECTORY TO INPUT FILE NAME
		chomp $input_file;
		my $full_input_file = $input_dir . $input_file;


               #RUN GENE MERGE ON CURRENT LISTS OF ESSENTIAL AND NON-ESSENTIAL GENES
               #./GeneMerge.pl gene-association.file description.file population.file study.file output.filename
               my $gm_out_file = $gene_merge_out_dir .  $input_file  . ".gm";

               `perl /home/esnitkin/packages/GeneMerge1.2/GeneMerge1.2.pl $gene_ass_file $gene_desc_file $population_file $full_input_file $gm_out_file`;


               #PARSE GENE MERGE OUTPUT AND WRITE DESIRED COLUMNS TO A FILE
               my $parsed_gm_out_file = $parsed_gene_merge_out_dir . $input_file . "_gm.parsed";

               &gm_parse_and_output($gm_out_file , $parsed_gm_out_file , \@gene_merge_cols , $rawE_col);


               #SORT OUTPUT FILE BY RAW SCORE ASCENDING
               `sort -n -k 3 $parsed_gm_out_file > temp_sort`;
               `mv temp_sort $parsed_gm_out_file`;

	}#end foreach


}#end gene_merge_batch



#sub gene_merge_DBbatch - Takes as input a directory
#                       - Performs gene merge analysis on the files in the directory
sub gene_merge_DBbatch
{
	my ($r_table_names , $r_percentiles, $exp_date , $gene_merge_out_dir , 
	    $parsed_gene_merge_out_dir , $gene_ass_file ,  $gene_desc_file, $population_file,
	    $r_gene_merge_cols , $rawE_col , $id2orf_table, $r_perc_sizes)			= @_;

	my @table_names = @$r_table_names;
	my @percentiles = @$r_percentiles;
	my @gene_merge_cols = @$r_gene_merge_cols;
	my @perc_sizes = @$r_perc_sizes;



	#GET LISTS OF GENES FROM DIFFERENT EXPERIMENTS AND PERCENTILES AND PERFORM ENRICHMENT ANALYSIS
	foreach my $table (@table_names)
	{
	
		for(my $i = 0;  $i < (@percentiles); $i++)
		{

	
			#GET GENES FOR CURRENT FEATURE SET AND PERCENTILE
			my @genes =  &get_gene_list($table , $id2orf_table , $percentiles[$i] , $perc_sizes[$i]);
	
		
			#WRITE GENES TO FILE
			my $gene_file = $gene_merge_out_dir . $table . "_" . $percentiles[$i] . ".genes";		
	
			&gene_list_to_file(\@genes , $gene_file);
	
			
			#RUN GENE MERGE ON CURRENT LISTS OF ESSENTIAL AND NON-ESSENTIAL GENES
			#./GeneMerge.pl gene-association.file description.file population.file study.file output.filename
			my $gm_out_file = $gene_merge_out_dir .  $table . "_" . $percentiles[$i] . ".gm";
	
			`perl GeneMerge1.2/GeneMerge1.2.pl $gene_ass_file $gene_desc_file $population_file $gene_file $gm_out_file`;

		

			#PARSE GENE MERGE OUTPUT AND WRITE DESIRED COLUMNS TO A FILE
			my $parsed_gm_out_file = $parsed_gene_merge_out_dir . $table . "_" . $percentiles[$i] . "_gm.parsed";
	
			&gm_parse_and_output($gm_out_file , $parsed_gm_out_file , \@gene_merge_cols , $rawE_col); 


			#SORT OUTPUT FILE BY RAW SCORE ASCENDING
			`sort -n -k 3 $parsed_gm_out_file > temp_sort`;
			`mv temp_sort $parsed_gm_out_file`;	
	

		}#end foreach


	}#end foreach


}#end gene_merge_DBbatch


#sub get_gene_list - Takes as input a table name and a column name in the lethality results datbase
#		     as well as a minimum nmber of occurences to call a gene essential or non-essential
#		   - Returns a list of ORF ids from the desired table whose value in the desired column
#		     are greater than the minimum number of occurences specified.
sub get_gene_list
{
	my ($table_name , $id2orf_table ,  $col_name , $num_genes) = @_;

	#CONNECT TO LETHALITY_RESULTS DATBASE
	my $dbh = &lethality_db_lib::leth_results_DB_connect();

	
	#GET THE LIST OF GENES FROM THE DESIRED TABLE WHOSE VALUE IN THE DESIRED COLUMN IS GREATER THAN MIN_OCCUR
	my @gene_list = &lethality_db_lib::get_gene_list($table_name , $id2orf_table , $col_name, $num_genes , $dbh);


	return @gene_list;

}#end get_gene_list


#sub gene_list_to_file - Takes as input a reference to a list of genes and the name of a file 
#		       - Outputs to the file the list of genes, one gene per line
sub gene_list_to_file
{
	my ($r_gene_list , $gene_file) = @_;

	open OUTPUT , ">$gene_file";


	foreach my $gene (@$r_gene_list)
	{

		print OUTPUT $gene . "\n";

	}#end foreach


	close OUTPUT;

}#end gene_list_to_file


#sub gm_parse_and_output - Takes as input the name of a gene merge output file, the columns in the file
#			   which are to be kept and the name of the file to output the columns to
#			 - Outputs columns to specified output file.
sub gm_parse_and_output
{
	my ($gm_output_file , $parsed_gm_output_file , $r_gm_cols , $rawE_col) = @_;
	my @gm_cols = @$r_gm_cols;


	#OPEN INPUT AND OUTPUT FILES
	open PARSED_GM_OUTPUT , ">$parsed_gm_output_file";


	#PUT THE CONTENTS OF THE GENE MERGE OUTPUT IN A VARIABLE AND GET THE MAIN CONTENTS
	my $file = `cat $gm_output_file`;
	my @file = split /\n\n+/ , $file;
	my @split_file = split /\n/ ,$file[2]; 


	#GO THROUGH GENE MERGE OUTPUT FILE WRITING ONLY DESIRED COLUMNS TO OUTPUT FILE
	foreach my $line (@split_file)
	{

		chomp $line;

		$line =~ s/\tNA/\t99999/g;
		
		my @line = split /\t/ , $line;

		if($line[$rawE_col] <= 1)
		{

			my $output_line = join "\t" , @line[@gm_cols];

			print PARSED_GM_OUTPUT $output_line . "\n";
		}

	}#end foreach


	#CLOSE INPUT AND OUTPUT FILES
	close PARSED_GM_OUTPUT;

}#end gm_parse_and_output
