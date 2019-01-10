package mummer;


#sub parse_tiling_file - Takes as input: 1) tiling_file	- A tiling file from from
#		       - Returns statistics from parsing the file
sub parse_tiling_file
{
	my ($tiling_file) = @_;

	
	#GO THROUGH FILE AND KEEP TRACK OF RELEVENT INFO
	my $num_contigs = 0;
	my $total_ref_seq = 0;
	my $total_query_seq = 0;

	#OPEN FILE AND DROP HEADER LINE
	open IN, $tiling_file;
	<IN>;

	foreach my $line (<IN>)
	{

		my ($start_r, $end_r, $gap, $contig_length, $align_cov, $perc_id, $orientation, $contig_id) = split /\t/, $line;

		$num_contigs++;

		$total_ref_seq += $end_r - $start_r;

		$total_query_seq += $contig_length * ($align_cov/100);

	}#end foreach


	close IN;


	return ($num_contigs, $total_ref_seq, $total_query_seq);


}#end parse_tiling_file


#sub get_contigs_from_tiling_file - Takes as input: 1) tiling_file - A tiling file from show-tiling
#				  - Returns an hash reference where the keys are the names of contigs used
#				    in the tiling and the values are the lengths of the corresponding contigs
sub get_contigs_from_tiling_file
{
	my ($tiling_file) = @_;

	my %contigs;


	#OPEN FILE AND DROP HEADER LINE
        open IN, $tiling_file;
        <IN>;

        foreach my $line (<IN>)
        {

                my ($start_r, $end_r, $gap, $contig_length, $align_cov, $perc_id, $orientation, $contig_id) = split /\t/, $line;


		$contigs{$contig_id} = $contig_length;

        }#end foreach


        close IN;


        return \%contigs;

}#end get_contigs_from_tiling_file


#sub parse_coords_file - Takes as input: 1) coords_file - A file produced by show-coords
#		       - Returns a hash where the keys are alignment IDs and the values are coordinates/identities/chromsomes
#			 alignment 
sub parse_coords_file
{
	my ($coords_file) = @_;

	my %coords_hash;


	#OPEN FILE AND REMOVE HEADER
	open COORDS, $coords_file;

	for(my $i = 0; $i < 5; $i++)
	{

		<COORDS>;

	}#end for	


	#GO THROUGH EACH LINE IN COORDS FILE AND PUT IN HASH
	my $aln_id = 1;

	foreach my $line (<COORDS>)
	{

		#GET RID OF NEW-LINE AND LEADING WHITE SPACE
		chomp $line;
		$line =~ s/^\s+//g;
		my @line = split /[\|\s]+/, $line;;
		
		#SPLIT UP LINE AND PUT IN HASH
		my ($start1, $end1, $start2, $end2, $length1, $length2, $perc_id, $total_length1,  $total_length2,  $coverge1,  $coverage2, $seq_id1, $seq_id2);

		if(@line == 9)
		{

			($start1, $end1, $start2, $end2, $length1, $length2, $perc_id, $seq_id1, $seq_id2) = @line;

		}elsif(@line == 13){


			($start1, $end1, $start2, $end2, $length1, $length2, $perc_id, $total_length1,  $total_length2,  $coverge1,  $coverage2, $seq_id1, $seq_id2) = @line;

		}else{

			print "Incorrect output format for show-coords!\n";
			exit;

		}#end if


		$coords_hash{$aln_id}{'S1'} = $start1; 
		$coords_hash{$aln_id}{'E1'} = $end1; 
		$coords_hash{$aln_id}{'S2'} = $start2; 
		$coords_hash{$aln_id}{'E2'} = $end2; 
		$coords_hash{$aln_id}{'L1'} = $length1; 
		$coords_hash{$aln_id}{'L2'} = $length2; 
		$coords_hash{$aln_id}{'ID'} = $perc_id; 
		$coords_hash{$aln_id}{'SEQ1'} = $seq_id1; 
		$coords_hash{$aln_id}{'SEQ2'} = $seq_id2; 


		$aln_id ++;

	}#end foreach
	

	#CLOSE FILE
	close COORDS;


	return \%coords_hash;

}#end parse_coords_file


#sub parse_tandem_file - Takes as input: 1) tandem_file - A file produced by exact-tandem
#		       - Returns a hash where the keys are repeat IDs and the values are start, extent, unit length and number of copies
#			 of the repeat
sub parse_tandem_file
{
	my ($tandem_file) = @_;

	my %tandem_hash;
	

	#OPEN FILE AND REMOVE HEADER
	open TANDEM, $tandem_file;

	for(my $i = 0; $i < 3; $i ++)
	{

		<TANDEM>;

	}#end for


	#GO THROUGH EACH LINE IN TANDEM FILE AND PUT IN HASH
	my $rep_id = 1;

	foreach my $line (<TANDEM>)
	{

		#REMOVE NEW-LINE AND REMOVE LEADING WHITE SPACE
		chomp $line;

		$line =~ s/^\s+//g;


		#SPLIT UP LINE AND STORE IN HASH
		my ($start, $total_length, $unit_length, $num_copies) = split /\s+/, $line;

		$tandem_hash{$rep_id}{'START'} = $start; 
		$tandem_hash{$rep_id}{'TOTAL_LEN'} = $total_length; 
		$tandem_hash{$rep_id}{'UNIT_LEN'} = $unit_length;
		$tandem_hash{$rep_id}{'NUM_COPIES'} = $num_copies;
	
		$rep_id ++;

	}#end foreach


	return \%tandem_hash;

}#end parse_tandem_file


1;
