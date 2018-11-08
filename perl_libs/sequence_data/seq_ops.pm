package seq_ops;


use Bio::Seq;
use POSIX;


#sub translate - Takes as input: 1) seq	- A nucleotide sequence
#	       - Returns a translated nucleotide sequence
sub translate
{
	my ($seq) = @_;

	my $seq_obj = Bio::Seq->new(-seq => $seq);

	my $trans_seq_obj = $seq_obj ->translate;

	my $trans_seq = $trans_seq_obj -> seq;

	return $trans_seq;

}#end translate


#sub translate_hash - Takes as input: 1) r_seq_hash	- A hash where the keys are sequence names and values
#							 are nucleotide sequences
#		    - Returns a hash with the same keys but with protein sequences as the values
sub translate_hash
{
	my ($r_seq_hash) = @_;

	my %seq_hash = %$r_seq_hash;
	my %trans_seq_hash;


	foreach my $seq_name (keys %seq_hash)
	{

		$trans_seq_hash{$seq_name} = &translate($seq_hash{$seq_name});


	}#end foreach	


	return \%trans_seq_hash;

}#end translate_hash 


#sub get_di_nucs - Takes as input:
#		 - Returns a list of di-nucleotides
sub get_di_nucs
{
	my @nucs = ('A', 'T', 'G', 'C');
	my @di_nucs = ();

	foreach my $nuc1 (@nucs)
	{

		foreach my $nuc2 (@nucs)
		{

			push @di_nucs, ($nuc1 . $nuc2);

		}#end foreach


	}#end foreach


	return @di_nucs;

}#end get_di_nucs


#sub get_kmers - Takes as input: 1) kmer_size - The size kmers desired
#	       - Returns a list of all possible nucleotide  kmers of desired size
sub get_kmers
{
	my ($kmer_size) = @_;

	my @bases = ('A', 'C', 'G', 'T');
	my @words = @bases;

     
	for (1 .. $kmer_size-1 ) 
	{
         
		my @newwords;

		foreach my $w (@words)
		{
             
			foreach my $b (@bases) 
			{
                 
				push (@newwords, $w.$b);

			}#end foreach

		}#end foreach

		@words = @newwords;

	}#end for

	
	return \@words;

}#end get_kmers


#sub rev_com - Takes as input: 1) seq	- A nucleotide sequence
#	     - Returns the reverse complement of the sequence
sub rev_com
{
	my ($seq) = @_;

	my $seq_obj = Bio::Seq->new(-seq => $seq);

	my $rev_seq_obj = $seq_obj->revcom;

	return $rev_seq_obj->seq;

}#end rev_com


#get_phylo_disc - Takes as input: seq_hash - A hash where values represent an alignment at a single position
#					     * all alignments must be the same length
#		- Returns a hash with the same keys, but the values are numbers indeicating the relationship
#		  between positions in the alignment (which base is the same as which)
sub get_phylo_disc
{
	my ($r_seq_hash) = @_;

	my %seq_hash = %$r_seq_hash;
	my @seq_hash_keys = keys %seq_hash;


	#CREATE HASH LINKING POSSIBLE ALIGNMENT STRUCTURES TO A NUMBER
	#GET AN ARBITRARY SEQUENCE TO SEE HOW LONG ALIGNMENT IS
	my $aln_length = length($seq_hash{$seq_hash_keys[0]});


	#CREATE A LIST OF ALL POSSIBLE ALIGNMENT STRUCTURES
	my @structs = (1);
     
	foreach my $pos (2 .. $aln_length) 
	{
         
		my @new_structs;

		foreach my $s (@structs)
		{
             
			foreach my $n (1 .. $pos) 
			{
                 
				my $n_minus1 = $n - 1;

				#ONLY ADD ONTO STRING IF PREVIOUS NUMBER IS PART OF IT 
				if(($n == 1) || ($s =~ /$n_minus1/))
				{

					push (@new_structs, $s.$n);

				}#end if

			}#end foreach

		}#end foreach

		@structs = @new_structs;

	}#end for

	#ASSIGN STRUCTURES TO NUMBER
	my $c = 1;
	my %struct_hash;

	foreach my $struct (@structs)
	{

		$struct_hash{$struct} = $c;
		#print "$c = $struct\n";			

		$c ++;

	}#end foreach
	

	#FOR EACH ALIGNMENT, DISCRETIZE AND ASSIGN TO PHYLO BIN
	my %disc_aln_hash;

	foreach my $key (keys %seq_hash)
	{

		#SPLIT ALIGNMENT
		$seq_hash{$key} =~ tr/a-z/A-Z/;
		my @aln = split //, $seq_hash{$key};

	
		#GO THROUGH ALIGNMENT AND DISCRETIZE
		my $nuc_count = 1;
		my %nuc_hash;
		my $disc_aln = '';

		foreach my $pos (@aln)
		{

			if(exists($nuc_hash{$pos}))
			{

				$disc_aln .= $nuc_hash{$pos};

			}else{

				$nuc_hash{$pos} = $nuc_count;
				$nuc_count ++;

				$disc_aln .= $nuc_hash{$pos};

			}#end if

		}#end foreach

		
		#SAVE DISCRETIZED ALINGMENT
		#print "$seq_hash{$key}, $disc_aln, $struct_hash{$disc_aln}\n";
		$disc_aln_hash{$key} = $struct_hash{$disc_aln};

	}#end foreach

	return (\%disc_aln_hash, \%struct_hash);

}#end get_phylo_disc


#get_phylo_disc2 - Takes as input: seq_hash - A hash where values represent an alignment at a single position
#					     * all alignments must be the same length
#		- Returns a hash with the same keys, but the values are numbers indeicating the relationship
#		  between positions in the alignment (which base is the same as which)
sub get_phylo_disc2
{
	my ($r_seq_hash) = @_;

	my %seq_hash = %$r_seq_hash;
	my @seq_hash_keys = keys %seq_hash;


	#FOR EACH ALIGNMENT, DISCRETIZE AND ASSIGN TO PHYLO BIN
	my %struct_hash;
	my $struct_count = 1;	

	my %disc_aln_hash;
	

	foreach my $key (keys %seq_hash)
	{

		#SPLIT ALIGNMENT
		$seq_hash{$key} =~ tr/a-z/A-Z/;
		my @aln = split //, $seq_hash{$key};

	
		#GO THROUGH ALIGNMENT AND DISCRETIZE
		my $nuc_count = 1;
		my %nuc_hash;
		my $disc_aln = '';

		foreach my $pos (@aln)
		{

			if(exists($nuc_hash{$pos}))
			{

				$disc_aln .= $nuc_hash{$pos};

			}else{

				$nuc_hash{$pos} = $nuc_count;
				$nuc_count ++;

				$disc_aln .= $nuc_hash{$pos};

			}#end if

		}#end foreach

		
		#CHECK IF CODE EXISTS FOR CURRENT ALIGNMENT STRUCTURE	
		if(!exists($struct_hash{$disc_aln}))
		{

			$struct_hash{$disc_aln} = $struct_count;
			$struct_count ++;

		}#end if


		#SAVE DISCRETIZED ALINGMENT
		$disc_aln_hash{$key} = $struct_hash{$disc_aln};

	}#end foreach

	return (\%disc_aln_hash, \%struct_hash);

}#end get_phylo_disc2


#sub cap_middle	- Takes as input: 1) seq - A sequence of odd length
#		- Returns the same sequence, except the entire sequence is lower case except the midle character
sub cap_middle
{
	my ($seq) = @_;

	#GET THE SIZE OF THE SURROUNDING SEQUENCE TO KEEP LOWERCASE
	my $seq_length = length($seq);
	my $surr_seq_len = floor($seq_length / 2);
	

	#GET THE MIDDLE BASE
	my $base = $seq;
        $base =~ s/^.{20}(.).{20}$/$1/;
        my $uc_base = uc $base;

	
	#MAKE THE SEQUENCE LOWER CASE, EXCEPT THE MIDDLE BASE
        $seq =~ tr/[A-Z]/[a-z]/;
        $seq =~ s/^(.{20}).(.{20})$/$1$uc_base$2/;


	return $seq;

}#end cap_middle


#sub get_num_seq_diffs - Takes as input: 1) seq1 - A sequence
#		         		 2) seq2 - A sequence the same length as seq1
#
#		       - Returns an array where 1's represent identities between two sequences and 0's differences
sub get_num_seq_diffs
{
	my ($seq1, $seq2) = @_;

	#GO THROUGH CODON AND DETERMINE NUMBER OF DIFFERENCES
	my @seq1 = split //, $seq1;
	my @seq2 = split //, $seq2;
	
	my @seq_diffs = ();
	
	for(my $p = 0; $p < @seq1; $p++)
	{
	
		if($seq1[$p] !~ /^$seq2[$p]$/i)
		{
	
			$seq_diffs[$p] = 1;    
	
		}else{
	
			$seq_diffs[$p] = 0;    
	
		}#end if
		
	}#end for

	return \@seq_diffs

}#end get_num_seq_diffs


#sub downsample_fastq - Takes as input: 1) fastq1 	- Forward reads
#					2) fastq2 	- Reverse reads
#					3) coverage	- The desired sequencing coverage
#					4) genome_size	- Estimated genome size for coverage calculations
#		      - Takes the first N lines of the fastq files, yielding desired coverage
sub downsample_fastq
{
	my ($fastq1, $fastq2, $coverage, $genome_size) = @_;

	#CREATE OUTPUT FILES
	my $fastq1_out = $fastq1; 
	$fastq1_out =~ s/.fastq$/_$coverage\X.fastq/;

	my $fastq2_out = $fastq2; 
	$fastq2_out =~ s/.fastq$/_$coverage\X.fastq/;
	
	#DETERMINE NUMBER OF LINES IN FASTQ FILE TO GET DESIRED COVERAGE
	my $read_length = length(`head -2 $fastq1 | tail -1`);
	
	my $num_seqs = int(($coverage*$genome_size)/$read_length);
	my $num_lines_sample = $num_seqs * 2; #x 4 for fastq format and /2 for paired reads


	#DETERMINE THE NUMBER OF LINES IN THE FASTQ FILE
	my $num_lines = `wc -l $fastq1`;
	chomp $num_lines;
	$num_lines =~ s/^(\d+) .*$/$1/;


	#DETERMINE ACTUAL COVERAGE FOR SINGLE OR PAIRED LIBRARIES
	my $total_coverage;

	if(-e $fastq2)
	{

		$total_coverage = int((($num_lines / 4) * 2 * $read_length) / $genome_size);

	}else{

		$total_coverage = int((($num_lines / 4) * $read_length) / $genome_size);

	}#end if


	#IF SUFFICIENT COVERAGE, THEN RETURN DOWNSAMPLES FASTQ FILES
	if($total_coverage > $coverage)
	{

		print "\n\nDownsampled from $total_coverage\x to $coverage\x\n\n";

		if(-e $fastq2)
		{

			`head -$num_lines_sample $fastq1 > $fastq1_out`;
			`head -$num_lines_sample $fastq2 > $fastq2_out`;

			return ($fastq1_out, $fastq2_out, $coverage);

		}else{

			$num_lines_sample = $num_lines_sample*2;
			`head -$num_lines_sample $fastq1 > $fastq1_out`;

			return ($fastq1_out, "NA", $coverage);

		}#end if

	}else{

		print "\n\nNot enough sequences for desired coverage, total coverage is $total_coverage\x\n\n";

		return ($fastq1, $fastq2, $total_coverage);

	}#end if

}#end downsample_fastq


#sub modify_fastq_read_names - Takes as input: 1) fastq		- The name of a fastq file
#					       2) pe_num	- Indicating which paired end read (1 or 2)
#			     - Creates streamlined read names
sub modify_fastq_read_names
{
	my ($fastq, $pe_num) = @_;

	my $fastq_temp = $fastq . "_temp";

	#GO THROUGH FIRST FILE                                                                                                                                 
        open FASTQ_IN, $fastq;                                                                                                                             
        open FASTQ_OUT, ">$fastq_temp";                                                                                                                        
                                                                                                                                                               
	my $count = 4;

        foreach my $line (<FASTQ_IN>)                                                                                                                          
        {                                                                                                                                                      
                                                                                                                                                               
		my $seq_count = int($count/4);

		if(($count % 2) == 0)
		{

	                $line =~ s/^([@+][^:\.]+)[:\.].*$/$1\.$seq_count\/$pe_num/;                                                                                                      
        
		}#end if
                                                                                                                                                       
                print FASTQ_OUT $line;                                                                                                                         
		
		$count ++;
                                                                                                                                                       
        }#end foreach                                                                                                                                          
        `mv $fastq_temp $fastq`;                                                                                                                                       
        close FASTQ_IN;                                                                                                                                        
        close FASTQ_OUT;


}#end modify_fastq_read_names


1;
