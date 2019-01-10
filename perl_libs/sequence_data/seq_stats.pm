package seq_stats;

use Bio::Tools::SeqWords;
use hash_lib;
use arith_ops;
use seq_ops;
use Bio::Restriction::Enzyme;
use Bio::Restriction::Analysis;


#sub init_monoNucHash - Takes as input: 
#		      - Returns a hash whose keys are the nucleotides and values are zero
sub init_monoNucHash
{
	my %monoNucHash = (A=>0, T=>0, G=>0, C=>0, N=>0);


	return \%monoNucHash;

}#end init_momNucHash


#sub init_diNucHash - Takes as input: 
#		    - Returns a hash whose keys are the di-nucleotides and values are zero
sub init_diNucHash
{
	my %diNucHash;
	my @nucs = ('A', 'T', 'G', 'C');

	
	foreach my $nuc1 (@nucs)
	{

		foreach my $nuc2 (@nucs)
		{

			my $di_nuc = $nuc1 . $nuc2;

			$diNucHash{$di_nuc} = 0;

		}#end foreach


	}#end foreach


	return \%diNucHash;

}#end init_momNucHash


#sub get_GC - Takes as input: 1) r_monoNucHash - A hash whose keys are the nucleotides and values are counts
#	    - Returns the GC percentage based on the counts in the hash
sub get_GC
{
	my ($r_monoNucHash) = @_;

	my %monoNucHash = %$r_monoNucHash;

	my $GC = 0;

	if ( ($monoNucHash{'G'} + $monoNucHash{'C'} + $monoNucHash{'A'} + $monoNucHash{'T'}) > 0 )
	{

		$GC = ( ($monoNucHash{'G'} + $monoNucHash{'C'})  / ($monoNucHash{'G'} + $monoNucHash{'C'} + $monoNucHash{'A'} + $monoNucHash{'T'}) );
	
	}#end if

	return &arith_ops::rnd($GC,2);

}#end get_GC


#sub get_N - Takes as input: 1) seq - A nucleotide sequence
#	   - Returns the fraction of ambiguous bases in a sequence
sub get_N
{
	my ($seq) = @_;

	my $seq_obj = Bio::Seq->new(-seq => $seq);

        #GET NUCLEOTIDE COUNTS
        my $r_monoNucHash = &init_monoNucHash();
        $r_monoNucHash = &seq_stats::get_monoNucCount($r_monoNucHash, $seq_obj);

        my $frac_N = $r_monoNucHash->{"N"} / &hash_lib::hash_sum($r_monoNucHash);


	return $frac_N;

}#end get_N


#sub get_monoNucCount - Takes as input: 1) r_mono_nuc_hash	- A hash whose keys are nucleotides and values are counts
#					2) seq			- A sequence to count nucleotides in
#		      - Returns a hash whose keys are nucleotides and values are uptdated counts based on frequencies in the sequence
sub get_monoNucCount
{
	my ($r_mono_nuc_hash, $seq_obj) = @_;

	my %mono_nuc_hash = %$r_mono_nuc_hash;

	
	my $r_temp_mono_nuc_hash = Bio::Tools::SeqWords->count_words($seq_obj, 1);
	my %temp_mono_nuc_hash = %$r_temp_mono_nuc_hash;


	foreach my $nuc (keys %temp_mono_nuc_hash)
	{

		$mono_nuc_hash{$nuc} += $temp_mono_nuc_hash{$nuc};

	}#end foreach


	return \%mono_nuc_hash;

}#end get_monoNucCount


#sub get_diNucCount - Takes as input: 1) r_di_nuc_hash      - A hash whose keys are di-nucleotides and values are counts
#                                     2) seq                - A sequence to count di-nucleotides in
#                     - Returns a hash whose keys are di-nucleotides and values are uptdated counts based on frequencies in the sequence
sub get_diNucCount
{
        my ($r_di_nuc_hash, $seq_obj) = @_;

	my %di_nuc_hash = %$r_di_nuc_hash;

	my $r_temp_di_nuc_hash = Bio::Tools::SeqWords->count_overlap_words($seq_obj, 2);
	my %temp_di_nuc_hash = %$r_temp_di_nuc_hash;

	my @di_nucs = &seq_ops::get_di_nucs();


	foreach my $di_nuc (@di_nucs)
	{

		$di_nuc_hash{$di_nuc} += $temp_di_nuc_hash{$di_nuc};

	}#end foreach


	return \%di_nuc_hash;

}#end get_diNucCount



#sub get_diNucSignature - Takes as input - 1) r_mono_nuc_hash	- A hash of mono-nucleotide counts 
#					   2) r_di_nuc_hash	- A hash of di-nucleotide counts
#			- Returns a hash of dinucleotide signatures	
sub get_diNucSignature
{
	my ($r_mono_nuc_hash, $r_di_nuc_hash) = @_;

	my %mono_nuc_hash = %$r_mono_nuc_hash;
	my %di_nuc_hash = %$r_di_nuc_hash;
	my %di_nuc_sig_hash;


	#CALCULATE DI-NUCLEOTIDE SIGNALS
	my $mono_nuc_count = &hash_lib::hash_sum($r_mono_nuc_hash);
	my $di_nuc_count = &hash_lib::hash_sum($r_di_nuc_hash);

	my @di_nucs = &seq_ops::get_di_nucs();

	foreach my $di_nuc (@di_nucs)
	{

		#GET MONO-NUCLEOTIDES
		my ($nuc1, $nuc2) = split //, $di_nuc;
		
		#CALCULATE DI-NUCLEOTIDE SIGNATURE
		$di_nuc_sig_hash{$di_nuc} = 			( $di_nuc_hash{$di_nuc} / $di_nuc_count ) 
					   						/ 
					( ( $mono_nuc_hash{$nuc1}  / $mono_nuc_count ) * ( $mono_nuc_hash{$nuc2}  / $mono_nuc_count ) );
	

	}#end foreach


	return \%di_nuc_sig_hash;

}#end get_diNucSignature


#sub get_diNucSignatureDiff - Takes as input - 1) r_di_nuc_sig_hash1	- A hash of di-nucleotide signatures
#                                              2) r_di_nuc_sig_hash12	- A second hash of di-nucleotide signatures
#                           - Returns the mean difference of signatures across all di-nucleotides
sub get_diNucSignatureDiff
{
        my ($r_di_nuc_sig_hash1, $r_di_nuc_sig_hash2) = @_;

	my %di_nuc_sig_hash1 = %$r_di_nuc_sig_hash1;
	my %di_nuc_sig_hash2 = %$r_di_nuc_sig_hash2;


        #MEAN DIFFERENCE IN DI-NUCLEOTIDE SIGNATURES
	my $di_nuc_sig_diff = 0;
	my @di_nucs = &seq_ops::get_di_nucs();


        foreach my $di_nuc (@di_nucs)
        {

		$di_nuc_sig_diff += abs( ($di_nuc_sig_hash1{$di_nuc} - $di_nuc_sig_hash2{$di_nuc}) );
		#print "$di_nuc($di_nuc_sig_hash1{$di_nuc} , $di_nuc_sig_hash2{$di_nuc}) | ";

        }#end foreach

	#print "\n";

	#NORMALIZE BY NUMBER OF DI-NUCLEOTIDES
	$di_nuc_sig_diff = $di_nuc_sig_diff / 16;


        return $di_nuc_sig_diff;

}#end get_diNucSignature


#sub kmer_count - Takes as input: 1) seq		- A sequence
#				  2) kmer_length	- The length of kmers to count
#		- Returns a hash where the keys are 
sub kmer_count
{
	my ($seq, $kmer_length) = @_;

	my $seq_obj = Bio::Seq->new(-seq => $seq);

	my $r_kmer_counts = Bio::Tools::SeqWords->count_overlap_words($seq_obj, $kmer_length);

	return $r_kmer_counts;

}#end kmer_count


#sub sliding_window_kmer_count - Takes as input: 1) seq		- A DNA sequence
#						 2) kmer_length	- The length of kmers to count
#						 3) window_size	- The window size in which kmers should be counted
#						 4) gap_size	- The gap between consecutive windows
#			       - Returns a 2D hash where the first key is the start position and the second key is the kmer. The
#				 values are the counts of kmers in the given windows 
sub sliding_window_kmer_count
{
	my ($seq, $kmer_length, $window_size, $gap_size) = @_;

	my %kmer_counts;


	for(my $start = 0; $start < length($seq); $start+=$gap_size)
	{

        	my $sub_seq = substr($seq, $start, &arith_ops::min($window_size, (length($seq) - $start)) );

        	$kmer_counts{$start} = &kmer_count($sub_seq, $kmer_length);

	}#end for


	return \%kmer_counts;

}#end sliding_window_kmer_count


#sub get_novelRE_cuts - Takes as input: 1) seq	- A bioperl nucleotide sequence object
#				   2) motif	- A motif recognized by the restriction enzyme
#				   3) position	- A position relative to the motif where the cut is made
#		      - Returns a list of the indicies in the sequence where a cut is made
sub get_novelRE_cuts
{
	my ($seq, $motif, $position) = @_;

	#DEFINE A BIOPERL RESTRICTION ENZYME
	my $re= Bio::Restriction::Enzyme->new(-enzyme=>"MY_RE", -seq=>$motif, -cut=>$position);


	#MAKE RESTRICTION ANALYSIS OBJECT
	my $ra = Bio::Restriction::Analysis->new(-seq=>$seq, -enzymes=>$re);


	my @cut_inds = $ra->positions("MY_RE");
	#my @fragments = $ra->fragments("MY_RE");

	#print $seq->seq() . "\n@cut_inds\n@fragments\n";

	return \@cut_inds;

}#end get_novelRE_cuts


#sub get_RE_cuts - Takes as input: 1) seq	- A bioperl nucleotide sequence object
#				   2) re_name	- The name of a RE in REbase
#		 - Returns a list of the indicies in the sequence where a cut is made
sub get_RE_cuts
{
	my ($seq, $re_name) = @_;

	#MAKE RESTRICTION ANALYSIS OBJECT
	my $ra = Bio::Restriction::Analysis->new(-seq=>$seq);


	my @cut_inds = $ra->positions($re_name);
	#my @fragments = $ra->fragments("MY_RE");

	#print $seq->seq() . "\n@cut_inds\n@fragments\n";

	return \@cut_inds;

}#end get_RE_cuts


#sub get_4fold_degenerate_positions - Takes as input: 1) seq - A nucleotide sequence
#				    - Returns the positions of 4-fold degenerate sites in the sequence
sub get_4fold_degenerate_positions
{
	my ($seq) = @_;

	my @four_fold_pos = ();
	my %four_fold_codons;

	$four_fold_codons{'CT'} = 0;
	$four_fold_codons{'GT'} = 0;
	$four_fold_codons{'TC'} = 0;
	$four_fold_codons{'CC'} = 0;
	$four_fold_codons{'AC'} = 0;
	$four_fold_codons{'GC'} = 0;
	$four_fold_codons{'CG'} = 0;
	$four_fold_codons{'GG'} = 0;


	for(my $p = 0;$p < length($seq);$p+=3)
	{
		$codon = substr($seq, $p, 2);

		if(exists($four_fold_codons{$codon}))
		{

			push @four_fold_pos, ($p+3);

		}#end if

	}#end for


	return \@four_fold_pos;

}#end get_4fold_degenerate_positions


1;
