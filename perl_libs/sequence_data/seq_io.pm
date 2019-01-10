package seq_io;

use Bio::SeqIO;


#sub fasta2hash - Takes as input: 1) fasta_file	- A fasta formatted file
#		- Returns a hash reference where the keys are the FASTA header
#		  and the values are the sequences
sub fasta2hash
{
	my ($fasta_file) = @_;

	my %fasta_hash;
	my @seq_names;

	my $seqio_object = Bio::SeqIO->new(-file => $fasta_file);

	while(my $seq_object = $seqio_object->next_seq)
	{

		$fasta_hash{$seq_object->display_id} = $seq_object->seq;
	
		push @seq_names, $seq_object->display_id;

	}#end while


	return (\%fasta_hash,\@seq_names);

}#end fasta2hash


#sub fastaQ2hash - Takes as input: 1) fasta_file        - A fasta formatted file of quality scores
#                - Returns a hash reference where the keys are the FASTA header
#                  and the values are space delimited quality scores
sub fastaQ2hash
{
        my ($fasta_file) = @_;

        my %fasta_hash;
        my @seq_names;

        my $seqio_object = Bio::SeqIO->new(-file => $fasta_file,
                                           -format => 'qual');
	print "$fasta_file\n\n";
        while(my $seq_object = $seqio_object->next_seq)
        {

                $fasta_hash{$seq_object->display_id} = join " ", @{$seq_object->qual};

		print $seq_object->display_id . "\n\n";
                push @seq_names, $seq_object->display_id;

        }#end while


        return (\%fasta_hash,\@seq_names);

}#end fastaQ2hash

#sub fastq2hash - Takes as input: 1) fastq_file        - A fastq formatted files
#                - Returns a hash reference where the keys are the FASTA header
#                  and the values are space delimited quality scores
sub fastq2hash
{
        my ($fastq_file) = @_;

        my %fastq_hash;
        my @seq_names;

        my $seqio_object = Bio::SeqIO->new(-file => $fastq_file,
                                           -format => 'fastq');

        while(my $seq_object = $seqio_object->next_seq)
        {

                $fastq_hash{$seq_object->display_id} = join " ", @{$seq_object->qual};

                push @seq_names, $seq_object->display_id;

        }#end while

        return (\%fastq_hash,\@seq_names);

}#end fastq2hash

#sub hash2fasta - Takes as input: 1) fasta_file 	- File to write to
#				  2) r_seq_hash		- A hash whose keys are sequence names and values are sequences
#		- Writes a fasta file
sub hash2fasta
{
	my ($fasta_file, $r_seq_hash) = @_;

	my %seq_hash = %$r_seq_hash;


	my $seq_out = Bio::SeqIO->new(-file   => ">$fasta_file",
                                     -format => "fasta");


	foreach my $header (keys %seq_hash)
	{

		my $seq_obj = Bio::Seq->new(-seq => $seq_hash{$header},                        
                                            -display_id => $header);
        
        	$seq_out->write_seq($seq_obj);

	}#end foreach


}#end hash2fasta


#sub fasta2bpobj - Takes as input: 1) fasta_file	- A fasta formatted file of sequences
#		 - Returns a reference to an array of bioperl sequence objects
sub fasta2bpobj
{
	my ($fasta_file) = @_;

	my @seq_objs;

	my $inseq = Bio::SeqIO->new(-file   => "<$fasta_file",
                                    -format => 'fasta' );

	 while (my $seq = $inseq->next_seq)
	{

		push @seq_objs, $seq;

	}#end while


	return \@seq_objs;

}#end fasta2bpobj


#sub fasta2bpobj_hash - Takes as input: 1) fasta_file	- A fasta formatted file of sequences
#		       - Returns a reference to a hash where the key is the sequence ID and the
#			 values is a corresponding  bioperl sequence object
sub fasta2bpobj_hash
{
	my ($fasta_file) = @_;

	my %seq_objs;

	my $inseq = Bio::SeqIO->new(-file   => "<$fasta_file",
                                    -format => 'fasta' );

	 while (my $seq = $inseq->next_seq)
	{

		$seq_objs{$seq->display_id} = $seq;

	}#end while


	return \%seq_objs;

}#end fasta2bpobj_hash


#sub write_fasta_file - Takes as input: 1) header	- The header for the fasta file
#					2) sequence	- The sequence to go in the file
#					3) fasta_file	- The name of the file
#		      - Outputs the desired fasta file
sub write_fasta_file
{
	my ($header, $sequence, $fasta_file) = @_;


	my $seq_out = Bio::SeqIO->new(-file   => ">$fasta_file",
                                     -format => "fasta");

	my $seq_obj = Bio::Seq->new(-seq => $sequence,                        
                          	-display_id => $header);
	
	$seq_out->write_seq($seq_obj);


}#end write_fasta_file


#sub append_fasta_file - Takes as input: 1) header       - The header for the fasta file
#                                       2) sequence     - The sequence to go in the file
#                                       3) fasta_file   - The name of the file
#                     - Outputs the desired fasta file
sub append_fasta_file
{
        my ($header, $sequence, $fasta_file) = @_;


        my $seq_out = Bio::SeqIO->new(-file   => ">>$fasta_file",
                                     -format => "fasta");

        my $seq_obj = Bio::Seq->new(-seq => $sequence,
                                -display_id => $header);

        $seq_out->write_seq($seq_obj);


}#end write_fasta_file


#sub translate_fasta - Takes as input: 1) fasta_n	- A fasta file of cds sequences
#				       2) fasta_p	- A file to hold translated sequences
#		     - Translates files in fasta_n and prints to fasta_p
sub translate_fasta
{
	my ($fasta_n, $fasta_p) = @_;

	#CREATE BIOPERL SEQIO OBJECTS
	my $seq_in = Bio::SeqIO->new('-file' => "<$fasta_n",
        	                     '-format' => 'fasta');

	my $seq_out = Bio::SeqIO->new('-file' => ">$fasta_p",
                	              '-format' => 'fasta');



	#GO THROUGH EACH SEQUENCE AND TRANSLATE
	while (my $inseq = $seq_in->next_seq)
	{

	        my $trans_seq = $inseq->translate;

	        $seq_out->write_seq($trans_seq);

	}#end while


}#end translate_fasta


#sub parse_genbank - Takes as input: 1) gb_file		- A genbank file
#		   - Returns a hash containing info for each gene in the file
sub parse_genbank
{
	my ($gb_file) = @_;

	my %annot_hash;
	my %seq_hash;


	#PARSE GENBANK FILE AND OUTPUT TO FILES
	my $seqio_object = Bio::SeqIO->new(-file => $gb_file );


	while(my $seq_object = $seqio_object->next_seq)
	{

	        #PRINT OUT FULL SEQUENCE
	        my $prim_seq = $seq_object->primary_seq->seq;
	        my $prim_seq_id = $seq_object->primary_seq->id;
	        my $prim_seq_desc = $seq_object->primary_seq->desc;

	
	        #PRINT OUT GENES AND ANNOTATIONS
	        for my $feat_object ($seq_object->get_SeqFeatures)
	        {

	                if ($feat_object->primary_tag eq "CDS")
	                {


	                        if($feat_object->has_tag("protein_id"))
	                        {


					#GET PERTINENT INFO FOR ENTRY
	                                my $seq = $feat_object->seq()->seq;
	
	                                my @protein_id = $feat_object->get_tag_values('protein_id');
					my $protein_id = $protein_id[0];
					$protein_id =~ s/^.*:(.*)$/$1/;
	
	                                my @annotation = $feat_object->get_tag_values('product');

	                                my $strand = $feat_object->location->strand;

	                                my $length = $feat_object->location->length;
				
					#DETERMINE STRAND, AND THEN NOTE START AND END
					my($start, $end);

					if($strand == 1)
					{

		                                $start = $feat_object->location->start;
		                                $end = $feat_object->location->end;

					}else{

		                                $end = $feat_object->location->start;
		                                $start = $feat_object->location->end;

					}#end if


					#STORE GENE INFO AND SEQUENCE
					$seq_hash{$protein_id} = $seq;
	
					if($feat_object->has_tag("note"))
					{
		
	                                	my @note = $feat_object->get_tag_values('note');

						$annot_hash{$protein_id}{'CONTIG'} = $prim_seq_id;
						$annot_hash{$protein_id}{'START'} = $start;
						$annot_hash{$protein_id}{'END'} = $end;
						$annot_hash{$protein_id}{'ANNOTATION'} = $annotation[0];

						$annot_hash{$protein_id}{'COG'} = $note[0];
						$annot_hash{$protein_id}{'COG'} =~ s/^(COG\d+) .*$/$1/;

						$annot_hash{$protein_id}{'COG_DESC'} = $note[0];
						$annot_hash{$protein_id}{'COG_DESC'} =~ s/^COG\d+ (.*)$/$1/;

					}else{

						$annot_hash{$protein_id}{'CONTIG'} = $prim_seq_id;
						$annot_hash{$protein_id}{'START'} = $start;
						$annot_hash{$protein_id}{'END'} = $end;
						$annot_hash{$protein_id}{'ANNOTATION'} = $annotation[0];
						$annot_hash{$protein_id}{'COG'} = "";
						$annot_hash{$protein_id}{'COG_DESC'} = "";

					}#end if

	                        }#end if
	
	                }#end if

	        }#end foreach

	}#end while


	return (\%seq_hash, \%annot_hash);

}#end parse_genbank


#sub fasta_file_subset - Takes as input: 1) fasta_file	- A fasta file
#					 2) start	- The first sequence to output
#					 3) end		- The last sequence to output
#		       - Creates a fasta file composed of the subset of sequences desired
sub fasta_file_subset
{
	my ($fasta_file, $start, $end) = @_;

	my $fasta_file_subset = $fasta_file . "_$start-$end";
	
	#CREATE INPUT AND OUTPUT STREAMS
	my $seqio_object = Bio::SeqIO->new(-file => $fasta_file);
	
	my $seq_out = Bio::SeqIO->new(-file   => ">$fasta_file_subset",
                                     -format => "fasta");

	#GO THROUGH FASTA AND OUTPUT SEQUENCES OF INTEREST
	my $seq_count = 0;

	while(my $seq_object = $seqio_object->next_seq)
	{

		#INCREMENT COUNTER AND DETERMINE IF SEQUENCE SHOULD BE OUTPUTED
		$seq_count ++;

		if($seq_count >= $start && $seq_count <= $end)
		{

                	$seq_out->write_seq($seq_object);

		}elsif($seq_count > $end){

			return $fasta_file_subset;

		}#end if

	}#end while

	return $fasta_file_subset;

}#end fasta_file_subset


1;
