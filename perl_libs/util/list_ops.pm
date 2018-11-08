package list_ops;


#sub sum - Takes as input a list reference
#	 - Returns the sum of the values in the list
sub sum
{
	my ($r_list) = @_;
	my $sum = 0;


	foreach my $value (@$r_list)
	{

		$sum += $value;

	}#end forach


	return $sum;

}#end sum


#sub mean - Takes as input a list reference
#	  - Returns the mean of the list
sub mean
{
	my ($r_list) = @_;
	my ($sum,$mean);
	

	#COMPUTE THE SUM OF THE LIST
	$sum = &sum($r_list);


	#COMPUTE THE MEAN AND RETURN IT
	$mean = $sum / scalar(@$r_list);


	return $mean;

}#end mean


#sub median - Takes as input a list reference
#	    - Returns the mean of the list
sub median
{
	my ($r_list) = @_;

	#SORT THE LIST
	@sorted_list = sort {$a<=>$b} @$r_list;


	#COMPUTE THE MEAN AND RETURN IT
	my $median;

	if((@sorted_list % 2) == 1)
	{

		$median = $sorted_list[(@sorted_list - 1)/ 2];

	}else{

		$median = ($sorted_list[(@sorted_list /2) -1] + $sorted_list[@sorted_list / 2]) / 2;

	}#end if

	return $median;

}#end median


#sub min - Takes as input a list reference
#	 - Returns the minimum value in the list
sub min
{
	my ($r_list) = @_;
	my $min = @$r_list[0];

	foreach $value (@$r_list)
	{

		if( $value < $min )
		{

			$min= $value;

		}#end if

	}#end foreach

	return $min;

}#end min


#sub max - Takes as input a list reference
#        - Returns the maximum value in the list
sub max
{
        my ($r_list) = @_;
	my @list = @$r_list;
        my $max = $list[0];

        foreach $value (@list)
        {

                if( $value > $max )
                {

                        $max = $value;

                }#end if

        }#end foreach

	return $max;

}#end max


#sub freq - Takes as input a list reference
#	  - Returns two list references. On contains the unique vallues in the list
#	    the second contains the number of time the corresponding value occurs
sub freq
{
	my ($r_list) = @_;
	my @freq;
	my @values;

	#SORT THE LIST
	my @s_list  = sort {$a <=> $b} @$r_list;


	#TRAVERSE THE LIST KEEPING TRACK OF VALUES AND THEIR FREQUENCIES
	my $old_value;
	my $value_index = -1;

	foreach my $value (@s_list)
	{


		if($value ne $old_value)
		{

			$value_index ++;	
			$freq[$value_index] = 1;
			$values[$value_index]  = $value;
			$old_value = $value;

		}else{
			
			$freq[$value_index] ++;

		}#end if		

	}#end foreach


	return (\@values , \@freq);

}#end freq



#sub cum_sum - Takes as input a list reference
#	     - Returns a reference to a list where each position contains the sum of all prior
#	       positions in the original list
sub cum_sum
{
	my ($r_list) = @_;
	my @cum_sum_list;
	my $cum_sum = 0;


	for(my $i = 0; $i < (@$r_list); $i++)
	{

		$cum_sum += @$r_list[$i];
		$cum_sum_list[$i] = $cum_sum;

	}#end for


	return \@cum_sum_list;

}#end cum_sum



#sub apply - Takes as input a list reference and a function
#	   - Returns a reference to a list with the function applied to each position
#	     on the list
sub apply
{


}#end apply


#sub abs - Takes as input a list reference
#	 - Returns the same list, with the absolute value of every element taken
sub abs
{
	my ($r_list) = @_;

	my @list = @$r_list;
	my @abs_list;

	foreach my $num ( @list ) 
	{


		if($num >= 0)
		{

			push @abs_list , $num;

		}else{

			my $abs_num = $num * -1;
			push @abs_list , $abs_num;

		}#end if


	}#end foreach


	return \@abs_list;

}#end abs


#sub diff - Takes as input two list references
#	  - Returns a reference to a new list which contains the differences of the numbers
#	    at corresponding positions in the two lists
sub diff
{
	my ($r_list1 , $r_list2) = @_;

	my @list_diff;
	my @list1 = @$r_list1;
	my @list2 = @$r_list2;


	for(my $i = 0; $i < (@list1); $i++)
	{

		$list_diff[$i] = $list2[$i] - $list1[$i];
		
	}#end for


	return \@list_diff;

}#end diff


#sub rnd - Takes as input a list reference and a desired precision
#	 - Returns the same list with all numbers rounded to the desired precision
sub rnd
{
	my ($r_list , $prec) = @_;
	
	my @rnd_list;
	my $prec_statement = "%." . $prec . "f";


	foreach my $ele (@$r_list)
	{

		my $rnd_ele = sprintf($prec_statement , $ele);

		push @rnd_list , $rnd_ele;


	}#end foreach


	return \@rnd_list;

}#end rnd


#sub list_perc_total - Takes as input: 1) r_list - A reference to a list of numbers
#		     - Returns a list of the same length as the input list, where the values
#		       are the percentage of the list sum that each value in the original list
#		       constitutes.
sub list_perc_total
{
	my ($r_list) = @_;

	my @list_perc;


	#GET THE TOTAL SUM OF VALUES IN THE LIST
	my $list_sum = &sum($r_list);


	#CREATE NEW LIST WHICH IS PERCENTAGE OF TOTAL
	foreach $ele (@$r_list)
	{

		push @list_perc, (($ele / $list_sum) * 100);

	}#end foreach


	return \@list_perc;

}#end 


1;
