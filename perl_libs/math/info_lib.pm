package info_lib;

#This file contains various subroutines relating to information theory


#$X = "1\t2\t3\t1\t2\t3\t1\t2\t3";
#$Y = "3\t2\t1\t3\t2\t1\t3\t2\t1";
#$Z = "1\t1\t1\t1\t1\t1\t1\t1\t1";
#&MI($X,$Y);
#&H($X);
#&H($Z);

########################################################################################

#sub MI - Takes as input two strings in the form of references to lists 
#	- Computes the mutual information between the two strings
#	 * MI(X;Y) = H(X) - H(X|Y) *
sub MI
{
	my ($X, $Y) = @_;

	my @X = @$X;
	my @Y = @$Y;

	my (%X_counts, %Y_counts, %XY_counts);

	#MAKE SURE X AND Y ARE THE SAME SIZE
	if(@X != @Y)
	{
		print "X : " . (@X) . "   Y : " . (@Y) . "\n"; 
		return "ERROR";
	}#end if

	#GO THROUGH THE TWO STRINGS AND GET SINGLE COUNTS AND PAIR COUNTS
	for(my $i=0; $i < @X; $i++)
	{
		#INCREMENT X COUNTS
		if(exists $X_counts{$X[$i]})
		{
			$X_counts{$X[$i]} ++;
		}else{
			$X_counts{$X[$i]} = 1;
		}

		#INCREMENT Y COUNTS
		if(exists $Y_counts{$Y[$i]})
		{
			$Y_counts{$Y[$i]} ++;
		}else{
			$Y_counts{$Y[$i]} = 1;
		} 
		 
		#INCREMENT XY COUNTS
		if(exists $XY_counts{"$X[$i],$Y[$i]"})
		{
			$XY_counts{"$X[$i],$Y[$i]"} ++;
		}else{
			$XY_counts{"$X[$i],$Y[$i]"} = 1;
		} 
	}#end for
	
	my $MI = 0;

	foreach my $XYkey (keys %XY_counts)
	{
		#RETRIEVE INDICIES
		$XYkey =~ /(\d+),(\d+)/;
		my $Xkey = $1;
		my $Ykey = $2;

		#COMPUTE MARGINAL AND JOINT PROBABILITIES
		my $pX = $X_counts{$Xkey} / @X;
		my $pY = $Y_counts{$Ykey} / @Y;
		my $pXY = $XY_counts{$XYkey} / @Y;

#		print "p(X) = $pX   p(Y) = $pY  p(XY) = $pXY\n";
#		print "XY:$XYkey\n";
#		print "X:$Xkey   Y:$Ykey\n\n";
		
		$MI += $pXY * (log($pXY/($pX*$pY)) / log(2));

	}#end foreach

	#print $MI . "\n";

	return $MI;

}#end MI


#sub multi_MI - Takes as input a reference to a list of N references to string lists
#	      - Returns the multi-information among the strings
sub multi_MI
{
	my ($r_strings) = @_;
	my @strings = @$r_strings;
	my %marginal_counts;
	my %joint_counts;
	my $string_length = scalar(@{$strings[0]});
	my $multi_I = 0;

        #GO THROUGH THE STRINGS AND GET SINGLE COUNTS AND MULTI COUNTS
        for(my $i=0; $i < $string_length; $i++)
        {

		#INCREMENT COUNTS FOR EACH INDIVIDUAL STRING AS WELL AS ACCUMATE THE JOIUNT STRING
		my @joint;
		for(my $j = 0; $j < (@strings); $j++)
		{
	
			if( exists($marginal_counts{$j}{ ((@{$strings[$j]})[$i]) }) )
			{	

				$marginal_counts{$j}{ ((@{$strings[$j]})[$i]) } ++;		

			}else{

				$marginal_counts{$j}{ ((@{$strings[$j]})[$i]) } = 1;

			}#end if

			push @joint , ((@{$strings[$j]})[$i]);
		}#end for


		#INCREMENT THE JOINT COUNTS
		my $joint = join "\t" , @joint;
		
		if( exists($joint_counts{$joint}) )
		{
			
			$joint_counts{$joint} ++;

		}else{

			$joint_counts{$joint} = 1;	
	
		}


        }#end for


	#CALCULATE THE JOINT AND MARGINAL PROBABILITIES
	foreach my $joint( keys %joint_counts)
	{

		#CALCULATE THE JOINT PROBABILITY
		my $P_joint = $joint_counts{$joint} / scalar(@{$strings[0]});

		
		#CALCULATE THE PRODUCT OF THE MARGINAL PROBABILTIES
		my $P_marginal = 1;

		my @joint = split /\t/ , $joint;

		for(my $i = 0; $i < (@joint); $i++)
		{

			$P_marginal = $P_marginal * ($marginal_counts{$i}{$joint[$i]} / $string_length);

		}#end for


		$multi_I += $P_joint * ( log( $P_joint/ $P_marginal) / log(2) );

	}#end foreach


	return $multi_I;

}#end multi_MI


#sub H - Takes as input a numerical tab delimited string
#      - Returns the entropy of the string
#	 * H(X) = sum X=x (p(x) log p(x)
sub H
{
	my ($X) = @_;
	my @X = @$X;	

	my %X_counts;

	#GO THROUGH THE TWO STRINGS AND GET SINGLE COUNTS AND PAIR COUNTS
	for(my $i=0; $i < @X; $i++)
	{
		#INCREMENT X COUNTS
		if(exists $X_counts{$X[$i]})
		{
			$X_counts{$X[$i]} ++;
		}else{
			$X_counts{$X[$i]} = 1;
		}
	}#end for

	my $H = 0;

	foreach my $Xkey (keys %X_counts)
	{
		my $pX = $X_counts{$Xkey} / @X;
#		print "X:$Xkey   p(X):$pX\n";
		$H += -($pX * (log($pX) / log(2)));
	}#end foreach

#	print "ENTROPY = $H\n\n";

	return $H;
}#end H


#sub profile_diff - Takes as input two nmerical tab delimited strings
#		  - Returns the sum of the absolute difference at each position
sub profile_diff
{
	my ($X , $Y) = @_;
	my $diff = 0;

	my @X = @$X;
	my @Y = @$Y;


	for(my $i=0; $i < (@X); $i++)
	{

		$diff += abs($X[$i] - $Y[$i]);

	}#end for


	return $diff;

}#end profile_diff

1;
