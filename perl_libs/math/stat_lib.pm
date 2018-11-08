package stat_lib;

#use Statistics::Basic::Correlation; 


#sub pearson - Takes as input to list references 
#	     - Returns the pearson's correlation between the two lists
sub pearson
{
	my ($r_list1 , $r_list2) = @_;


	#my $co = new Statistics::Basic::Correlation( $r_list1, $r_list2 );
	#my $corr = $co->query; 


	#return $corr;

}#end pearson


#sub sensitivity - Takes as input: 1) TP	- The number of true positives
#				   2) TN	- The number of true negatives
#				   3) FP	- THe number of false positives
#				   4) FN	- THe number of false negatives
#Returns the sensitivity calculation
sub sensitivity
{
	my ($TP, $TN, $FP, $FN) = @_;


	if( ($TP + $FN) > 0 )
	{

		return ($TP / ($TP + $FN));

	}else{

		return 0;

	}#end if

}#end sensitivity


#sub specificity - Takes as input: 1) TP        - The number of true positives
#                                  2) TN        - The number of true negatives
#                                  3) FP        - THe number of false positives
#                                  4) FN        - THe number of false negatives
#Returns the specificity calculation
sub specificity
{
	my ($TP, $TN, $FP, $FN) = @_;

        if( ($TN + $FP) > 0 )
	{

		return ($TN / ($TN + $FP));

        }else{

                return 0;

        }#end if


}#end specificity


#sub accuracy - Takes as input: 1) TP        - The number of true positives
#                               2) TN        - The number of true negatives
#                               3) FP        - THe number of false positives
#                               4) FN        - THe number of false negatives
#Returns the accuracy calculation
sub accuracy
{
	my ($TP, $TN, $FP, $FN) = @_;

        if( ($TP +$TN + $FP + $FN) > 0 )
	{

		return (($TP + $TN) / ($TP +$TN + $FP + $FN));

        }else{

                return 0;

        }#end if


}#end accuracy


#sub NPV - Takes as input: 1) TP        - The number of true positives
#                          2) TN        - The number of true negatives
#                          3) FP        - THe number of false positives
#                          4) FN        - THe number of false negatives
#Returns the negative predicitive value calculation
sub NPV
{
	my ($TP, $TN, $FP, $FN) = @_;

        if( ($TN + $FN) > 0 )
	{
        
		return ($TN / ($TN + $FN))

        }else{

                return 0;

        }#end if


}#end NPV



1;
