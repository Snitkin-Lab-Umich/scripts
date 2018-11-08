package arith_ops;

#sub min - Takes as input two numbers
#	 - Returns the minimum of the two numbers
sub min
{
	my($a , $b) = @_;

	if($a < $b)
	{

		return $a;

	}else{

		return $b;

	}#end if

}#end min


#sub max - Takes as input two numbers
#        - Returns the maximum of the two numbers
sub max
{
        my($a , $b) = @_;

        if($a > $b)
        {

                return $a;

        }else{

                return $b;

        }#end if

}#end max


#sub abs - Takes as input a number
#        - Returns the absolute value of that number
sub abs
{
        my ($num) = @_;


        if($num >= 0)
        {

        	return $num;

        }else{

                my $abs_num = $num * -1;
                return $abs_num;

        }#end if


}#end abs


#sub rnd - Takes as input a number and a desired precision
#        - Returns the same list with all numbers rounded to the desired precision
sub rnd
{
        my ($num , $prec) = @_;

        my $prec_statement = "%." . $prec . "f";


        my $rnd_num = sprintf($prec_statement , $num);


        return $rnd_num;

}#end rnd


1;
