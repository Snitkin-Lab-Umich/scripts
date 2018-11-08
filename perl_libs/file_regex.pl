#Takes as input a file and outputs a file with lines changed with regex

my ($in_file, $out_file) = @ARGV;

open INPUT , $in_file;
open OUTPUT , ">$out_file";


foreach $line (<INPUT>)
{

	$line =~ s/\'//g;

	print OUTPUT $line;

}

close INPUT;
close OUTPUT;
