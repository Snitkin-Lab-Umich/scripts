#Applies regex to names of file in directory
#
#Takes as input:
#		1) directory	- 
#		2) in_regex	-
#		3) out_regex	-
#
#Changes file names according to regex
#


#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl mv_file.pl directory in_regex out_regex\n\n";

	exit;

}#end if


#GET COMMAND LINE ARGUMENTS
my $directory = shift @ARGV;
my $in_regex = shift @ARGV;
my $out_regex = shift @ARGV;


#GET FILES
my @files = `ls $directory`;

#APPLY REGEX TO FILE NAMES
foreach my $file (@files)
{

	chomp $file;

	$new_file = $file;
	$new_file =~ s/$in_regex/$out_regex/;

	`mv $directory/$file $directory/$new_file`;
	print "mv $directory/$file $directory/$new_file\n";

}#end foreach
