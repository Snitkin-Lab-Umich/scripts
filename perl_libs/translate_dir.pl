#Takes as input a directory and applies a regular expression to all files in the directory


my $dir = shift @ARGV;

my @files = `ls $dir`;

foreach my $file (@files)
{

	chomp $file;	
	my $old_file = $dir . $file;


	$file =~ s/MOMA-LP/MOMALP/;
	my $new_file = $dir . $file;


	`mv $old_file $new_file`;


}
