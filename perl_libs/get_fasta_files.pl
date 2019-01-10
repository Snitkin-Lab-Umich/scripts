my @dirs = `ls -l | grep "^d"`;
`touch temp`; 

foreach my $dir (@dirs)
{
	chomp $dir;
	my @line = split /\s+/, $dir;

	#`cp $line[-1]/*.pfasta temp/`;
	print "cat $line[-1]/*.fa >> temp\n";
	`cat $line[-1]/*.fa >> temp\n`;
	`echo "" >> temp`;

}
