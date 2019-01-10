package bambus;


#sub parse_stats_file - Takes as inpput: 1) bambus_stats - A bambus file with stats of scafolding
#		      - Returns a hash where keys are stat names and values are the stats
sub parse_stats_file
{
	my ($bambus_stats) = @_;


	my %stats_hash;

	#OPEN FILE
	open IN, $bambus_stats;


	#GO THROUGH FILE EXTACTING INFO
	foreach my $line (<IN>)
	{

		chomp $line;

		if ($line =~ /\w+/)
		{

			my $stat_name = $line;
			$stat_name =~ s/^(.*)\:.*$/$1/;
			

			my $stat_value = $line;
			$stat_value =~ s/^.*\:\s+([\d\.\w]*)$/$1/;


			if($stat_value =~ /\d/)
			{

				$stats_hash{$stat_name} = $stat_value;

			}else{

				$stats_hash{$stat_name} = 0;

			}#end if


		}#emd if


	}#end foreach



	return %stats_hash;

}#end parse_stats_file


1;
