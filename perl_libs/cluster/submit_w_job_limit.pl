#BREAKS UP A SWARM FILE INTO A SET OF SWARM FILES OF DESIRED SIZ, AND SUBMITS IN SERIAL
#
#Takes as input:
#		1) swarm_file	- The name of a swarm file
#		2) num_jobs	- The number of jobs to run at a time

#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl submit_w_job_limit.pl swarm_file num_jobs\n\n";

	exit;

}#end if

use submit_jobs;

#GET COMMAND LINE VARIABLES
my $swarm_file = shift @ARGV;
my $num_jobs = shift @ARGV;

my $temp_swarm_file = $swarm_file . "_temp";


#OPEN SWARM FILE AND WRITE LINES TO NEW SWARM FILE, IN GROUPS OF 10
open SWARM_ORIG, $swarm_file;
open SWARM_SUB, ">$temp_swarm_file";

my $job_counter = 0;

foreach my $line (<SWARM_ORIG>)
{

	#IF FILE HAS MAX NUMBER OF ENTRIES THEN SUBMIT TO CLUSTER AND WAIT
	if($job_counter == $num_jobs)
	{

		&submit_jobs::submit_and_wait($temp_swarm_file, 1);

		close SWARM_SUB;
		open SWARM_SUB, ">$temp_swarm_file";

		$job_counter = 0;

	}#end if


	#WRITE CURRENT JOB TO FILE
	print SWARM_SUB $line;

	$job_counter ++;

}#end foreach

&submit_jobs::submit_and_wait($temp_swarm_file, 1);

close SWARM_SUB;

`rm $temp_swarm_file`;
