package submit_jobs;

#sub submit_and_wait - Takes as input: 1) swarm_file	- A file of commands to submit to the cluster
#				       2) bundle	- The number of jobs to run sequentially on a single core
# 		     - Submits jobs to cluster and returns only after all jobs are done
sub submit_and_wait
{
	my ($swarm_file, $bundle) = @_;

	#SUBMIT JOBS AND GET IDS
	#my $jobs = `swarm -f $swarm_file -V -b $bundle -R gpfs`;
	my $jobs = `swarm -f $swarm_file -V -g 4 -R gpfs`;
	#my $jobs = `swarm -f $swarm_file -V`;

	my @jobs = split /\n/, $jobs;


	#WAIT UNTIL ALL JOBS ARE FINISHED TO PROCEED
	my $not_done = 1;

	while($not_done > 0)
	{
		
		$not_done = 0;

		#GO THROUGH JOBS AND DETERMINE IF ANY ARE STILL RUNNING
		foreach my $job (@jobs)
		{

			my $status = `qstat $job`;

			if($status =~ /$job/)
			{

				$not_done ++;

			}#end if

		}#end foreach

		sleep(60);

	}#end while

}#end submit_and_wait


#sub submit_and_wait_mem - Takes as input: 1) swarm_file	- A file of commands to submit to the cluster
#				           2) bundle		- The number of jobs to run sequentially on a single core
#					   3) mem		- Amount of memory required per node
# 		     	 - Submits jobs to cluster and returns only after all jobs are done
sub submit_and_wait_mem
{
	my ($swarm_file, $bundle, $mem) = @_;

	#SUBMIT JOBS AND GET IDS
	my $jobs = `swarm -f $swarm_file -V -b $bundle -g $mem`;
	#my $jobs = `swarm -f $swarm_file -V`;

	my @jobs = split /\n/, $jobs;


	#WAIT UNTIL ALL JOBS ARE FINISHED TO PROCEED
	my $not_done = 1;

	while($not_done > 0)
	{
		
		$not_done = 0;

		#GO THROUGH JOBS AND DETERMINE IF ANY ARE STILL RUNNING
		foreach my $job (@jobs)
		{

			my $status = `qstat $job`;

			if($status =~ /$job/)
			{

				$not_done ++;

			}#end if

		}#end foreach

		sleep(60);

	}#end while

}#end submit_and_wait_mem
1;
