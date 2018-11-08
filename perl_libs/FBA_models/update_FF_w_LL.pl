
#ADD DIRECTORY CONTAINING DB_LIB TO PERL SEARCH PATH
BEGIN
{
        unshift (@INC, "/home/esnitkin/lib/perl");
}


use DBI;
use SBML_model_db_lib;
use IO_lib;
use file_lib;



#GET PARAMETERS FROM FILE
my $LL_param_file = shift @ARGV;
my $FF_param_file = shift @ARGV;

my $r_LL_params = &IO_lib::get_params_from_file($LL_param_file);
my %LL_params = %$r_LL_params;

my $r_FF_params = &IO_lib::get_params_from_file($FF_param_file);
my %FF_params = %$r_FF_params;


#DEFINE SOME PARAMETERS
my $LL_dbh = &SBML_model_db_lib::SBML_model_DB_connect( $LL_params{'data_source'} );
my $FF_dbh = &SBML_model_db_lib::SBML_model_DB_connect( $FF_params{'data_source'} );


#GET REACTION TO PROCESS RELATIONSHIPS FOR BOTH MODELS
my $r_LL_rxn2proc = &SBML_model_db_lib::get_rxns_procs_and_EC($LL_dbh);
my $r_FF_rxn2proc = &SBML_model_db_lib::get_rxns_procs_and_EC($FF_dbh);

my %LL_rxn2proc = %$r_LL_rxn2proc;
my %FF_rxn2proc = %$r_FF_rxn2proc;


#CHANGE IFF708 PROCESS ASSIGNMENTS WHERE APPLICABLE
foreach my $rxn (keys %FF_rxn2proc)
{

	if( exists( $LL_rxn2proc{$rxn} ) )
	{

		$FF_rxn2proc{$rxn}{'PROC'} = $LL_rxn2proc{$rxn}{'PROC'};

	}#end if

}#end foreach


#UPDATE IFF708 PROCESS ASSIGNMENTS
#&SBML_model_db_lib::update_rxn_eq_procs(\%FF_rxn2proc , $FF_dbh);
&SBML_model_db_lib::update_rxn_procs($FF_dbh);
