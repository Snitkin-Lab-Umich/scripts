package visML_coords_lib;

#This package contains subroutines relevant to the generation and parsing of VisML files

use DBI;
use SBML_model_db_lib;
use SBML_model_lib;
use list_ops;
use visML_coords_db_lib;


#sub parse_rxns - Takes as input a database handle for an SBML model DB
#               - Returns two 2D hashes. The first has key1 = rxn, key2 = metab and value = R/P/B. The
#                 second has key1 = metab, key2 = rxn, and value = R/P/B.
sub parse_rxns
{
        my ($dbh) = @_;

        my $r_rxn_hash;

        my $r_rxn2metab_hash;
        my $r_metab2rxn_hash;


        #GET THE REACTIONS FROM THE DATABASE
        $r_rxn_hash = &SBML_model_db_lib::get_rxn_equations_from_SBML_DB($dbh);


        #PARSE OUT REACTIONS TO EXTRACT WHICH METABOLITES PARTICIPATE IN WHICH REACTIONS AND VICA VERSA
        ($r_rxn2metab_hash , $r_metab2rxn_hash) = &SBML_model_lib::parse_rxn_participation($r_rxn_hash);


        return ($r_rxn2metab_hash , $r_metab2rxn_hash);

}#emd parse_rxns


#sub get_metab_info - Takes as input: 1) dbh    - An SBML model database handle
#                   - Returns a hash linking each metabolite with corresponding info, which will be
#                     displayed as the node description.
sub get_metab_info
{
        my ($dbh) = @_;

        my $r_metab_info;


        return $r_metab_info;

}#end get_metab_info


#sub get_rxn_info - Takes is input: 1) dbh      - An SBML model database handle
#                 - Returns a hash linking each reaction with corresponding info, which will be
#                   displayed as the node description.
sub get_rxn_info
{
        my ($dbh) = @_;

        my $r_rxn_info;


        #GET THE REACTIONS FROM THE DATABASE
        $r_rxn_info = &SBML_model_db_lib::get_rxn_equations_from_SBML_DB($dbh);


        return ($r_rxn_info);

}#end get_rxn_info


#sub get_rxn_proc_info - Takes as input: 1) dbh         -  a database handle
#                                        2) part_procs  - A list of processes which will be displayed in the network. Note
#                                                         that if list is empty, then all processes are assumed to participate
#                      - Returns a list of processes and a hash linking reactions to the processes they
#                        participate in.
sub get_rxn_proc_info
{
        my ($dbh , $r_part_procs) = @_;

        my %part_procs = %$r_part_procs;

        my $r_proc_list;
        my %proc_hash;

        my $r_rxn2proc;
        my %rxn2proc;

	my %procInd_hash;


        foreach (keys %part_procs)
        {

                print $_ . "\n";
        }


        #GET A LIST OF PROCESSES
        $r_proc_list = &SBML_model_db_lib::get_processes($dbh);


        #GET HASH WHICH LINKS REACTIONS TO THE PROCESSES THEY PARTICIPATE IN
        $r_rxn2proc = &SBML_model_db_lib::get_rxns_and_processes($dbh);
        %rxn2proc = %$r_rxn2proc;

        #ASSIGN VALUES TO HASH BASED ON WHETHER OR NOT PROCESS WILL BE DISPLAYED
        if ( scalar( keys(%part_procs) ) == 0 )
        {

                foreach my $proc (@$r_proc_list)
                {

                        $proc_hash{$proc} = "V";

                }#end foreach

        }else{


                foreach my $proc (@$r_proc_list)
                {
	

                        if( exists( $part_procs{$proc} ) )
                        {

                                $proc_hash{$proc} = "V";

                        }else{

                                $proc_hash{$proc} = "I";

                        }#end if


                }#end foreach


        }#end if


        return (\%proc_hash , \%rxn2proc);

}#end get_rxn_proc_info



#sub get_rxn_coords - Takes as input: 1) dbh            - A visML database handle
#                                     2) layout         - The name of the node layout desired
#                                     3) model          - THe name of the model for which to get the node layout
#                                     4) node_names     - THe names of all nodes in the network
#                   - Returns a 2D hash where the first key is the reaction and the second are X/Y, with
#                     the values being the corresponding X,Y coordinates.
sub get_rxn_coords
{
        my ($dbh, $model, $layout, $r_node_names) = @_;

        my $r_node_coords;
        my %node_coords;

        my @node_names = @$r_node_names;

        my $free_node = 1;


        #GET COORDINATES FROM THE DATABASE
        $r_node_coords = &visML_coords_db_lib::get_visML_coords($dbh , $layout, $model);
        %node_coords = %$r_node_coords;


        #SET COORDINATES FOR ANY NODE WHICH WAS NOT SPECIFIED IN THE DATABASE LAYOUT
        foreach my $node ( @node_names )
        {

                if( !( exists( $node_coords{$node} ) ) )
                {

                        $node_coords{$node} -> {"X"} = $free_node;
                        $node_coords{$node} -> {"Y"} = $free_node;


                        $free_node ++;

                }#end if

        }#end foreach


        return \%node_coords;

}#end get_rxn_coords

1;
