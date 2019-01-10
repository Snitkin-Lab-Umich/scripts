package FBA_model_lib;

#This package contains various methods for dealing with data from an SBML DB lib


###################################### FUNCTIONS TO CREATE GENERIC GPR FORMATED FILES DATA##################################
#sub GPR_to_file - Takes as input- 1) GPR_hash  - A hash containing reaction information
#                                  2) GPR_file  - A file to write the data in the hash to
sub GPR_to_file
{
        my ($r_GPR_hash , $GPR_file) = @_;
        my %GPR_hash = %$r_GPR_hash;

        open OUTPUT , ">$GPR_file";

        foreach my $rxn_code (keys %GPR_hash)
        {

                print OUTPUT $rxn_code . "\t" . $rxn_code . "\t" . $GPR_hash{$rxn_code} -> {'RXN'} . "\t" .
                             $GPR_hash{$rxn_code} -> {'GENE'} . "\t" . $GPR_hash{$rxn_code} -> {'PROTEIN'} . "\t" .
                             $GPR_hash{$rxn_code} -> {'CLASS'} . "\t" . $GPR_hash{$rxn_code} -> {'SYSTEM'} . "\n";

        }#end


}#end GPR_to_file


#sub PG_to_file - Takes as input- 1) PG_hash  - A hash containing protein to gene informations
#                                 2) PG_file  - A file to write the data in the hash to

sub PG_to_file
{
        my ($r_PG_hash , $PG_file) = @_;
        my %PG_hash = %$r_PG_hash;


        open OUTPUT , ">$PG_file";

        foreach my $protein (keys %PG_hash)
        {

                print OUTPUT $protein;

                foreach my $gene ( @{$PG_hash{$protein}} )
                {

                        print OUTPUT "\t" . $gene;

                }#end foreach


                print OUTPUT "\n";

        }#end


        close OUTPUT;


}#end PG_to_file


#sub RE_to_file - Takes as input- 1) RE_hash  - A hash containing reaction and protein complex information
#                                 2) RE_file  - A file to write the data in the hash to
sub RE_to_file
{
        my ($r_RE_hash , $RE_file) = @_;
        my %RE_hash = %$r_RE_hash;


        open OUTPUT , ">$RE_file";

        foreach my $rxn (keys %RE_hash)
        {


                print OUTPUT $rxn;


                foreach my $expr ( @{$RE_hash{$rxn}} )
                {

                        print OUTPUT "\t" . $expr;

                }#end foreach


                print OUTPUT "\n";

        }#end


        close OUTPUT;


}#end RE_to_file


#sub RP_to_file - Takes as input- 1) RP_hash  - A hash containing reactions and information on the proteins catalyzing them
#                                 2) RP_file  - A file to write the data in the hash to
sub RP_to_file
{
        my ($r_RP_hash , $RP_file) = @_;
        my %RP_hash = %$r_RP_hash;


        open OUTPUT , ">$RP_file";

        foreach my $rxn (keys %RP_hash)
        {

                print OUTPUT $rxn;

                foreach my $protein ( @{$RP_hash{$rxn}} )
                {

                        print OUTPUT "\t" . $protein;

                }#end foreach


                print OUTPUT "\n";

        }#end


        close OUTPUT;


}#end RP_to_file



################################################## FUNCTIONS TO PARSE DATABASE DATA#########################################

#sub get_metabolites - Takes as input: 1) rxn_hash - A hash where the keys are reactions codes and the values are reactions
#                    - Returns - A 2 dimensional hash in which the first key is a metbaolite, the second key is a reaction
#                      the metabolite participates in and the value is the signed coefficent of the metabolite in the
#                      reaction
sub get_metabolites
{
        my ($r_rxn_hash) = @_;

        my %rxn_hash = %$r_rxn_hash;

        my %metab2rxn_hash;


        foreach my $rxn_code ( keys %rxn_hash )
        {

                my $rxn = $rxn_hash{$rxn_code};
                chomp $rxn;


                my $reactants;
                my $products;


                #CHECK IF REVERSIBLE AND THEN SPLIT REACTION INTO REACTANCTS AND PRODUCTS
                if( $rxn =~ /<[-|=]+>/ )
                {

                        ($reactants , $products) = split /<[-|=]+>/ , $rxn;

                }else{


                        ($reactants , $products) = split /[-|=]+>/ , $rxn;

                }#end if

                #PARSE REACTANTS, GETTING ALL METABOLITES AND PLACING IN HASH
                my @reactants = split /\s*\+\s*/ , $reactants;

                foreach my $reactant (@reactants)
                {


                        #REMOVE TERMINAL WHITE SPACE
                        $reactant =~ s/\s+$//;
                        $reactant =~ s/^\s+//;


			#IF REACTANT IS PLACE HOLDER THEN IGNORE
			if( $reactant =~ /\*/ )
			{

				next;

			}#end if

                        #CHECK IF METABOLITE HAS COEFFICIENT, IF SO STORE, OTHERWISE STORE -1
                        if( $reactant =~ /^\(*([0-9|\.|e|-]+)\)*\s/ )
                        {

                                my $coeff = $reactant;

                                $coeff =~ s/^\(*([0-9|\.|e|-]+)\)*\s+[\w\_\-]+/$1/;

                                $reactant =~ s/^\(*([0-9|\.|e|-]+)\)*\s+([\w\_\-]+)$/$2/;

                                $metab2rxn_hash{$reactant} -> {$rxn_code} = -1 * $coeff;

                        }else{

                                $metab2rxn_hash{$reactant} -> {$rxn_code} = -1;


                        }#end if


                }#end foreach


                #PARSE PRODUCTS, GETTING ALL METABOLITES AND PLACING IN HASH
                my @products = split /\s+\+\s+/ , $products;

                foreach my $product (@products)
                {

                        #REMOVE TERMINAL WHITE SPACE
                        $product =~ s/\s+$//;
                        $product =~ s/^\s+//;


			#IF REACTANT IS PLACE HOLDER THEN IGNORE
			if( $product =~ /\*/ )
			{

				next;

			}#end if


                        #CHECK IF METABOLITE HAS COEFFICIENT, IF SO STORE, OTHERWISE STORE -1
                        if( $product =~ /^\(*([0-9|\.|e|-]+)\)*\s/ )
                        {


                                my $coeff = $product;
                                $coeff =~ s/^\(*([0-9|\.|e|-]+)\)*\s+[\w\_\-]+/$1/;

                                $product =~  s/^\(*([0-9|\.|e|-]+)\)*\s+([\w\_\-]+)$/$2/;

                                $metab2rxn_hash{$product} -> {$rxn_code} = $coeff;

                        }else{


                                $metab2rxn_hash{$product} -> {$rxn_code} = 1;


                        }#end if


                }#end foreach



        }#end foreach


        return \%metab2rxn_hash;


}#end get_metabolites


#sub parse_rxn_participation - Takes as input: 1) rxn_hash - A hash linking reaction codes to reaction equations
#			     - Returns two 2D hashes. The first hash reaction codes as the first key, metabolite
#			       codes as the second key and the value is R,P,B depending upon whether the metabolite
#			       is a reactant a product or both in the given reaction. THe second hash is the same
#			       format, except the first key is the metabolite and the second is the reaction
sub parse_rxn_participation
{
	my ($r_rxn_hash) = @_;

	my %rxn_hash = %$r_rxn_hash;

	my %parsed_rxn_hash;
	my %parsed_metab_hash;


	#FOR EACH REACTION PARSE AND SAVE RELEVANT INFORMATION
        foreach my $rxn_code ( keys %rxn_hash )
        {

                my $rxn = $rxn_hash{$rxn_code};
                chomp $rxn;


                my $reactants;
                my $products;

                my @reactants;
                my @products;

		my $is_rev = 0;


                #CHECK IF REVERSIBLE AND THEN SPLIT REACTION INTO REACTANCTS AND PRODUCTS
                if( $rxn =~ /<[-|=]+>/ )
                {

                        ($reactants , $products) = split /<[-|=]+>/ , $rxn;
			$is_rev = 1;
			
                }else{


                        ($reactants , $products) = split /[-|=]+>/ , $rxn;

                }#end if


                #PARSE REACTANTS, GETTING ALL METABOLITES AND PLACING IN HASH
                @reactants = split /\s*\+\s*/ , $reactants;

                foreach my $reactant (@reactants)
                {

                        #REMOVE TERMINAL WHITE SPACE
                        $reactant =~ s/\s+$//;
                        $reactant =~ s/^\s+//;


			#IF REACTANT IS PLACE HOLDER THEN IGNORE
			if( $reactant =~ /\*/ )
			{

				next;

			}#end if

                        #CHECK IF METABOLITE HAS COEFFICIENT, AND REMOVE IT
                        if( $reactant =~ /^\(*([0-9|\.|e|-]+)\)*\s/ )
                        {

                                $reactant =~ s/^\(*([0-9|\.|e|-]+)\)*\s+([\w\_\-]+)$/$2/;

                        }#end if


			#STORE ENTRY IN BOTH HASHES
			if( $is_rev == 1 )
			{

				#$parsed_rxn_hash{$rxn_code} -> {$reactant} = "B";
				#$parsed_rxn_hash{$rxn_code} -> {$reactant} = "B";
				$parsed_rxn_hash{$rxn_code} -> {$reactant} = "BR";
				$parsed_metab_hash{$reactant} -> {$rxn_code} = "BR";

			}else{

				$parsed_rxn_hash{$rxn_code} -> {$reactant} = "R";
				$parsed_metab_hash{$reactant} -> {$rxn_code} = "R";

			}#end if


                }#end foreach


                #PARSE PRODUCTS, GETTING ALL METABOLITES AND PLACING IN HASH
                @products = split /\s+\+\s+/ , $products;

                foreach my $product (@products)
                {

                        #REMOVE TERMINAL WHITE SPACE
                        $product =~ s/\s+$//;
                        $product =~ s/^\s+//;


			#IF REACTANT IS PLACE HOLDER THEN IGNORE
			if( $product =~ /\*/ )
			{

				next;

			}#end if


                        #CHECK IF METABOLITE HAS COEFFICIENT, IF SO REMOVE IT
                        if( $product =~ /^\(*([0-9|\.|e|-]+)\)*\s/ )
                        {

                                $product =~  s/^\(*([0-9|\.|e|-]+)\)*\s+([\w\_\-]+)$/$2/;

                        }#end if


			#STORE ENTRY IN BOTH HASHES
                        if( $is_rev == 1 )
                        {

                                #$parsed_rxn_hash{$rxn_code} -> {$product} = "B";
                                #$parsed_metab_hash{$product} -> {$rxn_code} = "B";
                                $parsed_rxn_hash{$rxn_code} -> {$product} = "BP";
                                $parsed_metab_hash{$product} -> {$rxn_code} = "BP";

                        }else{

                                $parsed_rxn_hash{$rxn_code} -> {$product} = "P";
                                $parsed_metab_hash{$product} -> {$rxn_code} = "P";

                        }#end if



                }#end foreach



        }#end foreach


	return (\%parsed_rxn_hash , \%parsed_metab_hash);


}#end parse_rxn_participation


1;
