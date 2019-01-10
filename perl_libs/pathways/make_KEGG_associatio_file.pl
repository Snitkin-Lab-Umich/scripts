#Takes as input a three column file (KEGGID PATHID JUNK)

#Produces a genemerge association file with the first column being
#the gene ID and the second column being a comma delimted list of pathways


#DEFINE SOME VARIABLES
my $KEGG_file = shift @ARGV;

my $gene_merge_ass = shift @ARGV;

my %path_hash;


#GO THROUGH KEGG FILE AND MAKE HASH LINKING GENES TO THEIR PATHS
open INPUT , $KEGG_file;

foreach $entry (<INPUT>)
{

	#PREPARE LINE
	chomp $entry;

	my ($gene , $path , $junk)  = split /\t/ , $entry;

	$gene =~ s/(.+):(.+)/$2/;
	$path =~ s/(.+):(.+)/$2/;

	
	#SEE IF GENE EXISTS IN PATH HASH, IF NOT CREATE NEW ENTRY, IF SO ADD TO LIST
	if( exists($path_hash{$gene}) )
	{

		push @{$path_hash{$gene}} , $path;

	}else{


		@{$path_hash{$gene}} = ($path);

	}#end if	


}#end foreach

close INPUT;



#GO THROUGH HASH OF PATHWAY ASSOCIATIONS AND WRITE TO OUTPUT FILE
open OUTPUT , ">$gene_merge_ass";

foreach $gene ( keys %path_hash )
{

	my $paths = join ";" , @{$path_hash{$gene}};

	print OUTPUT $gene . "\t" . $paths . "\n"; 

}#end foreach

close OUTPUT;

