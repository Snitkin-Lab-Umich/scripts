#!/usr/bin/perl

$t1 = $ARGV[0];

$outf = " > ".$t1.".phylip";

open (F, " < ".$t1);

$taxa = -1;

while($line = <F>)
  {
    if($line =~ />/)
      {    
	$taxa++;
	$name = $line;
	$name =~ s/[\s+|\(|\)|\,|;]/_/g;
	$name =~ s/,//g;
	$name =~ s/>//g;
	$taxonNames[$taxa] = $name;
      }
    else
      {
	$seq = $line;
	$seq =~ s/\s+//g;
	$sequences[$taxa] = $sequences[$taxa].$seq;
      }      
  }

close(F);

for($i = 0; $i <= $taxa; $i++)
  {
    print $taxonNames[$i]." ".(length($sequences[$i]))."\n";
  }

$s  = $taxa + 1;
$bp = length($sequences[0]);

open (F, $outf);

print F $s." ".$bp."\n";

for($i = 0; $i <= $taxa; $i++)
  {
    print F $taxonNames[$i]." ".$sequences[$i]."\n";
  }

