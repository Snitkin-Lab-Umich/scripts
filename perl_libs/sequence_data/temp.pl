use POSIX;

my $seq = "AAATAAA";

my $length = length($seq);

my $surr_seq_len = floor($length / 2);


print "$surr_seq_len\n";
