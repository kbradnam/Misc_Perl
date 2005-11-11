#!/usr/local/bin/perl

$/ = "";

open(IN,"csh_old_data.ace") || die "Couldn't open file\n";
open(OUT,">stripped_data.ace") || die "Couldn't create output file\n";

while(<IN>){
  print OUT if ($_ =~ /^PCR_product/ && $_ =~ /RNAi/);

}
close(IN);
close(OUT);


