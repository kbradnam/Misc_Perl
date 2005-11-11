#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;


my $tace = &tace;
my $dir = "/wormsrv1/geneace";
my @loci =`echo "table-maker -p ${dir}/wquery/loci_with_seq_but_no_pos_clone.def" | $tace /wormsrv1/geneace`;


my $command=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/loci_with_seq_but_no_pos_clone.def" quit
EOF

open (FH, "echo '$command' | $tace $dir | ") || die "Couldn't access geneace\n";
while (<FH>){
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){
    my $locus = $1;
    my $other = $2;
    $other =~ s/\..*//;
    print "\nLocus : \"$locus\"\nPositive_clone \"$other\" Inferred_automatically \"From sequence, transcript, pseudogene data\"\n";
  }
}
close(FH);



exit(0);





