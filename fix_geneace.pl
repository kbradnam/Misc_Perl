#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;


# read hash of locus2gene
my %locus2gene = &FetchData('cgc_name2gene');

#foreach (keys(%locus2gene)){
#  print "$_ = $locus2gene{$_}\n\n";
#}

# load data file
open(IN, "<links.txt") || die "Can't read file\n";

while(<IN>){
  s/\"//g;

  my ($multi,$locus,$polymorphism,$allele) = split /\t/;

  # skip gene clusters
  next if ($locus eq "act-123");

  my $gene = $locus2gene{$locus};
  print "Multi_pt_data : \"$multi\"\n";
  print "-D Locus $locus\n\n";

  print "Multi_pt_data : \"$multi\"\n";
  if(defined($allele)){
    print "Gene $gene $allele\n\n";
  }
  else{
    print "Gene $gene\n\n";
  }
}

close(IN);
exit;





