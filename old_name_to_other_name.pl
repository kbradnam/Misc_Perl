#!/usr/local/bin/perl5.6.1 -w

use Ace;

use strict;

$|=1;

# open a local database connection

my $path = "/wormsrv1/geneace/";
my $db = Ace->connect(-path  => $path);

my $query = "Find Locus; Old_name";
my @loci = $db->fetch(-query => $query);


foreach my $locus (@loci){
  my ($tag) = $locus->at('Name.Old_name');
  print "\nLocus : \"$locus\"\n";
  print "-D Old_name \"$tag\"\n";
  print "Other_name \"$tag\"\n";

}


