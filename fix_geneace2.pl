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
open(IN, "</wormsrv1/geneace/physical_map_data.ace") || die "Can't read file\n";

while(<IN>){

  if(m/Positive_locus/){

    m/Positive_locus \"(.*)\"/;
    my $locus = $1;

    my $gene = $locus2gene{$locus};
    if($gene){
      print "Positive_gene $gene\n";
    }
    else{
      print;
    }
  }
  else{
    print;
  }
}

close(IN);
exit;





