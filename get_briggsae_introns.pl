#!/usr/bin/perl -w
#
# map_thing.pl
# 
# test mapping script by Keith Bradnam aged 12 and a half

use strict;
use lib "/Korflab/lib/perl";
use Ace;

#############
# Paths etc #
#############
my $dbdir    = "/Korflab/Data_sources/WormBase/WS140";
my $giface   = "/Korflab/bin/giface";


#get list of genes using fetch many.  Note that we only want SMapped genes which will have a span
my $db = Ace->connect(-path  => $dbdir,
                      -program =>$giface) || do { print "Connection failure: ",Ace->error; die();}; 
my $i = $db->fetch_many(-query=>"Find briggsae_CDS");


# open output streams and files
open(GIFACE,"| $giface $dbdir") or die "failed to open giface connection to $dbdir\n";

# loop through genes
while (my $gene = $i->next){
    print "$gene\n";
  my $command = "gif seqget CDS:$gene; seqfeatures -version 2 +source curated -file briggsae.tmp\n";
  print GIFACE $command;
    `cat briggsae.tmp >> briggsae.gff`;
}

close(GIFACE);

$db->close;
exit(0);


