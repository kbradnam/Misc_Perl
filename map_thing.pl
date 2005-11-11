#!/usr/local/bin/perl5.8.0 -w
#
# map_thing.pl
# 
# test mapping script by Keith Bradnam aged 12 and a half

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;


#############
# Paths etc #
#############

my $giface   = &giface;
my $dbdir    = "/nfs/disk100/wormpub/DATABASES/current_DB";

#get list of genes using fetch many.  Note that we only want SMapped genes which will have a span
my $db = Ace->connect(-path  => $dbdir,
                      -program =>$giface) || do { print "Connection failure: ",Ace->error; die();}; 
my $i = $db->fetch_many(-query=>"Find Gene WHERE S_parent");


# open output streams and files
open(GIFACE,"| $giface $dbdir") or die "failed to open giface connection to $dbdir\n";
open(OUT, ">mappings.ace") or die "Cannot open mapping output file\n";
my $output   = "/nfs/team71/worm/krb/tmp_giface_output"; 


# loop through genes
while (my $gene = $i->next){

  my $command = "gif seqget Gene:$gene; seqfeatures -version 2 +source RNAi|Oligo_set|Orfeome|GenePair_STS -file $output\n";
  print GIFACE $command;
  print OUT "\n\nGene : $gene\n";

  open(TMP,"<$output") || die "Can't open file\n";
  while(<TMP>){
    if ($_ !~ m/^\#/){
      my @line = split(/\t/);
      print OUT "$line[8]";
    }
  }
  close(TMP);
}

  #gif seqget CDS:AH6.1; seqfeatures -version 2 +source RNAi|Oligo_set

close(OUT);
close(GIFACE);

$db->close;
exit(0);


