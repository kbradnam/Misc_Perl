#!/usr/local/bin/perl5.6.1 -w
 
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;
use Ace;

my $tace = &tace;

my $db = Ace->connect(-path=>"/nfs/disk100/wormpub/DATABASES/current_DB",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


my @genes= $db->fetch(-query=>'find worm_genes');

my %class;

foreach my $gene (@genes){
  if (defined $gene->at('Properties.Coding.CDS')){
    $class{$gene} = "CDS";
#    print "$gene CDS\n"; 
  }
  elsif (defined $gene->at('Properties.Transcript')){
    $class{$gene} = "RNA";
    print "$gene RNA\n";
  }
  elsif (defined $gene->at('Type.Coding_pseudogene') || defined $gene->at('Type.RNA_pseudogene')){
    $class{$gene} = "pseudo";
    print "$gene pseudo\n";
  }
}
$db->close;
