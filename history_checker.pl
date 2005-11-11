#!/usr/local/bin/perl5.8.0 -w

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use strict;
use Ace;

my $tace = &tace;




my $geneace = Ace->connect(-path=>"/wormsrv1/geneace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


# make hash of cds2gene from geneace
my %cds2gene;
my @genes = $geneace->fetch(-query=>"Find Gene_name WHERE Sequence_name_for");
foreach my $gene_name (@genes){
  my $gene = $gene_name->Sequence_name_for;
  $cds2gene{$gene_name} = $gene;
#  print "$gene_name -> $gene\n";
}

$geneace->close;

# open connection to DB
#my $db = Ace->connect(-path=>"/nfs/disk100/wormpub/DATABASES/current_DB",
my $db = Ace->connect(-path=>"/nfs/disk100/wormpub/camace_krb",
#my $db = Ace->connect(-path=>"/wormsrv2/autoace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


# grab set of problem history CDSs without Gene IDs
my @CDSs= $db->fetch(-query=>'find CDS WHERE Method = history AND !Gene_history');

my $lab;

# find out highest gene number in case new genes need to be created

foreach my $cds (@CDSs){

  # work out lab
  if ($cds->From_laboratory eq "RW"){
    $lab = "RW";
  }
  elsif ($cds->From_laboratory eq "HX"){
    $lab = "HX";
  }
   
  # get rid of wp suffix
  my $new_cds = $cds;
  $new_cds =~ s/:wp[0-9]*//;
  $new_cds =~ s/(\d+)[a-z]$/$1/;

  if ($new_cds=~ m/\.$/){
    print "// ERROR: $new_cds ($cds) has strange name\n";
    next;
  }

  if($lab eq "HX"){
    # lookup $new_cds in hash
    if($cds2gene{$new_cds}){
      my $gene = $cds2gene{$new_cds};
      print "// $cds does not have a gene connection but $new_cds -> $gene exists\n";
      print "CDS : \"$cds\"\n";
      print "Gene_history $gene\n\n";    
    }
    else{
      print "// $cds ($new_cds) does not have a gene connection.  Make new one?\n";
    }
  }
}


$db->close;


