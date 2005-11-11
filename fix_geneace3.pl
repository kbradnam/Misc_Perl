#!/usr/local/bin/perl5.8.0 -w

use strict;
use Ace;


my $geneace = Ace->connect(-path  => '/wormsrv1/geneace/');
#my $db = Ace->connect(-path  => '/nfs/disk100/wormpub/DATABASES/current_DB');
my $db = Ace->connect(-path  => '/wormsrv2/stlace');

my $class = "CDS";

my $query = "Find $class WHERE Method = history AND NOT (Gene || Gene_history)";
my @genes = $db->fetch(-query => $query);

foreach my $gene (@genes){
  my $history = $gene;
  $gene =~ s/:.*//;
  my $cds = $gene;

  # lose suffix if isoform
  chop ($cds) if ($gene =~ m/[a-z]$/);

  my $gene_name = $geneace->fetch(-class=>'Gene_name',-name=>"$cds");

  if (defined($gene_name)){
    my $gene_id = $gene_name->Public_name_for;
    if(defined($gene_id)){
      my $id = $gene_id->name;
      print "\n//$history has a gene name ($gene_name) connected to $id\n\n";
      print "CDS : \"$history\"\n";
      print "Gene_history $id\n\n";
      
    }
    else{
      print "\n//Panic!\n\n";
    }
  }
  else{
    print "\n//$history has no gene_name in geneace\n\n";
  }
}

$db->close;
$geneace->close;
exit;





