#!/usr/bin/perl -w
#
# duret.pl
#
# a script to test my PhD idea about optimal codon usage in protein domains
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use lib "/Korflab/lib/perl";
use Ace;
use strict;

# open a local database connection
my $db = Ace->connect(-path  =>  '/Korflab/Data_sources/WormBase/WS150') ;

my @genes = $db->fetch(-class => "elegans_CDS");

my $counter=0;

foreach my $gene (@genes){
  
  my $status = $gene->Prediction_status;
  
  my $lab = $gene->From_laboratory;

  my ($protein) = $gene->at("Visible.Corresponding_protein");
  next if !($protein);
  my $domain_space = 0;   # no of bases taken up by Pfam domains
  my $domain_counter = 0; # no of Pfam domains
  my @columns=();
 
  
  $protein = $db->fetch(Protein => "$protein") || die "Cannot fetch protein\n";
  my $length = $protein->at("Peptide")->right(2);

  @columns = $protein->at('Homol.Motif_homol');

  # want to look for overlapping pfam domains by comparing min and max coordinates
  # of each Pfam domain that matches protein
  my $query_min  = 0;
  my $query_max  = 0;
  my $target_min = 0;
  my $target_max = 0;
  my $domain_overlap = "";

  foreach my $col (@columns){

    # look for Pfam domain matches
    if ($col =~ m/PFAM/){
      $domain_counter++;

      my @pfam_row = $protein->at("Homol.Motif_homol.$col")->row;
      $query_min = $pfam_row[3];
      $query_max = $pfam_row[4];

      $domain_space += ($query_max - $query_min + 1);

      if($domain_counter == 1){
	# store coordinates of first Pfam domain as target to compare with
	$target_min = $query_min;
	$target_max = $query_max;
      }
      if($domain_counter >1){
	# test for overlapping domains
	# three potential scenarios for overlap
	if((($query_min >= $target_max) && ($query_max <= $target_max)) ||
	   (($query_min <= $target_max) && ($query_max >= $target_max)) ||
	   (($query_min <= $target_min) && ($query_max >= $target_min))){
	  $domain_overlap = "ERROR: Pfam overlap!";
	}
	# now set new $target_min, $target_max in case there are more than 2 domains
	$target_min = $query_min;
	$target_max = $query_max;

      }
    }
  }

  print "$gene,$lab,$status,$protein,$length,$domain_counter,$domain_space,$domain_overlap\n";
#  last if ($counter > 1000);
  $counter++;
  

  # kill AcePerl objects
  $gene->DESTROY();
  $protein->DESTROY();
}

print "\n";
$db->close;
exit;





