#!/usr/local/bin/perl5.6.1 -w

use Ace;
use strict;

$|=1;

# open a local database connection
my $db = Ace->connect(-path  =>  '/wormsrv2/current_DB') ;

my @genes = $db->fetch(-class => "Predicted_gene");

my $counter=0;

foreach my $gene (@genes){

	if($gene->at("Properties.Coding.Confirmed_by")){

	  # check coordinate system exons in relation to each other
	  my @exon_coord1 = sort by_number ($gene->get('Source_Exons',1));
	  my @exon_coord2 = sort by_number ($gene->get('Source_Exons',2));
	  my $intron_count = scalar(@exon_coord1)-1;
	  my $intron_lengths = 0;

	  for(my $i=1; $i<@exon_coord2;$i++){
	    $intron_lengths += ($exon_coord1[$i] - $exon_coord2[$i-1] +1);
		
	  }
	  
	  print "$gene $intron_count $intron_lengths\n"
	}	

}

print "\n";
$db->close;
exit;





sub by_number{ $a <=> $b;}
