#!/usr/local/bin/perl5.6.1 -w

use Ace;
use strict;

$|=1;

# open a local database connection
my $db = Ace->connect(-path  =>  '/wormsrv2/current_DB') ;

my @proteins = $db->fetch(-class => "Wormpep");

my $counter=0;

foreach my $protein (@proteins){
  my @columns=();
  
  @columns = $protein->at('Feature');
  foreach my $col (@columns){
    if ($col =~ m/tmhmm/){
      my @tm_row = $protein->at("Feature.$col")->row;
      my @tm_col = ();
      @tm_col = $protein->at("Feature.$col")->col;
      if((scalar(@tm_col) == 1) && ($tm_row[1] < 50)){
	print "$protein $tm_row[1] $tm_row[2]\n";
      }
    }
  }
}

print "\n";
$db->close;
exit;





