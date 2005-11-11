#!/usr/local/bin/perl5.8.0 -w

use strict;

my $weight = 0;
my $peptide = "AAAAAA";

my %weights = (
	    "A" => "1","C" => "2", "D" => "3","E" => "4", 
	    "F" => "5","G" => "6", "H" => "7","I" => "8", 
	    "K" => "9","L" => "10","M" => "1","N" => "2", 
	    "P" => "3","Q" => "4", "R" => "5", "S" => "6", 
	    "T" => "7","V" => "8","W" => "9","Y" => "10","U" => "2"
	    );


foreach (keys %weights){
  my $count = $peptide =~ s/$_//g;
  $weight += ( $weights{$_} * $count);
}

  print "weight is $weight\n";
