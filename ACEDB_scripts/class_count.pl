#!/usr/local/bin/perl -w


use Ace;

use strict;

$|=1;

my $path = shift;
# open a local database connection
my $db = Ace->connect(-path  =>  $path);


my %classes = $db->class_count();


foreach my $key (sort keys %classes){

   print "$key -> $classes{$key}\n";

}


