#!/usr/local/bin/perl -w

use Ace;

use strict;

$|=1;

# open a local database connection

my $path = shift;
my $db = Ace->connect(-path  => $path);

my $query = "Find Sequence; *.*;!yk*";
my @sequences = $db->fetch(-query => $query);

my $counter =0;
my $errors =0;
foreach my $sequence (@sequences){

    if ($sequence =~ m/a$/){
	my $previous = $sequences[$counter-1];
	my $copy = $sequence;
	$copy =~ s/[a-z]$//;
	if ($copy eq $previous){
	    print "WARNING: $sequence and $previous both exist in $path\n";
	    $errors++;
	}
    }
    $counter++;
}

print "\nFound $errors errors\n";
