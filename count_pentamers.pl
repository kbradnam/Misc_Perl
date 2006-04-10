#!/usr/bin/perl -w
#
# count_penatmers.pl
#
# A script to count and contrast pentamer frequencies of two sequence files
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use lib "/Korflab/lib/perl";
use FAlite;


my ($file1, $file2) = @ARGV;

open(OUT, ">pentamer_comparison.txt") or die;

my $merFH = {};
my $freqFH = {};
my $merNL = {};
my $freqNL = {};

# Assume we are looking at pentamer frequencies, but could do anything
my $window = 5;
my $ratio;


open(IN, $file1) or die "";
my $FA = new FAlite (\*IN);
my $counter = 0;
while (my $entry = $FA->nextEntry) {
    my $seq = $entry->seq;
    &mer($merFH, $seq, $window);
}	
close IN;
print "Counter is $counter\n";

$counter = 0;



open(IN, $file2) or die "";
my $FA2 = new FAlite (\*IN);
while (my $entry = $FA2->nextEntry) {
    my $seq = $entry->seq;
    &mer($merNL, $seq, $window);
}	
close IN;

print "Counter is $counter\n";

&freq($merFH);
&freq($merNL);

#------------- print -----------------------
print OUT "Pentamer Ratio File1 File2\n";

foreach my $mFH (keys %$merFH) {
    $ratio = ($merFH->{$mFH})/($merNL->{$mFH});
    print (OUT "$mFH $ratio $merFH->{$mFH} $merNL->{$mFH}\n");
#q    print (OUT "$mFH $merFH->{$mFH} $merNL->{$mFH}\n");
}

close OUT;


# -------------SUBROUTINE -------------------------------------- 

sub mer {
    my ($merclass, $seq, $window) = @_;
    for (my $i = 1; $i < length($seq) - $window + 1; $i++) {
	my $mer = substr($seq, $i, $window);
	$merclass->{$mer}++;
	$counter++;
    }	
}	 	

sub freq {
    my $total = 0;
    my ($merClass) = @_;
    foreach my $mer (sort keys %$merClass) {
	$total += $merClass->{$mer};
    }
    foreach my $mer (sort keys %$merClass) {
	$merClass->{$mer} /= $total;
    }
    print "Total is $total\n";
}
















