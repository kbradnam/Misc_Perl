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

# hash to store pentamer counts, key is pentamer value is an array, first element is counts of
# first subset of sequences, second element is count of second subset of sequences
my %pentamer_counts;

my ($file1, $file2) = @ARGV;


# Assume we are looking at pentamer frequencies, but could do anything
my $window = 5;
my $ratio;


#########################################
# loop through sequences in first file
#########################################

open(IN, $file1) or die "";
my $FA = new FAlite (\*IN);
my $pentamer_count1 = 0;
my $seq_count1 = 0;

while (my $entry = $FA->nextEntry) {
    my $seq = $entry->seq;
    $seq_count1++;
    
    # now need to count pentamers
    for (my $i = 1; $i < length($seq) - $window + 1; $i++) {
	$pentamer_count1++;
	my $pentamer = substr($seq, $i, $window);
	$pentamer_counts{$pentamer}[0]++;
    }	

    # also need to keep track of how many sequences contain each type of pentamer
    foreach my $key (keys %pentamer_counts){
	($pentamer_counts{$key}[2]++) if((defined($pentamer_counts{$key}[0])) && ($pentamer_counts{$key}[0] > 0));
    }

}	
close(IN);
print "$pentamer_count1 pentamers in first file ($seq_count1 sequences)\n";



#########################################
# loop through sequences in second file
#########################################

open(IN, $file2) or die "";
$FA = new FAlite (\*IN);
my $pentamer_count2 = 0;
my $seq_count2 = 0;

while (my $entry = $FA->nextEntry) {
    my $seq = $entry->seq;
    $seq_count2++;
    # now need to count pentamers
    for (my $i = 1; $i < length($seq) - $window + 1; $i++) {
	$pentamer_count2++;
	my $pentamer = substr($seq, $i, $window);
	$pentamer_counts{$pentamer}[1]++;
    }	

    # also need to keep track of how many sequences contain each type of pentamer
    foreach my $key (keys %pentamer_counts){
	($pentamer_counts{$key}[3]++) if((defined($pentamer_counts{$key}[1])) && ($pentamer_counts{$key}[1] > 0));
    }
}	
close(IN);
print "$pentamer_count2 pentamers in second file ($seq_count2 sequences)\n";



##############################################################
# Calculate frequencies of pentamers and ratio of frequencies
# between two files
##############################################################

open(OUT, ">pentamer_comparison.txt") or die;

# now need to compare ratios of pentamer frequencies
foreach my $key (sort keys %pentamer_counts){

	# can only calculate ratio if there is the same pentamer in both files
	next unless (defined($pentamer_counts{$key}[0]) && defined($pentamer_counts{$key}[1]));
	
    my $freq1 = $pentamer_counts{$key}[0]/$pentamer_count1;
    my $freq2 = $pentamer_counts{$key}[1]/$pentamer_count2;
	
#	next if ($pentamer_counts{$key}[0] < 20);
#	next if ($pentamer_counts{$key}[1] < 20);
	
    my $ratio = $freq1/$freq2;
    #print "$key $pentamer_counts{$key}[2] $pentamer_counts{$key}[3]\n";
    print OUT "$key $pentamer_counts{$key}[0] $pentamer_counts{$key}[1] $freq1 $freq2 $ratio\n";
}

close(OUT);

exit(0);



