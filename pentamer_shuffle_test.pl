#!/usr/bin/perl -w
#
# pentamer_shuffle_test.pl
#
# Take a file of sequences, make two new files with half of the sequences randomly split between them
# then compare pentamer frequencies, repeat 100 times
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use lib "/Korflab/lib/perl";
use FAlite;

# Need one array to store original sequences, and two subsequent arrays to hold random subsets
my (@all_seqs,@first,@second);

die "Specify number of iterations as second argument\n" if (!$ARGV[1]);

# how many iterations of random shuffling to perform?
my $n = $ARGV[1];

# what is lowest and highest ratios found in each comparison, default for min is high and max is low
my $min = 1;
my $max = 0;



# read file, read sequences and add to array

open(IN,"<$ARGV[0]") || die "Couldn't open input file\n";

my $FA = new FAlite (\*IN);
my $counter = 0;
while (my $entry = $FA->nextEntry) {
  my $seq = $entry->seq;
  $counter++;
  push(@all_seqs,$seq);
}
close(IN);
print "\n$counter sequences\n";


# main loop
# Make two equal sized subsets of @all_seqs which are chosen randomly, compare pentamer frequencies
# and repeat $n times

my $number;

# hash to store pentamer counts, key is pentamer value is an array, first element is counts of
# first subset of sequences, second element is count of second subset of sequences
my %pentamer_counts;

for(my $i=1; $i<=$n; $i++){
    print "Iteration $i\n";

    # clear hash
    %pentamer_counts = ();
    
    # copy @all_seqs to @first array and clear @second arrays
    my @first = @all_seqs;
    @second = ();

    # now randomly move half of @first sequences into @second
    for(my $j = 1; $j <= ($counter/2);$j++){
        $number = int(rand(scalar(@first)));

	# add to second array
        push(@second,$first[$number]);    

	# and remove from first array
        splice(@first,$number,1,@first[$number..-1]);
    }    

    # how many pentamers are there in all sequences?
    my $count1;
    my $count2;
    
    # now need to compare pentamer frequencies in @first and @second
    foreach my $seq (@first){
	for (my $i = 1; $i < length($seq) - 5 + 1; $i++) {
	    $count1++;
	    my $pentamer = substr($seq, $i, 5);
	    $pentamer_counts{$pentamer}[0]++;
	}	
    }
    foreach my $seq (@second){
	for (my $i = 1; $i < length($seq) - 5 + 1; $i++) {
	    $count2++;
	    my $pentamer = substr($seq, $i, 5);
	    $pentamer_counts{$pentamer}[1]++;
	}	
    }

    # now need to compare ratios of pentamer frequencies
    foreach my $key (sort keys %pentamer_counts){
	my $freq1 = $pentamer_counts{$key}[0]/$count1;
	my $freq2 = $pentamer_counts{$key}[1]/$count2;
	my $ratio = $freq1/$freq2;
#	print "$key $pentamer_counts{$key}[0] $pentamer_counts{$key}[1] $freq1 $freq2 $ratio\n";
	if ($ratio < $min){
	    $min = $ratio;
	    print "MIN: $ratio\n";
	}
	if ($ratio > $max){
	    $max = $ratio;
	    print "MAX: $ratio\n";
	}
    }
}

print "MIN = $min, MAX = $max\n";

exit(0);

















