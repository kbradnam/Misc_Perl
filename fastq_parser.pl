#!/usr/bin/perl
#
# fastq_parser.pl
#
# A script to check on how many sequences in a pair of FASTQ files might be too poor quality due to high Ns
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Getopt::Std;
use vars qw($opt_t $opt_n);
getopts('n:t:');

# defaults...look at 1 million sequences trying to find those with >50% Ns
my $THRESHOLD = 0.5;
my $N         = 1000000;

my $usage = "Usage: $0 <1st_fastq_file> <2nd_fastq_file> 
options:
  -t <float> threshold: maximum allowable fraction of N characters in a read [$THRESHOLD]
  -n <int> number of sequences to process in input files [$N]";

die "$usage\n" unless (@ARGV == 2);
my ($file1, $file2) = @ARGV;

$THRESHOLD = $opt_t if $opt_t;
$N         = $opt_n if $opt_n;

# need three counters, for when either sequence has higher %N than threshold, and when both exceed threshold 
my ($count1, $count2, $count_both) = (0, 0, 0);


# open files
open(my $in1, "<", $file1) or die "Can't open $file1\n";
open(my $in2, "<", $file2) or die "Can't open $file2\n";

# don't loop over entire file, just the first N sequences
for (my $i = 0; $i < $N; $i++){
	
	# read four lines at a time from both files, 2nd line is the sequence
	my $head1    = <$in1>;
	my $seq1     = <$in1>;
	my $id1      = <$in1>;
	my $quality1 = <$in1>;

	my $head2    = <$in2>;
	my $seq2     = <$in2>;
	my $id2      = <$in2>;
	my $quality2 = <$in2>;

	# calculate %N on two sequences
	my $n1 = calculate_percent_n($seq1);
	my $n2 = calculate_percent_n($seq2);
	
	# increase counters if necessary
	$count1++     if ($n1 > $THRESHOLD);
	$count2++     if ($n2 > $THRESHOLD);
	$count_both++ if ($n1 > $THRESHOLD && $n2 > $THRESHOLD);
}

close $in1;
close $in2;

# just report on final stats
print "Processed $N sequences in both files\n";
print "$count1 sequences in file 1 had %N above $THRESHOLD\n";
print "$count2 sequences in file 2 had %N above $THRESHOLD\n";
my $percent = sprintf("%.3f", $count_both / $N);
print "$count_both sequences ($percent%) in file 1 and file 2 had %N above $THRESHOLD\n";


# next steps
# 1) print out sequences to new file, but only if they are both good (low %N)
# 2) don't increment counter when either sequence is bad (high %N)

########################
# S u b r o u t i n e s 
########################

sub calculate_percent_n{
	my ($seq) = @_;
	$seq = uc($seq);
	chomp($seq);
	my $length = length($seq);
	my ($n) = $seq =~ tr/N/N/;
	my $percent = sprintf("%.1f", $n / $length);

	return($percent);
}