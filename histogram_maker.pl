#!/usr/bin/perl -w
#
# histogram_maker.pl
#
# A script to read in a file of numbers (one number per line) and sort into a table 
# ready to make a histogram from
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use Getopt::Long;

########################
# Command line options
########################

my $interval;      # what size bins to count things in?
my $min;           # force a start point for counting
my $max;           # force an end point for counting

GetOptions ("interval=i"      => \$interval,
	    	"min=i"           => \$min,
	    	"max=i"           => \$max);
	    

##################################################
# Sanity check on specified command line options
##################################################

# min should be less than max
if(($min && $max) && ($min > $max)){
   die "$min needs to be less than or equal to $max\n";
}

# array to store all numbers
my @numbers;

# input file should only contain numbers
open (IN, "<$ARGV[0]") || die "Failed to open input file\n\n"; 

while (<IN>) { 
	chomp;
	die "Non numeric character in input file\n" if (m/[^0-9\.]/);
	push(@numbers,$_);	
}	  	    
    
close(IN) || die "Couldn't close input file\n";  
 
# sort @numbers  
my @sorted = sort {$a <=> $b} @numbers;
my $records = scalar(@sorted);

# choose default interval,min, and max values if none specified 
if(!defined($interval)){
	# default to a fifth of the max size as interval
	$interval = int($sorted[-1]/5);
} 

if(!defined($min)){
	$min = $sorted[0];
}

if(!defined($max)){
	$max = $sorted[-1];
}

print "Using interval size of $interval\n";
print "Using minimum size of $min\n";
print "Using maximum size of $max\n";
print "\n";
print "Parsing $records records\n\n";

#################################################################
# now loop through bins
#################################################################

my $bin_start;
my $bin_end;
#my $number;
my $percent;
my $bin_counter;

BIN:for ($bin_start = $min; $bin_start < $max; $bin_start += $interval){
	
	# exit loop if we get here and there is no data left
	last BIN if (@sorted == 0);
	
	# reset bin counter for each bin
	$bin_counter = 0;

	# set bin end 
    $bin_end = $bin_start+$interval-1;
    
	#  if $min is greater than many records, we should first count everything less than $min
	if($sorted[0] <$min){	
		while ($sorted[0] < $min){
			$bin_counter++;
			shift(@sorted);
			last BIN if (@sorted == 0);
		}
		$percent = sprintf("%.4f",($bin_counter/$records));
		print "<${bin_start},$bin_counter,$percent\n";
		$bin_counter = 0;
	}

	#  move to next bin if the smallest value is greater than the current bin end size
	if ($sorted[0] > $bin_end){
		$percent = sprintf("%.4f",($bin_counter/$records));
		print "${bin_start}-${bin_end},$bin_counter,$percent\n";
		next BIN;
	}
	
	# if we are still here than the current record should be in this bin size and we can add 
	# to current total
	while($sorted[0] >= $bin_start && $sorted[0] <= $bin_end){
		$bin_counter++;
		shift(@sorted);
		last if (@sorted == 0);
	}
	$percent = sprintf("%.4f",($bin_counter/$records));
    print "${bin_start}-${bin_end},$bin_counter,$percent\n";

	# extra block to print out a final >x category if there are records above value of $max
	if(@sorted !=0 && $sorted[0] >$max){
		my $size = scalar(@sorted);
		$percent = sprintf("%.4f",($size/$records));
   		print ">${bin_end},$size,$percent\n";	
	}

}

