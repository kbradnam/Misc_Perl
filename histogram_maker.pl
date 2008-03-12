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

my $window;    # what size bins to count things in?
my $step;	   # what step size to use
my $min;       # force a start point for counting
my $max;       # force an end point for counting
my $integer;   # interger mode, easier to do the counting
my $float;	   # floating point mode, not so easy to do the counting

GetOptions ("window=f"      => \$window,
			"step=f"		=> \$step,
	    	"min=f"         => \$min,
	    	"max=f"         => \$max,
			"integer"		=> \$integer,
			"float"			=> \$float);

#################################################
# Sanity check on specified command line options
#################################################

# min should be less than max
if(($min && $max) && ($min > $max)){
   die "$min needs to be less than or equal to $max\n";
}

if(($float && $integer) || (!$float && !$integer)){
   die "Please specify -integer if only dealing with interers, *or* -float if only dealing with floating-point numbers\n";
}


# array to store all numbers
my @numbers;

# input file should only contain numbers
open (IN, "<$ARGV[0]") || die "Failed to open input file\n\n"; 

while (<IN>) { 
	chomp;
	die "Non numeric character in input file\n" if (m/[^0-9\.\-]/);
	push(@numbers,$_);	
}	  	    
    
close(IN) || die "Couldn't close input file\n";  
 
# sort @numbers  
my @sorted = sort {$a <=> $b} @numbers;
my $records = scalar(@sorted);

# choose default window, step, min, and max values if none specified 
if(!defined($window)){
	# default to a fifth of the max size as window
	$window = int($sorted[-1]/5);
} 

if(!defined($min)){
	$min = $sorted[0];
}

if(!defined($max)){
	$max = $sorted[-1];
}

# step size is equal to window size if not specified
if(!defined($step)){
	$step = $window;
}

print "Using window size of $window\n";
print "Using step size of $step\n";
print "Using minimum size of $min\n";
print "Using maximum size of $max\n";
print "\n";
print "Parsing $records records\n\n";

#################################################################
# now loop through bins
#################################################################

my $bin_start;
my $bin_end;
my $rounded_bin_end;
my $percent;
my $bin_counter;

# change $step size if dealing with floating point numbers
$step-- if ($float);
 
BIN:for ($bin_start = $min; $bin_start < $max; $bin_start += $step){
	# exit loop if we get here and there is no data left
	last BIN if (@sorted == 0);
	
	# reset bin counter for each bin
	$bin_counter = 0;

	# are we only bothering to count integers?
	if($integer){
		$bin_end = $bin_start+$window-1;
		$rounded_bin_end = $bin_start+$window-1;
	}
	# or do we have to deal with floating point numbers?
	elsif($float){
		# set bin end, we go down to 10 decimal places which might not be sensitive enough
	    # the rounded bin end figure is just for the final print out to look a bit neater
		$bin_end = $bin_start+$window-1.0000000001;
		$rounded_bin_end = $bin_start+$window-1.01;		
	}
	
#	print "# sorted[0] = *$sorted[0]*  bin_start = *$bin_start* bin_end = $bin_end\n";
	

	# if $min is greater than many records, we should first count everything less than $min
	if($sorted[0] <$min){	
		while ($sorted[0] < $min){
			$bin_counter++;
			shift(@sorted);
			last BIN if (@sorted == 0);
		}
		$percent = sprintf("%.4f",($bin_counter/$records));
		print "\"<${bin_start}\",$bin_counter,$percent\n";
		$bin_counter = 0;
	}

	#  move to next bin if the smallest value is greater than the current bin end size
	if ($sorted[0] > $bin_end){
		$percent = sprintf("%.4f",($bin_counter/$records));
		print "\"${bin_start} - ${rounded_bin_end}\",$bin_counter,$percent\n";
		next BIN;
	}

	# if we are still here than the current record should be in this bin size and we can add 
	# to current total
	# might have to keep some values which overlap with next bin, so need a counter to increment
	# through @sorted
	my $i = 0;
	while( (close_enough_or_bigger($sorted[$i],$bin_start)) && (close_enough_or_smaller($sorted[$i],$bin_end)) ){
#		print "$bin_counter\n";
		$bin_counter++;

		last if (@sorted == 0);

		# only get rid of number from array if it can't possibly be in the next bin
		if($sorted[$i] < ($bin_start + $step)){
			shift(@sorted);
		}
		else{
			$i++;
		}

	}
	
	$percent = sprintf("%.4f",($bin_counter/$records));
    print "\"${bin_start} - ${rounded_bin_end}\",$bin_counter,$percent\n";

	# extra block to print out a final >x category if there are records above value of $max 
	if(@sorted !=0 && (close_enough_or_bigger($sorted[0],$max))){
		my $size = scalar(@sorted);
		$percent = sprintf("%.4f",($size/$records));
   		print "\">${rounded_bin_end}\",$size,$percent\n";	
	}
}

# Need new functions to test for equality when dealing with floating point numbers
# E.g. 0.6 in perl might actually be stored as 0.60000000000000008882 which 
# means that it won't evaluate to true when you use the == operator
sub close_enough{
	my $a = shift;
	my $b = shift;
	
	# find the difference and compare to your desired theshold whether this is close
	# enough to be considered 'the same'. Differences have to be smaller than 10 decimal
	# places to be considered the same which should work in most common sense scenarios
	my $result = abs($a-$b);
	if ($result < 0.0000000001){
		return(1);
	}
	else{
		return(0);
	}
}

sub close_enough_or_bigger{
	my $a = shift;
	my $b = shift;
	
	# find the difference and compare to your desired theshold whether this is close
	# enough to be considered 'the same'. Differences have to be smaller than 10 decimal
	# places to be considered the same which should work in most common sense scenarios
	my $result = abs($a-$b);
	if ($result < 0.0000000001 || ($a>$b)){
		return(1);
	}
	else{
		return(0);
	}
}

sub close_enough_or_smaller{
	my $a = shift;
	my $b = shift;
	
	# find the difference and compare to your desired theshold whether this is close
	# enough to be considered 'the same'. Differences have to be smaller than 10 decimal
	# places to be considered the same which should work in most common sense scenarios
	my $result = abs($a-$b);
	if ($result < 0.0000000001 || ($a<$b)){
		return(1);
	}
	else{
		return(0);
	}
}
