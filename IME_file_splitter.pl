#!/usr/bin/perl
#
# IME_file_splitter.pl
#
# A script to separate a FASTA file into separate files based on matching a pattern in the FASTA header
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;

# script can be run in two main modes 
# i) SPLIT MODE - separate whole introns into separate files based on whether the *start* of the intron is within
#    a required range from the start of transcription (can be biased if just 1st base pair of intron 
#    is less than threshold)
# ii) PERCENT MODE - Separate just the parts of introns that fall within in the specified size category...this might mean
#     removing 5' and 3' ends. Can then calculate base composition on just those parts of an intron that are 
#     a defined range from the TSS. Only show calculated percentages, don't output sequences


# command line options
my $split; # see the two different descriptions above
my $percent; 
my $window; # how big a size window
my $max;    # how big to go up to
my $step;   # for sliding windows
my $prefix; # what prefix to give new file names (if using -split mode)

GetOptions ("split"    => \$split,
			"percent"  => \$percent,
			"window=i" => \$window,
			"max=i"    => \$max,
			"step=i"   => \$step,
			"prefix=s" => \$prefix);

die "Specify -split or -percent\n" if ($split && $percent);
die "Specify -split or -percent\n" if (!$split && !$percent);

# set some defaults
my $min = 1;
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);

# need to know end points of each window
my ($start,$stop);

# need to store merged sequences
my $new_seq;

# count how many sequences in each bin
my $counter;

for(my $start = $min;$start<$max; $start+= $step){
	$new_seq = "";
	$counter = 0;
	$stop = $start + $window;
	
	
	# only write output files if in split mode
	if($split){
		print "Processing $start - $stop\n";
		open(OUT, ">${prefix}_${start}_${stop}.fa") || die "Couldn't create output file\n";
	}
	
	open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0]\n";

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {	

		# process header
		my $header = $entry->def;
	   	$header =~ m/_i\d+_(\d+)_/;
		my $distance = $1;
		
		if(($distance >= $start) && ($distance <= $stop)){
			$counter++;
			my $seq = $entry->seq;
			if ($split){ 
				print OUT "$header\n$seq\n" 			
			}
			else{
			
				my $length = length($seq);

				# if sequence is longer than $stop, then just need part to extract
				# part of it
				if(($distance+$length) > $stop){
					my $tmp = substr($seq,1,$stop-$distance);
					$new_seq .= $tmp;
				}
				else{
					$new_seq .= $seq;
				}
			}
		}		
	}
	# can now calculate base composition
	if($percent){
		(my ($a,$c,$g,$t,$n,$o) = Keith::base_composition($new_seq,1));
		print "$start,$stop,$counter,$a,$c,$g,$t,$n\n";
	}
	close(FILE);
	close (OUT) if ($split);
}


exit(0);