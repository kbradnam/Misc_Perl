#!/usr/bin/perl
#
# IME_file_splitter.pl
#
# A script to separate a FASTA file into separate files based on coordinate information in the FASTA header
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use Getopt::Long;

# This script separates an intron or exon FASTA file into separate files based on whether the 
# *start* of the intron is within a required range from the start of transcription (can be  
# biased if just 1st base pairof intron is less than threshold)

# command line options
my $split; # see the two different descriptions above
my $percent; 
my $window; # how big a size window
my $max;    # how big to go up to
my $step;   # for sliding windows
my $prefix; # what prefix to give new file names (if using -split mode)

GetOptions ("window=i" => \$window,
			"max=i"    => \$max,
			"step=i"   => \$step,
			"prefix=s" => \$prefix);

# set some defaults
my $min = 1;
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);

# need to know end points of each window
my ($start,$end);

# need to store merged sequences
my $new_seq;

# count how many sequences in each bin
my $counter;

for(my $start = $min;$start<$max; $start+= $step){
	$new_seq = "";
	$counter = 0;
	$end = $start + $window -1;
	
	print "Processing $start - $end\n";
	open(OUT, ">${prefix}_${start}_${end}.fa") || die "Couldn't create output file\n";
	
	open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0]\n";
	
	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {
	    
		# process header
		my $header = $entry->def;

		# check whether we are working with intron or exon data
		my $distance;
		if($header =~ m/_i\d+_\d+_/){
			$header =~ m/_i\d+_(\d+)_/;
			$distance = $1;
		}
	   	elsif($header =~ m/_e\d+_\d+/){
			$header =~ m/_e\d+_(\d+)/;
			$distance = $1;
		}
		# check whether candidate intron/exon falls in size category
		if(($distance >= $start) && ($distance <= $end)){
			$counter++;
			my $seq = $entry->seq;
			if ($split){ 
				print OUT "$header\n$seq\n" 			
			}
			else{
			
				my $length = length($seq);

				# if sequence is longer than $end, then just need part to extract
				# part of it
				if(($distance+$length) > $end){
					my $tmp = substr($seq,1,$end-$distance);
					$new_seq .= $tmp;
				}
				else{
					$new_seq .= $seq;
				}
			}
		}		
	}
	
	close(FILE);
	close (OUT) if ($split);
}


exit(0);
