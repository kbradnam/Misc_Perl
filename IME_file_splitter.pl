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
use FAlite;
use Getopt::Long;

# script can be run in two main modes 
# i) separate whole introns into separate files based on whether the *start* of the intron is within
#    a required range from the start of transcription (can be biased if just 1st base pair of intron 
#    is less than threshold)
# ii) Separate just the parts of introns that fall within in the specified size category...this might mean
#     removing 5' and 3' ends. Can then calculate base composition on just those parts of an intron that are 
#     a defined range from the TSS


# command line options
my $mode1; # see the two different descriptions above
my $mode2; 
my $window; # how big a size window
my $max;    # how big to go up to
my $step;   # for sliding windows

GetOptions ("mode1" => \$mode1,
			"mode2" => \$mode2,
			"window=i" => \$window,
			"max=i" => \$max,
			"step=i" => \$step);

die "Specify -mode1 or -mode2\n" if ($mode1 && $mode2);
die "Specify -mode1 or -mode2\n" if (!$mode1 && !$mode2);

# set some defaults
my $min = 1;
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);
my ($start,$stop);

for(my $start = $min;$start<$max; $start+= $step){
	$stop = $start + $window;
	
	print "Processing $start - $stop\n";
	
	open(OUT, ">${ARGV[1]}_${start}_${stop}.fa") || die "Couldn't create output file\n";

	open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0]\n";

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {	

		# process header
		my $header = $entry->def;
	   	$header =~ m/_i\d+_(\d+)_/;
		my $distance = $1;
		
		if(($distance >= $start) && ($distance <= $stop)){
			my $seq = $entry->seq;
			print OUT "$header\n$seq\n";			
		}

	}
	close(FILE);
	close (OUT);
}


exit(0);