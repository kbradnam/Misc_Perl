#!/usr/bin/perl
#
# IME_seq_dump.pl
#
# A script to extract regions of intron sequences for further IME analysis
# a simpler version of IME_seq_extractor.pl
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


# this script takes a file of IME formatted intron sequences (FASTA header contains intron position and distance to TSS)
# it then extracts regions of intron sequence based on specified criteria (how large a window, and how far from the TSS)
# it also allows you to specify that x nt from any one intron must be present in window to be included. All sequences
# are dumped to a new output file.

# command line options
my $window;     # how big a size window
my $start;	    # start coordinate in transcript
my $min;        # minimum amount of sequence to come from any one intron sequence (avoid taking just 1 nt of an intron)

GetOptions ("window=i"     => \$window,
			"start=i"	   => \$min,
			"min=i"        => \$min,);


# check command line options 
die "Specify a name of of an IME-formatted FASTA file of intron sequences\n" if (@ARGV != 1);

# set some defaults
$start = 1    if (!$start);
$window = 500 if (!$window);
$min = 50    if (!$min);

my $file = $ARGV[0];


# Determine coordinate range (in transcript coordinates) for which we want to extract intron sequence
my $win_start = $start;
my $win_end = $start+$window-1;

open(FILE,"<$file") || die "Couldn't open $file\n";		

# open output file
open(OUT,">$file.$win_start-${win_end}") || die "couldn't create $file.$win_start-${win_end} file\n";



my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
    
	# grab basic details
	my $seq = $entry->seq;
	my $length = length($seq);
	my $header = $entry->def;
			
	# grab distance to TSS differently depending on what type of sequences we are dealing with
	my $distance;
	$header =~ m/_i\d+_(\d+)/;
	$distance = $1;
	
	# calculate end coordinate of intron
	my $end_coord = $distance + $length -1;


	###################################################################
	# check whether candidate sequence falls in various size categories
	###################################################################

	# CASE 1: intron sequence is wholly contained within window
	if(($distance >= $win_start) && ($end_coord <= $win_end)){
		my $tmp = substr($seq,0,$end_coord-$distance+1);
		
		# dump sequence if using -seqdump and there is at least $min nt
		print OUT "$header $length nt\n$tmp\n" if (length($tmp) >= $min);
	}

	# CASE 2: intron sequence is larger than  window
	elsif(($distance < $win_start) && ($end_coord > $win_end)){
		my $tmp = substr($seq,$win_start-$distance,$window);
		
		# dump sequence if using -seqdump and there is at least $min nt
		print OUT "$header $length nt\n$tmp\n" if (length($tmp) >= $min);
	}

	# CASE 3: intron sequence overlaps 5' edge of window
	elsif(($distance < $win_start) && ($end_coord >= $win_start)){
		my $tmp = substr($seq,$win_start - $distance,$end_coord - $win_start + 1);
		
		# dump sequence if using -seqdump and there is at least $min nt
		print OUT "$header $length nt\n$tmp\n" if (length($tmp) >= $min);
	}

	# CASE 4: intron sequence overlaps 3' edge of window
	elsif(($distance <= $win_end) && ($end_coord > $win_end)){
		my $tmp = substr($seq,0,$win_end - $distance + 1);
		
		# dump sequence if using -seqdump and there is at least $min nt
		print OUT "$header $length nt\n$tmp\n" if (length($tmp) >= $min);
	}		
			
}
close(FILE) || die "Can't close file\n";
close(OUT);


exit(0);
