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
my $max;		# maximum amount of sequence to include (starting from 5' end)
my $donor;		# how much to clip from 5' end of sequence
my $acceptor;	# how much to clip from 3' end of sequence
GetOptions ("window=i"   => \$window,
			"start=i"	 => \$start,
			"min=i"      => \$min,
			"max=i"		 => \$max,
			"donor=i"    => \$donor,  
			"acceptor=i" => \$acceptor,
			);


# check command line options 
die "Specify a name of of an IME-formatted FASTA file of intron sequences\n" if (@ARGV != 1);

# set some defaults
$start = 1     if (!$start);
$window = 500  if (!$window);
$min = 50      if (!$min);
# set $max to some very large value which will never be exceeded by default (hopefully!)
$max = 1000000 if (!$max);
$donor = 5     if (!$donor);
$acceptor = 10 if (!$acceptor);

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
	my $whole_seq = $entry->seq;

	# remove donor and acceptor;
	my $seq = substr($whole_seq,$donor,-$acceptor);
	
	# trim sequence at 3' end if it's length is greater than $max
	if(length($seq) > $max){
		$seq = substr($seq,0,$max)
	}
	
	# now calculate length of trimmed sequence
	my $length = length($seq);

	
	my $header = $entry->def;
			
	# grab distance to TSS differently depending on what type of sequences we are dealing with
	my $distance;
	$header =~ m/_i\d+_(\d+)/;
	$distance = $1;
	
	# modify distance to include the offset of the donor sequence that has been removed
	$distance += $donor;
	
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
