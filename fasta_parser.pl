#!/usr/bin/perl
#
# fasta_parser.pl
#
# A script to extract the same window of sequence from multiple sequences and write to new output files
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use Getopt::Long;

# This script windows through a set of sequences and adds just the specified region to a new file. The result
# is a set of files with DNA from 1-100 nt, 51-150, nt etc. Can work from 5' and/or 3' ends of the sequence

# command line options
my $window;      # how big a size window
my $min;	     # alternative start point
my $max;         # how big to go up to
my $step;        # for sliding windows
my $file;        # path to input file of sequences
my $five_prime;  # extract sequences from 5' end
my $three_prime; # extract sequences from 5' end
my $seqlength;	 # threshold for dumping sequences to file (don't want short sequences)

GetOptions ("window=i"    => \$window,
			"min=i"		  => \$min,
			"max=i"       => \$max,
			"step=i"      => \$step,
			"file=s"      => \$file,
			"five_prime"  => \$five_prime,
			"three_prime" => \$three_prime,
			"seqlength=i" => \$seqlength);


# check command line options 
if (!$file){
	die "Please specify a valid FASTA sequence file using the -file option\n" 	
}
if (!$five_prime && !$three_prime){
	die "Please specify -five_prime and/or -three_prime options\n" 	
}

# set some defaults
$min = 1 if (!$min);
$max = 1000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);
$seqlength = 25 if (!$seqlength);

print "\nStarting from position $min\n";
print "Only considering sequence up to position $max\n";
print "Using window size of $window nt\n";
print "Using step size of $step nt\n";
print "Will ignore any window with less than $seqlength nt\n\n";

# need to know end points of each window
my ($start,$end);


for(my $start = $min;$start+$window<$max; $start+= $step){
	$end = $start + $window -1;
	print "$start - $end\n";

	&process_sequence($file,$start,$end,$window);
}

# main subroutinee to loop through file of sequences and extract only sequence in certain range

sub process_sequence{

	my $file = shift;
	my $win_start = shift;
	my $win_end = shift;
	my $window = shift;
	
	# count how many sequences fall into each bin
	my $counter = 0;
	
	# store all sequences in range in one variable
	my $new_seq = "";
	
	open(FILE,"<$file") || die "Couldn't open $file\n";		

	my $fasta = new FAlite(\*FILE);

	if($five_prime){
		open(FOR,">${file}.${win_start}_${end}.fa") || die "couldn't create ${file}.${start}_${$end}.fa\n";		
	}
	if($three_prime){
		open(REV,">${file}.${win_start}_${end}.rev.fa") || die "couldn't create ${file}.${start}_${$end}.rev.fa\n";		
	}

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {
	    
		# grab basic details
		my $seq = $entry->seq;
		my $length = length($seq);
		my $header = $entry->def;
	
		# now extract relevant piece of sequence
		my $tmp_seq;	 

		# skip sequences that are shorter than start of current window position
		next if ($length <= $win_start);

		# treat 5' and 3' ends separately
		if($five_prime){
			# two conditions to consider
			if($window >= $length){
				$tmp_seq = substr($seq,$win_start-1,$length);
			}
			else{
				$tmp_seq = substr($seq,$win_start-1,$window);
			}
			# dump sequence if using -seqdump and there is at least 10 nt
			print FOR "${header}_${win_start}_${win_end}\n$tmp_seq\n" if (length($tmp_seq) >= $seqlength);
		}
		
		if($three_prime){
			if($window >= $length){
				$tmp_seq = substr($seq,$win_start,$length);
			}
			else{
				$tmp_seq = substr($seq,-$win_end,$window);
			}
			# dump sequence if using -seqdump and there is at least 10 nt
			print REV "${header}_${win_start}_${win_end}\n$tmp_seq\n" if (length($tmp_seq) >= $seqlength);
		}			
		
	}
	close(FILE) || die "Can't close file\n";
	close(FOR) || die "Can't close FOR output file handle\n" if ($five_prime);
	close(REV) || die "Can't close REV output file handle\n" if ($three_prime);
}



exit(0);
