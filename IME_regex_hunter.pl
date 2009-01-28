#!/usr/bin/perl
#
# IME_regex_hunter.pl
#
# A script to study regular expresison densities across a set of transcripts by using a window approach
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


# This script separates just the parts of sequences that fall within a specified size category...this might mean
# removing 5' and 3' ends. Script then take only those sequences that are a defined range from the TSS in order to
# calculate motif density (using a regular expression)

# command line options
my $window;     # how big a size window
my $min;	    # alternative start point
my $max;        # how big to go up to
my $step;       # for sliding windows
my $motif;	    # define a regular expression
my $intron;     # path to input file of intron sequences

GetOptions ("window=i"     => \$window,
			"min=i"		   => \$min,
			"max=i"        => \$max,
			"step=i"       => \$step,
			"motif=s"      => \$motif,
			"intron=s"     => \$intron);


# check command line options 
die "Please specify a FASTA file of introns which contains IME information in the header\n" if (!$intron); 	

die "Please specify a valid regular expression pattern with the -motif option\n" if (!$motif);

# set some defaults
$min = 1 if (!$min);
$max = 5000 if (!$max);
$window = 200 if (!$window);
$step = 50 if (!$step);

# need to know end points of each window
my ($start,$end);

# hashes to store all sequences that fall in size range
# key to hash are start coordinates, value is sequence
my %intron2seqs;

# hashes to count how many introns and CDSs fall in each size category
my %intron2count;

print "Start,End,Seq_count,Motif_count,Motif_density\n";

for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;
#	print "$start - $end\n";

	# count how many sequences fall into each bin
	my $counter = 0;

	# store all sequences in range in one variable
	my $new_seq = "";

	open(FILE,"<$intron") || die "Couldn't open $intron\n";		

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {

		# grab basic details
		my $seq = lc($entry->seq);
		my $length = length($seq);
		my $header = $entry->def;

		# grab distance to TSS from FASTA header
		my $distance;
		$header =~ m/_i\d+_(\d+)/;
		$distance = $1;
		

		# calculate end coordinate of sequence
		my $end_coord = $distance + $length -1;


		###################################################################
		# check whether candidate sequence falls in various size categories
		###################################################################

		# CASE 1: intron sequence is wholly contained within window
		if(($distance >= $start) && ($end_coord <= $end)){
			$counter++;
			my $tmp = substr($seq,0,$end_coord-$distance+1);
			$new_seq .= "$tmp";
		}

		# CASE 2: sequence is larger than  window
		elsif(($distance < $start) && ($end_coord > $end)){
			$counter++;
			my $tmp = substr($seq,$start-$distance,$window);
			$new_seq .= "$tmp";
		}

		# CASE 3: sequence overlaps 5' edge of window
		elsif(($distance < $start) && ($end_coord >= $start)){
			$counter++;
			my $tmp = substr($seq,$start - $distance,$end_coord - $start + 1);
			$new_seq .= "$tmp";
		}

		# CASE 4: sequence overlaps 3' edge of window
		elsif(($distance <= $end) && ($end_coord > $end)){
			$counter++;
			my $tmp = substr($seq,0,$end - $distance + 1);
			$new_seq .= "$tmp";
		}		
	}

	close(FILE) || die "Can't close file\n";
	

 	# use regexes to find motif
	$motif = lc($motif);

	my $motif_count = $new_seq =~ s/($motif)/ \U$1 /g;
	($motif_count = 0) if (!$motif_count);

			
	# count how many bases are in motif, and then calculate motif density
	my $motif_sequence = ($new_seq =~ tr /A-Z/A-Z/);
	my $percent_motif = sprintf("%.3f",($motif_sequence / length($new_seq)) * 100);
	
	print "$start,$end,$counter,$motif_count,$percent_motif\n";
}
exit(0);