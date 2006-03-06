#!/usr/bin/perl -w
#
# where_is_motif.pl
#
# A script to take motifs found by Nested MICA, convert the frequencies of
# each base in the motif to log likelihood scores and then find the locations
# (and scores) of those motifs in the original set of sequences used to
# generate the motifs.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use lib "/Korflab/lib/perl";
use Getopt::Long;
use FAlite;

########################
# Command line options
########################

my $motif;      # nested MICA *.xms file with (single) motif
my $target;     # file with sequences in which to find motif
my $threshold;  # limit for which to report motif hits
my $scores;     # Show scores
my $seqs;       # Show motif sequences in output
my $stats;      # report stats on log likelihood scores
my $min;        # Set minimum cut off for use with -stats option
my $max;        # Set maximum cut off for use with -stats option
my $interval;   # Set interval size for use with -stats option


GetOptions ("motif=s"     => \$motif,
	    "target=s"    => \$target,
	    "threshold=f" => \$threshold,
	    "scores"      => \$scores,
	    "seqs"        => \$seqs,	    
	    "stats"       => \$stats,
	    "min=f"       => \$min,
	    "max=f"       => \$max,
	    "interval=f"  => \$interval,
	    );

# check that both command line options are specified
die "Need to specify both -motif and -target options\n" if(!$motif || !$target);

# check that motif file looks like a valid file
die "-motif option must specify a valid Nested MICA *.xms output file\n" if($motif !~ m/\.xms$/);

# check files exist
die "$motif does not seem to exist\n" if (! -e $motif);
die "$target does not seem to exist\n" if (! -e $target);

# set threshold if none specified
$threshold = 0 if (!$threshold);

# keep track of scores if making stats
my @all_scores if ($stats);
    


##############################################################
#
#
# P A R T   I   - Extract data from Nested MICA output file
#
#
##############################################################


# log likelihood scores will be stored in an array of hashes each element of array
# corresponds to the base of the motif each value will be a key (A,C,G, or T) with
# the log likelihoods being the values
my @motif;
# will need to know motif length later on
my $motif_length;


# Need sets of expected nucleotide frequencies to compute log likelihood scores
# one set is based on (confirmed) intron base frequencies, the second is based on 
# the entire worm genome
my %expected_intronic = ("A" => "0.33339","C" => "0.16194",
			 "G" => "0.16021","T" => "0.34446");
my %expected_genomic  = ("A" => "0.32280","C" => "0.17733",
			 "G" => "0.17709","T" => "0.32279");


# track base position in motif
my $pos =0;
my $max_pos = 0;


# open motif file and read in one motif
open(MOTIF,"<$motif") || die "Could not open $motif file\n";

while(<MOTIF>){
    # keep track of motif position, need to stop if we get to the second motif
    if (m/<column pos=\"(\d+)\"/){
	$pos = $1;	
	($max_pos = $pos) if ($pos > $max_pos);
	last if ($pos < $max_pos);
	$motif_length++;
    }

    # get nucleotide frequencies from input filre
    if(m/weight symbol=\"([a-z])[a-z]+\">(0\.\d+)<\/weight/){
	my $base = uc($1);
	my $freq = $2;
	
	# take logarithm of observed over expected frequency and add to @motifs
	$motif[$pos]{$base} = log($freq/$expected_intronic{$base});
    }
}
close(MOTIF) || die "Couldn't close $motif\n";



##############################################################
#
#
# P A R T   II   - Find and score motifs in target sequence
#
#
##############################################################


open(TARGET,"<$target") || die "Couldn't open $target file\n";

my $fasta = new FAlite(\*TARGET);

# loop through each sequence in target file

while(my $entry = $fasta->nextEntry) {
    my $header = $entry->def;
    # trim header
    $header =~ s/ .*//;

    my $seq = $entry->seq;
    my $length = length($seq);

    # loop through sequence in windows equal to motif size
      for (my $i = 0; $i < length($seq)-$motif_length; $i++) {

	  # extract a window of sequence, split it, and place in array	
	  my $window = substr($seq, $i, $motif_length);
	  my @sequence = split(//,$window);

	  my $score =0;
	  for(my $j =0; $j<@sequence;$j++){
	      my $base = uc($sequence[$j]);
	      $score += $motif[$j]{$base};
	  }
	  # only want to print out scores above threshold, add scores to array if
	  # tracking stats
	  my $new_score = sprintf("%.2f",$score);
	  push(@all_scores,$new_score) if ($stats);

	  # print output if requested
	  if($scores || $seqs){
	      # want to show location of motif within original sequence
	      my $highlight = uc($window);
	      my $new_seq = $seq;
	      $new_seq =~ s/$window/$highlight/g;
	      if($score > $threshold){
		  my $start_coord = $i+1;
		  print "$header $new_score $start_coord/$length $window\n" if ($scores);
		  print "$new_seq\n" if ($seqs);	      
	      }
	  }
      }
}

close(TARGET) || die "Couldn't close $target\n";



# Process stats if required

if($stats){
    # structure to count log likelihood in different intervals
    my %counts;

    # set up what the bin sizes are going to be for counting if not specified on command line
    # need min, max, and interval settings, store details in %limits
    ($min = -40)    if (!$min);
    ($max = 10)     if (!$max);
    ($interval = 1) if (!$interval);

    my %limits;
    
    for(my $i=$min; $i<=$max; $i+=$interval){
	$limits{$i} = $i+$interval;
    }

    # work through all scores in @all_scores
    
    OUTER: while (@all_scores){
	# need to make sure that no score is outside min and max range
	my $flag = 0;
	my $score = shift(@all_scores);

	# now loop through all possible limit categories
	foreach my $key (sort {$limits{$a} <=> $limits{$b}}(keys(%limits))){
	    if(($score >= $key) && ($score < $limits{$key})){
		$counts{$key}++;
		$flag = 1;
		next OUTER;
	    }
	}
	# warn if any score is outside min and max
	print "ERROR! $score lies outside min ($min) and max ($max) boundaries\n" if ($flag ==0);
    }

    foreach my $key (sort {$limits{$a} <=> $limits{$b}}(keys(%limits))){

	# set upper limit to be slightly less
	my $upper = $limits{$key}-0.001;
	
	# print counts (if they exist), else print zero
	($counts{$key} = 0) if (!defined($counts{$key}));
	print "$key,$upper,$counts{$key}\n";
	   
    }

}
exit(0);
