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
my $score_only; # Only show DNA scores in output rather than sequence

GetOptions ("motif=s"     => \$motif,
	    "target=s"    => \$target,
	    "threshold=f" => \$threshold,
	    "scores"      => \$score_only);

# check that both command line options are specified
die "Need to specify both -motif and -target options\n" if(!$motif || !$target);

# check that motif file looks like a valid file
die "-motif option must specify a valid Nested MICA *.xms output file\n" if($motif !~ m/\.xms$/);

# check files exist
die "$motif does not seem to exist\n" if (! -e $motif);
die "$target does not seem to exist\n" if (! -e $target);

# set threshold if none specified
$threshold = 0 if (!$threshold);


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

    # need to 
    my $max_score=0;
    
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
	  # only want to print out scores above threshold
	  my $new_score = sprintf("%.2f",$score);

	  # want to show location of motif within original sequence
	  my $highlight = uc($window);
	  my $new_seq = $seq;
	  $new_seq =~ s/$window/$highlight/g;
	  if($score > $threshold){
	      my $start_coord = $i+1;
	      print "$header $new_score $start_coord/$length $window\n";
	      unless($score_only){
		  print "$new_seq\n";
	      }
	  }
      }
}

close(TARGET) || die "Couldn't close $target\n";

exit(0);
