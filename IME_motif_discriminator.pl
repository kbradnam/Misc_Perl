#!/usr/bin/perl
#
# IME_motif_discriminator.pl
#
# A script to take a set of motifs and test them against two sets of sequences, a set known to enhance, and a set known to not enhance
# Which motif gives the greatest difference in motif density between the two sets
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;
use List::Util qw(sum);

########################
# Command line options
########################

my $motif_dir;  # directory containing nested MICA *.xms files
my $plus;       # file with sequences that are known to enhance expression
my $minus;		# file with sequences that are known to not enhance expression
my $t_min;		# minimum threshold value to use
my $t_max;		# maximum threshold value to use
my $t_step;		# threshold step value to use
my $limit;		# only show matches with ratio values above a certain limit
my $sizes;		# calculate motif in different sized windows of sequence

GetOptions ("plus=s"       => \$plus,
			"minus=s"      => \$minus,
			"motif_dir=s"  => \$motif_dir,
			"t_min=f"      => \$t_min,
			"t_max=f"      => \$t_max,
			"t_step=f"     => \$t_step,
			"limit=f"	   => \$limit,
);

# are we using correct command-line options?
&pre_flight_checks;


# global variables

# log likelihood scores will be stored in an array of hashes each element of array
# corresponds to the base of the motif each value will be a key (A,C,G, or T) with
# the log likelihoods being the values
my @motif;

my $motif_length;

# limit for log-odds score for which to report motif hits
my $threshold;  


##############################################################
##############################################################
#
#    S T A R T 
#
##############################################################
##############################################################

# read directory containing motifs
opendir(DIR, "$motif_dir") || die "Can't read directory: $motif_dir\n";
my @files= readdir(DIR);

foreach my $motif (@files){

	next unless ($motif =~ m/\.xms$/);
	
	# first count how many motifs in the file
	my $num_motifs = &count_motifs($motif);

	# loop through each motif
	for (my $i = 0; $i < $num_motifs;$i++){
		
		# read motif file, fill @motif structure for motif $i
		&parse_motif($motif,$i);

		# loop through range of different thresholds
		for($threshold = $t_min; $threshold <= $t_max; $threshold += $t_step){

			# now caclulate motif count and density in plus and minus sequences
			my ($plus_motif_count, $plus_density, $plus_seq_count,$plus_zero_count) = &process_sequences($plus, $motif);
			my ($minus_motif_count, $minus_density, $minus_seq_count, $minus_zero_count) = &process_sequences($minus, $motif);

			my $plus_motifs_per_seq  = sprintf("%.2f",$plus_motif_count/$plus_seq_count);
			my $minus_motifs_per_seq = sprintf("%.2f",$minus_motif_count/$minus_seq_count);
			
			# fix possible divide by zero errors
			my $ratio;
			if ($minus_density == 0){
				$ratio = "NA";				
			}
			else{
				$ratio = sprintf("%.2f",$plus_density/$minus_density);				
			}
			
#			if((($ratio eq "NA") || ($ratio >= $limit)) && ($plus_motif_count > 0)){
			if(($ratio eq "NA") || ($ratio >= $limit)){
				print "$ratio\t${motif}_$i\t$threshold\t$plus_density\t$minus_density\t$plus_zero_count\t$plus_motifs_per_seq\t$minus_motifs_per_seq\n";		
			}
			
			
		}		
	}
}
close(DIR);
exit(0);


#############################################
#
#
#          S u b r o u t i n e s 
#
#
##############################################


sub pre_flight_checks{
	# check that valid command line options are specified
	die "Need to specify both -motif_dir and -plus and -minus options\n" if(!$motif_dir || !$plus || !$minus);

	# check files exist
	die "$motif_dir does not seem to exist\n" if (! -e $motif_dir);
	die "$plus does not seem to exist\n" if (! -e $plus);
	die "$minus does not seem to exist\n" if (! -e $minus);

	# set threshold values if none specified
	$t_min = 0   if (!$t_min);
	$t_max = 6   if (!$t_max);
	$t_step = 1  if (!$t_step);
	$limit = 0   if (!$limit);

}




##############################################################
#
# Extract data from Nested MICA output file
#
##############################################################

sub parse_motif{
	my $motif = shift;
	my $target_count = shift;
	
	my $motif_count = 0;
	
	# Need sets of expected nucleotide frequencies to compute log likelihood scores
	my %expected = ("a" => "0.2713","c" => "0.1534", "g" => "0.1701","t" => "0.4051");

	# track base position in motif
	my $pos = 0;
	my $max_pos = 0;

	$motif_length = 0;

	# open motif file and read in one motif
	open(MOTIF,"<$motif_dir/$motif") || die "Could not open $motif file\n";

	# skip to correct motif in file
	while(<MOTIF>){
		($motif_count++) if (m/<name>.*<\/name>/);		
		last if ($motif_count == $target_count);
	} 


	while(<MOTIF>){
		# keep track of motif position, need to stop if we get to the second motif
	    if (m/<column pos=\"(\d+)\"/){
			$pos = $1;	
			($max_pos = $pos) if ($pos > $max_pos);
			last if ($pos < $max_pos);
			$motif_length++;
	    }

	    # get nucleotide frequencies from input file
	    if(m/weight symbol=\"([a-z])[a-z]+\">(\-*0\.\d+)<\/weight/){
			my $base = lc($1);
			my $freq = $2;

			# if frequency is zero, set to be a very small positive value else can't take log
			($freq = 0.000001) if ($freq < 0.000001);
			
			# take logarithm of observed over expected frequency and add to @motifs
			$motif[$pos]{$base} = log($freq/$expected{$base});
	    }
	}
	close(MOTIF) || die "Couldn't close $motif\n";
	
}


##############################################################
#
# Find and score motifs in target sequence
#
##############################################################

sub process_sequences{
	 
	my $file = shift;
	my $motif = shift;
	
	# count all motifs and bases in a motif
	my $total_motif_count = 0;
	my $total_motif_base_count = 0;
	my $total_length = 0;
	my $total_seqs = 0;
	my $total_seq_length;
	# keep track of how many sequence have no motifs at all
	my $zero_count = 0;
	
	open(FILE,"<$file") || die "Couldn't open $file file\n";

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {
		
		$total_seqs++;

	    my $seq = lc($entry->seq);
		$total_seq_length += length($seq);
		
		# will count total amount of motif in each sequence (motifs may overlap so need to be sure we are not double counting)
		# to help this we will temporarily store a copy of $seq to mask out where any motifs are with a '-' character then
		# we can just count the dashes to know how many bases of a sequence are motif
		my $masked_seq = $seq;

		# will count how many motifs appear in each sequence	
		my $seq_motif_count = 0;

		# loop through sequence in windows equal to motif size
	    for (my $i = 0; $i < length($seq)-$motif_length; $i++) {

			# extract a window of sequence, split it, and place in array	
			my @sequence = split(//,substr($seq, $i, $motif_length));

			# Calculate motif score: only A,T,C,G bases counts towards score 
			my $score = 0;

			for(my $j = 0; $j<@sequence;$j++){
				($score += $motif[$j]{$sequence[$j]}) if ($sequence[$j] =~ m/[atcg]/);
			}
			$score = sprintf("%.2f",$score);	

			# if we've found some sequence that is above the threshold...
			if($score > $threshold){
				# count motif
				$total_motif_count++;
				$seq_motif_count++;
				# and mask the motif out of $masked_seq
				substr($masked_seq,$i,$motif_length) = ("-" x $motif_length);
			}
		}
		$zero_count++ if ($seq_motif_count == 0);
		
		# count motif in sequence, add to running total
		$total_motif_base_count += ($masked_seq =~ tr /-/-/);
	}

	close(FILE) || die "Couldn't close $file\n";

	my $total_motif_density = sprintf("%.2f",($total_motif_base_count / $total_seq_length) * 100);

	return($total_motif_count,$total_motif_density,$total_seqs,$zero_count);
}


###################################################################################
#
# Quick routine to just count motifs in the motif file
#
###################################################################################

sub count_motifs{
	my $motif = shift;
	my $motif_count = 0;
	
	# open motif file and read in one motif
	open(MOTIF,"<$motif_dir/$motif") || die "Could not open $motif file\n";

	while(<MOTIF>){
		($motif_count++) if (m/<name>.*<\/name>/);		
	}
	close(MOTIF) || die "Couldn't close $motif\n";	
	
	return($motif_count);
	
}

