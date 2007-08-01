#!/usr/bin/perl
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
use warnings;
use Getopt::Long;
use FAlite;

########################
# Command line options
########################

my $motif;      # nested MICA *.xms file with (single) motif
my $target;     # file with sequences in which to find motif
my $threshold;  # limit for log-odds score for which to report motif hits
my $scores;     # Show scores
my $seqs;       # Show motif sequences in output (one sequence per motif in each intron)
my $species;    # code to determine which species to use expected frequencies from
my $mdensity;	# count (and show) amount and percentage of motif in each sequence (one line per sequence)
my $mseqs;		# just show sequences of each intron that have motifs above threshold (one sequence per intron)
my $msummary;	# show motif count and percentage for all sequences combined

my $stats;      # report stats on all log likelihood scores found
my $min;        # Set minimum cut off for use with -stats option
my $max;        # Set maximum cut off for use with -stats option
my $interval;   # Set interval size for use with -stats option

GetOptions ("motif=s"    => \$motif,
	       "target=s"    => \$target,
	       "threshold=f" => \$threshold,
	       "scores"      => \$scores,
	       "seqs"        => \$seqs,	    	    
		   "species=s"   => \$species,
		   "mdensity"    => \$mdensity,  
		   "mseqs"		 => \$mseqs,
		   "msummary"    => \$msummary,
		   "stats"       => \$stats,
	       "min=f"       => \$min,
	       "max=f"       => \$max,
	       "interval=f"  => \$interval,
);

# are we using correct command-line options?
&pre_flight_checks;

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
# set of frequencies chosen by -species option
# following may have to be tidied up if I need to add more species
# Frequencies were calculated from chromosomes or various TAIR 7 downloads (apart from Introns which were from Tali)
my %expected;

if($species =~ m/cei/i){
	%expected = ("a" => "0.33339","c" => "0.16194", "g" => "0.16021","t" => "0.34446");
}
elsif($species =~ m/ceg/i){
	%expected = ("a" => "0.32280","c" => "0.17733","g" => "0.17709","t" => "0.32279");
}
elsif($species =~ m/^ati$/i){
	%expected = ("a" => "0.2769","c" => "0.1575", "g" => "0.1587","t" => "0.4069");
}
elsif($species =~ m/^atir$/i){
	%expected = ("a" => "0.4069","c" => "0.1587", "g" => "0.1575","t" => "0.2769");
}
elsif($species =~ m/atg/i){
	%expected = ("a" => "0.3195","c" => "0.1800", "g" => "0.1798","t" => "0.3192");
}
elsif($species =~ m/at5u/i){
	%expected = ("a" => "0.3016","c" => "0.2181", "g" => "0.1594","t" => "0.3210");
}
elsif($species =~ m/at3u/i){
	%expected = ("a" => "0.2960","c" => "0.1473", "g" => "0.1732","t" => "0.3835");
}
elsif($species =~ m/atig/i){
	%expected = ("a" => "0.3437","c" => "0.1568", "g" => "0.1562","t" => "0.3433");
}
elsif($species =~ m/atc/i){
	%expected = ("a" => "0.2868","c" => "0.2031", "g" => "0.2387","t" => "0.2713");
}
elsif($species =~ m/atu/i){
	%expected = ("a" => "0.3342","c" => "0.1685", "g" => "0.1651","t" => "0.3321");
}
elsif($species =~ m/atd/i){
	%expected = ("a" => "0.3211","c" => "0.1787", "g" => "0.1724","t" => "0.3278");
}
elsif($species =~ m/dmi/i){
	%expected = ("a" => "0.2942","c" => "0.2040", "g" => "0.1979","t" => "0.3040");
}
else{
	die "\'$species\' is not a valid species code.\n";
}


# track base position in motif
my $pos = 0;
my $max_pos = 0;
my $motif_size;

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
		my $base = lc($1);
		my $freq = $2;
	
		# take logarithm of observed over expected frequency and add to @motifs
		$motif[$pos]{$base} = log($freq/$expected{$base});
    }
}
close(MOTIF) || die "Couldn't close $motif\n";

# calculate size of motif
$motif_size = $max_pos+1;

##############################################################
#
#
# P A R T   II   - Find and score motifs in target sequence
#
#
##############################################################


open(TARGET,"<$target") || die "Couldn't open $target file\n";

my $fasta = new FAlite(\*TARGET);

# keep track of total length of sequence in each file and total length of motif, 
# total number of motifs and total number of sequences
my ($total_length,$total_motif,$seq_count,$total_motif_count);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	$seq_count++;
	
	# get and trim header
	my $header = $entry->def;
    $header =~ s/ .*//;

    my $seq = lc($entry->seq);
    my $length = length($seq);
	$total_length += $length;
	
	# will count total amount of motif in each sequence (motifs may overlap so need to be sure we are not double counting)
	# to help this we will temporarily store a copy of $seq to mask out where any motifs are with a '-' character then
	# we can just count the dashes to know how many bases of a sequence are motif
	my $masked_seq = $seq;
	
	# flag variable to know whether intron contained at least one motif above threshold
	my $above_threshold = 0;
	
    # loop through sequence in windows equal to motif size
    for (my $i = 0; $i < length($seq)-$motif_length; $i++) {
		# extract a window of sequence, split it, and place in array	
		my $window = substr($seq, $i, $motif_length);
		my @sequence = split(//,$window);
		
		# Calculate motif score: only A,T,C,G bases counts towards score 
		my $score = 0;
		for(my $j = 0; $j<@sequence;$j++){
			($score += $motif[$j]{$sequence[$j]}) if ($sequence[$j] =~ m/[atcg]/);
		}
		$score = sprintf("%.2f",$score);
		
		# Add score to @all_scores array if tracking stats
		push(@all_scores,$score) if ($stats);		
		
		
		# if we've found some sequence that is above the threshold...
		if($score > $threshold){
			# note that we are above the threshold
			$above_threshold = 1;
			
			# count motif
			$total_motif_count++;
			
			# and mask the motif out of $masked_seq
			substr($masked_seq,$i,$motif_size) = ("-" x $motif_length);
			
			# print out current motif details if -scores or -seqs specified
			if($scores){
				my $start_coord = $i+1;
			  	print "$header $score $start_coord/$length $window\n";
			}
			if($seqs){
				my $highlight = uc($window);
				my $new_seq = $seq;
				$new_seq =~ s/$window/ $highlight /g;
				print "$new_seq\n" if ($seqs);
			}
		}
	}
	
	# count motif in sequence, add to running total
	my $motif_count = ($masked_seq =~ tr /-/-/);
	$total_motif += $motif_count;
	my $percent_motif = sprintf("%.3f",($motif_count / $length) * 100);
	print "$header motif_density: $motif_count/$length $percent_motif%\n" if ($mdensity);
	
	# print out intron sequence if -mseqs specified and intron contains a motif above threshold
	print "$header\n$seq\n" if ($mseqs && $above_threshold);
}

close(TARGET) || die "Couldn't close $target\n";

# print motif summary if requested
if($msummary){
	my $percent_motif = sprintf("%.3f",($total_motif/$total_length) * 100);
	print "\nSUMMARY:\n"; 
	print "number_of_sequences: $seq_count total_sequence_length: $total_length\n";
	print "number_of_motifs: $total_motif_count total_motif_length: $total_motif motif_density: $percent_motif%\n\n\n";
}


#########################################
#
# Process stats if required
#
#########################################

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
		print "ERROR! $score lies outside min ($min) and max ($max) boundaries\n" if ($flag == 0);
    }

    foreach my $key (sort {$limits{$a} <=> $limits{$b}}(keys(%limits))){

		# set upper limit to be slightly less
		my $upper = $limits{$key} - 0.001;
	
		# print counts (if they exist), else print zero
		($counts{$key} = 0) if (!defined($counts{$key}));
		print "$key,$upper,$counts{$key}\n";   
    }

}
exit(0);


#############################################
#
#
#          S u b r o u t i n e s 
#
#
##############################################

sub pre_flight_checks{
	# check that both command line options are specified
	die "Need to specify both -motif and -target options\n" if(!$motif || !$target);

	# check that motif file looks like a valid file
	die "-motif option must specify a valid Nested MICA *.xms output file\n" if($motif !~ m/\.xms$/);

	# check files exist
	die "$motif does not seem to exist\n" if (! -e $motif);
	die "$target does not seem to exist\n" if (! -e $target);

	# check that species code has been chosen
	if(!$species){
		print "\nPlease specify a suitable species code using the -species option.\n";
		print "Species codes are required to determine the correct expected nucleotide\nfrequencies when scoring motifs\n\n";
		print "Current options (all case-insensitive) are:\n";
		print "AtI  - Arabidopsis thaliana introns\n";
		print "AtIR - Arabidopsis thaliana introns (reverse complemented)\n";
		print "AtG  - Arabidopsis thaliana genomic\n";
		print "AtC  - Arabidopsis thaliana CDS\n";
		print "AtIG - Arabidopsis thaliana intergenic\n";
		print "At5U - Arabidopsis thaliana 5' UTR (exons)\n";
		print "At3U - Arabidopsis thaliana 3' UTR (exons)\n";
		print "AtU  - Arabidopsis thaliana upstream region of genes (1000 bp 5' to transcript)\n";
		print "AtD  - Arabidopsis thaliana downstream region of genes (1000 bp 3' to transcript)\n";
		print "CeI - Caenorhabditis elegans introns\n";
		print "CeG - Caenorhabditis elegans genomic\n";
		print "DmI - Drosophila melanogaster introns\n\n";
		die "Choose one option only.\n\n";
	}

	# have we chosen an option to print some output?
	if (!$seqs && !$stats && !$scores && !$mdensity && !$msummary && !$mseqs){
		die "You have to choose at least one of the following options:\n-stats, -seqs, -scores, -mdensity, -msummary, -mseqs\n";	
	}
	# can't choose mdensity with other options
	if($mdensity && ($seqs || $stats || $mseqs)){
		die "Don't choose -mdensity with -seqs or -scores or -mseqs\n";	
	}  
	# can't choose mseqs with other options
	if($mseqs && ($seqs || $stats || $mdensity)){
		die "Don't choose -mseqs with -seqs or -scores or -mdensity\n";	
	}  	
	# don't specify a threshold if using -stats
	if($threshold && $stats){
		die "Don't specify -threshold option when using -stats option\n";
	}
}