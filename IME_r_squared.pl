#!/usr/bin/perl
#
# IME_r_squared.pl
#
# A script to take a set of motifs found by Nested MICA, and test them against introns that have
# known expression values in Arabidopsis
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
my $target;     # file with sequences in which to find motif
my $verbose;    # show details for motif count/density in each sequence
my $t_min;		# minimum threshold value to use
my $t_max;		# maximum threshold value to use
my $t_step;		# threshold step value to use
my $limit;		# only show matches with r2 values above a certain limit
my $counts;		# only show motif counts
my $density;    # only show motif density
my $sum;		# only show motif sum scores
my $all;		# show all measures
my $sizes;		# calculate motif in different sized windows of sequence

GetOptions ("target=s"     => \$target,
			"verbose"      => \$verbose,
			"motif_dir=s"  => \$motif_dir,
			"t_min=f"      => \$t_min,
			"t_max=f"      => \$t_max,
			"t_step=f"     => \$t_step,
			"limit=f"	   => \$limit,
			"density"      => \$density,
			"counts"       => \$counts,
			"sum"          => \$sum,
			"all"          => \$all,
			"sizes"		   => \$sizes
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

			# final loop which will use different lengths of sequence in the target introns 
			# if -sizes is specified will check motif in first 100 bp, first 200 bp, first 300 bp and all.
			my @sizes = (100000);
			@sizes = (100, 200, 300, 400, 500, 100000)  if ($sizes);
			foreach my $j (@sizes){

				# now score motif against set of verified experimental sequences to get three r2 values
				my ($count_r2, $density_r2, $sum_r2) = &process_sequence($j, $motif);

				# don't print out results that are below $limit or which can't be measured because of divide by zero values
				print "\n" if ($verbose);
				next if ($count_r2 eq "NA");
				next if ($density_r2 eq "NA");
				next if ($sum_r2 eq "NA");
				
				# for tidyness call $j 'All' if working with $j = 100000;
				($j = "All") if ($j == 100000);

				if($count_r2 > $limit || $verbose){
					print "$count_r2\t$motif\tm${i}\t$threshold\t${j}_nt\tmotif_count\n" if ($counts);					
				}
				if($density_r2 > $limit ||$verbose){
					print "$density_r2\t$motif\tm${i}\t$threshold\t${j}_nt\tmotif_density\n" if ($density);
				}
				if($sum_r2 > $limit ||$verbose){
					print "$sum_r2\t$motif\tm${i}\t$threshold\t${j}_nt\tmotif_sum\n" if ($sum);
				}
				print "\n" if ($verbose);
				
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
	die "Need to specify both -motif_dir and -target options\n" if(!$motif_dir || !$target);
	die "Need to specify either -counts, -density, -sum, or -all\n" if (!$counts && !$sum && !$density && !$all);
	die "Don't specify -all and -counts|-density|-sum\n" if ($all && ($counts || $sum || $density));

	# check files exist
	die "$motif_dir does not seem to exist\n" if (! -e $motif_dir);
	die "$target does not seem to exist\n" if (! -e $target);

	# set threshold values if none specified
	$t_min = 0   if (!$t_min);
	$t_max = 7   if (!$t_max);
	$t_step = 1  if (!$t_step);
	$limit = 0.7 if (!$limit);

	if($all){
		$counts = 1;
		$density = 1;
		$sum = 1;
	}
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
		($motif_count++) if (m/<name>motif\d+<\/name>/);		
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
	    if(m/weight symbol=\"([a-z])[a-z]+\">(0\.\d+)<\/weight/){
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

sub process_sequence{
	
	my $limit = shift;
	my $motif = shift;
	
	# will need to track quite a few stats for each sequence
	# values for increase in expression, counts of motif, motif density, sum of motif scores (above threshold)
	my @expression;
	my @counts;
	my @density;
	my @motif_sum;

	my $limit2 = $limit;
	($limit2 = "No length limit") if ($limit == 100000);
	print "\n\n$motif T = $threshold, Limit = $limit2\n\n" if ($verbose);


	open(TARGET,"<$target") || die "Couldn't open $target file\n";

	my $fasta = new FAlite(\*TARGET);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {

		# get and trim header
		my $header = $entry->def;
		$header =~ m/.* ([\d\.]+) /;
		my $expression = $1;
		$header =~ s/ .*//;

	    my $seq = lc($entry->seq);
	    my $length = length($seq);

		# will count total amount of motif in each sequence (motifs may overlap so need to be sure we are not double counting)
		# to help this we will temporarily store a copy of $seq to mask out where any motifs are with a '-' character then
		# we can just count the dashes to know how many bases of a sequence are motif
		my $masked_seq = $seq;

		# will count how many motifs appear in each sequence	
		my $motif_count = 0;

		# want to keep track of the sum of motif scores that are above the threshold
		my $sum_motif_score = 0;

	    # loop through sequence in windows equal to motif size
	    for (my $i = 0; $i < length($seq)-$motif_length; $i++) {

			# exit loop if we are past length limit
			last if ($i > $limit);
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
				$motif_count++;
				# and mask the motif out of $masked_seq
				substr($masked_seq,$i,$motif_length) = ("-" x $motif_length);
				
				# add to sum score
				$sum_motif_score += ($score - $threshold);
			}
		}

		# count motif in sequence, add to running total
		my $motif_base_count = ($masked_seq =~ tr /-/-/);

		# if we are not scanning the entire sequence. then we need to change $length before calculating
		# motif density. Set length to be whatever the limit is
		if($limit != 100000){
			$length = $limit;
		}
		my $percent_motif = sprintf("%.3f",($motif_base_count / $length) * 100);


		# print output for -verbose mode
		if ($verbose){
			$header =~ s/ /_/g;
			printf "%-25s", "$header"; 
			printf "%-9s", "${length}_nt  ";
			print "Expression: ";
			printf "%-4s","$expression";	
			print "  ";		
			print "motif_count: " if ($counts);
			printf "%-3s","$motif_count" if ($counts);
			print "motif_density: " if ($density);  
			printf "%-8s", "$percent_motif" if ($density);
			print "sum_motif_score: " if ($sum);
			printf "%.2f","$sum_motif_score" if ($sum);
			print "\n";
		}
		
		# add these values to array for later stats analysis
		push(@expression,$expression);
		push(@counts,$motif_count);
		push(@density,$percent_motif);
		push(@motif_sum,$sum_motif_score);

	}

	close(TARGET) || die "Couldn't close $target\n";

	# can always calculate three r2 values: 
	# 1) regression of expression against motif counts
	# 2) regression of expression against motif density
	# 3) regression of expression against sum of motif scores above the threshold

	my ($count_r2, $density_r2, $sum_r2) = (0,0,0);
	
	$count_r2 = &r2(\@expression,\@counts) if ($counts);
	$density_r2 = &r2(\@expression,\@density) if ($density);
	$sum_r2 = &r2(\@expression,\@motif_sum) if ($sum);
	return($count_r2,$density_r2,$sum_r2);
}



#######################################################
#
#          C a l c u la t e   R  - s q u a r e d
#
#######################################################

sub r2{
	# two arrays of data to regress against each other
	my $x = shift;
	my $y = shift;

	# need a bunch of (fairly standard) of statistics from both arrays
	my $n       = @{$x};
	my $sum_x   = sum(@{$x});
	my $sum_y   = sum(@{$y});
	my $sum_xy  = sum(map {${$x}[$_] * ${$y}[$_]} (0..$n-1));
	my $sum_x2  = sum(map {${$x}[$_]**2} (0..$n-1));
	my $sum_y2  = sum(map {${$y}[$_]**2} (0..$n-1));

	# can't calculate r2 if sum_x or sum_y = 0; so return NA instead
	if($sum_x == 0 || $sum_y == 0){
		return("NA");
	}
	# calculate r and return r2
	else{
		my $r = (($n * $sum_xy) - ($sum_x * $sum_y))/sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2)); 	
		return(sprintf("%.4f",$r**2));			
	}
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
		($motif_count++) if (m/<name>motif(\d+)<\/name>/);		
	}
	close(MOTIF) || die "Couldn't close $motif\n";	
	
	return($motif_count);
	
}

