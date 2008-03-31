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
use List::Util qw(sum);

########################
# Command line options
########################

my $target;     # file with sequences in which to find motif
my $verbose;    # show details for motif count/density in each sequence


GetOptions ("target=s"     => \$target,
			"verbose"      => \$verbose,
);

# are we using correct command-line options?
&pre_flight_checks;


# global variables

# motifs will be stored in an array of array of hashes
# 1st array is motif number, 2nd array is position within motif
# final hash has keys which are A,C,G, or T, and values which are nucleotide frequencies
my @motifs;

# how many cycles to go through
my $cycles = 1000;

# how many motifs to generate
my $num_motifs = 100;

# key is motif number, value is r2 score of motif
my %motif2score;

##############################################################
##############################################################
#
#    S T A R T 
#
##############################################################
##############################################################

my @seqs;
my @expression_scores;

&read_sequences;

&generate_motifs;

for(my $i = 0;$i<$cycles;$i++){
	
	&score_motifs;
	&cull_motifs;
	&reproduce;
	&mutate;
	&grow;
	&shrink;
	&recombine;
	
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
	die "Need to specify -target options\n" if(!$target);

	# check files exist
	die "$target does not seem to exist\n" if (! -e $target);
}



# simple routine to generate a small number of motifs
# default length is 6 nt, all base frequencies chosen randomly

sub generate_motifs{

	for(my $i = 0; $i < $num_motifs; $i++){

		# assume all motifs are 6 nt to start with
		for (my $j = 0; $j < 6; $j++){

			# set random frequencies for each base (they must sum to 1.0)
			my $a = rand(1);
			my $c = rand(1) * (1 - $a);
			my $g = rand(1) * (1 - $a - $c);
			my $t = 1 - $a - $c -$g;
			
			$motifs[$i][$j]{"A"} = $a;
			$motifs[$i][$j]{"C"} = $c;
			$motifs[$i][$j]{"G"} = $g;
			$motifs[$i][$j]{"T"} = $t;
		}
	}
}


# calculate r2 score for all motifs
sub score_motifs{

	# Need sets of expected nucleotide frequencies to compute log likelihood scores
	my %expected = ("a" => "0.2713","c" => "0.1534", "g" => "0.1701","t" => "0.4051");

	# will need to track quite a few stats for each sequence
	# values for increase in expression, counts of motif, motif density
	my @expression;
	my @counts;


	for(my $i = 0; $i < $num_motifs; $i++){

		# get motif length
		my $motif_length = @{$motifs[$i]};
		print "Length is $motif_length\n";
				
		my $motif_score;
		
		# loop through each of the 12 experimental introns, scoring them for motif		
		for (my $intron = 0; $intron < 12; $intron++){
			
			
		}
		


			    my $length = length($seq);

				# will count total amount of motif in each sequence (motifs may overlap so need to be sure we are not double counting)
				# to help this we will temporarily store a copy of $seq to mask out where any motifs are with a '-' character then
				# we can just count the dashes to know how many bases of a sequence are motif
				my $masked_seq = $seq;

				# will count how many motifs appear in each sequence	
				my $motif_count = 0;

			    # loop through sequence in windows equal to motif size
			    for (my $i = 0; $i < length($seq)-6; $i++) {

					# exit loop if we are past length limit
					last if ($i > $limit);

					# extract a window of sequence, split it, and place in array	
					my @sequence = split(//,substr($seq, $i, 6));

					# Calculate motif score: only A,T,C,G bases counts towards score 
					my $score = 0;
					for(my $j = 0; $j<@sequence;$j++){
		#				($score += $motif[$j]{$sequence[$j]}) if ($sequence[$j] =~ m/[atcg]/);
					}
					$score = sprintf("%.2f",$score);	

					# if we've found some sequence that is above the threshold...
					if($score > 3){
						# count motif
						$motif_count++;

						# and mask the motif out of $masked_seq
						substr($masked_seq,$i,6) = ("-" x 6);
					}
				}

				# count motif in sequence, add to running total
				my $motif_base_count = ($masked_seq =~ tr /-/-/);

				# if we are not scanning the entire sequence. then we need to change $length before calculating
				# motif density. Set length to be whatever the limit is
				if($limit != 10000){
					$length = $limit;
				}
				my $percent_motif = sprintf("%.3f",($motif_base_count / $length) * 100);
				print "$header\t${length}_nt\tExpression\t$expression\tmotif_count\t$motif_count\tmotif_density\t$percent_motif\n" if ($verbose);

				# add these values to array for later stats analysis
				push(@expression,$expression);
				push(@counts,$motif_count);
				push(@density,$percent_motif);

			}

			close(TARGET) || die "Couldn't close $target\n";

			# can always calculate two r2 values, one for regression of expression against motif counts
			# and one for regression of expression against motif density

			my $count_r2 = &r2(\@expression,\@counts);
			my $density_r2 = &r2(\@expression,\@density);

			return($count_r2,$density_r2);
	
	}

}


##############################################################
#
# Extract data from Nested MICA output file
#
##############################################################

sub parse_motif{
	
	my $motif = shift;
	my $count = shift;
	
	# Need sets of expected nucleotide frequencies to compute log likelihood scores
	my %expected = ("a" => "0.2713","c" => "0.1534", "g" => "0.1701","t" => "0.4051");

	# track base position in motif
	my $pos = 0;
	my $max_pos = 0;

	my $motif_length = 0;


		
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

		# take logarithm of observed over expected frequency and add to @motifs
#		$motif[$pos]{$base} = log($freq/$expected{$base});
    }
	
}


##############################################################
#
# Find and score motifs in target sequence
#
##############################################################

sub process_sequence{
	
	my $limit = shift;
	
	# will need to track quite a few stats for each sequence
	# values for increase in expression, counts of motif, motif density
	my @expression;
	my @counts;
	my @density;

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

	    # loop through sequence in windows equal to motif size
	    for (my $i = 0; $i < length($seq)-6; $i++) {

			# exit loop if we are past length limit
			last if ($i > $limit);
			
			# extract a window of sequence, split it, and place in array	
			my @sequence = split(//,substr($seq, $i, 6));

			# Calculate motif score: only A,T,C,G bases counts towards score 
			my $score = 0;
			for(my $j = 0; $j<@sequence;$j++){
#				($score += $motif[$j]{$sequence[$j]}) if ($sequence[$j] =~ m/[atcg]/);
			}
			$score = sprintf("%.2f",$score);	

			# if we've found some sequence that is above the threshold...
			if($score > 3){
				# count motif
				$motif_count++;

				# and mask the motif out of $masked_seq
				substr($masked_seq,$i,6) = ("-" x 6);
			}
		}

		# count motif in sequence, add to running total
		my $motif_base_count = ($masked_seq =~ tr /-/-/);

		# if we are not scanning the entire sequence. then we need to change $length before calculating
		# motif density. Set length to be whatever the limit is
		if($limit != 10000){
			$length = $limit;
		}
		my $percent_motif = sprintf("%.3f",($motif_base_count / $length) * 100);
		print "$header\t${length}_nt\tExpression\t$expression\tmotif_count\t$motif_count\tmotif_density\t$percent_motif\n" if ($verbose);

		# add these values to array for later stats analysis
		push(@expression,$expression);
		push(@counts,$motif_count);
		push(@density,$percent_motif);

	}

	close(TARGET) || die "Couldn't close $target\n";

	# can always calculate two r2 values, one for regression of expression against motif counts
	# and one for regression of expression against motif density

	my $count_r2 = &r2(\@expression,\@counts);
	my $density_r2 = &r2(\@expression,\@density);

	return($count_r2,$density_r2);
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
	else{
		my $r = (($n * $sum_xy) - ($sum_x * $sum_y))/sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2)); 	
		return(sprintf("%.4f",$r**2));			
	}
	# calculate r and return r2
}



sub read_sequences{
	@expression_scores = qw (12.3 12.3 6.9 4.3 4.1 4.0 1.9 1.8 1.4 1.2 1.1 0.6);

	$seqs[0] = "GTAAATTTCTGTGTTCCTTATTCTCTCAAAATCTTCGATTTTGTTTTCGTTCGATCCCAATTTCGTATATGTTCTT
	TGGTTTAGATTCTGTTAATCTTAGATCGAAGACGATTTTCTGGGTTTGATCGTTAGATATCATCTTAATTCTCGATTAGGGTTTCATA
	GATATCATCCGATTTGTTCAAATAATTTGAGTTTTGTCGAATAATTACTCTTCGATTTGTGATTTCTATCTAGATCTGGTGTTAGTTT
	CTAGTTTGTGCGATCGAATTTGTCGATTAATCTGAGTTTTTCTGATCTGCAG";

	$seqs[1] = "GTAATCAATTCTCCCTCTCTATCTATGTTTGTTTGAATTCTCTCTCGCATAGTTAAGATTCCTTTTTTCGTATTCT
	AGATCCATAGAATTTATCCAAAATTCATGAATTGTTTCTAAGACACGAAACGGTTTAAGTTCAGGTCATAGTTTTT
	CTGTAGATCTCGATTTACGTGAAAGTTTACTTACCAATAGATCTGAATTATCGAAATTGCAGTTCCTTTTCCTCGA
	GTGTCTCGTTCGGCTTCATGTCCTGTGGAATTTTTAATCTTTGTCCGATTCTGAATCCGGAAATTGTTAGGGATTT
	GGGTTTTATTCAAGATTTGTCATCGCTGTGAAAGTTTTCCCTTTTTTTATGTGGGTTCGAAGTTTTCTGAAAATCT
	CAATTAGTAAAAGGATCACTGGACTTGCTTATATTAATTCTACACCGTCTCTTCATAGATTTGGTCAATCCTGTAT
	CTCAGATCTTATTTGTTCATGTGATAGCCTTTAAATGTGTGAACTTTTCATGTACCTGCAG";

	$seqs[2] = "GTAAGTCTCGATGTGAATTATGCGATTGACTATCGATTTAGGGGATCCTCATGATCTAATTAGGTCGTACGGATTG
	ACTTAGAATTGGGAAATTTTGAATCTAGTTTTGGTAGATTATAGAATTTGTGACCTAGCTAGAACATAACAGAGGA
	TTATTCAAGAACCGCTTAGGAAAGTATCTCATAATGAACATCTTTAGATTGGTCTGAGATCTAAAATTGATTCATT
	TGCTCTTTTCTTGGTTGCTGCAG";

	$seqs[3] = "GTAAAGCCTCGATTTTTGGGTTTAGGTGTCTGCTTATTAGAGTAAAAACACATCCTTTGAAATTGTTTGTGGTCAT
	TTGATTGTGCTCTTGATCCATTGAATTGCTGCAG";

	$seqs[4] = "GTAAGATTATCTCTTCCCAAAATTGATTACTTTTATTATTGAACAATTATTAACCAATCATGGCTTAACGAACTGC
	AG";

	$seqs[5] = "GTTACTTTCACTCTGTCTCATCTCTGTGTTTGATTTGCTGATACTTCCTATTGCTTGTTTAAGCGTTTTATATAAT
	CTATACAACCAGAATTTGATCCTTGAAGTGTTTTGCGTGCTTTGTGTGGTAATTGATTGGAACAGAAGATCCAGTG
	AGATTTGAAAAAAAATTGAGACTTAAGTTGGCTTTAGCTAGTTTTGAAGTTTGAACTTGTGTGGTTTGTCTGCAG";

	$seqs[6] = "GTAAGTCTACATTCTTTCTTCTTTTAGTATCTTGCCTCATAAGTAAGGATCTTAGCAGGCAATGTTTATGGTATAC
	TATATTAGTATAGATTTTAGTGGAAATATGTTTGTTTTGAACTTATTTTATGATCATATTTGACTATTATCAAAGA
	TAAAGATTCATATACCGTACATTATATATCTCTATTTTTCTAGTTTACATGTATAGCTCTAAGTTTATTTGATGAT
	TCTGTTGACTACTTTTGGATATGTGTTTTGAAACCTTTGATAAATACTAAAATAAATTTTAATTTGAAAATGCTGC
	AG";

	$seqs[7] = "GTAAGTCATGCATCCACGGAGAAACTTTCTTTTATATATGTTTATATTTTTCATGTTTAATTTGTTTTAGCACACC
	CGTTTTTAACAATGATTTTAAGATGGAATGATAAGAGCTTTCATTCAAAAAAAAAAAATGAAAATGAAAACTTTAA
	CATGTAATTTATTAATGTCTTTCCAACTGTTTAGTAGAATTTTGTACCTTCCATGAATCATAAATTAATTATTTAT
	ATCTTTCTTGATCTTTCCTGCAG";

	$seqs[8] = "GTATGGTTCATCAACTCTTTCCATTTAATCGTAATGTTGGATCATGATCATCTTCAAAGCAATGAAATGACTAACA
	CAAGTCCTTGATTACTTTCTGCAG";

	$seqs[9] = "GTAATCTCTCTCTGTGTTGCACGTACATGGCTCCTTTGATATTATACGGAAAATCATATATAGCCTAAAGATATAT
	CTACAGCTGAAACCCAATATAGGTTCTCAAGAATCTCAACCAGTGGTAGATTCATTAATAACCATATAAAAATATT
	TTAATCATATAATCGAGTCTGATTGAGTAATTTCTGTTACAAGTTAAATATTAAAGTTTTATCTCACAGATAGTAA
	TGTAACAAATTATATCATTCAAAACTCAAAACCCTTAATTTGTATTATTTTTCTGCAG";

	$seqs[10] = "GTAATTTCACTTCAATATTATATTAGAAGTCACATATTTCCATATAGAAATGTGCAATCATATTCAAATCATAGTG
	GTGATTATATAAAAGATCTATGATGAACTAATAACGTTTAATTAATAAAACTAATACACATTTAAGGCTAGTACAA
	AATAAATACATGTAAATTAGTCCATGCAATGATGTTCTTGACTGTGGATCTCTAATTAACAAATATATACTGCAG";

	$seqs[11] = "GTTAGTTTCAATATCTCTAGTTTTCTTGAGAAACATATCTTTCTATAGAGCCATATCTTATATTAACTGGCTTGTA
	TAGTGTAACGTGCTAAACAAGTTAGGTTTAGCCGTCTCTCTTTTAGAGAGCATATTCCATATTAGTACACTTTTCG
	ATTTTAAATTTGCTTTAATTAGGAGTGTAAAACATAACACATTCAATGTTTTCGATTTTTCTGATAATCTATATCA
	TTGCCTCTTTGGTTCGGTTCATTGTGGTTTTGTGATTATGGTTTTTGTTGTTCTGCAG";
	
}


