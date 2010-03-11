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
my $background; # code to determine which background frequencies to use (from __DATA__ section)
my $background_file;     # read frequencies from file instead
my $mdensity;	# count (and show) amount and percentage of motif in each sequence (one line per sequence)
my $mseqs;		# just show sequences of each intron that have motifs above threshold (one sequence per intron)
my $mask;		# show all intron sequences (regardless of whether they have motifs) and mask out motif with -'s
my $msummary;	# show motif count and percentage for all sequences combined
my $mcount;     # just show count of motifs above threshold in each sequence
my $reverse;    # calculate frequencies for reverse strand
my $stats;      # report stats on all log likelihood scores found
my $min;        # Set minimum cut off for use with -stats option
my $max;        # Set maximum cut off for use with -stats option
my $interval;   # Set interval size for use with -stats option

GetOptions ("motif=s"    => \$motif,
	       "target=s"    => \$target,
	       "threshold=f" => \$threshold,
	       "scores"      => \$scores,
	       "seqs"        => \$seqs,	    	    
		   "bg:s"        => \$background,
		   "bg_file=s"   => \$background_file,
		   "mdensity"    => \$mdensity,  
		   "mseqs"		 => \$mseqs,
		   "mask"		 => \$mask,
		   "msummary"    => \$msummary,
		   "mcount"		 => \$mcount,
		   "stats"       => \$stats,
	       "min=f"       => \$min,
	       "max=f"       => \$max,
	       "interval=f"  => \$interval,
	       "reverse"     => \$reverse
) or die "\n";

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

# will need to know motif length for later on
my $motif_length;

# Need sets of expected nucleotide frequencies to compute log likelihood scores
# set of frequencies chosen by -bg_freqs or -bg_file options

my %expected;
read_background_frequencies();

# track base position in motif
my $pos = 0;
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

    # get nucleotide frequencies from input file
    if(m/weight symbol=\"([a-z])[a-z]+\">(\-*[10]\.\d+)<\/weight/){
		my $base = lc($1);
		my $freq = $2;
	
		# if frequency is zero, set to be a very small positive value else can't take log
		($freq = 0.00001) if ($freq < 0.00001);
			
		# take logarithm of observed over expected frequency and add to @motifs
		$motif[$pos]{$base} = log($freq/$expected{$base});
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

# keep track of:
# 1) total length of sequence in each file 
# 2) total length of motifs in each sequence 
# 3) number of sequences in input file
# 4) count of motifs in each sequence
# 5) total number of motifs in all sequences

my ($total_seq_length,$total_motif_length,$seq_count,$seq_motif_count,$total_motif_count);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	$seq_count++;
	$seq_motif_count = 0;
	# get and trim header
	my $header = $entry->def;
    $header =~ s/ .*//;

    my $seq = lc($entry->seq);
    my $length = length($seq);
	$total_seq_length += $length;
	
	# have one variable to hold all uppercased motifs in one variable
	my $all_motifs_seq = $seq;
	
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
			
			# count motifs
			$total_motif_count++;
			$seq_motif_count++;
			
			# and mask the motif out of $masked_seq
			substr($masked_seq,$i,$motif_length) = ("-" x $motif_length);
			
			# print out current motif details if -scores or -seqs specified
			if($scores){
				my $start_coord = $i+1;
			  	print "$header $score $start_coord/$length $window\n";
			}
			if($seqs || $mseqs){
				my $highlight = uc($window);
				my $new_seq = $seq;
				$new_seq =~ s/$window/ $highlight /g;
				print "$new_seq\n" if ($seqs);
				$all_motifs_seq =~ s/$window/$highlight/ig;
			}
		}
	}
	
	# count motif in sequence, add to running total
	my $motif_count = ($masked_seq =~ tr /-/-/);
	$total_motif_length += $motif_count;
	my $percent_motif = sprintf("%.3f",($motif_count / $length) * 100);
	print "$header motif_density: $motif_count/$length $percent_motif%\n" if ($mdensity);
	print "$header motif_count: $seq_motif_count\n" if ($mcount);
	
	# print out intron sequence if -mseqs specified 
	# if -mask is specified, will return all intron sequences, with just motifs masked out
	if ($mseqs && $mask){
		print "$header\n$masked_seq\n";
	}
	# if not in -mask mode, just show intron sequences that contain a motif above threshold
	elsif ($mseqs && !$mask && $above_threshold){
		# first have to add spaces to $all_motifs_seq
		$all_motifs_seq =~ s/([A-Z]+)/ $1 /g;
		print "$header\n$all_motifs_seq\n" 	
	}
}

close(TARGET) || die "Couldn't close $target\n";

# print motif summary if requested
if($msummary){
	my $percent_motif = sprintf("%.3f",($total_motif_length/$total_seq_length) * 100);
	print "\nSUMMARY:\n"; 
	print "motif_size: $motif_length\n";
	print "number_of_sequences: $seq_count total_sequence_length: $total_seq_length\n";
	print "number_of_motifs: $total_motif_count total_motif_length: $total_motif_length motif_density: $percent_motif%\n\n\n";
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

	# have we chosen an option to print some output?
	if (!$seqs && !$stats && !$scores && !$mdensity && !$msummary && !$mseqs && !$mcount){
		die "You have to choose at least one of the following options:\n-stats, -seqs, -scores, -mdensity, -msummary, -mseqs, -mcount\n";	
	}
	# can't choose mdensity with certain other options
	if($mdensity && ($seqs || $stats || $mseqs)){
		die "Don't choose -mdensity with -seqs or -scores or -mseqs\n";	
	}  
	# can't choose mdensity with certain other options
	if($mcount && ($seqs || $stats || $mseqs)){
		die "Don't choose -mcount with -seqs or -scores or -mseqs\n";	
	}
	# can't choose mseqs with other options
	if($mseqs && ($seqs || $stats || $mdensity || $mcount)){
		die "Don't choose -mseqs with -seqs or -scores or -mdensity or -mcount\n";	
	}  	
	# don't specify a threshold if using -stats
	if($threshold && $stats){
		die "Don't specify -threshold option when using -stats option\n";
	}
	# only specify -mask with -mseqs
	if($mask && !$mseqs){
		die "-mask option is only to be used with -mseqs option\n";
	}
	# always need to specify a background code, shouldn't just use -bg_file
	if($background_file && !$background){
#		die "Need to also specify -bg if using -bg_file\n";
	}
}


sub read_background_frequencies{
	my %background;

	# default filehandle
	my $fh = 'DATA';
	
	# use different filehandle if background file is supplied
	if($background_file){
		open(my $in, '<', $background_file) or die "Can't open file: \'$background_file\' $!";
		print "Reading from background file $background_file\n\n";
		$fh = $in;
	}

	# read <DATA> or $background_file
	while(<$fh>){
		chomp;
		next if m/^#/; # skip comment lines
		next if m/^$/; # skip blank lines

		# check we have enough fields
		my @line = split(/,/);
		if (@line < 6){
			print "The following background frequency information does not contain the minimum of 6 fields:\n$_\n\n";
			print "The correct format is a comma separated file with the following information:\n";
			print "<code>,<description>,A,C,G,T,optional notes\n\n";
			print "E.g.\n";
			print "Hsg,Homo sapiens genes,0.22,0.29,0.28,0.21,from my recent experiement\n\n";
			print "The <code> field should be used by the -bg option of this script.\n";
			exit;
		}
		my ($code,$details,$a,$c,$g,$t,$notes) = @line;

		# check frequencies, allow a little flexibility
		my $sum = $a + $c + $g + $t;
		die "The following background frequencies add up to $sum and not 1:\n$_\n" if ($sum < 0.99 or $sum > 1.01);
		
		# check $code field for whitespace
		if($code =~ m/\s+/){
			print "Removing spaces from \'$code\' field\n";
			$code =~ s/\s+//g;
		}			
		$background{$code}{'details'} = $details;
		$background{$code}{'a'} = $a;
		$background{$code}{'c'} = $c;
		$background{$code}{'g'} = $g;
		$background{$code}{'t'} = $t;
		$background{$code}{'notes'} = $notes;
	}

	# set an error status
	my $error = 0;
	if($background){
		print "You asked for background code \'$background\'. ";
		if($background{$background}){
			print "Using code ati: $background{$background}{'details'}\n";
			# load up %expected hash
			foreach my $base qw (a c g t){
				# remove white space if present from nucleotide frequency fields
				$background{$background}{$base} =~ s/\s+//g;

				# and now assign to final hash
				$expected{$base} = $background{$background}{$base}; 
			}	
		}
		else{
			$error = 1;
			print "This code does not exist\n";
		}	
	}
	else{
		$error = 1;
	}
	
	# print error message and summary of available error codes
	if($error){
		print "\nSpecify a background code (with the -bg option) from the following choices:\n\n";
		print "Code\tDescription\n----\t-----------\n";
		foreach my $key (sort keys %background){
			print "$key\t$background{$key}{'details'}\n";
		}		
		exit;
	}
}

__DATA__
# Code, Full name, Freq_A, Freq_C, Freq_G, Freq_T, Notes
cei,Caenorhabditis elegans introns,0.3333,0.1621,0.1597,0.3449,from 31044 WS180 confirmed introns
ceg,Caenorhabditis elegans genome,0.32280,0.17733,0.17709,0.32278,from WS180 chromosomes
ati,Arabidopsis thaliana introns,0.2713,0.1534,0.1701,0.4052,from 59260 high confidence TAIR7 introns
atg,Arabidopsis thaliana genome,0.3200,0.1802,0.1801,0.3197,from August 2006 TAIR set of chromosomes
at5u,Arabidopsis thaliana 5' UTRs,0.3017,0.2180,0.1593,0.3210,from 14862 high confidence TAIR7 5' UTRs
at3u,Arabidopsis thaliana 3' UTRs,0.2960,0.1473,0.1731,0.3836,from 13,459 high confidence TAIR7 3' UTRs
atig,Arabidopsis thaliana intergenic,0.3437,0.1568,0.1562,0.3433,from 304113 TAIR7 intergenic annotations
atc,Arabidopsis thaliana CDSs,0.2831,0.2067,0.2408,0.2694,from 70370 high confidence TAIR7 CDS exons
att,Arabidopsis thaliana transcripts,0.2815,0.1849,0.2078,0.3258,from 13163 high confidence TAIR7 transcripts
atu,Arabidopsis thaliana upstream regions,0.3342,0.1685,0.1651,0.3321,from 32041 TAIR7 1,000 bp upstream sequences
atd,Arabidopsis thaliana downstream regions,0.3211,0.1787,0.1724,0.3278,from 32041 TAIR7 1,000 bp downstream sequences
dmi,Drosophila melanogaster introns,0.2942,0.2040,0.1979,0.3040,
hs5u,Homo sapiens 5' UTRs,0.2134,0.2864,0.2970,0.2032,
osi,Oryza sativa introns,0.2769,0.1831,0.1875,0.3525,from 72066 high confidence japonica introns (TIGR 5.0 annotations)
os5u,Oryza sativa 5' UTRs,0.2083,0.3191,0.2491,0.2234,from 16267 high confidence japonica sequences (TIGR 5.0 annotations)
os3u,Oryza sativa 3' UTRs,0.2715,0.1874,0.2133,0.3278,from 14798 high confidence japonica sequences (TIGR 5.0 annotations)
osg,Oryza sativa genome,0.2822,0.2178,0.2178,0.2822,from 12 japonica genome sequences (TIGR 5.0 version)
osc,Oryza sativa CDSs,0.2421,0.2559,0.2774,0.2245,from 80515 high confidence japonica CDS exons (TIGR 5.0 annotations)
ost,Oryza sativa transcripts,0.2624,0.2125,0.2208,0.3043,from 13201 high confidence japonica transcripts (TIGR 5.0 annotations)
osig,Oryza sativa intergenic,0.2927,0.2073,0.2074,0.2926,from 56,208 japonica sequences (TIGR 5.0 annotations)
osu,Oryza sativa upstream regions,0.2960,0.2159,0.2020,0.2861,from 13201 japonica 1 kbp upstream regions
pti,Populus trichocarpa,0.2788,0.1581,0.1791,0.3839,from 26415 high confidence introns (v1.1 annotations)
pt5u,Populus trichocarpa 5' UTRs,0.3004,0.2249,0.1792,0.2954,26415 high confidence UTRs (v1.1 annotations)
pt3u,Populus trichocarpa 3' UTRs,0.2754,0.1658,0.1974,0.3614,from 7841 high confidence UTRs (v1.1 annotations)
ptg,Populus trichocarpa genome,0.3316,0.1687,0.1685,0.3311,from 22012 genome sequences (v1.0)
ptc,Populus trichocarpa CDSs,0.2816,0.2039,0.2437,0.2708,from 31232 high confidence poplar CDS exons (v1.1 annotations)
ptt,Populus trichocarpa transcripts,0.2803,0.1754,0.1996,0.3447,from 7089 high confidence poplar transcripts (v1.1 annotations)
ptu,Populus trichocarpa upstream regions,0.3473,0.1590,0.1480,0.3456,from 7089 1 kbp upstream sequences
vvi,Vitis vinifera introns,0.2963,0.1657,0.1822,0.3558,from 23091 high confidence introns
vv5u,Vitis vinifera 5' UTRs,0.2673,0.1987,0.1980,0.3360,from 5483 high confidence 5' UTRs 
vv3u,Vitis vinifera 3' UTRs,0.2823,0.1670,0.1937,0.3570,from 5232 high confidence 3' UTRs 
vvg,Vitis vinifera genome,0.3274,0.1727,0.1728,0.3272,from 35 genome sequences 
vvc,Vitis vinifera CDSs,0.2824,0.2023,0.2423,0.2729,from 27598 high confidence poplar CDS exons 
vvt,Vitis vinifera transcripts,0.2910,0.1740,0.1935,0.3416,from 5077 high confidence poplar transcripts 
vvu,Vitis vinifera upstream regions,0.3488,0.1513,0.1506,0.3492,from 50877 1 kbp upstream sequences
cri,Chlamydomonas reinhardtii introns,0.1774,0.2968,0.3261,0.1997,from 21979 high confidence introns
cru,Chlamydomonas reinhardtii upstream regions,0.1948,0.2934,0.3142,0.1976,from 4131 1 kbp upstream sequences
cr5u,Chlamydomonas reinhardtii 5' UTRs,0.2288,0.2869,0.2591,0.2252,from 4636 high confidence 5' UTRs 
cr3u,Chlamydomonas reinhardtii 3' UTRS,0.1979,0.2393,0.3451,0.2178,from 4430 high confidence 3' UTRs 
crg,Chlamydomonas reinhardtii genome,0.1802,0.3198,0.3198,0.1802,from 1266 genome sequences 
crc,Chlamydomonas reinhardtii CDSs,0.1812,0.3212,0.3361,0.1615,from 25286 high confidence CDS exons 
crt,Chlamydomonas reinhardtii transcripts,0.1878,0.2982,0.3211,0.1929,from 4131 high confidence poplar transcripts 
scg,Saccharomyces cerevisiae genome,0.3192,0.1903,0.1923,0.2982,from 5884 S. cerevisiae genes (including 1000 nt upstream and downstream plus introns) 
ppi,Physcomitrella patens introns,0.2788,0.1582,0.1791,0.3840,from 26415 P. patens introns 
sbi,Sorghum bicolor introns,0.2636,0.1965,0.2004,0.3395,from 57284 S. bicolor introns 
smi,Selaginella moellendorfii introns,0.2043,0.2115,0.2086,0.3755,from 11118 S. moellendorfii introns 
vci,Volvox carteri introns, 0.2071,0.2658,0.2724,0.2547,from 10960 V. carteri introns 