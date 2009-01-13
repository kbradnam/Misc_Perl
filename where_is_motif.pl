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
my $mask;		# show all intron sequences (regardless of whether they have motifs) and mask out motif with -'s
my $msummary;	# show motif count and percentage for all sequences combined
my $mcount;     # just show count of motifs above threshold in each sequence

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
		   "mask"		 => \$mask,
		   "msummary"    => \$msummary,
		   "mcount"		 => \$mcount,
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

# will need to know motif length for later on
my $motif_length;

# Need sets of expected nucleotide frequencies to compute log likelihood scores
# set of frequencies chosen by -species option
# following may have to be tidied up if I need to add more species

my %expected;

if($species =~ m/cei/i){
	# from 31,044 WS180 confirmed introns
	%expected = ("a" => "0.3333","c" => "0.1621", "g" => "0.1597","t" => "0.3449");
}
elsif($species =~ m/ceg/i){
	# from WS180 chromosomes
	%expected = ("a" => "0.32280","c" => "0.17733","g" => "0.17709","t" => "0.32279");
}
elsif($species =~ m/^ati$/i){
	# from 59,260 high confidence TAIR7 introns
	%expected = ("a" => "0.2713","c" => "0.1534", "g" => "0.1701","t" => "0.4051");
}
elsif($species =~ m/^atir$/i){
	# from 59,260 high confidence TAIR7 introns
	%expected = ("a" => "0.4051","c" => "0.1701", "g" => "0.1534","t" => "0.2713");
}
elsif($species =~ m/atg/i){
	# from August 2006 TAIR set of chromosomes
	%expected = ("a" => "0.3200","c" => "0.1802", "g" => "0.1801","t" => "0.3197");
}
elsif($species =~ m/at5u/i){
	# from 14,862 high confidence TAIR7 5' UTRs
	%expected = ("a" => "0.3017","c" => "0.2180", "g" => "0.1593","t" => "0.3210");
}
elsif($species =~ m/at3u/i){
	# from 13,459 high confidence TAIR7 3' UTRs
	%expected = ("a" => "0.2960","c" => "0.1473", "g" => "0.1731","t" => "0.3835");
}
elsif($species =~ m/atig/i){
	# from 30,4113 TAIR7 intergenic annotations
	%expected = ("a" => "0.3437","c" => "0.1568", "g" => "0.1562","t" => "0.3433");
}
elsif($species =~ m/atc/i){
	# from 70,370 high confidence TAIR7 CDS exons
	%expected = ("a" => "0.2831","c" => "0.2067", "g" => "0.2408","t" => "0.2694");
}
elsif($species =~ m/att/i){
	# from 13,163 high confidence TAIR7 transcripts
	%expected = ("a" => "0.2815","c" => "0.1849", "g" => "0.2078","t" => "0.3258");
}
elsif($species =~ m/atu/i){
	# from 32,041 TAIR7 1,000 bp upstream sequences
	%expected = ("a" => "0.3342","c" => "0.1685", "g" => "0.1651","t" => "0.3321");
}
elsif($species =~ m/atd/i){
	# from 32,041 TAIR7 1,000 bp downstream sequences
	%expected = ("a" => "0.3211","c" => "0.1787", "g" => "0.1724","t" => "0.3278");
}
elsif($species =~ m/^dmi$/i){
	%expected = ("a" => "0.2942","c" => "0.2040", "g" => "0.1979","t" => "0.3040");
}
elsif($species =~ m/^hs5u$/i){
	%expected = ("a" => "0.2134","c" => "0.2864", "g" => "0.2970","t" => "0.2032");
}
elsif($species =~ m/^osi$/i){
	# from 72,066 high confidence japonica introns (TIGR 5.0 annotations)
	%expected = ("a" => "0.2769","c" => "0.1831", "g" => "0.1875","t" => "0.3525");
}
elsif($species =~ m/os5u/i){
	# from 16,267 high confidence japonica sequences (TIGR 5.0 annotations)
	%expected = ("a" => "0.2083","c" => "0.3191", "g" => "0.2491","t" => "0.2234");
}
elsif($species =~ m/os3u/i){
	# from 14,798 high confidence japonica sequences (TIGR 5.0 annotations)
	%expected = ("a" => "0.2715","c" => "0.1874", "g" => "0.2133","t" => "0.3278");
}
elsif($species =~ m/osg/i){
	# from 12 japonica genome sequences (TIGR 5.0 version)
	%expected = ("a" => "0.2822","c" => "0.2178", "g" => "0.2178","t" => "0.2822");
}
elsif($species =~ m/osc/i){
	# from 80,515 high confidence japonica CDS exons (TIGR 5.0 annotations)
	%expected = ("a" => "0.2421","c" => "0.2559", "g" => "0.2774","t" => "0.2245");
}
elsif($species =~ m/ost/i){
	# from 13,201 high confidence japonica transcripts (TIGR 5.0 annotations)
	%expected = ("a" => "0.2624","c" => "0.2125", "g" => "0.2208","t" => "0.3043");
}
elsif($species =~ m/osig/i){
	# from 56,208 japonica sequences (TIGR 5.0 annotations)
	%expected = ("a" => "0.2927","c" => "0.2073", "g" => "0.2074","t" => "0.2926");
}
elsif($species =~ m/osu/i){
	# from 13,201 japonica 1 kbp upstream regions
	%expected = ("a" => "0.2960","c" => "0.2159", "g" => "0.2020","t" => "0.2861");
}
elsif($species =~ m/^pti$/i){
	# from 26,415 high confidence introns (v1.1 annotations)
	%expected = ("a" => "0.2788","c" => "0.1582", "g" => "0.1791","t" => "0.35384025");
}
elsif($species =~ m/pt5u/i){
	# from 26,415 high confidence UTRs (v1.1 annotations)
	%expected = ("a" => "0.3004","c" => "0.2249", "g" => "0.1792","t" => "0.2954");
}
elsif($species =~ m/pt3u/i){
	# from 7,841 high confidence UTRs (v1.1 annotations)
	%expected = ("a" => "0.2754","c" => "0.1658", "g" => "0.1974","t" => "0.3614");
}
elsif($species =~ m/ptg/i){
	# from 22,012 genome sequences (v1.0)
	%expected = ("a" => "0.3316","c" => "0.1687", "g" => "0.1685","t" => "0.3311");
}
elsif($species =~ m/ptc/i){
	# from 31,232 high confidence poplar CDS exons (v1.1 annotations)
	%expected = ("a" => "0.2816","c" => "0.2039", "g" => "0.2437","t" => "0.2708");
}
elsif($species =~ m/ptt/i){
	# from 7,089 high confidence poplar transcripts (v1.1 annotations)
	%expected = ("a" => "0.2803","c" => "0.1754", "g" => "0.1996","t" => "0.3447");
}
elsif($species =~ m/ptu/i){
	# from 7,089 1 kbp upstream sequences
	%expected = ("a" => "0.3473","c" => "0.1590", "g" => "0.1480","t" => "0.3456");
}
elsif($species =~ m/^vvi$/i){
	# from 23,091 high confidence introns
	%expected = ("a" => "0.2963","c" => "0.1657", "g" => "0.1822","t" => "0.3558");
}
elsif($species =~ m/vv5u/i){
	# from 5,483 high confidence 5' UTRs 
	%expected = ("a" => "0.2673","c" => "0.1987", "g" => "0.1980","t" => "0.3360");
}
elsif($species =~ m/vv3u/i){
	# from 5,232 high confidence 3' UTRs 
	%expected = ("a" => "0.2823","c" => "0.1670", "g" => "0.1937","t" => "0.3570");
}
elsif($species =~ m/vvg/i){
	# from 35 genome sequences 
	%expected = ("a" => "0.3274","c" => "0.1727", "g" => "0.1728","t" => "0.3272");
}
elsif($species =~ m/vvc/i){
	# from 27,598 high confidence poplar CDS exons 
	%expected = ("a" => "0.2754","c" => "0.2098", "g" => "0.2484","t" => "0.2464");
}
elsif($species =~ m/vvt/i){
	# from 5,077 high confidence poplar transcripts 
	%expected = ("a" => "0.2910","c" => "0.1740", "g" => "0.1935","t" => "0.3416");
}	
elsif($species =~ m/vvu/i){
	# from 5,0877 1 kbp upstream sequences
	%expected = ("a" => "0.3488","c" => "0.1513", "g" => "0.1506","t" => "0.3492");
}
elsif($species =~ m/^cri$/i){
	# from 21,979 high confidence introns
	%expected = ("a" => "0.1774","c" => "0.2968", "g" => "0.3261","t" => "0.1997");
}
elsif($species =~ m/cru/i){
	# from 4,131 1 kbp upstream sequences
	%expected = ("a" => "0.1948","c" => "0.2934", "g" => "0.3142","t" => "0.1976");
}
elsif($species =~ m/cr5u/i){
	# from 4,636 high confidence 5' UTRs 
	%expected = ("a" => "0.2288","c" => "0.2869", "g" => "0.2591","t" => "0.2252");
}
elsif($species =~ m/cr3u/i){
	# from 4,430 high confidence 3' UTRs 
	%expected = ("a" => "0.1979","c" => "0.2393", "g" => "0.3451","t" => "0.2178");
}
elsif($species =~ m/crg/i){
	# from 1,266 genome sequences 
	%expected = ("a" => "0.1802","c" => "0.3198", "g" => "0.3198","t" => "0.1802");
}
elsif($species =~ m/crc/i){
	# from 25,286 high confidence CDS exons 
	%expected = ("a" => "0.1812","c" => "0.3212", "g" => "0.3361","t" => "0.1615");
}
elsif($species =~ m/crt/i){
	# from 4,131 high confidence poplar transcripts 
	%expected = ("a" => "0.1878","c" => "0.2982", "g" => "0.3211","t" => "0.1929");
}
elsif($species =~ m/scg/i){
	# from 5,884 S. cerevisiae genes (including 1000 nt upstream and downstream plus introns) 
	%expected = ("a" => "0.3192","c" => "0.1903", "g" => "0.1923","t" => "0.2982");
}
else{
	die "\'$species\' is not a valid species code.\n";
}


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
		print "$header\n$seq\n" 	
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

	# check that species code has been chosen
	if(!$species){
		print "\nPlease specify a suitable species code using the -species option.\n";
		print "Species codes are required to determine the correct expected nucleotide\nfrequencies when scoring motifs\n\n";
		print "Current options (all case-insensitive) are:\n";
		print "AtI  - Arabidopsis thaliana introns\n";
		print "AtIR - Arabidopsis thaliana introns (reverse complemented)\n";
		print "AtG  - Arabidopsis thaliana genomic\n";
		print "AtC  - Arabidopsis thaliana CDSs\n";
		print "AtIG - Arabidopsis thaliana intergenic\n";
		print "At5U - Arabidopsis thaliana 5' UTR (exons)\n";
		print "At3U - Arabidopsis thaliana 3' UTR (exons)\n";
		print "AtT - Arabidopsis thaliana primary transcripts (unspliced)\n"; 
		print "AtU  - Arabidopsis thaliana upstream region of genes (1000 bp 5' to transcript)\n";
		print "AtD  - Arabidopsis thaliana downstream region of genes (1000 bp 3' to transcript)\n";
		print "CeI - Caenorhabditis elegans introns\n";
		print "CeG - Caenorhabditis elegans genomic\n";
		print "Hs5U - Homo sapiens 5' UTR (exons)\n";
		print "OsI  - Oryza sativa (japonica) introns\n";
		print "OsG  - Oryza sativa (japonica) genomic\n";
		print "OsC  - Oryza sativa (japonica) CDSs\n";
		print "OsIG - Oryza sativa (japonica) intergenic\n";
		print "Os5U - Oryza sativa (japonica) 5' UTR (exons)\n";
		print "Os3U - Oryza sativa (japonica) 3' UTR (exons)\n";
		print "OsT - Oryza sativa (japonica) primary transcripts (unspliced)\n";
		print "OsU - Oryza sativa (japonica) upstream regions (1 kbp)\n";
		print "PtI  - Populus trichocarpa introns (v1.1 annotations)\n";
		print "PtG  - Populus trichocarpa genomic (v1.0 sequence)\n";
		print "PtC  - Populus trichocarpa CDSs\n";
		print "Pt5U - Populus trichocarpa 5' UTR (exons)\n";
		print "Pt3U - Populus trichocarpa 3' UTR (exons)\n";
		print "PtT - Populus trichocarpa primary transcripts (unspliced)\n";
		print "PtU - Populus trichocarpa upstream regions (1 kbp)\n";
		print "VvI  - Vitis vinifera introns\n";
		print "VvG  - Vitis vinifera genomic\n";
		print "VvC  - Vitis vinifera CDSs\n";
		print "Vv5U - Vitis vinifera 5' UTR (exons)\n";
		print "Vv3U - Vitis vinifera 3' UTR (exons)\n";
		print "VvT - Vitis vinifera primary transcripts (unspliced)\n";
		print "VvU - Vitis vinifera upstream regions (1 kbp)\n";
		print "CrI  - Chlamydomonas reinhardtii introns\n";
		print "CrG  - Chlamydomonas reinhardtii genomic\n";
		print "CrC  - Chlamydomonas reinhardtii CDSs\n";
		print "Cr5U - Chlamydomonas reinhardtii 5' UTR (exons)\n";
		print "Cr3U - Chlamydomonas reinhardtii 3' UTR (exons)\n";
		print "CrT - Chlamydomonas reinhardtii primary transcripts (unspliced)\n";
		print "CrU - Chlamydomonas reinhardtii upstream regions (1 kbp)\n";
		print "DmI - Drosophila melanogaster introns\n";
		print "ScG - Saccharomyces cerevisiae genes (including 1,000 nt upstream and downstream)\n";
		die "Choose one option only.\n\n";
	}

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
}