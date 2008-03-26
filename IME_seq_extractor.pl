#!/usr/bin/perl
#
# IME_seq_extractor.pl
#
# A script to extract regions of sequences which appear to have the highest IME motif content
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


# This script windows through a set of 5' UTR or intron sequences and find the coordinates which gives the highest
# concentration of a specified motif. It then extracts all intron or UTR sequence in that range and writes it to
# a fasta file for further study (e.g. with nested mica)


# command line options
my $window;     # how big a size window
my $min;	    # alternative start point
my $max;        # how big to go up to
my $step;       # for sliding windows
my $motif;	    # path to a valid NestedMICA XMS motif file
my $intron;     # path to input file of intron sequences
my $five_utr;	# path to 5' UTR file
my $threshold;  # log-odds threshold value to use when scoring motifs
my $species;    # two-letter species abbreviation
my $seqdump;	# print sequences to a file, don't calculate density
my $seqlength;	# threshold for dumping sequences to file (don't want short sequences)

GetOptions ("window=i"     => \$window,
			"min=i"		   => \$min,
			"max=i"        => \$max,
			"step=i"       => \$step,
			"motif=s"      => \$motif,
			"intron=s"     => \$intron,
			"five_utr=s"   => \$five_utr,
			"threshold=i"  => \$threshold,
			"species=s"    => \$species,
			"seqdump"      => \$seqdump,
			"seqlength=i"  => \$seqlength);


# check command line options 
if (!$intron && !$five_utr){
	die "Please specify one of -intron or -five_utr\n" 	
}

if ($intron && $five_utr){
	die "Please specify just one of -intron or -five_utr\n" 	
}

die "Please specify a valid NestedMica motif file with the -motif option\n" if (!$motif);

die "Please specify a valid species code, e.g. At, Os, Pp\n" if (!$species);
# set some defaults
$min = 1 if (!$min);
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);
$threshold = 0 if (!$threshold);
$seqlength = 25 if (!$seqlength);

# need to know end points of each window
my ($start,$end);

# hashes to store all sequences that fall in size range
# key to hash are start coordinates, value is sequence
my %intron2seqs;
my %five_utr2seqs;

# hashes to count how many introns and CDSs fall in each size category
my %intron2count;
my %five_utr2count;
 
my $type;
($type = "intron") if ($intron);	
($type = "five_utr") if ($five_utr);	

# open output file if using -seqdump
if($seqdump){
	open(SEQ,">seqdump_${species}_${type}.fa") || die "couldn't write seqdump file\n";
}


for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;
	#print "$start - $end\n";

	# is there an intron file to process?
	if($intron){
		($intron2seqs{"$start"}, $intron2count{"$start"}) = &process_sequence($intron,$type,$start,$end,$window,"no");
	}

	if($five_utr){
		($five_utr2seqs{"$start"}, $five_utr2count{"$start"}) = &process_sequence($five_utr,$type,$start,$end,$window,"no");
	}

}

##################################################################
#
# MAIN OUTPUT
#
##################################################################

# keep track of best scores
my ($intron_max, $utr_max) = (0,0);

print "Start\tEnd\t";
print "5UTR count\t%motif in 5UTR region\t" if ($five_utr);
print "Intron count\t%motif in intron\n" if ($intron);


# One big loop to treat either 5' UTR or intron sequences
for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;

	my $data;
	my $percent;
	print "$start\t$end\t";

	if($five_utr2seqs{$start}){
		print "$five_utr2count{$start}\t";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$five_utr2seqs{$start}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species ${species}5u -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$3\t";
		$percent = $3;
		chop($percent);	
		
		# quit if we have not beaten best score so far
		if($percent >= $utr_max){
			$utr_max = $percent;
		}
		else{
			print "\n";
			my ($s2,$e2);
			$s2 = $start - $step;
			$e2 = $end - $step;
			print "Best window: $s2 - $e2\n";
			&process_sequence($five_utr,"five_utr",$s2,$e2,$window,"yes");
			close(SEQ);
			exit(0) if ($seqdump);
		}
	}


	if($intron2seqs{$start}){
		print "$intron2count{$start}\t";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";	
		print OUT "$intron2seqs{$start}\n";
		close(OUT);
		$data = `~keith/Work/bin/where_is_motif.pl -species ${species}i -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
	
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$3\t";
		$percent = $3;	
		chop($percent);
		
		# quit if we have not beaten best score so far
		if($percent >= $intron_max){
			$intron_max = $percent;
		}
		else{
			print "\n";
			my ($s2,$e2);
			$s2 = $start - $step;
			$e2 = $end - $step;
			print "Best window: $s2 - $e2\n";
			&process_sequence($intron,"intron",$s2,$e2,$window,"yes");
			close(SEQ);
			exit(0) if ($seqdump);
		}
	}
	print "\n";
}




# main subroutinee to loop through file of sequences and extract only sequence in certain range

sub process_sequence{

	my $file = shift;
	my $type = shift;
	my $win_start = shift;
	my $win_end = shift;
	my $window = shift;

	# dump option will just be yes or no, will only run seqdump if -seqdump is specified and $dump is 'yes'
	my $dump = shift;
	
	# count how many sequences fall into each bin
	my $counter = 0;
	
	# store all sequences in range in one variable
	my $new_seq = "";
	
	open(FILE,"<$file") || die "Couldn't open $file\n";		

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {
	    
		# grab basic details
		my $seq = $entry->seq;
		my $length = length($seq);
		my $header = $entry->def;
		
		
		# Can deal with introns and 5 UTRs in similar ways as these all have similar 
		# format FASTA headers which give the distance to the TSS 
			
		# grab distance to TSS differently depending on what type of sequences we are dealing with
		my $distance;

		if($type eq "intron"){
			$header =~ m/_i\d+_(\d+)/;
			$distance = $1;
		}
		elsif($type eq "five_utr"){
			$header =~ m/_5utr\d+_(\d+)/;
			$distance = $1;
		}

		# calculate end coordinate of sequence
		my $end_coord = $distance + $length -1;


		###################################################################
		# check whether candidate sequence falls in various size categories
		###################################################################

		# CASE 1: intron/CDS/utr/transcript sequence is wholly contained within window
		if(($distance >= $win_start) && ($end_coord <= $win_end)){
			$counter++;
			my $tmp = substr($seq,0,$end_coord-$distance+1);
			$new_seq .= "$tmp";
			
			# dump sequence if using -seqdump and there is at least 10 nt
			print SEQ "$header $length nt\n$tmp\n" if (($seqdump) && (length($tmp) >= $seqlength) && ($dump eq "yes"));
		}

		# CASE 2: sequence is larger than  window
		elsif(($distance < $win_start) && ($end_coord > $win_end)){
			$counter++;
			my $tmp = substr($seq,$win_start-$distance,$window);
			$new_seq .= "$tmp";
			
			# dump sequence if using -seqdump and there is at least 10 nt
			print SEQ "$header $length nt\n$tmp\n" if (($seqdump) && (length($tmp) >= $seqlength) && ($dump eq "yes"));
		}

		# CASE 3: sequence overlaps 5' edge of window
		elsif(($distance < $win_start) && ($end_coord >= $win_start)){
			$counter++;
			my $tmp = substr($seq,$win_start - $distance,$end_coord - $win_start + 1);
			$new_seq .= "$tmp";
			
			# dump sequence if using -seqdump and there is at least 10 nt
			print SEQ "$header $length nt\n$tmp\n" if (($seqdump) && (length($tmp) >= $seqlength) && ($dump eq "yes"));
		}

		# CASE 4: sequence overlaps 3' edge of window
		elsif(($distance <= $win_end) && ($end_coord > $win_end)){
			$counter++;
			my $tmp = substr($seq,0,$win_end - $distance + 1);
			$new_seq .= "$tmp";
			
			# dump sequence if using -seqdump and there is at least 10 nt
			print SEQ "$header $length nt\n$tmp\n" if (($seqdump) && (length($tmp) >= $seqlength) && ($dump eq "yes"));
		}		
				
	}
	close(FILE) || die "Can't close file\n";

	return($new_seq,$counter);
}




exit(0);
