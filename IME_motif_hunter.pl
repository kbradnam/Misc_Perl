#!/usr/bin/perl
#
# IME_motif_hunter.pl
#
# A script to study motif densities across a set of transcripts by using a windowing approach
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
# calculate motif density (using a NestedMICA motif)

# Can currently deal with introns, exons, 5' UTRs, upstream regions of transcripts, and whole transcripts


# command line options
my $window;     # how big a size window
my $min;	    # alternative start point
my $max;        # how big to go up to
my $step;       # for sliding windows
my $motif;	    # path to a valid NestedMICA XMS motif file
my $intron;     # path to input file of intron sequences
my $exon;       # path to input file of intron sequences
my $upstream;   # path to input file up promoter regions of genes
my $utr;	    # path to 5' UTR file
my $transcript; # path to a file of protein-coding transcripts

GetOptions ("window=i"     => \$window,
			"min=i"		   => \$min,
			"max=i"        => \$max,
			"step=i"       => \$step,
			"motif=s"      => \$motif,
			"intron=s"     => \$intron,
			"exon=s"       => \$exon,
			"upstream=s"   => \$upstream,
			"utr=s"        => \$utr,
			"transcript=s" => \$transcript);


# check command line options 
if (!$intron && !$exon && !$utr && !$upstream && !$transcript){
	die "Please specify at least one of -intron, -exon, -utr, -upstream, or -transcript\n" 	
}

die "Please specify a valid NestedMica motif file with the -motif option\n" if (!$motif);

# set some defaults
$min = 1 if (!$min);
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);

# need to know end points of each window
my ($start,$end);

# hashes to store all sequences that fall in size range
# key to hash are start coordinates, value is sequence
my %intron2seqs;
my %exon2seqs;
my %upstream2seqs;
my %utr2seqs;
my %transcript2seqs;

# hashes to count how many introns and exons fall in each size category
my %intron2count;
my %exon2count;
my %upstream2count;
my %utr2count;
my %transcript2count;

for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;
	#print "$start - $end\n";
	my $type;
	
	# is there an intron file to process?
	if($intron){
		$type = "intron";
		($intron2seqs{"$start"}, $intron2count{"$start"}) = &process_sequence($intron,$type,$start,$end,$window);
	}
	if($exon){
		$type = "exon";
		($exon2seqs{"$start"}, $exon2count{"$start"}) = &process_sequence($exon,$type,$start,$end,$window);
	}
	# only have 3000 bp of upstream sequence to work with
	if($upstream && ($end <= 3000)){
		$type = "upstream";
		($upstream2seqs{"$start"}, $upstream2count{"$start"}) = &process_sequence($upstream,$type,$start,$end,$window)
	}
	if($utr){
		$type = "utr";
		($utr2seqs{"$start"}, $utr2count{"$start"}) = &process_sequence($utr,$type,$start,$end,$window);
	}
	if($transcript){
		$type = "transcript";
		($transcript2seqs{"$start"}, $transcript2count{"$start"}) = &process_sequence($transcript,$type,$start,$end,$window);
	}
}

##################################################################
#
# MAIN OUTPUT
#
##################################################################

print "Start,End,";
print "Upstream_count,Upstream_bases_in_motif,Total_upstream_bases,%motif_in_upstream_region,";
print "5UTR_count,UTR_bases_in_motif,Total_UTR_bases,%motif_in_UTR_region,";
print "Intron_count,Intron_bases_in_motif,Total_intron_bases,%motif_in_intron,";
print "Exon_count,Exon_bases_in_motif,Total_Exon_bases,%motif_in_exon,";
print "Transcript_count,Transcript_bases_in_motif,Total_Transcript_bases,%motif_in_transcript\n";

# First deal with upstream sequence output (if present)
if($upstream){
	foreach my $key (sort {$a <=> $b}(keys %upstream2seqs)){
		$start = $key-3001;
		$end = $start + $window -1;
		print "$start,$end,$upstream2count{$key},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$upstream2seqs{$key}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species atu -mdensity -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,\n";
	}
}

# One big loop to treat either 5' UTR, intron, exon, or transcript sequences
for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;

	my $data;
	print "$start,$end,0,0,0,0,";

	if($utr2seqs{$start}){
		print "$utr2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$utr2seqs{$start}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species at5u -mdensity -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
	}


	if($intron2seqs{$start}){
		print "$intron2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$intron2seqs{$start}\n";
		close(OUT);
		$data = `~keith/Work/bin/where_is_motif.pl -species ati -mdensity -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
	}
	
	if($exon2seqs{$start}){
		print "$exon2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$exon2seqs{$start}\n";
		close(OUT);
		$data = `~keith/Work/bin/where_is_motif.pl -species atc -mdensity -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
	}
	
	if($transcript2seqs{$start}){
		print "$transcript2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$transcript2seqs{$start}\n";
		close(OUT);
		$data = `~keith/Work/bin/where_is_motif.pl -species att -mdensity -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
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
		
		# will treat upstream regions completely differently, for simplicity just deal
		# with upstream regions of forward strand genes		
		if(($type eq "upstream") && ($header =~ m/FORWARD/)){
			$counter++;
			my $tmp = substr($seq,$win_start-1,$window);
			$new_seq .= "$tmp";
		}
		
		# will treat transcript regions sort of like upstream (without the forward strand requirement)
		elsif($type eq "transcript"){
			# Can either have part of a transcript in a window
			if(($length >= $win_start) && ($length < $win_end)){
				$counter++;
				my $tmp = substr($seq,$win_start-1,$length-$win_start+1);
				$new_seq .= "$tmp";	
			}
			# or it should all be in the window, in which case we just take the entire window of sequence
			elsif($length >= $win_start){
				$counter++;
				my $tmp = substr($seq,$win_start-1,$window);
				$new_seq .= "$tmp";				
			}			
		}
		
		
		# Can deal with exons, introns, and 5' UTRs and transcripts in similar ways as these all have similar 
		# format FASTA headers which give the distance to the TSS (except transcripts where to the distance
		# to the TSS is alwasy the same as $win_start)
		elsif($type eq "intron" || $type eq "exon" || $type eq "utr"){
			
			# grab distance to TSS differently depending on what type of sequences we are dealing with
			my $distance;

			if($type eq "intron"){
				$header =~ m/_i\d+_(\d+)_/;
				$distance = $1;
			}
		   	elsif($type eq "exon"){
				$header =~ m/_e\d+_(\d+)/;
				$distance = $1;
			}
			elsif($type eq "utr"){
				$header =~ m/_5utr\d+_(\d+)/;
				$distance = $1;
			}
		

			# calculate end coordinate of sequence
			my $end_coord = $distance + $length -1;

			###################################################################
			# check whether candidate sequence falls in various size categories
			###################################################################

			# CASE 1: intron/exon/utr/transcript sequence is wholly contained within window
			if(($distance >= $win_start) && ($end_coord <= $win_end)){
				$counter++;
				my $tmp = substr($seq,0,$end_coord-$distance+1);
				$new_seq .= "$tmp";
			}

			# CASE 2: sequence is larger than  window
			elsif(($distance < $win_start) && ($end_coord > $win_end)){
				$counter++;
				my $tmp = substr($seq,$win_start-$distance,$window);
				$new_seq .= "$tmp";
			}

			# CASE 3: sequence overlaps 5' edge of window
			elsif(($distance < $win_start) && ($end_coord >= $win_start)){
				$counter++;
				my $tmp = substr($seq,$win_start - $distance,$end_coord - $win_start + 1);
				$new_seq .= "$tmp";
			}

			# CASE 4: sequence overlaps 3' edge of window
			elsif(($distance <= $win_end) && ($end_coord > $win_end)){
				$counter++;
				my $tmp = substr($seq,0,$win_end - $distance + 1);
				$new_seq .= "$tmp";
			}		
		}
				
	}
	close(FILE) || die "Can't close file\n";

	return($new_seq,$counter);
}




exit(0);
