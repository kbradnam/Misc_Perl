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

# Can currently deal with introns, CDSs, 5' UTRs, upstream regions of transcripts, and whole transcripts


# command line options
my $window;     # how big a size window
my $min;	    # alternative start point
my $max;        # how big to go up to
my $step;       # for sliding windows
my $motif;	    # path to a valid NestedMICA XMS motif file
my $intron;     # path to input file of intron sequences
my $rev_intron; # will reverse complement intron sequences
my $cds;        # path to input file of intron sequences
my $upstream;   # path to input file up promoter regions of genes
my $five_utr;	# path to 5' UTR file
my $three_utr;	# path to 3' UTR file
my $transcript; # path to a file of protein-coding transcripts
my $threshold;  # log-odds threshold value to use when scoring motifs
my $species;    # two-letter species abbreviation

GetOptions ("window=i"     => \$window,
			"min=i"		   => \$min,
			"max=i"        => \$max,
			"step=i"       => \$step,
			"motif=s"      => \$motif,
			"intron=s"     => \$intron,
			"rev_intron"   => \$rev_intron,
			"cds=s"        => \$cds,
			"upstream=s"   => \$upstream,
			"five_utr=s"   => \$five_utr,
			"three_utr=s"  => \$three_utr,
			"transcript=s" => \$transcript,
			"threshold=i"  => \$threshold,
			"species=s"    => \$species);


# check command line options 
if (!$intron && !$cds && !$five_utr && !$three_utr && !$upstream && !$transcript){
	die "Please specify at least one of -intron, -cds, -five_utr, -three_utr, -upstream, or -transcript\n" 	
}

die "Please specify a valid NestedMica motif file with the -motif option\n" if (!$motif);

die "Please specify a valid species code, e.g. At, Os, Pp\n" if (!$species);
# set some defaults
$min = 1 if (!$min);
$max = 5000 if (!$max);
$window = 250 if (!$window);
$step = 100 if (!$step);
$threshold = 0 if (!$threshold);

# need to know end points of each window
my ($start,$end);

# hashes to store all sequences that fall in size range
# key to hash are start coordinates, value is sequence
my %intron2seqs;
my %cds2seqs;
my %upstream2seqs;
my %five_utr2seqs;
my %three_utr2seqs;
my %transcript2seqs;

# hashes to count how many introns and CDSs fall in each size category
my %intron2count;
my %cds2count;
my %upstream2count;
my %five_utr2count;
my %three_utr2count;
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
	if($cds){
		$type = "cds";
		($cds2seqs{"$start"}, $cds2count{"$start"}) = &process_sequence($cds,$type,$start,$end,$window);
	}
	# only have 3000 bp of upstream sequence to work with in Arabidopsis
	if(($upstream) && ($species =~ m/at/i) && ($end <= 3000)){
		$type = "upstream";
		($upstream2seqs{"$start"}, $upstream2count{"$start"}) = &process_sequence($upstream,$type,$start,$end,$window)
	}
	# only have 1000 bp of upstream sequence to work with in rice
	if(($upstream) && ($species =~ m/os/i) && ($end <= 1000)){
		$type = "upstream";
		($upstream2seqs{"$start"}, $upstream2count{"$start"}) = &process_sequence($upstream,$type,$start,$end,$window)
	}
	if($five_utr){
		$type = "five_utr";
		($five_utr2seqs{"$start"}, $five_utr2count{"$start"}) = &process_sequence($five_utr,$type,$start,$end,$window);
	}
	if($three_utr){
		$type = "three_utr";
		($three_utr2seqs{"$start"}, $three_utr2count{"$start"}) = &process_sequence($three_utr,$type,$start,$end,$window);
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
print "Upstream count,Upstream bases in motif,Total upstream bases,%motif in upstream region,";
print "5UTR count,5UTR bases in motif,Total 5UTR bases,%motif in 5UTR region,";
print "Intron count,Intron bases in motif,Total intron bases,%motif in intron,";
print "CDS count,CDS bases in motif,Total CDS bases,%motif in CDS,";
print "Transcript count,Transcript bases in motif,Total Transcript bases,%motif in transcript,";
print "3UTR count,3UTR bases in motif,Total 3UTR bases,%motif in 3UTR region\n";

# First deal with upstream sequence output (if present)
if($upstream){
	foreach my $key (sort {$a <=> $b}(keys %upstream2seqs)){
		# different start points depending on species
		($start = $key-3001) if ($species =~ m/at/i);
		($start = $key-1001) if ($species =~ m/os/i);
		$end = $start + $window -1;
		print "$start,$end,$upstream2count{$key},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$upstream2seqs{$key}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species ${species}u -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,\n";
	}
}

# One big loop to treat either 5' UTR, 3'UTR, intron, CDS, or transcript sequences
for(my $start = $min;$start<$max; $start+= $step){
	$end = $start + $window -1;

	my $data;
	print "$start,$end,0,0,0,0,";

	if($five_utr2seqs{$start}){
		print "$five_utr2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$five_utr2seqs{$start}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species ${species}5u -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
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
	
		# use different background frequencies if looking on reverse strand
		if($rev_intron){
			# reverse complement sequence
			my $revcom = reverse $intron2seqs{$start};
		    $revcom =~ tr/ACGTacgt/TGCAtgca/;
			print OUT "$revcom\n";
			close(OUT);
			$data = `~keith/Work/bin/where_is_motif.pl -species ${species}ir -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;			
		}
		else{
			print OUT "$intron2seqs{$start}\n";
			close(OUT);
			$data = `~keith/Work/bin/where_is_motif.pl -species ${species}i -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
		}
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
	}
	
	if($cds2seqs{$start}){
		print "$cds2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$cds2seqs{$start}\n";
		close(OUT);
		$data = `~keith/Work/bin/where_is_motif.pl -species ${species}c -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;			
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
		$data = `~keith/Work/bin/where_is_motif.pl -species ${species}t -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
		$data =~ m/.*: (\d+)\/(\d+) (.*)/;
		print "$1,$2,$3,";
	}
	else{
		print "0,0,0,0,";
	}
	
	if($three_utr2seqs{$start}){
		print "$three_utr2count{$start},";
		open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">${start}_$end\n";
		print OUT "$three_utr2seqs{$start}\n";
		close(OUT);
		my $data = `~keith/Work/bin/where_is_motif.pl -species ${species}3u -mdensity -threshold $threshold -target /tmp/ime_seq -motif $motif`;
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
		# with upstream regions of forward strand genes	(in Arabidopsis)	
		if(($type eq "upstream") && ($species =~ m/at/i) && ($header =~ m/FORWARD/)){
			$counter++;
			my $tmp = substr($seq,$win_start-1,$window);
			$new_seq .= "$tmp";
		}
		# also have to capture rice upstream regions
		elsif(($type eq "upstream") && ($species =~ m/os/i)){
			# skip the one pesky rice sequence which is not 1,000 bp
			next if ($length < 1000);
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
		
		
		# Can deal with CDSs, introns, and 5' & 3' UTRs and transcripts in similar ways as these all have similar 
		# format FASTA headers which give the distance to the TSS (except transcripts where to the distance
		# to the TSS is alwasy the same as $win_start)
		elsif($type eq "intron" || $type eq "cds" || $type eq "five_utr" || $type eq "three_utr"){
			
			# grab distance to TSS differently depending on what type of sequences we are dealing with
			my $distance;

			if($type eq "intron"){
				$header =~ m/_i\d+_(\d+)/;
				$distance = $1;
			}
		   	elsif($type eq "cds"){
				$header =~ m/_cds\d+_(\d+)/;
				$distance = $1;
			}
			elsif($type eq "five_utr"){
				$header =~ m/_5utr\d+_(\d+)/;
				$distance = $1;
			}
			elsif($type eq "three_utr"){
				$header =~ m/_3utr\d+_(\d+)/;
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
