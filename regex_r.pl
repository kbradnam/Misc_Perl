#!/usr/bin/perl
#
# regex_breeder.pl
#
# A script to calculate the r2 value for any regex motif against a set of introns
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use List::Util qw(sum);

die "Usage:
regex_r.pl <motif regex> <number of strands to search> <distance threshold> <file name(s) containing introns>\n" unless (@ARGV == 4);
my ($motif, $strands, $distance, @files) = @ARGV;

# turn on autoflush
$| = 1;

# want to keep track of sequences from input file, plus their expression level, plus counts of motifs
my @seqs;
my @expression;
my @counts;

# first convert any IUPAC characters to regex equivalents
my %bases_to_regex = (
	R => '(A|G)',
	Y => '(C|T)',
	S => '(C|G)',
	W => '(A|T)',
	K => '(T|G)',
	M => '(C|A)',
	B => '(C|G|T)',
	D => '(A|G|T)',
	H => '(A|C|T)',
	V => '(A|C|G)',
	N => '(A|C|G|T)'
);

$motif = uc($motif);
foreach my $code (keys %bases_to_regex){
	$motif =~ s/$code/$bases_to_regex{$code}/g;
}



foreach my $file (@files){
	# reset arrays
	@counts = ();
	@seqs = ();
	@expression = ();

	my ($plus_count, $rev_count) = (0, 0);
	# read sequence files and also extract expression values
	read_file($file);

	# now loop through each sequence
	for (my $i=0; $i < @seqs; $i++){
		my $seq = $seqs[$i];

		# count motif in sequence
		my ($count, $revcount) = (0,0);
		$count = $seq =~ s/($motif)/$1/g;
		($count = 0) if (!$count);

		# now check reverse strand if needed
		if($strands == 2){
			my $revcomp = Keith::revcomp($seq);
			$revcount = $revcomp =~ s/($motif)/$1/g;
			($revcount = 0) if (!$revcount);
		}
		$counts[$i]    = ($count + $revcount);
		$plus_count += $count;
		$rev_count += $revcount;
	}

	# can now calculate correlation (r) for that motifset
	my ($r,$n) = r2(\@expression,\@counts);

	print "$file: r = $r, $n introns, +/- strand motfs = $plus_count/$rev_count\n";
	
}

exit;

sub read_file{
	my $file = shift;
	open(IN,"<$file") || die "Couldn't open $file file\n";

    my $fasta = new FAlite(\*IN);

    # loop through each sequence in target file
    while(my $entry = $fasta->nextEntry) {

 		# get and trim header
    	my $header = $entry->def;
    	$header =~ m/.* x(\d{1,2}\.\d)/;
    	my $expression = $1;

		die "$header\n" if (!$expression);

        my $seq = uc($entry->seq);
		# add sequence and expression values to arrays
		#trim sequences as necessary
		$seq = substr $seq, 0, $distance;
		
		push(@seqs,$seq);
		push(@expression,$expression);
	}
	close(IN);
	
}

#######################################################
#
#          C a l c u la t e   R  - s q u a r e d
#
#######################################################

sub r2{
	# three arrays of data, need to regress against each other
    my $x = shift;
    my $y = shift;

    # need a bunch of (fairly standard) of statistics from both arrays
    my $n       = @{$x};
    my $sum_x   = sum(@{$x});
    my $sum_y   = sum(@{$y});
    my $sum_xy  = sum(map {${$x}[$_] * ${$y}[$_]} (0..$n-1));
    my $sum_x2  = sum(map {${$x}[$_]**2} (0..$n-1));
    my $sum_y2  = sum(map {${$y}[$_]**2} (0..$n-1));

    # can't calculate r2 if sum_x or sum_y = 0; so return 0 instead
    if($sum_x == 0 || $sum_y == 0){
    	return("0");
    }
    # calculate r and return r2
    else{
        my $r = (($n * $sum_xy) - ($sum_x * $sum_y))/sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2));        
    	return(sprintf("%.4f",$r),$n);                  
    }
}