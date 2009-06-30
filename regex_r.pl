#!/usr/bin/perl
#
# regex_breeder.pl
#
# A script to 
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use List::Util qw(sum);

die "Specify a regular expression for a motif, and at least one file name containing introns\n" if (!$ARGV[1]);

# turn on autoflush
$| = 1;

# want to keep track of sequences from input file, plus their expression level, plus counts of motifs
my @seqs;
my @expression;
my @counts;

my ($motif,@files) = @ARGV;

foreach my $file (@files){
	# reset arrays
	@counts = ();
	@seqs = ();
	@expression = ();

	# read sequence files and also extract expression values
	read_file($file);

	# now loop through each sequence
	for (my $i=0; $i < @seqs; $i++){
		my $seq = $seqs[$i];

		# count motif in sequence
		my $count = 0;
		$count = $seq =~ s/($motif)/$1/g;
		($count = 0) if (!$count);

		$counts[$i] = $count;			
	}

	# can now calculate correlation (r) for that motifset
	my ($r,$n) = r2(\@expression,\@counts);

	print "$file: n = $n r = $r\n";
	
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