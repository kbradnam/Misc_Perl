#!/usr/bin/perl 
#
# IME_sequence_scanner.pl
#
# A script to look for IME like regions in genomic DNA
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


# command line options
my $window;     # how big a size window
my $step;       # for sliding windows
my $motif;	    # path to a valid NestedMICA XMS motif file
my $min;		# min starting coordinate
my $max;        # max starting coordinate
my $sequence;   # path to input file of FASTA sequences

GetOptions ("window=i"     => \$window,
			"min=i"		   => \$min,
			"max=i"        => \$max,
			"step=i"       => \$step,
			"motif=s"      => \$motif,
			"sequence=s"   => \$sequence);

# check command-line options
die "Please specify a valid NestedMica .xms motilf file\n" if (!$motif);

# set some defaults
$min = 0 if (!$min);
$window = 400 if (!$window);
$step = 200 if (!$step);



# open sequence file and process
open(IN,"<$sequence") || die "Couldn't open $sequence\n\n";

my $fasta = new FAlite(\*IN);

# grab first entry from file (assumes one sequence per file, i.e. chromosomes)
my $entry = $fasta->nextEntry;
my $seq = $entry->seq;
my $length = length($seq);

close(IN);

for(my $i = $min; $i < $length; $i+=$step){
	my $end = $i+$window;
    # get window of sequence
    my $temp_seq = substr($seq,$i,$window);
#	print "$temp_seq\n";

	# write temp file with sequence
    open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
	print OUT ">$i\n";
	print OUT "$temp_seq\n";
	close(OUT);
	
	# need to use upstream motif background frequencies and then transcript motif frequencies
	my $data = `~keith/Work/bin/where_is_motif.pl -species atg -mdensity -target /tmp/ime_seq -motif $motif`;
	$data =~ m/.*: (\d+)\/(\d+) (.*)/;
	print "$i - $end: $1,$2,$3\n";

	last if ($i > $max);
}

exit(0);







