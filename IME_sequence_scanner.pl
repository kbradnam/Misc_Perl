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
$min = 1 if (!$min);
$max = 10000000 if (!$max);
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


# now loop through sequence grabbing the window of sequence

for(my $i = $min; $i < $length; $i+=$step){
	my $end = $i+$window-1;
    # get window of sequence
    my $temp_seq = substr($seq,$i-1,$window);

	print "\n";

	my ($start,$split,$split_seq);
	$split = $window / 4;
	$start = 1;

	# now split big window into three equal sized regions and analyze those separately	
	foreach my $j (1..4){
		my $start2 = $start+$i;
		my $end = $i+ $start + $split-1;
		print "$start2 - $end : ";
		
		$split_seq = substr($temp_seq,$start-1,$split);

		# write temp file with sequence
	    open(OUT,">/tmp/ime_seq") || die "Can't write to output file\n";
		print OUT ">tmp_seq\n";
		print OUT "$split_seq\n";
		close(OUT);
		
		# first motif test will use upstream motif frequencies, other tests will use transcript motif frequencies
		my $script_output;
		if($j == 1){
			$script_output = `~keith/Work/bin/where_is_motif.pl -species atu -mdensity -target /tmp/ime_seq -motif $motif`;
			$script_output =~ m/.*: (\d+)\/(\d+) (.*)/;
			print "$1,$2,$3\n";

		}
		elsif($j > 2){
			$script_output = `~keith/Work/bin/where_is_motif.pl -species att -mdensity -target /tmp/ime_seq -motif $motif`;
			$script_output =~ m/.*: (\d+)\/(\d+) (.*)/;
			print "$1,$2,$3\n";
		}
		else{
			print "\n";
		}
		
		$start += $split;
	}

	

	


	last if ($i > $max);
}

exit(0);







