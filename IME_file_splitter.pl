#!/usr/bin/perl
#
# split_FASTA.pl
#
# A script to separatea FASTA file into separate files based on matching a pattern in the FASTA header
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use Getopt::Long;

my $min = 1;
foreach my $max qw(250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 3750 4000 4250 4500 4750 5000){
	
	print "Processing $min - $max\n";
	
	open(OUT, ">${ARGV[1]}_${min}_${max}.fa") || die "Couldn't create output file\n";

	open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0]\n";

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {	

		# process header
		my $header = $entry->def;
	   	$header =~ m/_i\d+_(\d+)_/;
		my $distance = $1;
		
		if(($distance >= $min) && ($distance <= $max)){
			my $seq = $entry->seq;
			print OUT "$header\n$seq\n";			
		}

	}
	close(FILE);
	close (OUT);

	$min += 250;
}


exit(0);