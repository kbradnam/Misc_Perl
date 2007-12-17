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



########################
# Command line options
########################

my $file; 	  # input FASTA file
my $pattern;  # pattern to match in FASTA headers
my $out;	  # optional name for an output file

GetOptions("file=s"    => \$file,
		   "out=s"     => \$out,
           "pattern=s" => \$pattern);

die "Use -file option to specify a FASTA sequence file\n" if (!defined($file));
die "Use -pattern option to specify a pattern in FASTA header\n" if (!defined($pattern));

##############################################################
#
# Main loop
#
##############################################################

open(FILE,"<$file") || die "Couldn't open $file\n";

my $output = "${file}.new";
($output = $out) if ($out);

open(OUT, ">$output") || die "Couldn't create output file\n";

my $fasta = new FAlite(\*FILE);

# count 5' UTR introns for each gene
my %gene2utr_count;

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {	

	# process header
	my $header = $entry->def;
   	if($header =~ m/$pattern/){
		my $seq = $entry->seq;
		print OUT "$header\n$seq\n";
	}

}
close(FILE);
close (OUT);

exit(0);