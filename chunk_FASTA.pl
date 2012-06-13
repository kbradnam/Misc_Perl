#!/usr/bin/perl
#
# <name>.pl
#
# A script to 
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;

# parse info from FASTA file
open(FASTA,"$seqs") || die "Can't open $seqs\n";
my $fasta = new FAlite(\*FASTA);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry){
    my $seq = $entry->seq;
	my $header = $entry->def;
}
close(FASTA);
