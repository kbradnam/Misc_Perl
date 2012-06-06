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

use strict; use warnings;
use FAlite;

die "Usage: $0 <subread length> <filename>" unless @ARGV == 2;

my ($subread_length, $file) = @ARGV;

# parse info from FASTA file
open(my $fasta_file,"$file") or die "Can't open $file\n";
my $fasta = new FAlite(\*$fasta_file);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry){
    my $seq = $entry->seq;
	my $header = $entry->def;
}
close($fasta_file);