#!/usr/bin/perl
#
# chunk_FASTA.pl
#
# A script to split a FASTA file into separate files of approximately equal size
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "Usage: $0 <input fasta> <files-to-split-into>" unless (@ARGV == 2);

my ($input_file, $n) = @ARGV;

# count sequences in input file to determine number of seqs per split file
my $total_seq_count = `grep -c \">\" $input_file`;
chomp($total_seq_count);

my $target;
# will we have any remainder, if so increase size of first files
# try to ensure that there is always an even number in each file (just in case input file comes in pairs)
if (($total_seq_count / $n) % 2  == 0){
	$target = int($total_seq_count / $n);
} else{
	$target = int($total_seq_count / $n) + 1;
}

print "Will split $total_seq_count sequences in $input_file into $n files containing approximately $target sequences each\n";



# get rid of trailing fasta, fa etc file extensions from file name
my $file_prefix = $input_file;
$file_prefix =~ s/\.(dna|fasta|fa)$//;

# keep track of how many sequences we've seen
my $seq_count = 0;

# want to increment a file suffix to be part of output files
my $file_count = 1;
my $output_file = $file_prefix . ".chunk$file_count.fa";
my $out;
open ($out, ">", $output_file) or die "Can't create $output_file\n";

#### MAIN LOOP ####

# parse info from FASTA file
open(my $in,"$input_file") or die "Can't read from $input_file\n";
my $fasta = new FAlite($in);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry){
    my $seq = $entry->seq;
	my $header = $entry->def;

	$seq_count++;
	
	# do we need to split to a new file?
	# last file may end up with slightly more sequence, so don't go beyond $n
	if($seq_count > $target and $file_count != $n){
		close($out);
		$seq_count = 1;
		$file_count++;
		$output_file = $file_prefix . ".chunk$file_count.fa";
		open ($out, ">", $output_file) or die "Can't create $output_file\n";
	}
	print $out "$header\n$seq\n";

}
close($in);
close($out);