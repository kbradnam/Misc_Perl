#!/usr/bin/perl
#
# find_homopolymers.pl 
#
# A script to find homopolymer runs in FASTA sequences
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;
use Keith;


die "
find_homopolmers.pl - a program to find homopolymer runs (A,C,G or T)
in sequences in a FASTA file

Usage: find_homopolymers.pl <fasta file of worm genome> <min homopolymer length> <max homopolymer length>\n" unless @ARGV == 3;



my ($file,$min,$max) = @ARGV;

############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$file") || die "Can't open $file, please specify a valid FASTA sequence\n\n";

my $fasta = new FAlite(\*FILE);

# store counts of homopolymers in hash array structure:
# key is base (A,C,G, or T)
# value is an array, and the index position is a length of a homopolymer run
my %seqs;

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
                        
		my $seq = uc($entry->seq);
		foreach my $base qw(A C G T){
			for(my $i = $min; $i<=$max; $i++){
				my $count;
				# homopolymer runs must be flanked by a different character or be at the start or end
				# of the sequence. Homopolymers are removed from sequence to speed up subsequent look ups
				if ($base eq 'A'){
					$count = $seq =~ s/([CGTN]|\b)(${base}{$i})([CGTN]|\b)/${1}${3}/g;
				}
				elsif ($base eq 'C'){
					$count = $seq =~ s/([AGTN]|\b)(${base}{$i})([AGTN]|\b)/${1}${3}/g;
				}
				elsif($base eq 'G'){
					$count = $seq =~ s/([ACTN]|\b)(${base}{$i})([ACTN]|\b)/${1}${3}/g;	
				}
				else{
					$count = $seq =~ s/([ACGN]|\b)(${base}{$i})([ACGN]|\b)/${1}${3}/g;
				}		
			
				$seqs{$base}[$i] += $count;			
			}
		}
}
close(FILE);

# loop through each length, and then each base for output
for(my $i = $min; $i<=$max; $i++){
	foreach my $base qw(A C G T){
		my $percent = sprintf("%.1f",$seqs{$base}[$i] / ($seqs{'A'}[$i] + $seqs{'C'}[$i] + $seqs{'G'}[$i] + $seqs{'T'}[$i]) * 100);
		print "$file,$i,$base,$seqs{$base}[$i],$percent\n";			
	}
}

exit;