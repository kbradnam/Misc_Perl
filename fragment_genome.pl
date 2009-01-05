#!/usr/bin/perl
#
# fragment_genome.pl 
#
# a script to fragment a FASTA sequence file into smaller FASTA sequences of a specified size
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
fragment_genome.pl - a program shredding a genome sequence (or any fasta file) into randomly
selected oligomers of a desired size.

Usage: fragment_genome.pl <fasta file of worm genome> <word size>\n" unless @ARGV == 2;


my ($FASTA,$WORD) = @ARGV;


open(FILE,"$FASTA") || die "Can't open $FASTA\n";
my $fasta = new FAlite(\*FILE);

my $genome;

# loop through each sequence in input file and add to one big string 
while(my $entry = $fasta->nextEntry) {
	my $seq = uc($entry->seq);  
	$genome .= $seq;

}
close(FILE);

my $length = length($genome); 

for(my $i = 0; $i <= ($length/$WORD);$i++){
	
	# choose a random position (must be less than size of sequence minus the word size)
	my $rand = int(rand(1) * ($length-$WORD-1));

	# extract fragment of sequence from genome and convert back to string
	my $fragment = substr($genome,$rand,$WORD);
	
	print ">random_read_$i\n$fragment\n";
}



__END__

