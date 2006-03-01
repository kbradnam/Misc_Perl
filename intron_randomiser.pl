#!/usr/bin/perl -w
#
# intron_randomiser.pl
#
# A script to take a file of real worm introns and then write a new file of random sequences
# with the same length as the real introns, but with nucleotide frequencies based on real
# intron nucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use lib "/Korflab/lib/perl";
use FAlite;



# Need sets of expected nucleotide frequencies to compute log likelihood scores
# one set is based on (confirmed) intron base frequencies, the second is based on 
# the entire worm genome
my %expected_intronic = ("A" => "0.33339","C" => "0.16194",
			 "G" => "0.16021","T" => "0.34446");
my %expected_genomic  = ("A" => "0.32280","C" => "0.17733",
			 "G" => "0.17709","T" => "0.32279");


open(IN,"<$ARGV[0]") || die "Couldn't open $ARGV[0] file\n";
open(OUT,">${ARGV[0]}.random") || die "Couldn't write output file\n";

my $fasta = new FAlite(\*IN);

# loop through each sequence in input file

while(my $entry = $fasta->nextEntry) {
    my $header = $entry->def;
    my $seq = $entry->seq;
    print OUT "\n$header\n";
    my $length = length($seq);

    my $counter = 0;
    # loop through sequence in windows equal to motif size
      for (my $i = 0; $i < length($seq); $i++) {
	  if($counter == 60){
	      print OUT "\n";
	      $counter = 0;
	  }
	  $counter++;

	  my $rand = rand(1);
	  my $base;
	  if($rand < 0.33340){$base = "A";}
	  elsif($rand < 0.49534){$base = "C";}
	  elsif($rand < 0.65555){$base = "G";}
	  else{$base = "T";}
	  print OUT "$base";
      }
    print OUT "\n";
}

close(IN) || die "Couldn't close $ARGV[0]\n";
close(OUT) || die "Couldn't close output file\n";

exit(0);
