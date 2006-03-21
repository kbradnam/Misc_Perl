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
use Getopt::Long;
use lib "/Korflab/lib/perl";
use FAlite;


########################
# Command line options
########################

my $mono;  # create random sequence based on mononucleotide frequencies
my $dinuc; # create random sequence based on dinucleotide frequencies

GetOptions ("mono"     => \$mono,
            "dinuc"    => \$dinuc);

# sanity check
die "Need to specify -mono OR -dinuc option\n" if (($mono && $dinuc) || (!$mono && !$dinuc));


# Need set of expected nucleotide frequencies to compute log likelihood scores
# two hashes.  One is mono nucleotide frequencies of confirmed introns (used for
# first base of random intron). The second is dinucleotide frequencies (used for
# rest of intron)

my %mono = ("A" => "0.33339","C" => "0.16194",
	    "G" => "0.16021","T" => "0.34446");

my %dinucs = ("A" => ["0.457743938","0.58825315","0.729569921"],
	      "C" => ["0.351293566","0.542848977","0.711407452"],
	      "G" => ["0.332089564","0.51511319","0.704411625"],
	      "T" => ["0.205221746","0.373847922","0.534890808"]);

open(IN,"<$ARGV[0]") || die "Couldn't open $ARGV[0] file\n";
open(OUT,">${ARGV[0]}.random") || die "Couldn't write output file\n";

my $fasta = new FAlite(\*IN);

# loop through each sequence in input file

while(my $entry = $fasta->nextEntry) {
    my $header = $entry->def;
    my $seq = $entry->seq;
    print OUT "\n$header\n";
    my $length = length($seq);
    
    #  make new random intron corresponding to real intron length
    # first base of random intron is always based on momonucleotide frequencies
    
    my $counter = 1;
    my $base;
    my $rand = rand(1);
    if   ($rand <= 0.33339){$base = "A";}
    elsif($rand <= 0.49533){$base = "C";}
    elsif($rand <= 0.65554){$base = "G";}
    else{$base = "T";}
    
    print OUT "$base";

    my $last_base = $base;
    
    for (my $i = 1; $i < length($seq); $i++) {
	if($counter == 60){
	    print OUT "\n";
	    $counter = 0;
	}
	$counter++;
	$rand = rand(1);

	# this bit varies depending on whether we are using mono- or di-nucleotide frequencies
	if($mono){
	    if   ($rand <= $dinucs{$last_base}[0]){$base = "A";}
	    elsif($rand <= $dinucs{$last_base}[1]){$base = "C";}
	    elsif($rand <= $dinucs{$last_base}[2]){$base = "G";}
	    else{$base = "T";}
	}
	elsif($dinuc){
	    if   ($rand <= 0.33339){$base = "A";}
	    elsif($rand <= 0.49533){$base = "C";}
	    elsif($rand <= 0.65554){$base = "G";}
	    else{$base = "T";}
	}	    
	
	print OUT "$base";
	$last_base = $base;
    }
#    print OUT "\n";
}


     
close(IN) || die "Couldn't close $ARGV[0]\n";
close(OUT) || die "Couldn't close output file\n";

exit(0);
