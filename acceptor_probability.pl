#!/usr/bin/perl -w

use strict;

# the consensus splice acceptor sequence is TTTTCAG/R

my @matrix;

# fill position matrix for splice acceptor with probabilities of each
# base (probabilities from Kent & Zahler 2000 Genome Research paper)

&fill_matrix;

# read input file
open(IN, "<$ARGV[0]") || die "Can't open input file\n";
while(<IN>){
    chomp;
    # get rid of slash
    s/\///;
    my @acceptor = split //;
    my $prob;
    for(my $i=0; $i<@acceptor;$i++){
	if($i == 0){
	    $prob = $matrix[$i]{$acceptor[$i]};
	}
	else{
	    $prob *= $matrix[$i]{$acceptor[$i]};
	}
    }
    print "$_) $prob\n";
}


exit(0);

sub fill_matrix{

    $matrix[0]{"A"} = 0.28;
    $matrix[0]{"C"} = 0.08;
    $matrix[0]{"G"} = 0.06;
    $matrix[0]{"T"} = 0.58;
    
    $matrix[1]{"A"} = 0.05;
    $matrix[1]{"C"} = 0.03;
    $matrix[1]{"G"} = 0.02;
    $matrix[1]{"T"} = 0.90;
    
    $matrix[2]{"A"} = 0.01;
    $matrix[2]{"C"} = 0.01;
    $matrix[2]{"G"} = 0.003;
    $matrix[2]{"T"} = 0.97;
    
    $matrix[3]{"A"} = 0.09;
    $matrix[3]{"C"} = 0.16;
    $matrix[3]{"G"} = 0.08;
    $matrix[3]{"T"} = 0.97;
    
    $matrix[4]{"A"} = 0.03;
    $matrix[4]{"C"} = 0.84;
    $matrix[4]{"G"} = 0.001;
    $matrix[4]{"T"} = 0.13;
    
    $matrix[5]{"A"} = 1;
    $matrix[5]{"C"} = 0;
    $matrix[5]{"G"} = 0;
    $matrix[5]{"T"} = 0;
    
    $matrix[6]{"A"} = 0;
    $matrix[6]{"C"} = 0;
    $matrix[6]{"G"} = 1; 
    $matrix[6]{"T"} = 0;
    
    $matrix[7]{"A"} = 0.41;
    $matrix[7]{"C"} = 0.16;
    $matrix[7]{"G"} = 0.31;
    $matrix[7]{"T"} = 0.13;
}
