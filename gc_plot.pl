#!/usr/bin/perl -w

# Keith Bradnam, 1998

use strict;
use Getopt::Long;


# set window and step size
my $window = 3000;
my $step_size = $window;
#my $step_size = 50000;


# open sequence file and process
open(IN,"$ARGV[0]") || die "Couldn't open file\n\n";

# store chromosome sequence
my $seq;

while(my $temp=<IN>){
    if($temp !~ />/){
	chomp($temp);
	$temp =~  tr/a-z/A-Z/;
	$seq .= $temp;
    }
}
close(IN);


# now need to loop through sequence in fixed window sizes
my $total_length = length($seq);

# to store % A, %T etc.
my ($A,$T,$C,$G,$N);

# expected values for bases
my ($e_A, $e_T, $e_C, $e_G);

$e_A = 0.32085;
$e_T = 0.32169;
$e_C = 0.17885;
$e_G = 0.17861;

# gonna track cumulative deviations from expected whole chromosome values
my $A_dev;
my $T_dev;
my $C_dev;
my $G_dev;


my $A_sum;
my $T_sum;
my $C_sum;
my $G_sum;


for(my $i = 0; $i < $total_length; $i+=$step_size){
    # get window of sequence
    my $temp_seq = substr($seq,$i,$window);

    # Calculate A,T,C,G, N
    $A=$T=$C=$G=$N=0;
    
    $A = $temp_seq =~ tr/A/A/; 
    $T = $temp_seq =~ tr/T/T/;
    $C = $temp_seq =~ tr/C/C/;
    $G = $temp_seq =~ tr/G/G/;
    $N = $temp_seq =~ tr/N/N/;

    
    my $GC_percent  = (($G+$C)/$window)*100;
    my $A_percent   = ($A/$window)*100;
    my $T_percent   = ($T/$window)*100;
    my $C_percent   = ($C/$window)*100;
    my $G_percent   = ($G/$window)*100;
    my $N_percent   = ($N/$window)*100;

    $A_dev = $A - ($e_A*$window);
    $T_dev = $T - ($e_T*$window);
    $C_dev = $C - ($e_C*$window);
    $G_dev = $G - ($e_G*$window);

    $A_sum += ($A_dev*$window)/$total_length;
    $T_sum += ($T_dev*$window)/$total_length;
    $C_sum += ($C_dev*$window)/$total_length;
    $G_sum += ($G_dev*$window)/$total_length;
    
    my $start = $i+1;
    my $stop  = $start+$step_size;

    
    print "$start,$stop,$A_sum,$T_sum,$C_sum,$G_sum\n";
}

exit(0);







