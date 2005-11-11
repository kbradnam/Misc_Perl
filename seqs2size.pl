#!/usr/bin/perl -w
#
# seqs2size.pl
#
# a quick script to calculate the size of sequences in a
# fasta file
#
# by Keith Bradnam, 14th September 2005
##############################################

use strict;

my $flag;
my $seq;
my $length;

 SEQ: while(<>){
     if(/^>/){
	 if($flag){
	     $length = length($seq);
	     print "$length\n";
	 }
	 $seq = "";
	 $flag = 1;
        next SEQ;
     }
     $seq .= $_;
}

# get last sequence in file
$length = length($seq);
print "$length\n";

