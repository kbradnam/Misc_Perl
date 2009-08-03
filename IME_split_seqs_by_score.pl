#!/usr/bin/perl
#
# IME_split_seqs_by_score.pl
#
# A script to take a FASTA file of Arabidopsis introns (with IME info in header), plus a CSV file containing IMEter scores
# and create new output files containing introns of fixed scores
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;

die "Specify a CSV file and a FASTA file on the command line\n" if (@ARGV < 2);

my ($csv,$seqs) = @ARGV;

my %intron_data;
my %gene_data;

# parse info from CSV file
open(CSV,"<$csv") or die "Can't open $csv file\n";
while(<CSV>){
	next if m/^Key/; # skip header line
	my ($id, $type, $position, $length, $distance, $ime_v1, $ime_v2, $motifs_t0, $motifs_t3, $motifs_t6) = split(/,/);
	# load data into hash, don't need everything for now
	$intron_data{$id}{type} = $type;
	$intron_data{$id}{position} = $position;
	$intron_data{$id}{score} = $ime_v2;
}
close(CSV);

# parse info from FASTA file
open(FASTA,"$seqs") || die "Can't open $seqs\n";
my $fasta = new FAlite(\*FASTA);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry){
    my $seq = $entry->seq;
	my $header = $entry->def;
	my ($id) = $header =~ m/>(\w+\.\d+_i\d+)_/;
	$intron_data{$id}{seq} = $seq;
	$intron_data{$id}{header} = $header;
}
close(FASTA);

open(A,">at_ime_0.fa")       or die "Can't open at_ime_0.fa\n";
open(B,">at_ime_0.1_4.9.fa")  or die "Can't open at_ime_0.1_4.9.fa\n";
open(C,">at_ime_5_9.9.fa")   or die "Can't open at_ime_5_9.9.fa\n";
open(D,">at_ime_10_19.9.fa") or die "Can't open at_ime_10_19.9.fa\n";
open(E,">at_ime_20_plus.fa") or die "Can't open at_ime_20_plus.fa\n";

foreach my $key (keys(%intron_data)){
	my $score = $intron_data{$key}{score};
	if($score == 0){
		print A "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score > 0 && $score < 5){
		print B "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 5 && $score < 10){
		print C "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 10 && $score < 20){
		print D "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	else{
		print E "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
}

close(A);
close(B);
close(C);
close(D);
close(E);