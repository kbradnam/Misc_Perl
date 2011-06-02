#!/usr/bin/perl
#
# IME_split_seqs_by_score.pl
#
# A script to take a FASTA file of Arabidopsis introns (with IME info in header), plus a TSV file containing IMEter scores
# and create new output files containing introns of fixed scores
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;

die "Specify a TSV file (produced by the_imeter.pl) and a FASTA file on the command line\n" if (@ARGV < 2);

my ($tsv,$seqs) = @ARGV;

my %intron_data;
my %gene_data;

# parse info from CSV file
open(TSV,"<$tsv") or die "Can't open $tsv file\n";
while(<TSV>){
	chomp;
	my ($id, $ime_v2) = split(/\t/);
	
	my ($key, $intron, $position, $type) = $id =~ m/([A-Z0-9\.]+)_i(\d+)_(\d+)_(CDS|5UTR|3UTR)/;

	# trim $id
	$id =~ s/_\d+_(CDS|5UTR|3UTR)//;

	# load data into hash, don't need everything for now
	$intron_data{$id}{type} = $type;
	$intron_data{$id}{position} = $position;
	$intron_data{$id}{score} = $ime_v2;
}
close(TSV);


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
#	print "$id\n$seq\n$header\n\n";
}
close(FASTA);


open(A,">at_ime_0.fa")       or die "Can't open at_ime_0.fa\n";
open(B,">at_ime_0.1_4.9.fa")  or die "Can't open at_ime_0.1_4.9.fa\n";
open(C,">at_ime_5_9.9.fa")   or die "Can't open at_ime_5_9.9.fa\n";
open(D,">at_ime_10_19.9.fa") or die "Can't open at_ime_10_19.9.fa\n";
open(E,">at_ime_20_plus.fa") or die "Can't open at_ime_20_plus.fa\n";

__END__
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

__END__
# alternative scoring scheme to use (splits into equal sized datasets, at least for TAIR10 data)
open(A,">at2_ime_0_1.68.fa")       or die "Can't open at_ime_0.fa\n";
open(B,">at2_ime_1.69_3.53.fa")  or die "Can't open at_ime_0.1_4.9.fa\n";
open(C,">at2_ime_3.54_5.64.fa")   or die "Can't open at_ime_5_9.9.fa\n";
open(D,">at2_ime_5_65_9.25.fa") or die "Can't open at_ime_10_19.9.fa\n";
open(E,">at2_ime_9.26_77.9.fa") or die "Can't open at_ime_20_plus.fa\n";


foreach my $key (keys(%intron_data)){
	my $score = $intron_data{$key}{score};
	if($score >= 0 && $score < 1.69){
		print A "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 1.69 && $score < 3.54){
		print B "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 3.54 && $score < 5.65){
		print C "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 5.65 && $score < 9.26){
		print D "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
	elsif($score >= 9.26){
		print E "$intron_data{$key}{header} IME_v2 = $score\n",Keith::tidy_seq($intron_data{$key}{seq}),"\n";
	}
}
