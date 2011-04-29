#!/usr/bin/perl
# A script to quickly extract CDS features from a pair of GFF & FASTA files
#
# Written by Keith Bradnam
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Getopt::Std;
use FAlite;
use Keith;
use vars qw($opt_f $opt_g);
getopts('f:g:');

my $usage = " 
Usage: $0 -f <FASTA file - can be gzipped> -g <GFF file>

";

die "$usage" if (!$opt_f or !$opt_g);

my ($fasta_file, $gff_file) = ($opt_f, $opt_g);

# will store input sequences in a hash: key is chromosome number (0,1,2), value is sequence
my %seqs;

# Read sequences from FASTA file
read_file($fasta_file);	

# now process GFF file to extract required sequecnes
parse_gff($gff_file);


exit();
	
	

###############
#			  #
# Subroutines #
#			  #
###############

# read sequences from FASTA file to a hash
sub read_file {
	
	my ($file) = @_;
	
	# if dealing with gzipped FASTA file, treat differently
	my $input;
    if($file =~ m/\.gz$/){
            open($input, "gunzip -c $file |") or die "Can't open a pipe to $file\n";
    } else{
            open($input, "<", "$file") or die "Can't open $file\n";
    }

	my $fasta_file = new FAlite(\*$input);

	while(my $entry = $fasta_file->nextEntry){
		my ($chr) = $entry->def =~ m/\.chr(\d+)/;
		$seqs{$chr} = uc $entry->seq; #everything after '.' on fasta header
	}
	close $input;
}	
	
sub parse_gff {

	my ($gff_file) = @_;
		
	# Process GFF file
	open my $gff, "<", $gff_file or die "Can't open $gff_file $!\n";

	# keep track of CDS sequence and it's strand, plus how many exons we've seen
	my %cds_to_strand;
	my %cds_to_seq;
	my %cds_to_exons;
	
	while (<$gff>) {
		# only want to consider CDS entries that belong to each gene
		next unless ($_ =~ m/CDS/);
	
		my @gff_line = split(/\t/, $_);

		# extract gene number to form ID
		my ($id) = $gff_line[8] =~ m/gene_index (\d+)/;
				
		# add chromosome, strand, length, and coordinate info to hash for each CDS
		my ($chr) = $gff_line[0] =~ m/chr(\d+)/;
		my $strand = $gff_line[6];
		my ($start, $end) = ($gff_line[3], $gff_line[4]);
		my $length = $end - $start + 1;
		
		# extract sequence of exon
		my $seq = substr ($seqs{$chr}, $start - 1, $length);
		$cds_to_seq{$id}   .= $seq;
		$cds_to_strand{$id} = $strand;
		$cds_to_exons{$id}++;
	}	
	close $gff;

	# loop over CDSs to print them out
	foreach my $key (sort {$a <=> $b} keys(%cds_to_seq)){
		my $seq = $cds_to_seq{$key};
		my $length = length($seq);
		
		# reverse complement sequence?
		$seq = Keith::revcomp($seq) if ($cds_to_strand{$key} eq '-');

		# tidy sequence
		$seq = Keith::tidy_seq($seq);

		print ">$key strand=$cds_to_strand{$key} exons=$cds_to_exons{$key} length=$length\n$seq\n";
	}
	
}

