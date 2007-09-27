#!/usr/bin/perl
#
# find_deletions_in_alleles.pl
#
# quick script to process some gff files to find deletion alleles overlapping with first introns
# by Keith Bradnam
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;


my %transcript2strand;
my %transcript2start;
my %transcript2stop;
my %transcript2chr;

open (GFF, "<transcript.gff") || die "Failed to open transcripts gff file\n\n";

while (<GFF>) {
	chomp;

	my @gff_line = split /\t/;
	my $chr        = $gff_line[0];
	my $start      = $gff_line[3];
	my $stop       = $gff_line[4];
	my $strand     = $gff_line[6];
	my $transcript = $gff_line[8];

	# clean up $transcript
	$transcript =~ s/Transcript //;
	$transcript =~ s/"//g;
#	print "$transcript\n";
	$transcript2strand{$transcript} = $strand;
	$transcript2chr{$transcript} = $chr;
	
	if($strand eq "+"){
		$transcript2start{$transcript} = $start;
		$transcript2stop{$transcript}  = $stop;
	}
	else{
		$transcript2start{$transcript} = $stop;
		$transcript2stop{$transcript}  = $start;
	}
}
close(GFF);

# for each transcript, keep track of which intron is the closest to the TSS
# i.e. which is the first intron
# might also want to track all introns that are 'near' (within 1000 bp?) to the TSS
my %transcript2closest_intron;
my %transcript2closest_distance;

open (GFF, "<intron.gff") || die "Failed to open intron gff file\n\n";

while (<GFF>) {
	chomp;

	my @gff_line = split /\t/;
	my $chr        = $gff_line[0];
	my $start      = $gff_line[3];
	my $stop       = $gff_line[4];
	my $strand     = $gff_line[6];
	my $transcript = $gff_line[8];

	# clean up $transcript
	$transcript =~ s/Transcript "//;
	$transcript =~ s/".*//g;
	
	# how far is intron from TSS?
	my $distance;
	if($strand eq "+"){
		$distance = $start - $transcript2start{$transcript};
	}
	else{
		$distance = $transcript2start{$transcript} - $stop;
	}
	
	# is this the first intron seen for this transcript?
	# if so set the distance from the TSS as the closest distance
	if(!defined($transcript2closest_distance{$transcript})){
		$transcript2closest_distance{$transcript} = $distance;
		$transcript2closest_intron{$transcript} = "$start-$stop";
	}
	# if not, check against the closest distance and change if it is closer
	elsif($distance < $transcript2closest_distance{$transcript}){
		$transcript2closest_distance{$transcript} = $distance;
		$transcript2closest_intron{$transcript} = "$start-$stop";
	}
}
close(GFF);


my $counter = 0; 

foreach my $key (keys (%transcript2closest_intron)){
#	print "TRANSCRIPT: $key $transcript2chr{$key} $transcript2strand{$key} $transcript2start{$key}-$transcript2stop{$key}\n";
#	print "FIRST INTRON: $transcript2closest_intron{$key} $transcript2closest_distance{$key} bp from TSS\n\n";

	open (GFF, "<deletion_alleles.gff") || die "Failed to open alleles gff file\n\n";
	while (<GFF>) {
		chomp;
		my @gff_line = split /\t/;
		my $chr    = $gff_line[0];
		next if ($chr ne $transcript2chr{$key});
	
		my $allele_start  = $gff_line[3];
		my $allele_stop   = $gff_line[4];
		my $allele        = $gff_line[8];
		$allele =~ s/Variation //;
		$allele =~ s/"//g;
		
		my ($intron_start,$intron_stop) = split(/-/,$transcript2closest_intron{$key});
		
		#does intron contain deletion???	
		if(($allele_start > $intron_start) && ($allele_stop < $intron_stop)){
			$counter++;

			print "$counter) TRANSCRIPT: $key $transcript2chr{$key} $transcript2strand{$key} $transcript2start{$key}-$transcript2stop{$key}\n";
			print "FIRST INTRON: $transcript2closest_intron{$key} ($transcript2closest_distance{$key} bp from TSS) contains ";			
			print "$allele ($allele_start-$allele_stop)\n\n";
			print "http://www.wormbase.org/db/seq/gbrowse/elegans/?name=$allele;class=Variation\n\n";
		}
	}	
	close(GFF);

}

