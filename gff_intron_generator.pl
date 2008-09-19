#!/usr/bin/perl
#
# fragment_genome.pl 
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;

my $exon_count = 0;
my ($prev_exon_end_coord,$next_exon_start_coord,$intron_start,$intron_end);
my $mrna;
while(<>){
	my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split /\t/;	
	next unless ($feature eq "mRNA" || $feature eq "exon");
		
	# get mRNA id from comment field
	if($feature eq "mRNA"){
		$mrna = $comment;
		chomp($mrna);
		$mrna =~ s/^ID=//;
		$mrna =~ s/\;.*//;
		$exon_count = 0;
	}
	
	# count each exon
	if($feature eq "exon"){
		$exon_count++;
		
		# is this an odd numbered exon, if so record the end coordinate
		if($exon_count >1){
			$intron_start = $prev_exon_end_coord +1;
			$intron_end = $start -1;
			print "$location\t.\tintron\t$intron_start\t$intron_end\t.\t$strand\t.\tParent=$mrna\n";
			$prev_exon_end_coord = $stop;
		}
		else{
			$prev_exon_end_coord = $stop;
			
		}
	}
	print;
	
}               


__END__