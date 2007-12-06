#!/usr/bin/perl
#
# parse_IME_intron_data.pl
#
# A script to work through the A. thaliana intron files made by Genis and extract intron sizes 
# at each intron position, but to do this separately for 5' UTR intron data
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use List::Util qw(sum);

# has a file been specified?
die "Please specify an input file (in fasta format): parse_IME_intron_data.pl <filename>\n" if (!$ARGV[0]);

##############################################################
#
# Main loop
#
##############################################################

open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0] file\n";

my $fasta = new FAlite(\*FILE);

# count 5' UTR introns for each gene
my %gene2utr_count;

# hashes to keep track of intron sizes for introns that are just in UTRs, just in CDSs, or in both
my %utr_data;
my %cds_data;
my %all_data;

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {	
	# process header
	my $header = $entry->def;
   	my ($gene,$intron,$tss_distance,$type) = split(/_/,$header);
	$gene =~ s/>//;
	$intron =~ s/i//;
	
	# skip 3' UTR introns
	next if ($type eq "3UTR");
	
	# increment 5' UTR intron counter if necessary
	$gene2utr_count{$gene}++ if ($type eq "5UTR");
	
	# 1) Add size of any 5' UTR introns to %utr_data hash
	push(@{$utr_data{$intron}},length($entry->seq)) if ($type eq "5UTR");
	
	# 2) For any intron add to %all_data hash
	push(@{$all_data{$intron}},length($entry->seq));
	
	# 3) now look out for CDS introns for genes for which we have already seen at least one 5' UTR intron
	# for these introns we need to offset the value of $intron by however many 5'UTR introns come before it
	if(($type eq "CDS") && defined($gene2utr_count{$gene})){
		push(@{$cds_data{$intron-$gene2utr_count{$gene}}},length($entry->seq));
	}
}
close(FILE) || die "Couldn't close $ARGV[0]\n";

# now process data to get averages at each position
foreach my $position (sort {$a <=> $b} (keys(%all_data))){
	
	# for each intron position, want 3 sets of stats
	my @stats1 = &calc_stats(\%all_data,$position);
	my @stats2 = &calc_stats(\%cds_data,$position);
	my @stats3;
	if (defined(${$utr_data{$position}}[0])){
		@stats3 = &calc_stats(\%utr_data,$position); 
	}
	else{
		@stats3 = "";
	}
	# join stats together in one big happy list
	my $output = join (',',$position,@stats1,@stats2,@stats3);
	print "$output\n";
	last if ($position >14);
}



# subroutine to calculate basic statistics
sub calc_stats{
	my $ref = shift;
	my $position = shift;
	print "$position\n\n";
	my $n = scalar(@{$$ref{$position}});
	my $mean = sprintf("%.2f",sum(@{$$ref{$position}})/$n);
	my $stdev = sprintf("%.2f",sqrt(sum(map {($_ - $mean) ** 2} @{$$ref{$position}}) / ($n-1)));

	# standard error
	my $se = sprintf("%.2f",$stdev / sqrt($n));

	# 95% confidence limits of the mean
	my $ci = $se * 1.96;
	return($n,$mean,$stdev,$se,$ci);
}

	
exit(0);