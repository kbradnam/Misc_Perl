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

# hashes to keep track of intron sizes for introns that are:
# 1) just in UTRs
# 2) just in CDSs
# 3) in both UTRs and in CDSs
# 4) in CDSs that have 5' UTR introns
# 5) in CDSs that don't have 5' UTR introns 
my %utr_data;
my %cds_data;
my %all_data;
my %cds_data2;
my %cds_data3;

# first loop through file once, just to work out which genes have 5' UTR introns
while(my $entry = $fasta->nextEntry) {	
	# process header
	my $header = $entry->def;
	$header =~ s/ .*//;
   	
	my ($gene,$intron,$tss_distance,$type) = split(/_/,$header);
	$gene =~ s/>//;

	# skip 3' UTR introns
	next if ($type eq "3UTR");
	
	# increment 5' UTR intron counter if necessary
	$gene2utr_count{$gene}++ if ($type eq "5UTR");
}
close(FILE);

open(FILE,"<$ARGV[0]") || die "Couldn't open $ARGV[0] file\n";
$fasta = new FAlite(\*FILE);

# now loop through each sequence again
while(my $entry = $fasta->nextEntry) {	
	# process header
	my $header = $entry->def;
	$header =~ s/ .*//;
   	
	my ($gene,$intron,$tss_distance,$type) = split(/_/,$header);
	$gene =~ s/>//;
	$intron =~ s/i//;
	
	# skip 3' UTR introns
	next if ($type eq "3UTR");
	
	# 1) Add size of any 5' UTR introns to %utr_data hash
	push(@{$utr_data{$intron}},length($entry->seq)) if ($type eq "5UTR");
	
	# 2) For any intron add to %all_data hash
	push(@{$all_data{$intron}},length($entry->seq));
	
	# 3) now look out for CDS introns 
	if($type eq "CDS"){
		# If the genes already has at least one 5' UTR intron then we need to offset the 
		# value of $intron by however many 5'UTR introns come before it
		if(defined($gene2utr_count{$gene})){
			push(@{$cds_data {$intron-$gene2utr_count{$gene}}},length($entry->seq));
			push(@{$cds_data2{$intron-$gene2utr_count{$gene}}},length($entry->seq));
		}
		else{
			push(@{$cds_data {$intron}},length($entry->seq));
			push(@{$cds_data3{$intron}},length($entry->seq));
		}	
		
	}
		

}
close(FILE) || die "Couldn't close $ARGV[0]\n";


# print header line for stats output
my $header = "Intron position,";
$header .= "N-all,Mean_all,STDEV,SE,95% CI,";
$header .= "N-cds,Mean_cds,STDEV,SE,95% CI,";
$header .= "N-UTR,Mean_UTR,STDEV,SE,95% CI,";
$header .= "N-cds2,Mean_cds2,STDEV,SE,95% CI,";
$header .= "N-cds3,Mean_cds3,STDEV,SE,95% CI";
print "$header\n";

# now process data to get averages at each position
foreach my $position (sort {$a <=> $b} (keys(%all_data))){
	
	# for each intron position, want 5 sets of stats
	# Only calculate if there is intron data at that position

	my @stats1;
	if (defined(${$all_data{$position}}[0])){
		@stats1 = &calc_stats(\%all_data,$position);
	}
	else{
		@stats1 = qw (0 0 0 0 0);
	}

	my @stats2;
	if (defined(${$cds_data{$position}}[0])){
		@stats2 = &calc_stats(\%cds_data,$position);
	}
	else{
		@stats2 = qw (0 0 0 0 0);
	}
	
	my @stats3;	
	if (defined(${$utr_data{$position}}[0])){
		@stats3 = &calc_stats(\%utr_data,$position); 
	}
	else{
		@stats3 = qw (0 0 0 0 0);
	}

	my @stats4;	
	if (defined(${$cds_data2{$position}}[0])){
		@stats4 = &calc_stats(\%cds_data2,$position); 
	}
	else{
		@stats4 = qw (0 0 0 0 0);
	}

	my @stats5;	
	if (defined(${$cds_data3{$position}}[0])){
		@stats5 = &calc_stats(\%cds_data3,$position); 
	}
	else{
		@stats5 = qw (0 0 0 0 0);
	}
		 
	my $output .= join (',',$position,@stats1,@stats2,@stats3,@stats4,@stats5);

	print "$output\n";

	# don't need to look at all intron positions
	last if ($position >14);
}



# subroutine to calculate basic statistics
sub calc_stats{
	my $ref = shift;
	my $position = shift;
	my $n = scalar(@{$$ref{$position}});
	my $mean = sprintf("%.2f",sum(@{$$ref{$position}})/$n);

	# need to bail out if $n = 1 because can't count stdev with sample size of 1
	return($n,$mean,"0","0","0") if ($n == 1);
 	
	my $stdev = sprintf("%.2f",sqrt(sum(map {($_ - $mean) ** 2} @{$$ref{$position}}) / ($n-1)));

	# standard error
	my $se = sprintf("%.2f",$stdev / sqrt($n));

	# 95% confidence limits of the mean
	my $ci = $se * 1.96;
	return($n,$mean,$stdev,$se,$ci);
}

	
exit(0);