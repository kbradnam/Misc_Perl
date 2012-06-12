#!/usr/bin/perl
#
# find_breakpoints.pl
#
# A script to use mapping information of subreads (from BLAT) to identify possible breakpoints
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;

die "Usage: $0 <number of subreads that were BLATted> <BLAT psl file>" unless @ARGV == 2;

my ($subread_count, $file) = @ARGV;
# keep track of how many matches for each sequence there are (ideally should only be one)
my %seq_data;


my $min_indel_size = 3;

open(my $in, "<", $file) or die "Can't read from $file\n";

while(<$in>){
	# skip PSL header lines
	next unless m/^\d+/;
	chomp;
#	print;

	my @fields = split;
	my $seq_id = $fields[9];

	$seq_data{$seq_id}{'hits'}++;

	# following data will be overwritten if there is more than 1 hit per query sequence
	$seq_data{$seq_id}{'matches'}     = $fields[0];
	$seq_data{$seq_id}{'mismatches'}  = $fields[1];
	$seq_data{$seq_id}{'q_gap_count'} = $fields[4];
	$seq_data{$seq_id}{'q_gap_bases'} = $fields[5];
	$seq_data{$seq_id}{'s_gap_count'} = $fields[6];
	$seq_data{$seq_id}{'s_gap_bases'} = $fields[7];
	$seq_data{$seq_id}{'s_id'}        = $fields[13];
	$seq_data{$seq_id}{'s_start'}     = $fields[15];
	$seq_data{$seq_id}{'s_end'}       = $fields[16];
	$seq_data{$seq_id}{'blocks'}      = $fields[17];
	
	
#	print "INS\t$_\n"   if ($fields[4]  > 0 and $fields[6] == 0 and $fields[5] > $min_indel_size);
#	print "DEL\t$_\n"   if ($fields[4] == 0 and $fields[6]  > 0 and $fields[7] > 2);
#	print "INDEL\t$_\n" if ($fields[4]  > 0 and $fields[6]  > 0 and $fields[5] > 2 and $fields[7] > 2);
	 
	# calculate pair of subread
	my $subread_pair;
	if ($seq_id =~ m/(\S+)5$/){
		$subread_pair = $1 . "3";
	} else{
		$seq_id =~ m/(\S+)3$/;
		$subread_pair = $1 . "5";
	}
	$seq_data{$seq_id}{'subread_pair'} = $subread_pair;
}
close($in);


# track how many hits different subreads have (expect 1 hit per subread in an ideal world)
my %counts;

# track distance between subreads when they match on the same chromosome
my %distances;

# track how many subread pairs only:
# 	have 1 match (pair is unique)
#	are on different chromosomes
#	are on the same chromosome
# 	have a consistent read pairs (other read pair sequence maps to same chromosome)
#	might represent insertions in the shattered chromosome
#	might represent deletions in the shattered chromosome
#	might represent other complex indels in the shattered chromosome
# 	might represent a regular, straightforward 1-to-1 match (with < $min_indel_size insertions/mismatches

my $unique_pair_count = 0;
my $same_chr = 0;
my $diff_chr = 0;
my $consistent_read_pair = 0;
my $insertions = 0;
my $deletions = 0;
my $indels = 0;
my $normal_match = 0;
# track IDs of pairs on different chromosomes
my @different;


foreach my $key (sort {$seq_data{$b}{'hits'} <=> $seq_data{$a}{'hits'}} keys %seq_data){

	my $count = $seq_data{$key}{'hits'};

	if ($count == 1){
		$counts{1}++;
		
		# now check whether matching pair is unique
		my $subread_pair = $seq_data{$key}{'subread_pair'};
		if (exists $seq_data{$subread_pair}{'hits'} and ($seq_data{$subread_pair}{'hits'} == 1)){
			$unique_pair_count++;
			
			# now check whether they are on same or different chromosomes
			if ($seq_data{$key}{'s_id'} ne $seq_data{$subread_pair}{'s_id'}){
				$diff_chr++;
				push(@different, $key);
			} else {
				$same_chr++;

				# check that read pair is also on the same chromosome
				my $read_pair_1 = $key;
				my $read_pair_2 = $read_pair_1;
				$read_pair_2 = $seq_data{$key}{'subread_pair'};
				
				
				if ($key =~ m/_1/){
					$read_pair_1 =~ s/_1/_2/;
					$read_pair_2 =~ s/_1/_2/;

				} else{
					$read_pair_1 =~ s/_2/_1/;
					$read_pair_2 =~ s/_2/_1/;

				}
				# do the two paired subreads also have a unique hit on the same chromosome?
				next unless (exists $seq_data{$read_pair_1}{'hits'} and 
				    $seq_data{$read_pair_1}{'hits'} == 1   and
				    exists $seq_data{$read_pair_2}{'hits'} and 
				    $seq_data{$read_pair_2}{'hits'} == 1   and
				    $seq_data{$read_pair_1}{'s_id'} eq $seq_data{$key}{'s_id'} and
				    $seq_data{$read_pair_2}{'s_id'} eq $seq_data{$key}{'s_id'});

				$consistent_read_pair++;

				# check for insertions/deletions/indels
				if ($seq_data{$key}{'q_gap_count'} > 0 and 
				    $seq_data{$key}{'s_gap_count'} == 0 and 
				    $seq_data{$key}{'s_gap_bases'} > $min_indel_size){
					$insertions++;
				} 
				elsif ($seq_data{$key}{'q_gap_count'} == 0 and 
				    $seq_data{$key}{'s_gap_count'} > 0 and 
				    $seq_data{$key}{'s_gap_bases'} > $min_indel_size){
					$deletions++;
				} elsif ($seq_data{$key}{'q_gap_count'} > 0 and 
				    $seq_data{$key}{'s_gap_count'} > 0 and 
				    $seq_data{$key}{'q_gap_bases'} > $min_indel_size and 
				    $seq_data{$key}{'s_gap_bases'} > $min_indel_size){
					$indels++;
				} else{
					$normal_match++;
				}

				# check how far apart matches are
				my $distance = abs($seq_data{$key}{'s_start'} - $seq_data{$subread_pair}{'s_start'});
				if ($distance < 100){
					$distances{'< 100'}++;
				} elsif($distance < 1000){
					$distances{'< 1,000'}++;
				} elsif($distance < 10000){
					$distances{'< 10,000'}++;
				} elsif($distance < 100000){
					$distances{'< 100,000'}++;
				} elsif($distance <= 1000000){
					$distances{'< 1,000,000'}++;
				} else{
					$distances{'>= 1,000,000'}++;
				}
			}
		}
	} 
	elsif ($count == 2){
		$counts{'2'}++;
	}
	elsif ($count == 3){
		$counts{'3'}++;
	}
	elsif ($count < 11){
		$counts{'4-10'}++;
	}
	elsif ($count < 26){
		$counts{'11-25'}++;
	}
	elsif ($count < 51){
		$counts{'26-50'}++;
	}
	elsif ($count < 101){
		$counts{'51-100'}++;
	}
	else {
		$counts{'100+'}++;
	}

}
my $matching_subreads = keys(%seq_data);
my $unmatched_subreads = $subread_count - $matching_subreads;

my $matched_percent   = sprintf("%.2f", ($matching_subreads  / $subread_count) * 100);
my $unmatched_percent = sprintf("%.2f", ($unmatched_subreads / $subread_count) * 100);

print "\n";
print "There were $matching_subreads subreads that matched genome ($matched_percent%)\n";
print "There were $unmatched_subreads subreads that *didn't* match genome ($unmatched_percent%)\n\n";

print "Matches\tN\t%\n";

no warnings;
foreach my $key (sort {$a <=> $b} keys %counts){
	my $percent = sprintf("%.2f", ($counts{$key}/$matching_subreads) * 100);
	print "$key\t$counts{$key}\t$percent%\n";
}
use warnings;


my $percent;
print "\n";
print "$counts{'1'} (100.0%)\tNumber of subreads with 1 BLAT match\n";

$percent = sprintf("%.2f", ($unique_pair_count / $counts{'1'}) * 100);
print "$unique_pair_count ($percent%)\tNumber of subread pairs that have 1 BLAT match each\n";

$percent = sprintf("%.2f", ($same_chr / $unique_pair_count) * 100);
print "$same_chr ($percent%)\tNumber of subread pairs that have 1 BLAT match each on *same* chromosome\n";

$percent = sprintf("%.2f", ($consistent_read_pair / $unique_pair_count) * 100);
print "$consistent_read_pair ($percent%)\tNumber of subread pairs that have matching read pairs with 1 BLAT match each on *same* chromosome\n";

$percent = sprintf("%.2f", ($normal_match / $unique_pair_count) * 100);
print "$normal_match ($percent%)\tNumber of subread pairs that are suggestive of no indels in shattered chromosome sequence\n";

$percent = sprintf("%.2f", ($deletions / $unique_pair_count) * 100);
print "$deletions ($percent%)\tNumber of subread pairs that are suggestive of deletions in shattered chromosome sequence\n";

$percent = sprintf("%.2f", ($indels / $unique_pair_count) * 100);
print "$indels ($percent%)\tNumber of subread pairs that are suggestive of complex indels in shattered chromosome sequence\n";

$percent = sprintf("%.2f", ($insertions / $unique_pair_count) * 100);
print "$insertions ($percent%)\tNumber of subread pairs that are suggestive of insertions in shattered chromosome sequence\n";



no warnings;
foreach my $key ('< 100', '< 1,000', '< 10,000', '< 100,000', '< 1,000,000', '>= 1,000,000'){
	$percent = sprintf("%.2f", ($distances{$key}/$unique_pair_count) * 100);
	print "\t\t", sprintf("%12s", $key), "\t$distances{$key}\t$percent%\n";
}
use warnings;




# print out different IDs
#@different = sort @different;
#foreach my $id (@different){
#	print "$id\n";
#}