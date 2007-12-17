#!/usr/bin/perl
#
# chidi.pl 
#
# A script to perform a chi-squared test of the dinucleotide frequencies of two FASTA files
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;



# sanity checks
die "Usage: chidi.pl <file1> <file2>\n" if (!$ARGV[1]);

my @dinucs = qw (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);

# hashes for obersered and expected dinucleotide frequencies of both files

my %file1_ob;
my %file2_ob;
my %file1_ex;
my %file2_ex;
								
############################################################
# Read sequence file 1
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0]\n";
my $fasta = new FAlite(\*FILE);

# loop through each sequence in file 1
while(my $entry = $fasta->nextEntry) {	
	my $seq = uc($entry->seq);
	# to count dinucleotides, loop through sequence, take 2 bp and increment the hash counter
	foreach my $i (0..length($seq)){
	    my $tmp = substr($seq,$i,2);		
		$file1_ob{$tmp}++;
	}
}
close(FILE);


############################################################
# Read sequence file 2
############################################################

open(FILE,"$ARGV[1]") || die "Can't open $ARGV[1]\n";
$fasta = new FAlite(\*FILE);

# loop through each sequence in file 1
while(my $entry = $fasta->nextEntry) {	
	my $seq = uc($entry->seq);
	# to count dinucleotides, loop through sequence, take 2 bp and increment the hash counter
	foreach my $i (0..length($seq)){
	    my $tmp = substr($seq,$i,2);		
		$file2_ob{$tmp}++;
	}
}
close(FILE);


############################################################
# Perform chi-squared test
############################################################

# need total of all counts in both sequences, plus totals of 'rows' in chi-square table

my $total;
my $row1;
my $row2;

foreach my $di (@dinucs){
	$row1  += $file1_ob{$di};
	$row2  += $file2_ob{$di};
	$total += ($file1_ob{$di} + $file2_ob{$di});
}


# now calculate expected values

foreach my $di (@dinucs){
	# calculate (column total * row total) / $total
	$file1_ex{$di} = sprintf("%.2f",(($file1_ob{$di}+$file2_ob{$di}) * $row1) / $total);
	$file2_ex{$di} = sprintf("%.2f",(($file1_ob{$di}+$file2_ob{$di}) * $row2) / $total);	
}

# now calculate chi-squared values
my ($chi1,$chi2);
my $chi_total;
print "\tObs1\tExp2\t\tChi1\tObs2\tExp2\t\tChi2\n";
foreach my $di (@dinucs){
	$chi1 = sprintf("%.2f",(($file1_ob{$di} - $file1_ex{$di})**2)/$file1_ex{$di});
	$chi2 = sprintf("%.2f",(($file2_ob{$di} - $file2_ex{$di})**2)/$file2_ex{$di});	
	print "$di\t$file1_ob{$di}\t$file1_ex{$di}\t$chi1\t$file2_ob{$di}\t$file2_ex{$di}\t$chi2\n";

	$chi_total += ($chi1+$chi2);
}

printf  "Chi squared value = %6.2f\n", $chi_total;
				
print "Significance level at 5% = 25.00\n";
print "Significance level at 1% = 30.58\n";


exit(0);
