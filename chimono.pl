#!/usr/bin/perl
#
# chimono.pl 
#
# A script to perform a chi-squared test of the mononucleotide frequencies of two FASTA files
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;



# sanity checks
die "Usage: chidi.pl <file1> <file2>\n" if (!$ARGV[1]);

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
	# to count mononucleotides, loop through sequence, take 2 bp and increment the hash counter
    $file1_ob{"A"} += ($seq =~ tr/A/A/); 
    $file1_ob{"C"} += ($seq =~ tr/C/C/); 
    $file1_ob{"G"} += ($seq =~ tr/G/G/); 
    $file1_ob{"T"} += ($seq =~ tr/T/T/); 
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
	# to count mononucleotides, loop through sequence, take 2 bp and increment the hash counter
    $file2_ob{"A"} += ($seq =~ tr/A/A/); 
    $file2_ob{"C"} += ($seq =~ tr/C/C/); 
    $file2_ob{"G"} += ($seq =~ tr/G/G/); 
    $file2_ob{"T"} += ($seq =~ tr/T/T/);
}
close(FILE);


############################################################
# Perform chi-squared test
############################################################

# need total of all counts in both sequences, plus totals of 'rows' in chi-square table

my $total;
my $row1;
my $row2;

foreach my $nt qw (A C G T){
	$row1  += $file1_ob{$nt};
	$row2  += $file2_ob{$nt};
	$total += ($file1_ob{$nt} + $file2_ob{$nt});
}


# now calculate expected values
foreach my $nt qw (A C G T){
	# calculate (column total * row total) / $total
	$file1_ex{$nt} = sprintf("%.2f",(($file1_ob{$nt}+$file2_ob{$nt}) * $row1) / $total);
	$file2_ex{$nt} = sprintf("%.2f",(($file1_ob{$nt}+$file2_ob{$nt}) * $row2) / $total);	
}

# now calculate chi-squared values
my ($chi1,$chi2);
my $chi_total;
print "\tObs1\tExp2\t\tChi1\tObs2\tExp2\t\tChi2\n";
foreach my $nt qw (A C G T){
	$chi1 = sprintf("%.2f",(($file1_ob{$nt} - $file1_ex{$nt})**2)/$file1_ex{$nt});
	$chi2 = sprintf("%.2f",(($file2_ob{$nt} - $file2_ex{$nt})**2)/$file2_ex{$nt});	
	print "$nt\t$file1_ob{$nt}\t$file1_ex{$nt}\t$chi1\t$file2_ob{$nt}\t$file2_ex{$nt}\t$chi2\n";

	$chi_total += ($chi1+$chi2);
}

printf  "Chi squared value = %6.2f\n", $chi_total;
				
print "Significance level at 5% = 7.81\n";
print "Significance level at 1% = 11.34\n";


exit(0);
