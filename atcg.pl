#!/usr/bin/perl
#
# atcg.pl 
#
# A script to count GC content plus mono- and dinucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;

my $help;    # print help
my $dinuc;   # calculate dinucleotide frequencies as well as mononucleotide frequencies
my $chi;     # perform chi-squared analysis on dinucleotide frequencies

GetOptions ("dinuc"   => \$dinuc,
	    "chi"     => \$chi,
	    "help"    => \$help);

# sanity checks
die "Can only use -chi option if -dinuc option is chosen\n" if ($chi && !$dinuc);



if($help){
	print "atcg.pl - a program to calculate nucleotide and dinucleotide composition\n\n";
    print "1) atcg.pl can analyse multiple sequences, but will concatenate the results together.\n";
    print "2) use -dinuc option to also calculate dinucleotide frequencies\n";
    print "3) use -chi option to perform chi-squared analysis of dinucleotide frequencies\n";
    print "4) 'Other' counter is a check for any non ATCG, or N characters in sequence\n";
    print "5) GC content is calculated irrespective of Ns in input sequence, Masked GC% is\n";
    print "   calculated using length of sequence MINUS any N characters\n\n";
    exit(0);
}


my %mono = ('A' => 0,'C' => 0,'G' => 0,'T' => 0,'N' => 0);

# hashes for dinucleotide frequencies and percentages
my %di;
my %di_percent;

								
############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0], please specify the FASTA sequence file you wish atcg.pl\nto analyse on the command line, e.g. atcg.pl filename.\n\n";


my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {	

	my $seq = $entry->seq;

	# count mononucleotide frequencies
	$mono{"A"} += ($seq =~ tr/A/A/); 
	$mono{"C"} += ($seq =~ tr/C/C/); 
	$mono{"G"} += ($seq =~ tr/G/G/); 
	$mono{"T"} += ($seq =~ tr/T/T/); 
	$mono{"N"} += ($seq =~ tr/N/N/); 

	# do we also want to calculate dinucleotides? Skip if we don't
	next unless ($dinuc);
	
	# to count dinucleotides, loop through sequence, take 2 bp and increment the hash counter
	foreach my $i (0..length($seq)){
	    my $tmp = substr($seq,$i,2);		
		$di{$tmp}++;
	}
}
close(FILE);

# Calculate final statistics
# need to determine number of Ns to work out masked length and GC content of masked sequence
my $total_length = $mono{'A'} + $mono{'C'} + $mono{'G'} + $mono{'T'} + $mono{'N'};
my $masked_length     = $total_length - $mono{'N'};
my $total_GC_percent  = sprintf("%3.2f",(($mono{'G'}+$mono{'C'})/$total_length)*100);
my $masked_GC_percent = sprintf("%3.2f",(($mono{'G'}+$mono{'C'})/$masked_length)*100);
my $A_percent         = ($mono{'A'}/$total_length)*100;
my $T_percent         = ($mono{'T'}/$total_length)*100;
my $C_percent         = ($mono{'C'}/$total_length)*100;
my $G_percent         = ($mono{'G'}/$total_length)*100;
my $N_percent         = ($mono{'N'}/$total_length)*100;


print "\nLength\tGC%\tMasked length\tMasked GC%:\n";
print "$total_length\t$total_GC_percent\t$masked_length\t$masked_GC_percent\n\n";

print "A\tC\tG\tT\tN\n";
print "$mono{'A'}\t$mono{'C'}\t$mono{'G'}\t$mono{'T'}\t$mono{'N'}\n\n";

print "A%\tC%\tG%\tT%\tN%\n";
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n\n", $A_percent, $C_percent, $G_percent, $T_percent, $N_percent;

exit(0) unless ($dinuc);



######################################
# calculate dinucleotide percentages
######################################

				

my $dinuc_total = $di{'AA'}+$di{'AC'}+$di{'AG'}+$di{'AT'} + 
				  $di{'CA'}+$di{'CC'}+$di{'CG'}+$di{'CT'} +
				  $di{'GA'}+$di{'GC'}+$di{'GG'}+$di{'GT'} +
				  $di{'TA'}+$di{'TC'}+$di{'TG'}+$di{'TT'};

# will make four lines of output text
my ($out1,$out2,$out3,$out4);				
foreach my $i qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT){
	
	# set to 0 if no values exist for that dinucleotide
	($di{$i} = 0) if (!defined($di{$i}));
	
	$di_percent{$i} = (($di{$i}/$dinuc_total)*100);
	$out1 .= "$i\t";
	$out2 .= "$di{$i}\t";
	$out3 .= "$i%\t";
	$out4 .= sprintf("%3.2f\t",$di_percent{$i});
}

print "$out1\n$out2\n\n$out3\n$out4\n\n";
exit;


# only proceed if we are doing chi-squared analysis
exit(0) unless ($chi);

my $tot_obs = $dinuc_total;
#my $AA_exp = (($AA+$AT+$AC+$AG)*($AA+$TA+$CA+$GA))/$tot_obs;
#my $AT_exp = (($AA+$AT+$AC+$AG)*($TA+$TT+$TC+$TG))/$tot_obs;
#my $AC_exp = (($AA+$AT+$AC+$AG)*($CA+$CT+$CC+$CG))/$tot_obs;
#my $AG_exp = (($AA+$AT+$AC+$AG)*($GA+$GT+$GC+$GG))/$tot_obs;

#my $TA_exp = (($TA+$TT+$TC+$TG)*($AA+$AT+$AC+$AG))/$tot_obs;
#my $TT_exp = (($TA+$TT+$TC+$TG)*($TA+$TT+$TC+$TG))/$tot_obs;
#my $TC_exp = (($TA+$TT+$TC+$TG)*($CA+$CT+$CC+$CG))/$tot_obs;
#my $TG_exp = (($TA+$TT+$TC+$TG)*($GA+$GT+$GC+$GG))/$tot_obs;

#my $CA_exp = (($CA+$CT+$CC+$CG)*($AA+$AT+$AC+$AG))/$tot_obs;
#my $CT_exp = (($CA+$CT+$CC+$CG)*($TA+$TT+$TC+$TG))/$tot_obs;
#my $CC_exp = (($CA+$CT+$CC+$CG)*($CA+$CT+$CC+$CG))/$tot_obs;
#my $CG_exp = (($CA+$CT+$CC+$CG)*($GA+$GT+$GC+$GG))/$tot_obs;

#my $GA_exp = (($GA+$GT+$GC+$GG)*($AA+$AT+$AC+$AG))/$tot_obs;
#my $GT_exp = (($GA+$GT+$GC+$GG)*($TA+$TT+$TC+$TG))/$tot_obs;
#my $GC_exp = (($GA+$GT+$GC+$GG)*($CA+$CT+$CC+$CG))/$tot_obs;
#my $GG_exp = (($GA+$GT+$GC+$GG)*($GA+$GT+$GC+$GG))/$tot_obs;


#my $chi_squared = ((($AA-$AA_exp)**2)/$AA_exp)+((($AT-$AT_exp)**2)/$AT_exp)+
#                  ((($AC-$AC_exp)**2)/$AC_exp)+((($AG-$AG_exp)**2)/$AG_exp)+
#                  ((($TA-$TA_exp)**2)/$TA_exp)+((($TT-$TT_exp)**2)/$TT_exp)+
#                  ((($TC-$TC_exp)**2)/$TC_exp)+((($TG-$TG_exp)**2)/$TG_exp)+
#                  ((($CA-$CA_exp)**2)/$CA_exp)+((($CT-$CT_exp)**2)/$CT_exp)+
#                  ((($CC-$CC_exp)**2)/$CC_exp)+((($CG-$CG_exp)**2)/$CG_exp)+
#                  ((($GA-$GA_exp)**2)/$GA_exp)+((($GT-$GT_exp)**2)/$GT_exp)+
#                  ((($GC-$GC_exp)**2)/$GC_exp)+((($GG-$GG_exp)**2)/$GG_exp);

#printf  "Chi squared value = %6.2f\n", $chi_squared;
				
#print "Significance level at 5% = 16.92\n";
#print "Significance level at 1% = 21.67\n";


exit(0);
