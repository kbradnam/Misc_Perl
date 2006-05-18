#!/usr/bin/perl -w
#
# atcg.pl 
#
# A script to count GC content plus mono- and dinucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use Getopt::Long;

my $help;    # print help
my $dinuc;   # calculate dinucleotide frequencies as well as mononucleotide frequencies
my $chi;     # perform chi-squared analysis on dinucleotide frequencies
my $verbose; # extra output information

GetOptions ("dinuc"   => \$dinuc,
	    "chi"     => \$chi,
	    "verbose" => \$verbose,
	    "help"    => \$help);

# sanity checks
die "Can only use -chi option if -dinuc option is chosen\n" if ($chi && !$dinuc);


print "atcg.pl - a program to calculate nucleotide and dinucleotide composition\n" if ($verbose);

if($help){
    print "1) atcg.pl can analyse multiple sequences, but will concatenate the results together.\n";
    print "2) use -dinuc option to also calculate dinucleotide frequencies\n";
    print "3) use -chi option to perform chi-squared analysis of dinucleotide frequencies\n";
    print "4) use -verbose option for some extra (minor) output\n";
    print "5) 'Other' counter is a check for any non ATCG, or N characters in sequence\n";
    print "6) GC content is calculated irrespective of Ns in input sequence, Masked GC% is\n";
    print "   calculated using length of sequence MINUS any N characters\n\n";
    exit(0);
}


############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(IN,"$ARGV[0]") || die "Can't open $ARGV[0], please specify the FASTA sequence file you wish atcg.pl\nto analyse on the command line, e.g. atcg.pl filename.\n\n";

my $seq = "";

while(my $temp=<IN>){
    if($temp !~ />/){
	chomp($temp);
	$temp =~  tr/a-z/A-Z/;
	$seq .= $temp;
    }
}
close(IN);
    
  

# Calculate A,T,C,G, N
my ($A,$T,$C,$G,$N);
$A=$T=$C=$G=$N=0;

$A = $seq =~ tr/A/A/; 
$T = $seq =~ tr/T/T/;
$C = $seq =~ tr/C/C/;
$G = $seq =~ tr/G/G/;
$N = $seq =~ tr/N/N/;



# Calculate final statistics
# need to determine number of Ns to work out masked length and GC content of masked sequence
my $total_length      = length($seq);
my $masked_length     = $total_length - $N;
my $total_GC_percent  = sprintf("%3.2f",(($G+$C)/$total_length)*100);
my $masked_GC_percent = sprintf("%3.2f",(($G+$C)/$masked_length)*100);
my $other             = $total_length - $A - $T - $C - $G - $N;
my $A_percent         = ($A/$total_length)*100;
my $T_percent         = ($T/$total_length)*100;
my $C_percent         = ($C/$total_length)*100;
my $G_percent         = ($G/$total_length)*100;
my $N_percent         = ($N/$total_length)*100;



print "\nOverall statistics:\n-------------------" if ($verbose);
print "\n";

print "Length\tGC%\tMasked length\tMasked GC%:\n";
print "$total_length\t$total_GC_percent\t$masked_length\t$masked_GC_percent\n";



print "\nNucleotide counts:\n------------------" if ($verbose);
print "\n";

print "A\tC\tG\tT\tN\tOther\n";
print "$A\t$C\t$G\t$T\t$N\t$other\n\n";

print "A%\tC%\tG%\tT%\tN%\n";
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n\n", $A_percent, $C_percent, $G_percent, $T_percent, $N_percent;



exit(0) unless ($dinuc);



######################################
# calculate dinucleotide frequencies
######################################

print "Dinucleotide counts:\n--------------------\n" if ($verbose);

my ($AA,$AT,$AC,$AG,$TA,$TT,$TC,$TG,$CA,$CT,$CC,$CG,$GA,$GT,$GC,$GG);
$AA=$AT=$AC=$AG=$TA=$TT=$TC=$TG=$CA=$CT=$CC=$CG=$GA=$GT=$GC=$GG=0;

my ($AA_percent,$AC_percent,$AG_percent,$AT_percent,
    $CA_percent,$CC_percent,$CG_percent,$CT_percent,
    $GA_percent,$GC_percent,$GG_percent,$GT_percent,
    $TA_percent,$TC_percent,$TG_percent,$TT_percent);

$AA_percent=$AC_percent=$AG_percent=$AT_percent = 0;
$CA_percent=$CC_percent=$CG_percent=$CT_percent = 0;
$GA_percent=$GC_percent=$GG_percent=$GT_percent = 0;
$TA_percent=$TC_percent=$TG_percent=$TT_percent = 0;

    
foreach my $i(0..length($seq)){
    my $tmp = substr($seq,$i,2);	
  SWITCH: {		       
      if ($tmp =~ /AA/){$AA++; last SWITCH;}
      if ($tmp =~ /AT/){$AT++; last SWITCH;}
      if ($tmp =~ /AC/){$AC++; last SWITCH;}
      if ($tmp =~ /AG/){$AG++; last SWITCH;}
      if ($tmp =~ /TA/){$TA++; last SWITCH;}
      if ($tmp =~ /TT/){$TT++; last SWITCH;}
      if ($tmp =~ /TC/){$TC++; last SWITCH;}
      if ($tmp =~ /TG/){$TG++; last SWITCH;}
      if ($tmp =~ /CA/){$CA++; last SWITCH;}
      if ($tmp =~ /CT/){$CT++; last SWITCH;}
      if ($tmp =~ /CC/){$CC++; last SWITCH;}
      if ($tmp =~ /CG/){$CG++; last SWITCH;}
      if ($tmp =~ /GA/){$GA++; last SWITCH;}
      if ($tmp =~ /GT/){$GT++; last SWITCH;}
      if ($tmp =~ /GC/){$GC++; last SWITCH;}
      if ($tmp =~ /GG/){$GG++; last SWITCH;}
  }
}

my $dinuc_total = $AA+$AC+$AG+$AT+$CA+$CC+$CG+$CT+$GA+$GC+$GG+$GT+$TA+$TC+$TG+$TT;
$AA_percent = ($AA/$dinuc_total)*100;
$AC_percent = ($AC/$dinuc_total)*100;
$AG_percent = ($AG/$dinuc_total)*100;
$AT_percent = ($AT/$dinuc_total)*100;
$CA_percent = ($CA/$dinuc_total)*100;
$CC_percent = ($CC/$dinuc_total)*100;
$CG_percent = ($CG/$dinuc_total)*100;
$CT_percent = ($CT/$dinuc_total)*100;
$GA_percent = ($GA/$dinuc_total)*100;
$GC_percent = ($GC/$dinuc_total)*100;
$GG_percent = ($GG/$dinuc_total)*100;
$GT_percent = ($GT/$dinuc_total)*100;
$TA_percent = ($TA/$dinuc_total)*100;
$TC_percent = ($TC/$dinuc_total)*100;
$TG_percent = ($TG/$dinuc_total)*100;
$TT_percent = ($TT/$dinuc_total)*100;

# print dinuc output (counts and percentages)
print "AA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n";
print "$AA\t$AC\t$AG\t$AT\t$CA\t$CC\t$CG\t$CT\t$GA\t$GC\t$GG\t$GT\t$TA\t$TC\t$TG\t$TT\n\n";

print "AA%\tAC%\tAG%\tAT%\tCA%\tCC%\tCG%\tCT%\tGA%\tGC%\tGG%\tGT%\tTA%\tTC%\tTG%\tTT%\n";
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t", $AA_percent, $AC_percent, $AG_percent, $AT_percent;
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t", $CA_percent, $CC_percent, $CG_percent, $CT_percent;
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t", $GA_percent, $GC_percent, $GG_percent, $GT_percent;
printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\n\n", $TA_percent, $TC_percent, $TG_percent, $TT_percent;


# only proceed if we are doing chi-squared analysis
exit(0) unless ($chi);

my $tot_obs= $AA+$AT+$AC+$AG+$TA+$TT+$TC+$TG+$CA+$CT+$CC+$CG+$GA+$GT+$GC+$GG;
my $AA_exp = (($AA+$AT+$AC+$AG)*($AA+$TA+$CA+$GA))/$tot_obs;
my $AT_exp = (($AA+$AT+$AC+$AG)*($TA+$TT+$TC+$TG))/$tot_obs;
my $AC_exp = (($AA+$AT+$AC+$AG)*($CA+$CT+$CC+$CG))/$tot_obs;
my $AG_exp = (($AA+$AT+$AC+$AG)*($GA+$GT+$GC+$GG))/$tot_obs;

my $TA_exp = (($TA+$TT+$TC+$TG)*($AA+$AT+$AC+$AG))/$tot_obs;
my $TT_exp = (($TA+$TT+$TC+$TG)*($TA+$TT+$TC+$TG))/$tot_obs;
my $TC_exp = (($TA+$TT+$TC+$TG)*($CA+$CT+$CC+$CG))/$tot_obs;
my $TG_exp = (($TA+$TT+$TC+$TG)*($GA+$GT+$GC+$GG))/$tot_obs;

my $CA_exp = (($CA+$CT+$CC+$CG)*($AA+$AT+$AC+$AG))/$tot_obs;
my $CT_exp = (($CA+$CT+$CC+$CG)*($TA+$TT+$TC+$TG))/$tot_obs;
my $CC_exp = (($CA+$CT+$CC+$CG)*($CA+$CT+$CC+$CG))/$tot_obs;
my $CG_exp = (($CA+$CT+$CC+$CG)*($GA+$GT+$GC+$GG))/$tot_obs;

my $GA_exp = (($GA+$GT+$GC+$GG)*($AA+$AT+$AC+$AG))/$tot_obs;
my $GT_exp = (($GA+$GT+$GC+$GG)*($TA+$TT+$TC+$TG))/$tot_obs;
my $GC_exp = (($GA+$GT+$GC+$GG)*($CA+$CT+$CC+$CG))/$tot_obs;
my $GG_exp = (($GA+$GT+$GC+$GG)*($GA+$GT+$GC+$GG))/$tot_obs;


my $chi_squared = ((($AA-$AA_exp)**2)/$AA_exp)+((($AT-$AT_exp)**2)/$AT_exp)+
                  ((($AC-$AC_exp)**2)/$AC_exp)+((($AG-$AG_exp)**2)/$AG_exp)+
                  ((($TA-$TA_exp)**2)/$TA_exp)+((($TT-$TT_exp)**2)/$TT_exp)+
                  ((($TC-$TC_exp)**2)/$TC_exp)+((($TG-$TG_exp)**2)/$TG_exp)+
                  ((($CA-$CA_exp)**2)/$CA_exp)+((($CT-$CT_exp)**2)/$CT_exp)+
                  ((($CC-$CC_exp)**2)/$CC_exp)+((($CG-$CG_exp)**2)/$CG_exp)+
                  ((($GA-$GA_exp)**2)/$GA_exp)+((($GT-$GT_exp)**2)/$GT_exp)+
                  ((($GC-$GC_exp)**2)/$GC_exp)+((($GG-$GG_exp)**2)/$GG_exp);

printf  "Chi squared value = %6.2f\n", $chi_squared;
				
print "Significance level at 5% = 16.92\n";
print "Significance level at 1% = 21.67\n";


exit(0);
