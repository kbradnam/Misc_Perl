#!/usr/bin/perl -w

# Keith Bradnam, 1998

use strict;
use Getopt::Long;

my $dinuc;
my $help;
GetOptions ("dinuc"  => \$dinuc,
	    "help"   => \$help);

print "atcg - a program to calculate nucleotide and dinucleotide composition\n";

if($help){
    print "1) atcg can analyse multiple sequences, but will concatenate the results together.\n";
    print "2) use -dinuc option to also calculate dinucleotide frequencies (slower)\n";
    print "3) 'Other' counter is a check for any non ATCG, or N characters in sequence\n";
    print "4) GC content is calculated irrespective of Ns in input sequence, Masked GC% is\n";
    print "   calculated using length of sequence MINUS any N characters\n\n";
    exit(0);
}


# open sequence file and process
open(IN,"$ARGV[0]") || die "Couldn't open file, please specify the FASTA sequence file you wish atcg\nto analyse on the command line, e.g. atcg filename.\n\n";

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
my $GC_percent        = (($G+$C)/$total_length)*100;
my $masked_GC_percent = (($G+$C)/$masked_length)*100;
my $other             = $total_length - $A - $T - $C - $G - $N;
my $A_percent         = ($A/$total_length)*100;
my $T_percent         = ($T/$total_length)*100;
my $C_percent         = ($C/$total_length)*100;
my $G_percent         = ($G/$total_length)*100;
my $N_percent         = ($N/$total_length)*100;



print "Nucleotide count:\n";
print "A  = $A\tC  = $C\tG  = $G\tT  = $T\tN  = $N\tOther = $other\n";
printf "A%% = %3.3f\tC%% = %3.3f\tG%% = %3.3f\tT%% = %3.3f\tN%% = %3.3f\n\n", $A_percent, $C_percent, $G_percent, $T_percent, $N_percent;


print "Length (bp)\tChromosome GC %\t\tMasked length\tMasked GC %:\n";

print "$total_length\t$GC_percent\t$masked_length\t$masked_GC_percent\n\n";





exit(0) unless ($dinuc);

my ($AA,$AT,$AC,$AG,$TA,$TT,$TC,$TG,$CA,$CT,$CC,$CG,$GA,$GT,$GC,$GG);
$AA=$AT=$AC=$AG=$TA=$TT=$TC=$TG=$CA=$CT=$CC=$CG=$GA=$GT=$GC=$GG=0;


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
print "Dinucleotide count: (rows show 5' nucleotide, columns show 3' nucleotide\n";
print "\tA\tC\tG\tT\t\n";
print "\t--------------------------\n";
print "A\t$AA\t$AC\t$AG\t$AT\n";
print "\n";

print "C\t$CA\t$CC\t$CG\t$CT\n";
print "\n";

print "G\t$GA\t$GC\t$GG\t$GT\n";				
print "\n";

print "T\t$TA\t$TC\t$TG\t$TT\n";				#
print "\n";



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


my $chi=((($AA-$AA_exp)**2)/$AA_exp)+((($AT-$AT_exp)**2)/$AT_exp)+
    ((($AC-$AC_exp)**2)/$AC_exp)+((($AG-$AG_exp)**2)/$AG_exp)+
    ((($TA-$TA_exp)**2)/$TA_exp)+((($TT-$TT_exp)**2)/$TT_exp)+
    ((($TC-$TC_exp)**2)/$TC_exp)+((($TG-$TG_exp)**2)/$TG_exp)+
    ((($CA-$CA_exp)**2)/$CA_exp)+((($CT-$CT_exp)**2)/$CT_exp)+
    ((($CC-$CC_exp)**2)/$CC_exp)+((($CG-$CG_exp)**2)/$CG_exp)+
    ((($GA-$GA_exp)**2)/$GA_exp)+((($GT-$GT_exp)**2)/$GT_exp)+
    ((($GC-$GC_exp)**2)/$GC_exp)+((($GG-$GG_exp)**2)/$GG_exp);

printf  "Chi squared value = %6.2f\n", $chi;
				
print "Significance level at 5% = 16.92\n";
print "Significance level at 1% = 21.67\n";
