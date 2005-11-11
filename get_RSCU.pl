#!/usr/bin/perl -w
#
# get_RSCU.pl
#
# simple program to count codons and calculate RSCU for a FASTA file
#
# by Keith Bradnam
#
# 3/28/2005
#
#############################################################################################

use strict;

my %count_codons;

open (INFILE, "$ARGV[0]") || die "could not open input file";

my $seq = "";;

while(my $temp = <INFILE>){
    # stop at each header line to process previous gene
    if ($temp =~ m/^>/ && defined($seq)){	

	# process that gene's sequence and add count of codons to hash     
	while($seq =~ s/(^[AUCG]{3})//){
	    $count_codons{$1}++;
	}
    }
    else{
	chomp($temp);
	$temp =~ tr/a-z/A-Z/;
	die "non alphanumeric character in $temp\n" if ($temp =~ m/\W/);
	die "not ATCG character in $temp\n" if ($temp =~ m/[BDEFHIJKLMNOPQRSUVWXYZ]/);

	# now convert Thymine (T) to Uracil (U) to save time later
	$temp =~ tr/T/U/;

	# now append to $seq
	$seq .= $temp;
    }    
}    
close (INFILE);

# process last sequence in file
while($seq =~ s/(^[AUCG]{3})//){
    $count_codons{$1}++;
}


# Now output data
my @types = ("Phe UUU-UUC", "Leu UUA-UUG-CUU-CUC-CUA-CUG","Iso AUU-AUC-AUA", "Met AUG",
	     "Val GUU-GUC-GUA-GUG","Ser UCU-UCC-UCA-UCG-AGU-AGC", "Pro CCU-CCC-CCA-CCG", "Thr ACU-ACC-ACA-ACG",
	     "Ala GCU-GCC-GCA-GCG", "Tyr UAU-UAC", "Ter UAA-UAG-UGA", "His CAU-CAC", "Gln CAA-CAG",
	     "Asn AAU-AAC", "Lys AAA-AAG", "Asp GAU-GAC", "Glu GAA-GAG", "Cys UGU-UGC",
	     "Trp UGG", "Arg CGU-CGC-CGA-CGG-AGA-AGG", "Gly GGU-GGC-GGA-GGG");


foreach my $type (@types){
    my ($name,$codon_list) = split(/ /,$type);

    printf("%3s %23s",$name,$codon_list);

    my @codons = split(/\-/,$codon_list);

    # need a codon counter for each amino acid to be able to count RSCU    
    my $total = 0;
    
    # now just loop through codons to print codon counts and calculate total number of codons
    foreach my $codon (@codons){
	$total += $count_codons{$codon};
	printf(" %7d",$count_codons{$codon});
    }

    # now loop through again calculating RSCU values
    # RSCU = ration of observed codons compared with what would be expected if distributed equally
    foreach my $codon (@codons){       
	my $RSCU = $count_codons{$codon}/($total/scalar(@codons));
	printf(" %.4f",$RSCU);
    }
    print "\n";       
}
