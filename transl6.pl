#!/usr/bin/perl -w

if ($ARGV[0] eq "" || $ARGV[1] eq ""){
    print "\n\t\tTRANSL6 version 1.0\n";
    print "TRANSL6 translate a nucleotide sequence in all six open\n";
    print "reading frames and outputs them to a single file.\n";
    print "Can be run from a single command line, e.g.\n\n";
    print "transl6 \"input_file\" \"output_file\"\n\n";

    print "Enter input file:  ";
    chomp($file1 = <STDIN>);
    print "Enter output file:  ";
    chomp($file2 =<STDIN>);
}
else{
    $file1 = $ARGV[0];
    $file2 = $ARGV[1];
}

`echo "\n" >> $file1`;
open (IN, "$file1") || die "Could not open file";
open (OUT, ">$file2") || die "couldn't create output file";
&read_codons;

while(<IN>){
    chomp;
    if(m/^>/){
	$name =$_;
    }
    else{                       
        $temp.= $_;
    }
}

@seq = split(//,$temp);
$tmp3 = uc($temp);
$tmp3 =~ s/A/Q/g;
$tmp3 =~ s/T/W/g;
$tmp3 =~ s/C/E/g;
$tmp3 =~ s/G/R/g;
$tmp3 =~ s/Q/T/g;
$tmp3 =~ s/W/A/g;
$tmp3 =~ s/E/G/g;
$tmp3 =~ s/R/C/g;
	
@rev_seq = reverse(split(//,$tmp3));


########################################

for ($a=0;$a<@seq;$a+=3){
    if(defined($seq[$a]) && defined($seq[$a+1]) && defined($seq[$a+2])){
	$codon = uc($seq[$a].$seq[$a+1].$seq[$a+2]);
	$pep1 .= "$codons{$codon}";
    }
    if(defined($rev_seq[$a]) && defined($rev_seq[$a+1]) && defined($rev_seq[$a+2])){
	$codon = uc($rev_seq[$a].$rev_seq[$a+1].$rev_seq[$a+2]);
	$rev_pep1 .= "$codons{$codon}";
    }

}
print OUT "$name 1 +\n$pep1\n";
print OUT "$name 2 -\n$rev_pep1\n";


for ($a=1;$a<@seq;$a+=3){
    if(defined($seq[$a]) && defined($seq[$a+1]) && defined($seq[$a+2])){
	$codon = uc($seq[$a].$seq[$a+1].$seq[$a+2]);
	$pep2 .= "$codons{$codon}";
    }
    if(defined($rev_seq[$a]) && defined($rev_seq[$a+1]) && defined($rev_seq[$a+2])){
	$codon = uc($rev_seq[$a].$rev_seq[$a+1].$rev_seq[$a+2]);
	$rev_pep2 .= "$codons{$codon}";
    }
}

print OUT "$name 3 +\n$pep2\n";
print OUT "$name 4 -\n$rev_pep2\n";


for ($a=2;$a<@seq;$a+=3){
    if(defined($seq[$a]) && defined($seq[$a+1]) && defined($seq[$a+2])){	
	$codon = uc($seq[$a].$seq[$a+1].$seq[$a+2]);
	$pep3 .= "$codons{$codon}";
    }
    if(defined($rev_seq[$a]) && defined($rev_seq[$a+1]) && defined($rev_seq[$a+2])){
	$codon = uc($rev_seq[$a].$rev_seq[$a+1].$rev_seq[$a+2]);
	$rev_pep3 .= "$codons{$codon}";
    }
}

print OUT "$name 5 +\n$pep3\n";
print OUT "$name 6 -\n$rev_pep3\n";

close(IN);
close(OUT);
       
####################################################
sub read_codons{

    %codons=("TTT","F","TTC","F","TTA","L","TTG","L","CTT","L","CTC","L","CTA","L","CTG","L",
	     "ATT","I","ATC","I","ATA","I","ATG","M","GTT","V","GTC","V","GTA","V","GTG","V",
	     "TCT","S","TCC","S","TCA","S","TCG","S","CCT","P","CCC","P","CCA","P","CCG","P",
	     "ACT","T","ACC","T","ACA","T","ACG","T","GCT","A","GCC","A","GCA","A","GCG","A",
	     "TAT","Y","TAC","Y","TAA","X","TAG","X","CAT","H","CAC","H","CAA","Q","CAG","Q",
	     "AAT","N","AAC","N","AAA","K","AAG","K","GAT","D","GAC","D","GAA","E","GAG","E",
	     "TGT","C","TGC","C","TGA","X","TGG","W","CGT","R","CGC","R","CGA","R","CGG","R",
	     "AGT","S","AGC","S","AGA","R","AGG","R","GGT","G","GGC","G","GGA","G","GGG","G");

}

