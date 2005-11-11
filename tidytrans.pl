#!/usr/local/bin/perl 

$answer=1;
unless ($ARGV[0]){
    print "\n\ttidytrans v1.5 by Keith Bradnam and Roy Chaudhuri\n\n";
    print "tidytrans translates a file containing nucleotide sequences into\n";
    print "protein sequence.  It also tidies them at the same time.\n";
    print "tidytrans can also cope with unspecified nucleotide sequence (N)\n";
    print "ONLY if N's occur at the 3rd position of quartet codons, e.g.\n";
    print "CCN will be translated into Proline (P)\n\n";
    print "Version 1.5 is adapted for use with gaps3. If a codon is not\n";
    print "recognised then an X is inserted into the protein sequence, to\n";
    print "maintain sequence length. This does not affect the stop codon \n";
    print "count.\nStop codons are now represented by asterisks (*) rather"; 
    print "than an X\n";
    print "In addition, tidytrans will now translate codons containing\n";
    print "Y and R nucleotides, as long as the amino acid is not ambiguous.\n";
    print "For example, YUG will be translated into Leucine (L).\n\n";
    print "tidytrans can be also be run in command line mode, e.g.\n\n";
    print "tidytrans \"input filename\" \"output filename\"\n\n";


    print"please enter the input file name: ";
    chop ($ARGV[0]=<STDIN>);
}
until (open (INFILE,$ARGV[0])){
    print "sorry, $ARGV[0] was not found, ";
    print"please re-enter the file name: ";
    chop ($ARGV[0]=<STDIN>);
}

unless ($ARGV[1]){
    print"please enter the output file name: ";
    chop ($ARGV[1]=<STDIN>);

    print "\nPlease select the genetic code you wish to use when translating:\n\n";
    print "1) Universal Genetic Code \t\t\t\tTGA=* TAA=* TAG=*\n";
    print "2) Vertebrate Mitochondrial code \t\t\tAGR=* ATA=M TGA=W\n";
    print "3) Yeast Mitochondrial code \t\t\t\tCTN=* ATA=M TGA=W\n";
    print "4) Filamentous fungi Mitochondrial code \t\tTGA=W\n";
    print "5) Insects and Plathyhelminthes Mitochondrial code \tATA=M TGA=W AGR=S\n";
    print "6) Nuclear code of Cilitia \t\t\t\tUAA=Q UAG=Q\n";
    print "7) Nuclear code of Euplotes \t\t\t\tUGA=C\n";
    print "8) Mitochondrial code of Echinoderms \t\t\tUGA=W AGR=S AAA=N\n\n";
    
    $answer=0;
    until($answer>0 && $answer<9){
	print "Enter choice (1-8): ";
	chop($answer = <STDIN>);
	print "your choice was $answer\n";
	print "Please enter a value beween 1 and 8\n" if($answer <1 || $answer >8);
    }
}

open (OUT,">$ARGV[1]")|| die "$!";
&read_code1 if ($answer == 1);
&read_code2 if ($answer == 2);
&read_code3 if ($answer == 3);
&read_code4 if ($answer == 4);
&read_code5 if ($answer == 5);
&read_code6 if ($answer == 6);
&read_code7 if ($answer == 7);
&read_code8 if ($answer == 8);



##################################################
SEQ:while(<INFILE>){
    if(/^>/){
	if($flag){
	    print OUT &tidy_seq($seq),"\n";
	}
	print OUT;
	$seq = "";
	$flag = 1;
	next SEQ;
    }
    $seq .= $_;
}
if($flag){
    print OUT &tidy_seq($seq),"\n";
}

print "\nSequences have been tidied and translated.\n";
if($error1>0){
    print "\n$error1 sequence\(s\) contained incomplete codons.\n";
}
if($error2>0){
    print "$error2 sequence\(s\) contained multiple stop codons.\n";
}
print "\n";				
###########################################################
sub tidy_seq{
#adds a new line character every 60 bases  
    my ($x_counter) = 0;
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;
    my ($tmp) = "";
    my ($output_seq) = "";
    my ($translated_seq) = "";
    my (@seq2) = "";
    my ($start,$end);
    my ($untranslated_length) = length($seq);


    $error1++ if(($untranslated_length % 3) >0);

    @seq2 = split(//,$seq);
    for ($a=0;$a<@seq2;$a++){
	$tmp = $seq2[$a].$seq2[$a+1].$seq2[$a+2];
	$translated_seq .= $codons{$tmp};
	if ($codons{$tmp} eq "") {$translated_seq .= "X"}
	$x_counter++ if($codons{$tmp} eq "X");
	last if (($a+3>$untranslated_length));

	$a+=2;
    }
    my ($translated_length) = length($translated_seq);
    my ($to_add) = int($translated_length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($translated_seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $error2++ if ($x_counter>1);
    $output_seq .= substr($translated_seq,$start);
    return ($output_seq);
}
########################################
sub read_code1{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","I","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","*","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","R","AGG","R",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G",
             "TTY","F","TTR","L","CTY","L","CTR","L",
             "ATY","I","GTY","V","GYR","V","TCY","S",
             "TCR","S","CCY","P","CCR","P","ACY","T",
             "ACR","T","GCY","A","GCR","A","TAY","Y",
             "TAR","*","CAY","H","CAR","Q","AAY","N",
             "AAR","K","GAY","D","GAR","E","TGY","C",
             "CGY","R","CGR","R","AGY","S","AGR","R",
             "GGY","G","GGR","G","YTR","L","YTA","L",
             "YTG","L");

}
######################################
sub read_code2{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","M","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","W","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","*","AGG","*",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code3{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","*","CTC","*","CTA","*","CTG","*",
             "ATT","I","ATC","I","ATA","M","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","W","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","R","AGG","R",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","*","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code4{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","I","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","W","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","R","AGG","R",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code5{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","M","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","W","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","S","AGG","S",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code6{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","I","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","Q","TAG","Q",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","*","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","R","AGG","R",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code7{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","I","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","K","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","C","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","R","AGG","R",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
sub read_code8{
    %codons=("TTT","F","TTC","F","TTA","L","TTG","L",
	     "CTT","L","CTC","L","CTA","L","CTG","L",
             "ATT","I","ATC","I","ATA","I","ATG","M",
	     "GTT","V","GTC","V","GTA","V","GTG","V",
             "TCT","S","TCC","S","TCA","S","TCG","S",
	     "CCT","P","CCC","P","CCA","P","CCG","P",
             "ACT","T","ACC","T","ACA","T","ACG","T",
	     "GCT","A","GCC","A","GCA","A","GCG","A",
             "TAT","Y","TAC","Y","TAA","*","TAG","*",
	     "CAT","H","CAC","H","CAA","Q","CAG","Q",
             "AAT","N","AAC","N","AAA","N","AAG","K",
	     "GAT","D","GAC","D","GAA","E","GAG","E",
             "TGT","C","TGC","C","TGA","W","TGG","W",
	     "CGT","R","CGC","R","CGA","R","CGG","R",
             "AGT","S","AGC","S","AGA","S","AGG","S",
	     "GGT","G","GGC","G","GGA","G","GGG","G",
	     "TCN","S","CTN","L","CCN","P","CGN","R",
	     "ACN","T","GTN","V","GCN","A","GGN","G");
}
######################################
