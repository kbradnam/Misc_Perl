#!/usr/local/bin/perl5.6.1 -w

use Ace;

$db = Ace->connect(-path  =>  '/wormsrv2/current_DB') || die "Can't connect\n";

my @seqs = $db->fetch(Sequence => "AH6.*");

foreach $seq (@seqs){
    $fasta_seq = $seq->asDNA();
    print "$fasta_seq\n";
}






