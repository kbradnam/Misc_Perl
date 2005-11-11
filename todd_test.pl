#!/usr/bin/perl -w

# This program will pull out N base pairs upstream of each
# gene in C. elegans

use lib "/Korflab/lib/perl";
use Ace;
use strict;
use Ace::Sequence;

my $upstream = 10;
my $length   = 10;

my $db = Ace->connect('/Korflab/Data_sources/WormBase/WS140') || die "Connection
failure: ",Ace->error;

my @genes = $db->fetch(Coding_transcripts => '2*');
for my $cd (@genes) {

 if (my $seq = Ace::Sequence->new($cd)) {
   print "$cd $seq \n";
   print $seq->features(),"\n";
    my @first_gene = sort {$a->start <=> $b->start} $seq->features('gene');

    if (my $g = $first_gene[0]) {
        my $up = Ace::Sequence->new(-seq=>$g,
                                    -offset=>(-$upstream),
                                   -length => ($length));
        my $current_dna = $up->dna;
        chomp $current_dna;
        if ($up !~ /^$/) {
        print "\>",$seq->name,"\n",$up->dna,"\n";
        }
      }
    }

}

