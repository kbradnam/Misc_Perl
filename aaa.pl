#!/usr/local/bin/perl5.6.1 -w


use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;
use Ace;


  # Find out Locus->Sequence connections
  my $tace = &tace;
  my $db = Ace->connect(-path  => "/wormsrv2/autoace",
                        -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  my $query = "Find Locus WHERE (Genomic_sequence OR Transcript) AND CGC_approved";
  my $locus_seq_count = $db->fetch(-query=> "$query");
  $db->close;

print "There are $locus_seq_count locus->sequence connections\n";


