#!/usr/local/bin/perl5.6.0 -w
# 
# current_DB_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the current_DB database
# to look for bogus sequence entries
#
# Last updated by: $Author$
# Last updated on: $Date$

use Ace;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use strict;

 
##############################
# command-line options       #
##############################

my $database = "/wormsrv2/current_DB/";



############################################
# Initialise variables
############################################

open(IN,"/wormsrv2/autoace/fixfiles/new_connections.ace") || die "Couldn't open the input file darling\n";



our $tace = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || do { print "Whoops. Connection failure: ",Ace->error; die();};

open(OUT1,">stlouis_EST_CDS_fix.ace") || die "What a darned shame, couldn't create the infernal output thingy\n";
open(OUT2,">sanger_EST_CDS_fix.ace") || die "What a darned shame, couldn't create the infernal output thingy\n";

print OUT1 "//Sequences from St. Louis\n";
print OUT2 "//Sequences from Sanger\n";

my $seq_counter;
my $hx_counter;
my $rw_counter;
my $other_counter;

while(<IN>){
  chomp;
  if (/Sequence/){
    $seq_counter++;
    my $seq = $_;
    $seq =~ s/Sequence : //;
    $seq =~ s/\"//g;
    last if ($seq eq "");
#    print "$seq\n";
    my $matching_cDNA = <IN>;
    chomp($matching_cDNA);
    $matching_cDNA =~ s/Matching_cDNA //;

    my $lab = &lookup_lab($seq);
    if ($lab eq "RW"){
#      print "Here! $lab $seq $matching_cDNA\n";
      print OUT1 "Sequence : $seq\nMatching_cDNA $matching_cDNA\n\n";
      $rw_counter++;
    }
    elsif($lab eq "HX"){
      print OUT2 "Sequence : $seq\nMatching_cDNA $matching_cDNA\n\n";
      $hx_counter++;
    }
    else {print "PANIC!!! $seq -> $matching_cDNA -> $lab\n";$other_counter++;}
  }
  
}
print "Found $seq_counter sequences\n";
print "Found $rw_counter sequences for St. Louis\n";
print "Found $hx_counter sequences for Sanger\n";
print "Found $other_counter other sequences\n";

sub lookup_lab{
  my $seq = shift;
  my $db_seq = $db->fetch(-class=>'Sequence',-name=>"$seq");
  my ($lab) = $db_seq->at('Origin.From_Laboratory');

  if(defined($lab)){
    return($lab);
  }
  else{
    $lab = "OTHER";
    return($lab);
  }
}




$db->close;
close(IN);
close(OUT1);
close(OUT2);
exit(0);
