#!/usr/local/bin/perl5.6.1 -w


$|=1;
use IO::Handle;
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;



my $dbpath = "/wormsrv1/geneace";

&test;

sub test{

  my $answer = check_write_access($dbpath);

  print "Answer is $answer\n";
  
  exit;
  print "now I'm here\n";

}
sub check_write_access{

  my $database = shift;
  my $answer = "yes";

  $answer = "no" if (-e "${database}/database/lock.wrm");
  return($answer);

}
    

