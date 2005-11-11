#!/usr/local/bin/perl5.6.0 -w

use Ace;

$db = Ace->connect(-path  =>  '/wormsrv2/current_DB') || die "Can't connect\n";

my @papers = $db->fetch(-class => 'Paper');

foreach $paper (@papers){
  if($paper->Abstract){
    $longtext = $paper->Abstract;
    $longtext_details = ($db->fetch(Longtext => "$longtext"))->asAce;
    print "$longtext_details\n";
  }
}






