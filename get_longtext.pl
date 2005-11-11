#!/usr/local/bin/perl5.6.1 -w
 
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;
use Ace;

my $tace = &tace;

my $db = Ace->connect(-path=>"/nfs/disk100/wormpub/DATABASES/current_DB",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


my $text1= $db->fetch(-query=>'find LongText WBPaper00023205');
my $text2= $db->fetch(-query=>'find LongText WBPaper00017513');

print "Longtext 1 is:\n$text1";
print "Longtext 2 is:\n$text2";

$db->close;
