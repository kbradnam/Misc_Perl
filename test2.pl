#!/usr/local/bin/perl5.6.0 -w

use Ace;

my $db = Ace->connect (-path => "/wormsrv2/autoace", -program => "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_5/tace");

print "Hello\n";   
#$db->close;

