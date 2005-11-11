#!/usr/local/bin/perl5.8.0 -w

use strict;

my $path = "/nfs/disk100/wormpub";
my @users = ("ar2","dl1","pad","krb");

# is today monday?
my $age = 1;
$age = 3 if (`date +%u` == 1);

foreach my $user (@users){
  my $dbpath = $path."/camace_".$user."/database";

  my (@date) =`ls -ltr $dbpath/block*.wrm`;
  
  # get the last block file in @date and split to find date
  my @line = split(/\s+/,$date[$#date]);

  # Was file modified in last day?
  if (-M $line[8] >$age){
   print "camace_${user} not modified in last day, no need to re-run camcheck\n";
  }
}
