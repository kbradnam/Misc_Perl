#!/usr/local/bin/perl5.8.0 -w

my $dbpath = "/nfs/disk100/wormpub/TEST_BUILD/autoace";
my $command = "mv $dbpath/database/bog.wrm $dbpath/database/bog.old";

my $status = system($command);
if($status != 0){
  print "ERROR: $command failed\n";
}

print "Status is $status\n";
