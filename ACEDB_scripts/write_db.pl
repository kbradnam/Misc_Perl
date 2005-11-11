#!/usr/local/bin/perl5.6.0 -w

my $tace="/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";
my $db = glob("~krb/testace");

my $command=<<END;
pparse aaa.ace
save
quit
END

&DbWrite($command,"$tace -tsuser \"big boy blue\"",$db);


sub DbWrite {
  my ($command,$exec,$dir)=@_;
  open (WRITEDB,"| $exec $dir");
  print WRITEDB $command;
  close WRITEDB;
}
