#!/usr/local/bin/perl5.6.1 -w
#
# db_backup_and_compare.pl
#
# backup database and compare to last backed up database to look for lost data
#
# Last updated by: $Author$     
# Last updated on: $Date$      

use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use Carp;

######################################
# variables and command-line options # 
######################################

our (@backups);

our $backup_dir = "/nfs/disk100/wormpub/DATABASES/BACKUPS";
our $date = `date +%y%m%d`; chomp $date;


# Find what backups have been made for the database


  my $db = "geneace";
#  my $backup_dbs = "/${backup_dir}/${db}_backup.*";
#  open (BACKUP_DB, "/bin/ls -d -1 -t $backup_dbs |")  || croak "cannot open $backup_dbs\n";
#  while (<BACKUP_DB>) { 
#    chomp;
#    (/$db\_backup\.(\d+)$/); 
#    # All database dates get added to @backups array, first element is last backup
#    push(@backups,$1);
#  }
#  close(BACKUP_DB);


  # keep TransferDB logs in backup directory
  chdir("$backup_dir") || print "Couldn't cd to $backup_dir\n";

  my $return = system("TransferDB.pl -start /wormsrv1/$db -end ${backup_dir}/${db}_backup\.test -database -wspec -name ${db}"); 
  if($return != 0){
    print "ERROR: Couldn't run TransferDB.pl correctly. Exit = $return  Check TransferDB log\n";
    exit;
  }

exit(0);
