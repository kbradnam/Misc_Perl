#!/usr/local/bin/perl5.6.1 -w

use Ace;

use strict;

$|=1;

# open a local database connection

my $path = "/wormsrv1/geneace/";
my $db = Ace->connect(-path  => $path);

my $query = "Find Laboratory Representative != \"WBPerson*\"";
my @labs = $db->fetch(-query => $query);


foreach my $lab (@labs){
#  print "Laboratory : $lab\n";
  my ($representative) = $lab->at('CGC.Representative');
#  print "Representative \"$representative\"\n";
  my @people = $lab->at('Person');
#  print "Person \"$person\"\n";
#  print "\n";

  $representative =~ s/ [A-Z]+$//;
#  print "Representative \"$representative\"\n";

  my @names = $db->fetch(-class => "Person_name", -name => "$representative");
  foreach my $person_name (@names){
#    print "Person_name : \"$person_name\"\n";
    my @lastnames = $person_name->at('Name.Last_name_of');
    foreach my $lastname (@lastnames){
#      print "Last_name_of $lastname\n";
      foreach my $bloke (@people){
	if ($lastname eq $bloke){
	  print "Laboratory : $lab\n";
	  print "Representative \"$bloke\"\n\n";
	}
      }
    }
  }
  print "\n\n";
}


$db->close;
