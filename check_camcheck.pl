#!/usr/local/bin/perl


use Ace;




my $dbpath = "/wormsrv2/camace"; 

my $db = Ace->connect(-path=>$dbpath);
 
#$clone = "ZK637";
$clone = "ZK643";

$obj = $db->fetch(Sequence=>$clone);
if (!defined ($obj)) {
  print "Could not fetch sequence $clone\n";
}
else{
  print "Fetched $obj\n\n";
  $dna = $db->fetch(DNA=>$clone);
  print "DNA is:\n$DNA\n\n";
}

#$seq_ace=$obj->DNA;
#if (!$seq_ace) {
#  print"$clone NOT_IN_ACEDB $clone\n" ;
#  
#}
#else{
#  print "$seq_ace\n";
#}
