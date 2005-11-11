#!/usr/bin/perl -w

use lib "/Korflab/lib/perl";
use Ace;
use strict;

# open a local database connection
my $db = Ace->connect(-path  =>  '/Korflab/Data_sources/WormBase/WS140');


my $subseq = $db->fetch(CDS => "AH6.1");

print "Subseq is $subseq\n";

exit(0);

my $tag = $subseq->at('Properties.Coding.CDS');
  print "Tag is $tag\n";

if(defined($subseq->at('Properties.Coding.CDS'))){
  print "CDS tag is present in $subseq\n";
}
else{
  print "CDS tag NOT present\n";
}

$db->close;

exit(0);


if($subseq->Transposon){
  print "Tag present\n";
}


my $subseq2 = $db->fetch(Sequence => "AC3.7");

print "$subseq2\n";

my $tag2 = $subseq2->Properties('CDS');
  print "Tag2 is $tag2\n";

exit;

# grab object
my $obj = $db->fetch(Predicted_gene => 'ZK154.3');
my $sequence  = $db->fetch(Sequence => 'ZK637');
my $count     = $db->count(Sequence => 'D*');

my @remarks = $db->fetch(-query=>'find Predicted_gene; DB_remark');

print "Object is $obj\nSequence is $sequence\nCount is $count\n\n";

print "@remarks\n";


