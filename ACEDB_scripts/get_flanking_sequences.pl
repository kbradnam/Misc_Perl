#!/usr/local/bin/perl5.6.1 -w
#
# get_flanking_sequence
#
# dl1

use Ace;
use strict;

my $sequence     = shift;
my $x_coordinate = shift;
my $y_coordinate = shift;
my $flanking_seq_length    = 30;

# if you only provide one coordinate, use the x again (i.e. 1 bp feature)
if (!defined $y_coordinate) {$y_coordinate = $x_coordinate;}

print "Looking for flanking sequences to $x_coordinate - $y_coordinate in $sequence\n";

my $dbpath = "/wormsrv2/current_DB";                        # path of ACEDB database
my $tace   = "/nfs/disk100/wormpub/ACEDB/bin_ALPHA/tace";   # path of tace executable

my %right;
my %rightOffset;
my %left;
my $seq;
my $source;
my $source_start;
my $source_end;
my $clone = 0;
my $cds = 0;

my $db = Ace->connect(-path=>$dbpath,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};


my $obj = $db->fetch(Sequence=>$sequence);
if (!defined ($obj)) {
  print LOG "Could not fetch sequence $sequence\n";
  next;
}

# Grab sequence
$seq =$obj->asDNA();
$seq =~ s/\>[\w\.]+//mg;
$seq =~ tr/a-z/A-Z/;
$seq =~ s/\W+//mg;
 
# munge coordinates to mark flanking sequence
my $left_boundary = $x_coordinate - $flanking_seq_length - 1;
my $flank_left  = $y_coordinate;
my $sequence_length   = length ($seq);


if(defined($obj->at('Properties.Coding.CDS'))){
  $cds = 1;
  # grab coordinates of gene in parent
  $source = $obj->Source;
  my $new_obj = $obj;
  $new_obj =~ s/\./\\\./; # Need to escape the dot else AcePerl has kittens
  ($source_start) = $source->at("Structure.Subsequence.$new_obj");
  ($source_end) = $source->at("Structure.Subsequence.$new_obj")->right(2);
}


if(defined($obj->at('Properties.Genomic_canonical'))){
  $clone =1; 
  $right{$obj}            = $obj->Overlap_right ;
  $rightOffset{$obj}      = $obj->Overlap_right(2) ;
  $left{$obj}             = $obj->Overlap_left ;
}



# checks
if (($x_coordinate - $flanking_seq_length - 1) < 1) {
    print "You have specified a origin which lies to the left of the sequence [You want $left{$obj}]\n";

    if($clone){
      my $obj_left = $db->fetch(Sequence=>$left{$obj});
      if (!defined ($obj_left)) {
	print LOG "Could not fetch sequence $sequence\n";
	next;
      }
      my $left_seq = $obj_left->asDNA();
      $left_seq =~ s/\>\w+//mg;
      $left_seq =~ tr/a-z/A-Z/;
      $left_seq =~ s/\W+//mg;
      
      $rightOffset{$obj_left}  = $obj_left->Overlap_right(2);
      
      my $overlap = ($rightOffset{$obj_left} - 51 );      
      my $newseq = substr($left_seq,$overlap,50) . $seq;
      $seq = $newseq;
      
      $left_boundary = $left_boundary + 50;
      $flank_left  = $flank_left  + 50;
      
      $obj_left->DESTROY();
    }


    if($cds){
      my $left_seq = $source->asDNA();
      $left_seq =~ s/\>\w+//mg;
      $left_seq =~ tr/a-z/A-Z/;
      $left_seq =~ s/\W+//mg;
      my $diff = $flanking_seq_length - $x_coordinate;
      my $newseq = substr($left_seq,$source_start,$diff);
      print "diff is $diff\n";
      print "$source_start\n";
      print "Newseq\n$newseq\n";
      $seq = $newseq;
    }

}


if ($y_coordinate > ($sequence_length - $flanking_seq_length)) {
    print "You have specified a origin which lies to the right of the clone [You want $right{$obj}]\n";

    my $obj_right = $db->fetch(Sequence=>$right{$obj});
    if (!defined ($obj_right)) {
	print LOG "Could not fetch sequence $sequence\n";
	next;
    }
    my $right_seq = $obj_right->asDNA();
    $right_seq =~ s/\>\w+//mg;
    $right_seq =~ tr/a-z/A-Z/;
    $right_seq =~ s/\W+//mg;
    
    my $overlap = ($sequence_length - $rightOffset{$obj} +1);
    $seq .= substr($right_seq,$overlap,50);

    $obj_right->DESTROY();

}

# extract flanking sequence
my $seq_right = substr($seq,$left_boundary,$flanking_seq_length);
my $seq_left  = substr($seq,$flank_left,$flanking_seq_length);

# report
print "$seq_right\t$seq_left\n";




$obj->DESTROY();
$db->close;
exit(0);
