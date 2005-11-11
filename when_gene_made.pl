#!/usr/local/bin/perl5.8.0 -w

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use strict;
use Ace;
use Date::Calc qw(Day_of_Week);

my $tace = &tace;



# open connection to DB
my $db = Ace->connect(-path=>"/nfs/disk100/wormpub/camace_krb",
#my $db = Ace->connect(-path=>"/wormsrv2/stlace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


# grab set of problem history CDSs without Gene IDs
my @CDSs= $db->fetch(-query=>'find CDS WHERE Method = history AND NOT Remark = *auto-generated*');
#my @CDSs= $db->fetch(-query=>'find CDS WHERE Method = curated');
#my @CDSs= $db->fetch(-query=>"find CDS *wp*");

my $lab;

foreach my $history (@CDSs){

  # process name to get to stem and release info, and determine whether it is an isoform
  my $name = $history;
  $name =~ s/[a-z]*:wp(\d+)//;
  my $release = $1;
  my $isoform;
  if($history =~ m/\.[0-9]+([a-z])/){
    $isoform = "$1";
  }
  else{
    $isoform = "NO";
  }

  my $tmp = $history;
  $tmp =~ m/\.([0-9]+)[a-z]*/;
  my $suffix = $1;

  # what day was history gene made?
  my ($born,$bornyear) = &day_created($history);  
  print "History CDS - $history, made on $born\n";

  # now find if there are curated objects for history CDSs
  my @variants= $db->fetch(-query=>"find CDS ${name}*");

  my $counter = 0;

  for (my $i=0; $i<@variants; $i++){
#    print "\tMATCHES: $variants[$i]\n";
    next if ($history eq $variants[$i]);

    my ($day,$year) = &day_created($variants[$i]);   
    # ignore data if timestamp year is not from >2000
    if($year !~ m/^200/){
      next;
    }
    
    # what to do if you match a history gene
    if ($variants[$i] =~ m/wp/){
      next;
      my $release2 = $variants[$i];
      $release2 =~ s/.*:wp//;

      # must skip if match is from older release than source
      if ($release2 < $release){
#	print "\tWARNING: matches history object from older release ($release2)\n";
	next;
      }
      # must skip if match is auto generated
      my $remark = $variants[$i];
      if($remark =~ m/auto-generated/){
#	print "\tWARNING: auto-generated history object\n";
	next;
      }
    }

    # skip .a .b .c variants etc.
    if($variants[$i] =~ m/\.[a-z]/){
      #print "\tWARNING: non numbered suffix gene\n";
      next;
    }

    # check that history CDS isn't an isoform but matching CDS is
    my $tmp = $variants[$i];
    $tmp =~ s/:wp.*//;
    if(($isoform eq "NO") && ($tmp =~ m/[a-z]$/)){
#      print "\tWARNING: history CDS is not isoform, but match is\n";
      next;
    }

    # check that both are not isoforms with different letters
    if($variants[$i] =~ m/\.[0-9]+([a-z])/){
      my $variant_isoform = "$1";
      if($isoform ne "NO"){
	unless($isoform eq $variant_isoform){
#	  print "\tWARNING: isoform matches different isoform\n";
	  next;
	}
	
      }
    }
    


    # check that .6a or .6b doesn't match .6
    my $tmp = $variants[$i];
    $tmp =~ s/:wp.*//;
    if($isoform ne "NO" && $tmp =~ m/\.[0-9]+$/){
#      print "\tWARNING: isoform matches a non-isoform\n";
      next;

    }

    # check that .1 isn't matching .10, .11 etc.
    if($variants[$i] =~ m/\.([0-9]+)[a-z]*/){
      my $new_suffix = $1;
      if($suffix != $new_suffix){
#	print "\tWARNING: different numbered suffixes\n";
	next;
      }
    }



    $counter++;   
    print "\tRESULT: $counter) $history (born $born)-> $name* matches $variants[$i] made on $day\n";

  }
  
  # if no hits, try checking Pseudogene and transcript?
#  print "$day: CHANGED\n";
  print "\n";

}





sub day_created{

  my $obj = shift;


  # now get timestamp from curated object not history
  my $aql_query = "select s,s->Properties.node_session from s in object(\"CDS\",\"$obj\")";
  my @aql = $db->aql($aql_query);
  my $timestamp = $aql[0]->[1];

  # if original then use Structure tag
  if ($timestamp =~ m/original/){
    $aql_query = "select s,s->Structure.node_session from s in object(\"CDS\",\"$obj\")";
    @aql = $db->aql($aql_query);
    $timestamp = $aql[0]->[1];
    
    # if still original use DB_info tag
    if ($timestamp =~ m/original/){
      $aql_query = "select s,s->DB_info.node_session from s in object(\"CDS\",\"$obj\")";
      @aql = $db->aql($aql_query);
      $timestamp = $aql[0]->[1];
    }
  }
  
  my ($year, $month, $day, $time,$user) = $timestamp =~ m/(\d{4})\-(\d{2})\-(\d{2})_(\d{2}:\d{2}:\d{2}\.*\d*)_(.*)$/;
  
  my $weekday = Day_of_Week($year, $month, $day);

  return($weekday,$year);

}
$db->close;


