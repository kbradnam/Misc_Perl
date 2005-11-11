#!/usr/local/bin/perl5.6.1 -w
#
# works on a acefile/dump of transposons from the database
#
# Kerstin Jekosch

use strict;
use lib '/wormsrv2/scripts/';
use Wormbase;

$|=1;

my (%feat,$cds,$timestamp);
my $tace    = &tace;
my $ts = "-O \".*-*_*_*\"";
my $lab;
my $log;

my $db = Ace->connect(-path=>"/wormsrv2/stlace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};


open(SANGER,">sanger_transcripts.ace") || die "Couldn't create Sanger output file\n";
open(STLOUIS,">stlouis_transcripts.ace") || die "Couldn't create St. Louis output file\n";



while (<>) {
 
 next if (/\/\//);
	
  # generate Transcript class object
  if (/Sequence\s+:\s+\"(\S+)\"\s+$ts/) {

    $cds = $1;

    print SANGER "$log" if ($lab eq "HX");
    print STLOUIS "$log" if ($lab eq "RW");

    $lab = ""; # reset at each new object


    my $obj = $db->fetch(Sequence=>$cds);

    my $source = $obj->Source;
    my $new_obj = $obj;
    $new_obj =~ s/\./\\\./; # Need to escape the dot else AcePerl has kittens
    my ($source_start) = $source->at("Structure.Subsequence.$new_obj");
    my ($source_end) = $source->at("Structure.Subsequence.$new_obj")->right(2);
    
    # Now remove old subsequence from parent
    $log  = "Sequence : \"$source\"\n";
    $log .= "-D Subsequence \"$cds\"\n\n";

    # Remove subsequence altogether
    $log .= "-D Sequence \"$cds\"\n\n";

    # Add link to Transcript
    $log .= "Sequence : \"$source\"\n";
    $log .= "S_Child Transcript_child $cds $source_start $source_end\n\n";

    # Add details of new transcript
    $log .= "Transcript : \"$cds\"\n";
  }
  elsif (/^Origin\s+$ts From_Laboratory\s+$ts \"(\w+)\"/) {
    $lab = $1;
    s/From_Laboratory/From_laboratory/;
    $log .= $_;
  }
  elsif (/^Structure\s+-O \".*\" From -O \".*\" Source_Exons/) {
    s/Structure\s+$ts\s+From\s+$ts\s+Source_Exons/Source_exons/;
    $log .=  $_;
  }
  elsif (/^Visible\s+$ts\s+Locus_genomic_seq/) {
    s/Visible\s+$ts\s+Locus_genomic_seq/Locus/;
    $log .= $_;
  }
 # tags we don't want anymore
 elsif (/CDS_predicted_by/) {
   next;
 }
 elsif (/^Properties\s+$ts\s+RNA/) {
#    print "Transcript\n";
   next;
 }
 elsif (/^Structure\s+-O \".*\" From $ts Source $ts \"(\S+)\" $ts/) {
#    print "SMap S_parent Canonical_parent $1\n";
#    $log .= "SMap S_parent Sequence $1\n";
  # don't need this as the Sequence object will create the SMap relationship
   next;
  }
 
 # keep the rest
 else {
   $log .= $_;
 }
}


print SANGER "\n\n\n$log" if ($lab eq "HX");
print STLOUIS "\n\n\n$log" if ($lab eq "RW");

$db->close;


close(SANGER);
close(STLOUIS);

exit(0);
