#!/usr/local/bin/perl5.6.1 -w
#
# Convert old allele data to new allele model (to cope for KO deletion alleles)
#
# Keith Bradnam 17/12/2002

use strict;
use lib '/wormsrv2/scripts/';
use Wormbase;

$|=1;


my $ts = "-O \".*-*_*_*\"";



while (<>) {
	
  # No Source tag in new model
  if (/^Source\s+$ts/) {
    s/^Source\s+$ts Gene/Gene/;
    print "$_";
  }
  # What to do with sequences
  elsif (/^Sequence_details\s+$ts Sequence\s+$ts \"([\w\.]+)\"/) {
    my $seq = $1;
    # if sequences is a CDS, use new Predicted_gene tag
    if($seq =~ /\./){
      s/^Sequence_details\s+$ts Sequence\s+$ts \"([\w\.]+)\"/Predicted_gene \"$seq\"/;
    }
    # else use existing (but moved) Sequence tag
    else{
       s/^Sequence_details\s+$ts Sequence\s+$ts \"([\w\.]+)\"/Sequence \"$seq\"/;
    }
    print "$_";
  }
 # Allelic_difference + [A/G] style info  
  elsif (/^Sequence_details\s+$ts Allelic_difference\s+$ts \"(\[\w+\\\/\w+\])\"/) {
    my $diff = $1;
    s/^Sequence_details\s+$ts Allelic_difference\s+$ts \"(\[\w+\\\/\w+\])\"/Substitution \"$diff\"/;    
    print "$_";
  }

  # Allelic_difference + Deletion tag
  elsif (/^Sequence_details\s+$ts Allelic_difference\s+$ts \"Deletion\"/) { 
    s/^Sequence_details\s+$ts Allelic_difference\s+$ts \"Deletion\"/Deletion/;
    print "$_";							  
  }

 # Insertion tags
  elsif (/^Description\s+$ts\s+Insertion/) {
    s/Description\s+$ts\s+Insertion/Insertion/;
    print "$_";
  }

 # Insertion tags
  elsif (/^Description\s+$ts\s+Deletion/) {
    s/Description\s+$ts\s+Deletion/Deletion/;
    print "$_";
  }

 # keep the rest
 else {
   print "$_"
 }
}



exit(0);
