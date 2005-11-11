#!/usr/bin/perl -w

use lib "/Korflab/lib/perl";
use Ace;
use strict;

# open a local database connection
my $db = Ace->connect(-path  =>  '/Korflab/Data_sources/WormBase/WS140') ;

# fetch C.elegans CDSs
my @genes = $db->fetch(-class => "elegans_CDS");

my $counter=0;

foreach my $gene (@genes){
    $counter++;
    print "$counter) Gene is $gene\n";  

    # get some gene info
    my $lab = $gene->From_laboratory;
    my $status = $gene->Prediction_status;
    print "Lab is $lab, status is $status\n\n";
    
    # kill AcePerl objects
    $gene->DESTROY();
}


print "\n";
$db->close;
exit;





