#!/usr/bin/perl
#
# meme2xms.pl 
#
# A script to take position weight matrix information (copied from a meme output file) and convert to a nested mica compatible format
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Motif;

die "Specify a file containing meme motifs\n" if (@ARGV != 1);

my $file = $ARGV[0];

my $motifset = Motifset->new;

$motifset->meme_import("$file");       

# output logo file of motif using Motif.pm (this generates forward and reverse strand motif images)
$motifset->nestedmica_export("$file");


exit(0);

