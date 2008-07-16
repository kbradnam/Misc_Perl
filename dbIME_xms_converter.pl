#!/usr/bin/perl
#
# dbIME_motif_draw.pl
#
# A script to generate images of motifs (in PNG and SVG formats) from an xms file
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Motif;
use Cwd;

die "$ARGV[0] does not appear to be a valid XMS motif file\n" if ($ARGV[0] !~ m/\.xms$/);

# keep main part of file name for later
my $file = $ARGV[0];
$file =~ s/\.xms$//;

my $dir = getcwd;

my $motifset = Motifset->new();
$motifset->nestedmica_import("$ARGV[0]");

# output logo file of motif using Motif.pm (this generates forward and reverse strand motif images)
$motifset->logo("$dir");

# remove reverse complement file tas we don't need this
unlink("Rev".$file.".svg") || print "Couldn't remove reverse complement SVG file\n";

# convert to PNG (if convert tool is available)
if(-e "/opt/local/bin/convert"){
	system("/opt/local/bin/convert $file.svg $file.png");
}
else{
	print "ImageMagick 'convert' tool was not available to convert file to PNG format\n\n";
}
exit(0);




