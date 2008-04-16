#!/usr/bin/perl
#
# IME_tree_prep.pl
#
# A script to calculate distances between a set of Nested Mica motifs, and prepare data to make a tree from those motifs
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use Motif;

########################
# Command line options
########################

my $motif_dir;  # directory containing nested MICA *.xms files

GetOptions ("motif_dir=s"  => \$motif_dir);


# read directory containing motifs
opendir(DIR, "$motif_dir") || die "Can't read directory: $motif_dir\n";
my @files = grep { /\.xms/ && -f "$motif_dir/$_" } readdir(DIR);
closedir(DIR);

my $motifset = Motifset->new();

foreach my $file (@files){
	$motifset->nestedmica_import("$motif_dir/$file");	
} 

foreach my $method qw(NM KLD SW PCC ED){
	$motifset->score_motifset("score"=>$method,"align"=>"global");
	$motifset->splitstree("treefile.global.$method");
	$motifset->score_motifset("score"=>$method,"align"=>"gapped");
	$motifset->splitstree("treefile.gapped.$method");
}
exit(0);




