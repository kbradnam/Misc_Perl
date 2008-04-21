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
my $nm;
my $kld;
my $sw;
my $pcc;
my $ed;
my $gapped;
my $global;

GetOptions ("motif_dir=s"  => \$motif_dir,
			"nm" => \$nm,
			"kld" => \$kld,
			"sw" => \$sw,
			"pcc" => \$pcc,
			"ed" => \$ed,
			"gapped" => \$gapped,
			"global" => \$global);


die "must specify at least distance calculating method from -nm, -kld, -sw, -pcc, or -ed\n" if (!$nm && !$kld && !$sw && !$pcc && !$ed);
die "must specify -gapped, -global or both options\n" if (!$gapped && !$global);

# load methods as specified by command line options
my @methods;
push(@methods,"NM")  if ($nm);
push(@methods,"KLD") if ($kld);
push(@methods,"SW")  if ($sw);
push(@methods,"PCC") if ($pcc);
push(@methods,"ED")  if ($ed);



# read directory containing motifs
opendir(DIR, "$motif_dir") || die "Can't read directory: $motif_dir\n";
my @files = grep { /\.xms/ && -f "$motif_dir/$_" } readdir(DIR);
closedir(DIR);

my $motifset = Motifset->new();

foreach my $file (@files){
	$motifset->nestedmica_import("$motif_dir/$file");	
} 

foreach my $method (@methods){
	if($global){
		$motifset->score_motifset("score"=>$method,"align"=>"global");
		$motifset->splitstree("treefile.global.$method");		
	}
	if($gapped){
		$motifset->score_motifset("score"=>$method,"align"=>"gapped");
		$motifset->splitstree("treefile.gapped.$method");		
	}
}
exit(0);




