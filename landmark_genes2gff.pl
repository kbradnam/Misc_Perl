#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

#my $database = "/wormsrv1/autoace";
my $database = "/nfs/disk100/wormpub/DATABASES/current_DB";

# get list of landmark genes from database, and add to hash
my %landmarks;
&get_landmark_genes;

# now loop through GFF files to look for existing gene spans

my @chromosomes = qw( I II III IV V X );
 
foreach (@chromosomes) {
  my $gff = "$database/CHROMOSOMES/CHROMOSOME_$_.gff";
  open (GFF,"<$gff") || die "Can't open $gff\n";
  my $chromosome = "CHROMOSOME_$_";
  
  while (<GFF>) {
    
    # only want to match the following lines
    # CHROMOSOME_II   gene    gene    23347   24428   .       +       .       Gene "WBGene00005017"
    next unless (/gene/ && /WBGene/);

    # more exact check by splitting line and checking fields
    my @data = split(/\t/);
    next unless ($data[1] eq "gene"&& $data[2] eq "gene");

    # modify 9th GFF column to look up in hash to get CGC name (or Public name)
    my $gene = $data[8];
    $gene =~ s/Gene //;
    $gene =~ s/\"//g;
    chomp($gene);
    if($landmarks{$gene}){
      print "$data[0]\tlandmark\tgene\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\tLocus $landmarks{$gene}\n";
    }
  }
  close GFF;

}

# II framework gene  5651408 5652390 .       -       .       Locus bli-2




##############################
# get list of landmark genes
##############################

sub get_landmark_genes{

  my $tace = &tace;
  my $database = "/wormsrv1/geneace";
  
  my $db = Ace->connect(-path  => $database,
			-program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
  
  # only want landmark genes which will be in GFF files, i.e. those with Sequence name 
  # this is virtually all of them
  my @landmarks = $db->fetch(-query=>"Find Gene WHERE Landmark_gene AND Sequence_name");

  foreach my $gene (@landmarks){
    my $public_name = $gene->Public_name->name;
    $landmarks{$gene->name} = $public_name;
  }  
  $db->close;

}
