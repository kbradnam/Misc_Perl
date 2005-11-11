#!/usr/bin/perl -w
#
# find_overlapping_genes.pl
#
# by Keith Bradnam
# 3/1/2005
#
#############################################################################################

use strict;

#############
# Paths etc #
#############

my $tace        = "/Korflab/bin/tace";                       # tace executable path
my $dbdir       = "/Korflab/Data_sources/WormBase/WS140";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 

print "Using $dbdir for database directory\n";
print "Using $tace for tace binary\n";


##############################################################
# first get list of RNA genes and Pseudogenes from database
# these will potentially be filtered out of analysis
##############################################################

# rna genes
my $rna_data = glob("~keith/Work/Overlapping_genes/WS140_RNA_genes.txt");
open(RNA, "<$rna_data") || die "Could't open $rna_data\n";

# store Gene IDs of RNA genes in hash, value is not really used
my %rna_genes;
while(<RNA>){
    chomp;
    if(m/^Gene/){
	s/Gene : \"(WBGene\d+)\"/$1/;
	$rna_genes{$_} = 1;
    }
}
close(RNA);

# pseudogenes
my $pseudogene_data = glob("~keith/Work/Overlapping_genes/WS140_pseudogenes.txt");
open(PSEUDO, "<$pseudogene_data") || die "Could't open $pseudogene_data\n";

# store Gene IDs of pseudogenes in hash, value is not really used
my %pseudogenes;
while(<PSEUDO>){
    chomp;
    if(m/^Gene/){
	s/Gene : \"(WBGene\d+)\"/$1/;
	$pseudogenes{$_} = 1;
    }
}
close(PSEUDO);



###########################################
# get gene spans from GFF file
###########################################        

# keep count of overlapping genes
my $count = 0; 

foreach my $chromosome (@chromosomes) {
    print "\nProcessing Chromosome $chromosome\n";

    my @genes;     # gene name
    my %genecount; # count of genes of each name
    my %genespan;  # store start and stop coordinates
    
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;
	
	# ignore non gene lines
	my @gff_line = split /\t/;
	next unless ($gff_line[1] eq "gene");
	
	# get gene name from last gff column
	my $name = $gff_line[8];
	$name =~ s/Gene \"(WBGene\d+)\"/$1/;
	
	$genecount{$name}++;
	$genespan{$name} = [$gff_line[3],$gff_line[4]];
	
	# Check that gene names don't appear more than once in each GFF file
	if($genecount{$name}>1){
	    die "ERROR: $name appears $genecount{$name} times\n"; 
	}
	else{
	    # load up array with all genes on that chromosome
	    push(@genes,$name);
	}
    }
    close(GFF_in);


    
    # now cycle through genes looking for overlaps
    foreach my $gene (@genes){

	# set query gene status
	my $query_status = "C";
	$query_status    = "R" if($rna_genes{$gene});
	$query_status    = "P" if($pseudogenes{$gene});

	my $query_start  = $genespan{$gene}->[0];
	my $query_end    = $genespan{$gene}->[1];
	my $query_length = $query_end - $query_start +1;

	foreach my $newgene (@genes){

	    # set target gene status
	    my $target_status = "C";
	    $target_status    = "R" if($rna_genes{$newgene});
	    $target_status    = "P" if($pseudogenes{$newgene});

	    my $final_status = "${target_status}in${query_status}";
	    
	    my $target_start  = $genespan{$newgene}->[0];
	    my $target_end    = $genespan{$newgene}->[1];
	    my $target_length = $target_end - $target_start +1;

	    # only looking for genes wholly contained within other genes
	    if (($target_start > $query_start) && ($target_end < $query_end)){
		$count++;
		my $out1 = sprintf("%4d %1s %4s %12s %8d %8d %6d",$count,$chromosome,$final_status,$gene,$query_start,$query_end,$query_length);
		my $out2   .= sprintf("%12s %8d %8d %6d",$newgene,$target_start,$target_end,$target_length);
		print "$out1  ->  $out2\n";
	    }
	    
	}
    }
  
}

