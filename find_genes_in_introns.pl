#!/usr/bin/perl
#
# find_genes_in_introns.pl
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;
 
#############
# Paths etc #
#############

my $dbdir       = "/Korflab/Data_sources/WormBase/WS150";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 

print "Using $dbdir for database directory\n";


##############################################################
# first get list of RNA genes and Pseudogenes from database
# these will potentially be filtered out of analysis
##############################################################

# rna genes
my $rna_data = glob("~keith/Work/Overlapping_genes/WS150_RNA_genes.ace");
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
my $pseudogene_data = glob("~keith/Work/Overlapping_genes/WS150_pseudogenes.ace");
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



##############################################################
# get gene and Coding_transcript intron spans from GFF files
##############################################################

# keep count of genes in introns
my $count = 0; 

foreach my $chromosome (@chromosomes) {
    print "\nProcessing Chromosome $chromosome\n";

    my @genes;      # list of gene name
    my %genecount;  # count of genes of each name
    my %genespan;   # store start and stop coordinates
    my %intronspan; # span of each intron
	
	############################################################
	#
	# PART 1 - Process GFF file to extract relevant information
	#
	#############################################################
	
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";

    while (<GFF_in>) {
		chomp;
	
		# skip header info
		next if m/^\#/;
	
		# ignore non gene lines
		my @gff_line = split /\t/;

		# Extract gene data
		if($gff_line[1] eq "gene"){

		    # get gene name from last gff column
		    my $name = $gff_line[8];
		    $name =~ s/Gene \"(WBGene\d+)\"/$1/;
	    
		    $genecount{$name}++;
		    $genespan{$name} = [$gff_line[3],$gff_line[4]];
	    
		    # Check that gene names don't appear more than once in each GFF file
			die "ERROR: $name appears $genecount{$name} times\n" if($genecount{$name}>1); 

			# load up array with all genes on that chromosome
			push(@genes,$name);	    
		}

		# Extract intron data
		elsif(($gff_line[1] eq "Coding_transcript") && ($gff_line[2] eq "intron")){

		    # use intron coords as name in hash, this means multiple transcript objects with the same intron
		    # will just overwrite each other and this will remove redundancy	   
		    my $name = "$gff_line[3]-$gff_line[4]";
		    $intronspan{$name} = [$gff_line[3],$gff_line[4]];
		}
		else{
		    # skip everything else
		    next;
		}
	}
	
    close(GFF_in);


	############################################################
	#
	# PART 2 - Process information to determine overlaps
	#
	#############################################################

    # now cycle through introns looking for those containing a gene
    foreach my $intron (keys %intronspan){

		my $intron_start  = $intronspan{$intron}->[0];
		my $intron_end    = $intronspan{$intron}->[1];
		my $intron_length = $intron_end - $intron_start +1;

		foreach my $newgene (@genes){

		    # determine gene status (CDS, Transcript, Pseudogene)
		    my $status = "CDS";
		    $status    = "RNA"    if ($rna_genes{$newgene});
		    $status    = "Pseudo" if ($pseudogenes{$newgene}); 

		    # skip non-CDS genes
	#	    next if ($status ne "CDS");
    
		    my $gene_start  = $genespan{$newgene}->[0];
		    my $gene_end    = $genespan{$newgene}->[1];
		    my $gene_length = $gene_end - $gene_start +1;

		    # only looking for genes wholly contained within other genes
		    if (($gene_start >= $intron_start) && ($gene_end <= $intron_end)){
				$count++;
				my $out1  = sprintf("%4d %1s %12s %8d %8d %6d",$count,$chromosome,$intron,$intron_start,$intron_end,$intron_length);
				my $out2 .= sprintf("%12s %6s %8d %8d %6d",$newgene,$status,$gene_start,$gene_end,$gene_length);
				print "$out1  ->  $out2\n";

				# print out in CSV format
		#		my $out1  = sprintf("%4d,%1s,%8d,%8d,%6d",$count,$chromosome,$intron_start,$intron_end,$intron_length);
	    #		my $out2 .= sprintf("%12s,%8d,%8d,%6d",$newgene,$gene_start,$gene_end,$gene_length);
	    #		print "$out1,$out2\n";

		    }
		}
    }
}

exit(0);

