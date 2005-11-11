#!/usr/bin/perl -w
#
# find_genes_in_introns.pl
#
# by Keith Bradnam
# 7/1/2005
#
#############################################################################################

use strict;

#############
# Paths etc #
#############

my $dbdir       = "/Korflab/Data_sources/WormBase/WS150";    # Database path
#my $dbdir       = glob("~Nod/Work/WS140");    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 


# load hash of CDS name 2 peptide length
#my $file = glob("~Nod/cds2peptide_length.txt");
my $file = glob("~keith/Work/Genes_in_introns/cds2peptide_length.WS150.txt");

open(CDS, "<$file") || die "Could't open $file\n";

# store CDS name as key, peptide length as value
my %cds2peptide;
while(<CDS>){
    chomp;
    s/\"//g;
    my @a = split(/\s+/);
    $cds2peptide{$a[0]} = $a[1];
}
close(CDS);

##############################################################
# get gene and Coding_transcript intron spans from GFF files
##############################################################

# keep count of genes in introns
my $count = 0; 

foreach my $chromosome (@chromosomes) {
#    print "\nProcessing Chromosome $chromosome\n";

#    exit if ($chromosome ne "I");
    my @genes;     # gene name
    my %genecount; # count of genes of each name
    my %genespan;  # store start and stop coordinates

    my %intronspan;# span of each intron
	
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;
	
	# ignore non gene lines
	my @gff_line = split /\t/;

	# get coords of potential nested CDSs
	if(($gff_line[1] eq "curated") && ($gff_line[2] eq "CDS")){

	    # get gene name from last gff column
	    my $name = $gff_line[8];	  
	    $name =~ s/CDS \"([\w\d]+\.[\d\w]+)\" ;.*/$1/;

	    # what strand is gene on?
	    my $nested_strand = $gff_line[6];
	    
	    $genecount{$name}++;
	    $genespan{$name} = [$gff_line[3],$gff_line[4],$nested_strand];
	    
	    # Check that gene names don't appear more than once in each GFF file
	    if($genecount{$name}>1){
		die "ERROR: $name appears $genecount{$name} times\n"; 
	    }
	    else{
		# load up array with all genes on that chromosome
		push(@genes,$name);
	    }	    
	}
	# now look for possible flanking introns, this will include some that are in UTRs
	elsif(($gff_line[1] eq "Coding_transcript") && ($gff_line[2] eq "intron")){

	    # ignore unconfirmed introns
	    next if ($gff_line[8] !~ m/Confirmed_/);

#	    # ignore confirmed introns
#	    next if ($gff_line[8] =~ m/Confirmed_/);

	    # get name of intron containing transcript from last gff column
	    my $transcript = $gff_line[8];
	    $transcript =~ s/Transcript \"([\w\d\.]+)\".*/$1/;

	    # what strand is this on?
	    my $flanking_strand = $gff_line[6];

	    
	    # convert to CDS name if transcript has a dot number suffix
	    # i.e. remove second dot suffix
	    ($transcript =~ s/\.\d+$//) if($transcript =~ m/[\w\d]+\.[\w\d]+\./);

	    # use intron coords as name in hash, this means multiple transcript objects with the same intron
	    # will just overwrite each other and this will remove redundancy	   
	    my $name = "$gff_line[3]-$gff_line[4]";
	    $intronspan{$name} = [$gff_line[3],$gff_line[4],$transcript,$flanking_strand];
	}
	else{
	    # skip everything else
	    next;
	}
    }

    close(GFF_in);

    # now cycle through introns looking for those containing a gene
    foreach my $intron (keys %intronspan){
	
	my $intron_start  = $intronspan{$intron}->[0];
	my $intron_end    = $intronspan{$intron}->[1];
	my $intron_length = $intron_end - $intron_start +1;
	my $parent        = $intronspan{$intron}->[2];
	my $intron_strand = $intronspan{$intron}->[3];
	
	foreach my $newgene (@genes){

	    my $gene_start  = $genespan{$newgene}->[0];
	    my $gene_end    = $genespan{$newgene}->[1];
	    my $gene_strand = $genespan{$newgene}->[2];
	    my $gene_length = $gene_end - $gene_start +1;
	    
	    # only looking for genes wholly contained within other genes
	    if (($gene_start >= $intron_start) && ($gene_end <= $intron_end)){
		$count++;
#		my $out1  = sprintf("%4d %1s %12s %8d %4d %8d %6d",$count,$chromosome,$parent,$cds2peptide{$parent},$intron_start,$intron_end,$intron_length);
#		my $out2 .= sprintf("%12s %8d %8d %6d",$newgene,$gene_start,$gene_end,$gene_length);
#		my $out3 .= sprintf("%4d",$cds2peptide{$newgene});

		my $result = "";
		$result = "@" if($cds2peptide{$newgene} < $cds2peptide{$parent});
#		print "$out1  ->  $out2 $out3 (aa) $result\n";
		print "$count,$chromosome,$parent,$cds2peptide{$parent},$intron_start,$intron_end,$intron_length,$newgene,$gene_start,$gene_end,$gene_length,$cds2peptide{$newgene},${gene_strand}${intron_strand}\n";
	    }	    
	}
    } 
}

