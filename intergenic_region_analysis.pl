#!/usr/bin/perl -w
#
# intergenic_region_analysis.pl
#
# by Keith Bradnam
# 3/14/2005
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

# store sizes of chromosomes in hash (calculated from WS140 dna files)
my %chr_size = ("I"   => "15080553",
		"II"  => "15279312",
		"III" => "13783318",
		"IV"  => "17493786",
		"V"   => "20922232",
		"X"   => "17718851");


# need array to store representation of chromosome
my @dna;

# hash to keep track of size of intergenic regions (key = size, value = count)
my %intergenic;


###########################################
# get gene spans from GFF file
###########################################        


# fill array for chromosome with 'I' characters



foreach my $chromosome (@chromosomes) {

    # just use chromosome I for now.
    exit if ($chromosome ne "I");
    print "Chromosome $chromosome\n";
    
    # reset dna array and recreate to size of chromosome
    # all elements are just 'I' at this stage
    @dna =();
    @dna = ("I") x $chr_size{$chromosome};
    
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";

    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;
	
	# ignore non gene lines
	my @gff_line = split /\t/;

#	next unless ($gff_line[1] eq "curated" && $gff_line[2] eq "CDS");
	next unless ($gff_line[1] eq "gene" && $gff_line[2] eq "gene");

	my $start     = $gff_line[3];
	my $stop      = $gff_line[4];
	my $length    = $stop-$start+1;
	my $direction = $gff_line[6];
	my $name      = $gff_line[8];
	$name =~ s/CDS \"([\w.]+)\" ;.*/$1/;

	my @replace = ("$direction") x $length;
	
	splice(@dna,$start,$length, @replace);
    }
    close(GFF_in);

    # now to loop through array to count sizes of intergenic regions    
    # need to track endpoints to count size of each intergenic region
    my $intergenic_start;
    my $intergenic_end;
    my $size;
       
    # need to start at 1 as DNA coordinates from GFF start at 1
    # whereas @dna array starts at 0.  I.e. we effectively always ignore $dna[0]

    for (my $i = 1; $i<$chr_size{$chromosome};$i++){

	# wait until you enter an intergenic region
	if($dna[$i] eq "I"){

	    $intergenic_start = $i;
	    # now start counting until you leave intergenic region
	    $size = 0;
	    while($dna[$i] eq "I"){
		$size++;
		$i++;
	    }
	    # add size details to hash;
	    $intergenic_end = $i-1;
	    print "$intergenic_start - $intergenic_end) size is $size\n";
	    $intergenic{$size}++;
	}
    }
    
    
}

print "End\n";
exit(0);

