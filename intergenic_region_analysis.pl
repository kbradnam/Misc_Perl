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
my $dbdir       = "/Korflab/Data_sources/WormBase/WS150";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 

# store sizes of chromosomes in hash (calculated from WS150 dna files)
my %chr_size = ("I"   => "15072418",
		"II"  => "15279313",
		"III" => "13783317",
		"IV"  => "17493785",
		"V"   => "20922233",
		"X"   => "17718851");


# need array to store representation of chromosome
my @dna;

# hash to keep track of size of intergenic regions (key = size, value = count)
my %intergenic;
# separate hashes to count intergenic sizes when between Forward (F) and Reverse (R) genes (and all combinations thereof)
my %intergenic_FF;
my %intergenic_FR;
my %intergenic_RF;
my %intergenic_RR;


###########################################
# get gene spans from GFF file
###########################################        


foreach my $chromosome (@chromosomes) {
    
    # reset dna array and recreate to size of chromosome
    # all elements are just 'I' at this stage to make a virtual chromosome with just intergenic states (I) for now
    @dna =();
    @dna = ("I") x $chr_size{$chromosome};
    
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";

    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;
	
	# ignore non gene lines
	my @gff_line = split /\t/;

	next unless ($gff_line[1] eq "gene" && $gff_line[2] eq "gene");

	my $start     = $gff_line[3];
	my $stop      = $gff_line[4];
	my $length    = $stop-$start+1;
	my $direction = $gff_line[6];

	# substitute the appropriate part of our virtual chromosome with pluses (+) for forward genes
	# and minuses (-) for reverse genes
	my @replace = ("$direction") x $length;	
	splice(@dna,$start,$length, @replace);
	
    }
    close(GFF_in);

    # now to loop through array to count sizes of intergenic regions    
    # need to track endpoints to count size of each intergenic region
    my $intergenic_start;
    my $intergenic_end;
    my $size;

    # need a marker to be able to move past first non-coding (telomeric) sequence
    my $flag = 0;
    
    # need to start at 1 as DNA coordinates from GFF start at 1
    # whereas @dna array starts at 0.  I.e. we effectively always ignore $dna[0]   
    SEQ: for (my $i = 1; $i<$chr_size{$chromosome};$i++){

	# want to ignore the first stretch of non-coding sequence at the telomer
	# as this is not strictly intergenic sequence
	if($flag == 0){
	    while($dna[$i] eq "I"){
		$i++;
	    }
	    # now must be in the first coding region so can change $flag value
	    $flag = 1;
	}
	
	
	# wait until you enter an intergenic region
	if($dna[$i] eq "I"){

	    $intergenic_start = $i;

	    # what direction was previous gene?
	    my $previous = $dna[$i-1];
	    
            # now start counting until you leave intergenic region
	    $size = 0;

	    while($dna[$i] eq "I"){
		$size++;
		$i++;
		# Need to quit if we are at the end of the chromosome so
		# last telomeric region doesnt't get counted as intergenic
		last SEQ if ($i == $chr_size{$chromosome});
	    }

	    # add size details to hash;
	    $intergenic_end = $i-1;
	    
	    # what is direction of next gene?
	    my $next = $dna[$i];
		
	    print "$chromosome,$intergenic_start,$intergenic_end,$size,${previous}${next}\n";
	    $intergenic{$size}++;
	    $intergenic_FF{$size}++ if ($previous eq "+" && $next eq "+");
	    $intergenic_FR{$size}++ if ($previous eq "+" && $next eq "-");
	    $intergenic_RF{$size}++ if ($previous eq "-" && $next eq "+");
	    $intergenic_RR{$size}++ if ($previous eq "-" && $next eq "-");
	}
    }    
}

#foreach my $key (sort {$intergenic{$a} <=> $intergenic{$b}} (keys %intergenic)){
#    print "$key - $intergenic{$key}\n";
#}

exit(0);

