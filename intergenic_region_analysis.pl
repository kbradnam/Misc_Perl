#!/usr/bin/perl -w
#
# intergenic_region_analysis.pl
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use Getopt::Long;

########################
# Command line options
########################

my $seqs;  # write intergenic output to files
my $stats; # write details of intergenic regions to screen
my $min;   # restrict a subset of regions to those with a minimum size
my $max;   # restrict a subset of regions to those with a maximum size

GetOptions ("seqs"     => \$seqs,
            "stats"    => \$stats,
	    "min=i"    => \$min,
	    "max=i"    => \$max);


die "Specify -seqs and/or -stats\n" if (!$seqs && !$stats);
if($seqs){
    die "Specify min + max\n" if (($min && !$max) || (!$min && $max));    
    die "Max needs to be bigger than min\n" if ($max < $min);
}
# Set $subset flag to be true if using min and max
my $subset = 0;
$subset = 1 if ($min && $max);


#############
# Paths etc 
#############

my $tace        = "/Korflab/bin/tace";                       # tace executable path
my $dbdir       = "/Korflab/Data_sources/WormBase/WS150";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 



###################
# Misc. variables
###################

my $seq;              # will hold chromosome sequence as a string...
my @dna;              # ...which will then get split into an array...
my @operons;          # ...and then copied for a separate array to process operon data
my $chr_length;       # obvious really
my $intergenic;       # will store sequence of an intergenic region
my $intergenic_start; # start coord
my $intergenic_end;   # stop coord
my $intergenic_size;  # length
my $previous;         # strandedness of gene (+ or -) 5' to intergenic region
my $next;             # strandedness of gene (+ or -) 3' to intergenic region
my $in_operon;        # true if intergenic region is in operon, otherwise false


if($seqs){
    # will be writing various output files for intergenic sequences
    open (FF, ">intergenic_FF.dna")   || die "Failed to open intergenic_FF file\n\n";
    open (FR, ">intergenic_FR.dna")   || die "Failed to open intergenic_FR file\n\n";
    open (RF, ">intergenic_RF.dna")   || die "Failed to open intergenic_RF file\n\n";
    open (FFO, ">intergenic_FFO.dna") || die "Failed to open intergenic_FFO file\n\n";
    open (FRO, ">intergenic_FRO.dna") || die "Failed to open intergenic_FRO file\n\n";
    open (RFO, ">intergenic_RFO.dna") || die "Failed to open intergenic_RFO file\n\n";
    
    # second set of files for intergenic regions between 50-1000 bp
    open (FF2, ">intergenic_FF_subset.dna")   || die "Failed to open intergenic_FF2 file\n\n";
    open (FR2, ">intergenic_FR_subset.dna")   || die "Failed to open intergenic_FR2 file\n\n";
    open (RF2, ">intergenic_RF_subset.dna")   || die "Failed to open intergenic_RF2 file\n\n";
    open (FFO2, ">intergenic_FFO_subset.dna") || die "Failed to open intergenic_FFO2 file\n\n";
    open (FRO2, ">intergenic_FRO_subset.dna") || die "Failed to open intergenic_FRO2 file\n\n";
    open (RFO2, ">intergenic_RFO_subset.dna") || die "Failed to open intergenic_RFO2 file\n\n";
}

####################################
# Main loop through each chromosome
####################################

foreach my $chromosome (@chromosomes) {
    
    # get chromosome sequence, load into $seq string
    open (DNA, "<$gffdir/CHROMOSOME_${chromosome}.dna") || die "Failed to open dna file\n\n";
    $seq = "";
    while(my $tmp =<DNA>){
        chomp($tmp);
        # skip header line
        next if ($tmp =~ m/^>/);
        $seq .= $tmp;
    }
    close(DNA);


   
    # to make things easier to calculate, will add a character to left end of $seq and @dna, such that position 1 of dna
    # sequence becomes array element 1 (rather than zero).  Add 'S' for 'Start'
    $seq = "S".$seq;

    # split to an array and calculate size
    @dna=();
    @operons=();
    @dna = split(//,$seq);
    unshift(@dna,"S");
    $chr_length = scalar(@dna)-1;

    # now want a copy of @dna to treat the operons separately (bit of a waste of memory to do it this way)
    @operons = @dna;


    
    #################################################
    # GFF step: extract details of genes and operons
    #################################################
    
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";

    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;
	
	# ignore any line which isn't gene or operon
	my @gff_line = split /\t/;
	next unless (($gff_line[1] eq "gene") || ($gff_line[1] eq "operon"));
	next unless (($gff_line[2] eq "gene") || ($gff_line[2] eq "operon"));

	# extract details
	my $start     = $gff_line[3];
	my $stop      = $gff_line[4];
	my $length    = $stop-$start+1;
	my $direction = $gff_line[6];

	# substitute the appropriate part of @dna and @operons with pluses (+) for forward genes
	# and minuses (-) for reverse genes.	
	my @replace = ("$direction") x $length;	

	splice(@dna,$start,$length, @replace)     if ($gff_line[1] eq "gene");
	splice(@operons,$start,$length, @replace) if ($gff_line[1] eq "operon");
	
    }
    close(GFF_in);

    
    ########################################
    # Main loop through chromosome sequence
    ########################################
    
    
    # need a marker to be able to move past first non-coding (telomeric) sequence
    my $flag = 0;

    SEQ: for (my $i = 1; $i<$chr_length;$i++){

	# want to ignore the first stretch of non-coding sequence at the telomer
	# as this is not strictly intergenic sequence
	if($flag == 0){
	    while(($dna[$i] ne "+") && ($dna[$i] ne "-")){
		$i++;
	    }
	    # now must be in the first coding region so can change $flag value
	    $flag = 1;
	}
	
	
	# wait until you enter an intergenic region, i.e. not a plus or minus
	if(($dna[$i] ne "+") && ($dna[$i] ne "-")){

	    $intergenic_start = $i;

	    # what direction was previous gene?
	    $previous = $dna[$i-1];
	    
            # now start counting until you leave intergenic region
	    $intergenic_size = 0;

	    while(($dna[$i] ne "+") && ($dna[$i] ne "-")){
		$intergenic_size++;
		$i++;
		# Need to quit if we are at the end of the chromosome so
		# last telomeric region doesnt't get counted as intergenic
		last SEQ if ($i == $chr_length);
	    }

	    # add size details to hash;
	    $intergenic_end = $i-1;
	    
	    # what is direction of next gene?
	    $next = $dna[$i];


	    # are we in an operon?
	    if(($operons[$i-1] eq "+") || ($operons[$i-1] eq "-")){
		$in_operon = 1;
		print "$chromosome,$intergenic_start,$intergenic_end,$intergenic_size,${previous}${next},1\n" if ($stats);
		
	    }
	    else{
		$in_operon = 0;
		print "$chromosome,$intergenic_start,$intergenic_end,$intergenic_size,${previous}${next},0\n" if ($stats);
	    }
	    
	    # now get sequence
	    $intergenic = substr($seq,$intergenic_start,$intergenic_size);

	    # write sequences to file?
	    if($seqs){

		if(($previous eq $next) && ($in_operon == 0)){
		    print FF ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}
		if(($previous eq "+") && ($next eq "-") && ($in_operon == 0)){
		    print FR ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}
		if(($previous eq "-") && ($next eq "+") && ($in_operon == 0)){
		    print RF ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}
		if(($previous eq $next) && ($in_operon == 1)){
		    print FFO ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}
		if(($previous eq "+") && ($next eq "-") && ($in_operon == 1)){
		    print FRO ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}
		if(($previous eq "-") && ($next eq "+") && ($in_operon == 1)){
		    print RFO ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		}

		# also print subset of sequences in certain size range?
		if($subset){
		    if(($previous eq $next) && ($in_operon == 0) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print FF2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		    if(($previous eq "+") && ($next eq "-") && ($in_operon == 0) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print FR2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		    if(($previous eq "-") && ($next eq "+") && ($in_operon == 0) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print RF2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		    if(($previous eq $next) && ($in_operon == 1) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print FFO2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		    if(($previous eq "+") && ($next eq "-") && ($in_operon == 1) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print FRO2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		    if(($previous eq "-") && ($next eq "+") && ($in_operon == 1) && ($intergenic_size >= $min) && ($intergenic_size <= $max)){
			print RFO2 ">${chromosome}_${intergenic_start}_${intergenic_end}\n$intergenic\n";
		    }
		}
	    }
	}
    }    
}

if($seqs){
    close(FF);
    close(FR);
    close(RF);
    close(FFO);
    close(FRO);
    close(RFO);
    close(FF2);
    close(FR2);
    close(RF2);
    close(FFO2);
    close(FRO2);
    close(RFO2);
}



exit(0);

