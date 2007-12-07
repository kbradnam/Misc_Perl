#!/usr/bin/perl
#
# get_prashanth_worm_data.pl
#
# a script to extract some data from C.elegans GFF files needed by Prashanth
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;
use Getopt::Long;

########################
# Command line options
########################

my $test;  # test mode, only test one chromosome for speed

GetOptions (
	    "test"   => \$test);

#############
# Paths etc 
#############

my $dbdir       = "/Korflab/Data_sources/WormBase/WS180";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 
@chromosomes = qw (I) if ($test);

###################
# Misc. variables
###################

my $seq;              # will hold chromosome sequence as a string...
my @dna;              # ...which will then get split into an array...
my $chr_length;       # obvious really
my $intergenic;       # will store sequence of an intergenic region
my $intergenic_start; # start coord
my $intergenic_end;   # stop coord
my $intergenic_size;  # length
my $previous;         # strandedness of gene (+ or -) 5' to intergenic region
my $next;             # strandedness of gene (+ or -) 3' to intergenic region

####################################
# Main loop through each chromosome
####################################

foreach my $chromosome (@chromosomes) {
    
    # get chromosome sequence, load into $seq string
	if($test){
		open (DNA, "<CHROMOSOME_I.dna") || die "Failed to open dna file\n\n";
	}
    else{
		open (DNA, "<$gffdir/CHROMOSOME_${chromosome}.dna") || die "Failed to open dna file\n\n";
    }
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
    @dna = split(//,$seq);

    $chr_length = scalar(@dna)-1;

    #################################################
    # GFF step: extract details of genes
    #################################################

    if($test){
		open (GFF_in, "<CHROMOSOME_I.gff") || die "Failed to open gff file\n\n";
	}
    else{
		open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
	}
	
    while (<GFF_in>) {
		chomp;
	
		# skip header info
		next if m/^\#/;
	
		# ignore any line which isn't gene a gene
		my @gff_line = split /\t/;

		#  want gene lines
		next unless (($gff_line[1] eq "gene") && ($gff_line[2] eq "gene"));
	
		# extract details
		my $start     = $gff_line[3];
		my $stop      = $gff_line[4];
		my $length    = $stop-$start+1;
		my $direction = $gff_line[6];

		# substitute the appropriate part of @dna and with pluses (+) for forward genes
		# and minuses (-) for reverse genes.	
		my @replace = ("$direction") x $length;	

		splice(@dna,$start,$length, @replace);
    }
    close(GFF_in);
	    
    ########################################
    # Main loop through chromosome sequence
    ########################################
    
    # need a marker to be able to move past first non-coding (telomeric) sequence
    my $flag = 0;
	
	# for the 3' telomere calculation, need to keep track of the end of hte last gene
	my $previous_gene_end;
	
    SEQ: for (my $i = 1; $i<$chr_length;$i++){

		# want to ignore the first stretch of non-coding sequence at the telomere
		# as this is not strictly intergenic sequence
		if($flag == 0){
	    	while(($dna[$i] ne "+") && ($dna[$i] ne "-")){
				$i++;
	    	}
	    	# now must be in the first coding region so can change $flag value and print out
			# telomere coordinates
	    	$flag = 1;
			my $j = $i-1;
			print "Telomere 1 $j\n";
		}	

		# keep skipping forwards until you enter an intergenic region, i.e. not a plus or minus
		
		if (($dna[$i] eq "+") || ($dna[$i] eq "-")){
			$previous_gene_end = $i;
			next SEQ;	
		}

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
			if ($i == $chr_length){
				$previous_gene_end++;
				print "Telomere $previous_gene_end $i\n";
				last SEQ;
			}
	    }

	    # add size details to hash;
	    $intergenic_end = $i-1;
	    
	    # what is direction of next gene?
	    $next = $dna[$i];
		print "$chromosome $intergenic_start $intergenic_end $intergenic_size ${previous}${next}\n";
		}
	    
}


exit(0);

