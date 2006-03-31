#!/usr/bin/perl -w
#
# intergenic_region_analysis.pl
#
# Last updated by: $Author$
# Last updated on: $Date$
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


# will get chromosome sequence as a string and then split into an array
my $seq;
my @dna;
my $chr_length;
my $intergenic;

foreach my $chromosome (@chromosomes) {

    print "Processing chromosome $chromosome\n";
    
    # get chromosome sequence
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
    @dna = split(//,$seq);
    unshift(@dna,"S");
    
    $chr_length = scalar(@dna)-1;
    print "Chromosome $chromosome - $chr_length bp\n";


    

    # now scan GFF file for details of genes
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

	# substitute the appropriate part of the chromosome sequence with pluses (+) for forward genes
	# and minuses (-) for reverse genes.  Need to subtract 1 from $start as base 1 of a sequence
	# will be position 0 in the @dna array
	my @replace = ("$direction") x $length;	
	splice(@dna,$start,$length, @replace);
	
    }
    close(GFF_in);

    print "Here\n";

    
    # now to loop through array to count sizes of intergenic regions    
    # need to track endpoints to count size of each intergenic region
    my $intergenic_start;
    my $intergenic_end;
    my $intergenic_size;

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
	    my $previous = $dna[$i-1];
	    
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
	    my $next = $dna[$i];
		
	    print "$chromosome,$intergenic_start,$intergenic_end,$intergenic_size,${previous}${next}\n";

	    # now get sequence, need t
	    $intergenic = substr($seq,$intergenic_start,$intergenic_size);
	    print "$intergenic\n\n";
	    
	}
    }    
}


exit(0);

