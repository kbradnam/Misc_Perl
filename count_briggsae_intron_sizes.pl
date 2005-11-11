#!/usr/bin/perl -w
#
# count_briggsae_intron_sizes.pl
#
# by Keith Bradnam
# 6/6/2005
#
#############################################################################################

use strict;
use Getopt::Long;


########################
# Command line options
########################

my $count;       # separate mode for just counting introns in different size classes
my $interval;    # what size bins to count introns in?
my $limit;       # what size to stop counting introns at (use if you only want to count short introns for example)

GetOptions ("count"       => \$count,
	    "interval=i"  => \$interval,
	    "limit=i"     => \$limit);


##################################################
# Sanity check on specified command line options
##################################################

# set default limit if none specified, this will pretty much cover all real introns
$limit = 35000 if (!$limit);



#############
# Paths etc #
#############

my %all_introns;

# fourth hash for keeping track of just unique introns (different transcripts generate multiple GFF intron lines for the same intron)
my %unique_introns;

open (GFF_in, "<briggsae.gff") || die "Failed to open gff file\n\n";
while (<GFF_in>) {
    chomp;
    
    # skip header info
    next if m/^\#/;
    
    # get basic details from GFF line
    my @gff_line = split /\t/;
    
    # only want curated intron lines from GFF file, this ignores possible confirmed introns from UTRs
    next unless (($gff_line[1] eq "curated") && ($gff_line[2] eq "intron"));

    my $cds = $gff_line[8];
    $cds =~ s/Sequence \"//;
    $cds =~ s/\"//g;
    
    # need to ignore the same intron in different GFF lines
    # use chromosome number and coords as a unique key
    my $key = $cds.":".$gff_line[3]."-".$gff_line[4];

    
    if(exists $unique_introns{$key}){
	# can skip intron
	next;
    }
    else{
	# add to hash of unique introns so as not to count it again
	$unique_introns{$key} = "1";				   
	
	# Calculate length and increment counter in respective hash
	my $length = $gff_line[4] - $gff_line[3]+1;
	if($length >= 50000){
	    print "$cds $length bp\n";
	}
	$all_introns{$length}++;
    }
}	  	    
    
close(GFF_in);  
print "\n";       


    
#################################################################
# need to work out intron counts for specified bin size
#################################################################
my $bin_start;
my $bin_end;

my $max_intron_size = $limit;

for ($bin_start = 1; $bin_start < $max_intron_size; $bin_start += $interval){

    $bin_end = $bin_start+$interval-1;

    # quit if limit exceeded
    if($bin_start > $limit){
	exit(0);
    }
    
    # now work out total of intron in current bin
    my $all_bin_total         = 0;

    for(my $count = $bin_start; $count <= $bin_end; $count++){

	# can only sum up introns if they exist at that size
	$all_bin_total         += $all_introns{$count}         if ($all_introns{$count});
    }
    print "${bin_start}-${bin_end}\t$all_bin_total\n";

}

exit(0);
