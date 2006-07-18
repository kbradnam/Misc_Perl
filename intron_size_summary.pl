#!/usr/bin/perl -w
#
# intron_size_summary.pl
#
# A script to process and print intron information from a GFF file
#
# by Keith Bradnam
# 3/9/2005
#
#############################################################################################

use strict;
use Getopt::Long;


########################
# Command line options
########################

my $count;         # separate mode for just counting introns in different size classes
my $interval;      # what size bins to count introns in?
my $limit;         # what size to stop counting introns at (use if you only want to count short introns for example)
my $min;           # print GFF lines matching above minimum intron length
my $max;           # print GFF lines matching below maximum intron length
my $confirmed;     # only show introns confirmed by transcript data
my $unconfirmed;   # only show introns not confirmed by transcript data
my $all;           # show all introns (this is also the default)
my $release;       # which WormBase release to use (will default to latest)
my $print;         # send output to a file
my $ignore_n;      # ignore any intron which contains N's (repeat masked sequence)
my $flank;         # how many bases of flanking exon sequences to extract
my $donor;         # Only extract $end_size bp from start of intron 
my $donor_size;    # size of intron donor sequence to extract
my $acceptor;      # Only extract $end_size bp from start of intron 
my $acceptor_size; # size of intron donor sequence to extract
my $gc;            # only print/count introns that start with a GC
my $test;          # run in test mode

GetOptions ("count"           => \$count,
	    "interval=i"      => \$interval,
	    "limit=i"         => \$limit,
	    "min=i"           => \$min,
	    "max=i"           => \$max,
	    "confirmed"       => \$confirmed,
	    "unconfirmed"     => \$unconfirmed,
	    "all"             => \$all,
	    "release=s"       => \$release,
	    "print"           => \$print, 
	    "ignore_n"        => \$ignore_n,
	    "flank=i"         => \$flank,
	    "donor"           => \$donor,
	    "donor_size=i"    => \$donor_size,
	    "acceptor"        => \$acceptor,
	    "acceptor_size=i" => \$acceptor_size,
	    "gc"              => \$gc,
	    "test"            => \$test);


##################################################
# Sanity check on specified command line options
##################################################

# set default limit if none specified, this will pretty much cover all real introns
$limit = 25000 if (!$limit);

# which WormBase release to use
if($release && $release !~ m/WS\d{3}/){
    die "please specify WormBase release as 'WSXXX' (where XXX is an integer)\n";
}
# set a default release   
$release = "WS160" if (!$release);

# Check confirmed/unconfirmed/all flags
if(($all && $confirmed) || ($all && $unconfirmed) || ($confirmed && $unconfirmed)){
    die "Only specify one of -confirmed, -unconfirmed and -all\n";
}

#default to -all if not specified
($all = 1) if(!$all && !$confirmed && !$unconfirmed);

# min should be less than max
if(($min && $max) && ($min > $max)){
   die "$min needs to be less than or equal to $max\n";
}

# set default flanking sequence to be zero
$flank = 0 if(!$flank);

# check that min and max have been set
if(!$min || !$max){
    die "You must specify a minimum and maximum intron size (they may be the same)\n";
}

# check that if -donor is used that -donor_size is also used
if($donor && !$donor_size){
    die "You must use the -donor_size option with the -donor option\n";
}

# check that if -acceptor is used that -acceptor_size is also used
if($acceptor && !$acceptor_size){
    die "You must use the -acceptor_size option with the -acceptor option\n";
}

# check that both -donor and -acceptor are not specified at the same time
if($acceptor && $donor){
    die "You can only use -donor OR -acceptor in one run of this script, not both\n";
}


#############
# Paths etc #
#############
my $dbdir       = "/Korflab/Data_sources/WormBase/$release";    
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 
@chromosomes = qw(I) if ($test);                             # just use one chromosome when in test mode


# store data in simple hash, key = intron size, value = count
my %confirmed_introns;

# second hash for nasty unconfirmed introns
my %unconfirmed_introns;

# third hash for both good and bad introns
my %all_introns;

# fourth hash for keeping track of just unique introns (different transcripts generate multiple GFF intron lines for the same intron)
my %unique_introns;


if($print){
    open(OUT,">introns.${min}-${max}bp") || die "Couldn't open output file\n";
}

my @dna;
my $seq;
my $counter =0;

foreach my $chromosome (@chromosomes) {

    print "Processing chromosome $chromosome\n";
    
    # need to get DNA sequence for each chromosome unless just counting introns
    unless($count){
	$seq = &fetch_seq($chromosome);

	# turn into array
	@dna = "";
	@dna = split(//,$seq);
	$seq = "";
    }
    
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
	chomp;
	
	# skip header info
	next if m/^\#/;

	# get basic details from GFF line
	my @gff_line = split /\t/;

	# only want curated intron lines from GFF file, this ignores possible confirmed introns from UTRs
	next unless (($gff_line[1] eq "curated") && ($gff_line[2] eq "intron"));

	# skip unconfirmed introns if -confirmed tag has been set
	next if ($confirmed && ($gff_line[8] !~ m/Confirmed/));

	# skip confirmed introns if -unconfirmed tag has been set
	next if ($unconfirmed && ($gff_line[8] =~ m/Confirmed/));

	my $chr = $gff_line[0];
	$chr =~ s/CHROMOSOME_//;
	my $strand = $gff_line[6];
	my $cds = $gff_line[8];
	$cds =~ s/CDS \"//;
	$cds =~ s/\" .*//g;
	
	# need to ignore the same intron in different GFF lines
	# use chromosome number and coords as a unique key
	my $key = $chromosome.":".$gff_line[3].$gff_line[4];
	
	if(exists $unique_introns{$key}){
	    # can skip intron
	    next;
	}
	else{
	    # add to hash of unique introns so as not to count it again
	    $unique_introns{$key} = "1";				   

	    # Calculate length and increment counter in respective hash
	    my $length = $gff_line[4] - $gff_line[3]+1;

	    # add to relevant hashes (depending on whether it is confirmed intron or not)
	    if($gff_line[8] =~ m/Confirmed/){
		$confirmed_introns{$length}++;		
	    }
	    else{
		$unconfirmed_introns{$length}++;
	    }		
	    $all_introns{$length}++;

	   
	    # as long as we are not in -count mode, check intron size is in desired range and print details
	    if(!$count && ($length >= $min) && ($length <= $max)){

		# increment intron counter
		$counter++;
		
		# set counter to make new intron names for FASTA header file
		my $name = "${counter}_${cds} ${chr} ${gff_line[3]}-${gff_line[4]} $strand ${length} bp";
		
		# if -flank specified, additionally extract extra bases of flanking exons
		# need to substract 1 from start coordinate to work in array coordinatges
		my $upstream_flank   = $flank+1;
		my $downstream_flank = $flank-1;
		
		my @intron;

		# Just get specified amounts from 5' end of intron if -donor is specified
		if($donor){
		    # want first part of intron sequence as specified by $donor_size
		    if($strand eq "+"){
			@intron = @dna[$gff_line[3]-$upstream_flank..$gff_line[3]+$donor_size-2];
		    }
		    else{
			@intron = @dna[$gff_line[4]-$donor_size..$gff_line[4]+$downstream_flank];
		    }
		}

		# Or just get specified amounts from 3' end of intron if -acceptor is specified
		elsif($acceptor){
		    # want last part of intron sequence as specified by $acceptor_size
		    if($strand eq "+"){
			@intron = @dna[$gff_line[4]-$acceptor_size..$gff_line[4]+$downstream_flank];
		    }
		    else{
			@intron = @dna[$gff_line[3]-$upstream_flank..$gff_line[3]+$acceptor_size-2];
		    }
		}

		# or finally just get entire intron
		else{
		    @intron = @dna[$gff_line[3]-$upstream_flank..$gff_line[4]+$downstream_flank];
		}
		


		my $intron;		
		
		# reverse complement intron sequence if on negative strand
		# if flanking sequence is requested, put in upper case though have to change the behaviour a bit
		# if -donor or -acceptor is also specified as there is only one flanking sequence not two 
		if($strand eq "-"){
		    my @flipped = reverse(@intron);
		    # upper case flanking regions
		    if($flank){
			my $i = 0;
			while($i<$flank){
			    $flipped[$i] = uc($flipped[$i]) unless ($acceptor);
			    $flipped[-$i-1] = uc($flipped[-$i-1]) unless ($donor);
			    $i++;
			}
		    }

		    $intron =  (join "",@flipped);
		    $intron =~ tr/atcg/tagc/;
		    # also change flanking regions if necessary
		    ($intron =~ tr/ATCG/TAGC/) if ($flank);
		}
		else{
		    if($flank){
			my $i = 0;
			while($i<$flank){
			    $intron[$i] = uc($intron[$i]) unless ($acceptor);
			    $intron[-$i-1] = uc($intron[-$i-1]) unless ($donor);
			    $i++;
			}
		    }

		    $intron =  (join "",@intron);
		}
		# change thymine to uracil
#		$intron =~ tr/t/u/;

		# warn if intron may have a repeat sequence in it
		if($intron =~ m/n/){
		    print "WARNING: intron contains a repeat\n";			
		}

		# skip non GC.. introns if -gc option is used
		if($gc){
		    # need to ignore flanking sequences though else sequence will never start with gc
		    my $test_intron = $intron;
		    $test_intron =~ s/[A-Z]//g if($flank);

		    next if ($test_intron !~ /^gc/);
		}

		# print out details
		if($print){
		    print OUT ">$name\n";
		    print OUT "$intron\n" unless (($intron =~ m/n/) && ($ignore_n));
		}
		else{
		    print ">$name\n";
		    print "$intron\n\n" unless (($intron =~ m/n/) && ($ignore_n));
		}
	    }	    
	}
    }	  	    
    
    close(GFF_in);  
    print "\n";       
}

# no need to go any further unless we are counting introns
if (!$count){
    close(OUT) if ($print);
    exit(0);
}
    
#################################################################
# need to work out intron counts for specified bin size
#################################################################
my $bin_start;
my $bin_end;

# hard code intron size...assume that nothing will go over this?
my $max_intron_size = 25000;

for ($bin_start = 1; $bin_start < $max_intron_size; $bin_start += $interval){

    $bin_end = $bin_start+$interval-1;

    # quit if limit exceeded
    if($bin_start > $limit){
	exit(0);
    }
    
    # now work out total of intron in current bin
    my $confirmed_bin_total   = 0;
    my $unconfirmed_bin_total = 0;
    my $all_bin_total         = 0;

    for(my $count = $bin_start; $count <= $bin_end; $count++){

	# can only sum up introns if they exist at that size
	$confirmed_bin_total   += $confirmed_introns{$count}   if ($confirmed_introns{$count});
	$unconfirmed_bin_total += $unconfirmed_introns{$count} if ($unconfirmed_introns{$count});
	$all_bin_total         += $all_introns{$count}         if ($all_introns{$count});
    }
    print "${bin_start}-${bin_end}\t$all_bin_total\t$confirmed_bin_total\t$unconfirmed_bin_total\n";

}



#######################################################################################
#
# Simple subroutine to fetch a chromosome sequence (as DNA) and store it as a string
#
#######################################################################################

sub fetch_seq{

    my $chromosome = shift;
    # open sequence file (use masked sequences to help see if repeat sequence is in intron)
#    open (DNA, "<$gffdir/CHROMOSOME_${chromosome}_masked.dna") || die "Failed to open dna file\n\n";
    open (DNA, "<$gffdir/CHROMOSOME_${chromosome}.dna") || die "Failed to open dna file\n\n";

    my $seq;
    while(my $tmp =<DNA>){
	chomp($tmp);
	# skip header line
	next if ($tmp =~ m/^>/);
	$seq .= $tmp;
    }
    close(DNA);
    return($seq);
}
