#!/usr/bin/perl -w
#
# find_transposons_in_introns.pl
#
# by Keith Bradnam
# 
## created 18th November  2005
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

my $release;     # which WormBase release to use (will default to latest)                 

GetOptions ("release=i"   => \$release);

# which WormBase release to use
if($release && $release !~ m/\d{3}/){
    die "please specify WormBase release as 'NNN' (where NNN is an integer)\n";
}
# set a default release   
$release = "150" if (!$release);


#############
# Paths etc #
#############

my $tace        = "/Korflab/bin/tace";                           # tace executable path
my $dbdir       = "/Korflab/Data_sources/WormBase/WS${release}"; # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                          # GFF splits directory
my @chromosomes = qw( I II III IV V X );                         # chromosomes to parse 

print "Using $dbdir for database directory\n";
print "Using $tace for tace binary\n";



##############################################################
# get transpson and Coding_transcript intron spans from GFF files
##############################################################

# keep count of transposons in introns
my $count = 0; 

foreach my $chromosome (@chromosomes) {
    print "\n# Processing Chromosome $chromosome\n";

    my @transposons;     # transposon name
    my %transposon_count; # count of transposons of each name
    my %transposon_span;  # store start and stop coordinates
    my %intron_span;# span of each intron
	
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
	chomp;
	
	# skip header info and then split line into array
	next if m/^\#/;
	my @gff_line = split /\t/;

	
	# want to match all transposons and transposon CDSs (both have the same GFF method field)
	if($gff_line[2] eq "transposable_element"){

	    # get transposon name from last gff column
	    my $name = $gff_line[8];

	    # Have to remove one of two things depending if it is a Transposon or Transposon_CDS
	    $name =~ s/(CDS|Transposon) //;
	    $name =~ s/\"//g;	    
	    
	    # load up details to hashes (need to store start coord, stop coord, and strand)
	    $transposon_count{$name}++;
	    $transposon_span{$name} = [$gff_line[3],$gff_line[4],$gff_line[6]];
	    
	    # Check that transposon names don't appear more than once in each GFF file
	    if($transposon_count{$name}>1){
		die "ERROR: $name appears $transposon_count{$name} times\n"; 
	    }
	    else{
		# load up array with all transposons on that chromosome
		push(@transposons,$name);
	    }	    
	}

	# now also want to match introns from coding transcripts
	elsif(($gff_line[1] eq "Coding_transcript") && ($gff_line[2] eq "intron")){
	    # get intron name from last gff column
	    my $transcript = $gff_line[8];
	    $transcript =~ s/Transcript \"([\w\d\.]+)\".*/$1/;

	    # use intron coords as name in hash, this means multiple transcript objects with the same intron
	    # will just overwrite each other and this will remove redundancy	   
	    my $name = "$gff_line[3]-$gff_line[4]";
	    # store coords with name and strand in hash array
	    $intron_span{$name} = [$gff_line[3],$gff_line[4],$transcript,$gff_line[6]];
	}
	else{
	    # skip everything else
	    next;
	}
    }

    close(GFF_in);

    
    # now cycle through introns looking for those containing a transposon
    foreach my $intron (keys %intron_span){
	
	my $intron_start  = $intron_span{$intron}->[0];
	my $intron_end    = $intron_span{$intron}->[1];
	my $intron_name   = $intron_span{$intron}->[2];
	my $intron_strand = $intron_span{$intron}->[3];
	my $intron_length = $intron_end - $intron_start +1;

	foreach my $new_transposon (@transposons){

	    
	    my $transposon_start  = $transposon_span{$new_transposon}->[0];
	    my $transposon_end    = $transposon_span{$new_transposon}->[1];
	    my $transposon_strand = $transposon_span{$new_transposon}->[2];
	    my $transposon_length = $transposon_end - $transposon_start +1;

	    # only looking for transposons wholly contained within other genes
	    if (($transposon_start >= $intron_start) && ($transposon_end <= $intron_end)){
		$count++;
#		my $out1  = sprintf("%4d %1s %12s %8d %8d %6d",$count,$chromosome,$intron,$intron_start,$intron_end,$intron_length);
#		my $out2 .= sprintf("%12s %8d %8d %6d",$new_transposon,$transposon_start,$transposon_end,$transposon_length);
#		print "$out1  ->  $out2\n";

		# print out in CSV format
		print "$count,$chromosome,$intron_name,$intron_start,$intron_end,$intron_length,$intron_strand,";
		print "$new_transposon,$transposon_start,$transposon_end,$transposon_length,$transposon_strand\n";

	    }
	    
	}
    }
  
}

