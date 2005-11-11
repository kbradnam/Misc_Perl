#!/usr/bin/perl -w

use strict;

# species name is key, list of intron or exon sizes are values
my %species2intron;
my %species2exon;

# read input file
open(IN, "<$ARGV[0]") || die "Can't open input file\n";
while(<IN>){

    # only want header lines
    next unless m/^>/;

    # strip line to start with species name
    s/.*?; +//;

    # Fix some obvious species name errors/contractions
    s/A\.thaliana/Arabidopsis thaliana/;
    s/C\.elegans/Caenorhabditis elegans/;
    s/C\. elegans/Caenorhabditis elegans/;
    s/D\. melanogaster/Drosophila melanogaster/;
    s/D\.melanogaster/Drosophila melanogaster/;
    s/D\.simulans/Drosophila simulans/;
    s/D\.pseudoobscura/Drosophila pseudoobscura/;
    s/D\.discoideum/Dictyostelium discoideum/;
    s/Fugu rubripes/Takifugu rubripes/;
    s/Human DNA /Homo sapiens /;
    s/H\.sapiens/Homo sapiens /;
    s/M\.musculus/Mus musculus/;
    s/R\.norvegicus/Rattus norvegicus/;
    s/S\.cerevisiae/Saccharomyces cerevisiae/;
    s/S\.pombe/Schizosaccharomyces pombe/;
    s/Z\.mays/Zea mays/;
    s/Zebrafish DNA /Danio rerio /;
    s/sativa,/sativa/;

    
    # ignore cases where species name is truncated, e.g. A.astacus
    next unless m/^[A-Za-z]+ [A-Za-z]/;

    # Remove the 'Genomic sequence for' and 'Sequence of BAC ... from' text which messes up
    #some rice and Arabidopsis entries
    s/Genomic sequence for //i;
    s/Sequence of BAC .*? from //;

    # cut out gene descriptions
    # remove phase information and intron/exon sum stats
    # remove splice info
    s/(\S+) (\S+) .*?; (.*) /$1 $2 $3/;
    s/phase:[012]+,//;
    s/intr_sum:[-0-9]+//;
    s/ex_sum:[-0-9]+//;
    s/[; ]+\{splice.*//;
    s/size://g;

    # match up species name and exon/intron coordinates
    m/(\S+ \S+) intron\((.*?)\); exon\((.*?)\)/; 
    
    my $name   = $1;
    my $intron = $2;
    my $exon   = $3;

#    print "*$name* *$intron* *$exon*\n";

    # cycle through introns and add to hash
    my @introns = split(/,/,$intron);
    while(my $int = pop(@introns)){
	#skip unclassified introns
	next if ($int eq "u");
	push(@{$species2intron{$name}},$int);
    }

    # cycle through exons and add to hash
    my @exons = split(/,/,$exon);
    while(my $ex = pop(@exons)){
	push(@{$species2exon{$name}},$ex);
    }
}



# new hash, species name is key, value is array of intron size counts (position 1 = size 1 etc.)
my %species2intronsize;

# loop through each species to order intron sizes into new hash
foreach my $key (keys(%species2intron)){    
    # count total number of introns for each species
    my $count =@{$species2intron{$key}};
    #only want species with at least 250 introns
    if ($count > 250){
	# increment counter of introns of size N by using incrementing integer at position N in array
	foreach my $size (@{$species2intron{$key}}){
	    ${$species2intronsize{$key}}[$size]++;
	}
    }
}

# Only look at introns between a minimum and maximum size
my $min_intron = 15;
my $max_intron = 200;

# print header for csv file
print "Species,";
my $start = $min_intron;
while($start < $max_intron){
    my $end = $start+2;
    print "$start-$end,";   
    $start+=3;
}
print "\n";


foreach my $key (keys(%species2intronsize)){
    print "$key,";
    for(my $i =$min_intron; $i <$max_intron; $i+=3){
	my $count=0;
	# count introns in n..n+2 size intervals
	if(defined(${$species2intronsize{$key}}[$i])){
	    $count += "${$species2intronsize{$key}}[$i]";
	}
	if(defined(${$species2intronsize{$key}}[$i+1])){
	    $count += "${$species2intronsize{$key}}[$i+1]";
	}
	if(defined(${$species2intronsize{$key}}[$i+2])){
	    $count += "${$species2intronsize{$key}}[$i+2]";
	}
	print "$count,";
    }
    print "\n";
}
exit(0);

