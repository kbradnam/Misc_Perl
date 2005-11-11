#!/usr/bin/perl -w

use strict;

# species name is key, list of intron sizes are values
my %species;

# read input file
open(IN, "<$ARGV[0]") || die "Can't open input file\n";
while(<IN>){
    chomp;
    # get rid of first characte
    s/>//;
    # skip lines with no intron data (chloroplast etc.)
    next if($_ !~ m/(\S+) (\S+) .*intronLen=([0-9,]+)/);
    m/(\S+) (\S+) .*intronLen=([0-9,]+) /;
    my $name = ($1." ".$2);
 #   print "$name\n";

    my $intron = $3;
    my @introns = split(/,/,$intron);
    while(my $int = pop(@introns)){
	push(@{$species{$name}},$int);
    }
}

# new hash, species name is key, value is array of intron size counts (position 1 = size 1 etc.)
my %species2size;

# loop through each species to order intron sizes into new hash
foreach my $key (keys(%species)){    
    # count total number of introns for each species
    my $count =@{$species{$key}};
    if ($count > 250){
	foreach my $size (@{$species{$key}}){
	    ${$species2size{$key}}[$size]++;
	}
    }
}

foreach my $key (keys(%species2size)){
    print "$key,";
    for(my $i =10; $i <150; $i+=3){
	#print "i = $i\n";
	my $count=0;
	if(defined(${$species2size{$key}}[$i])){
	    $count += "${$species2size{$key}}[$i]";
	}
	if(defined(${$species2size{$key}}[$i+1])){
	    $count += "${$species2size{$key}}[$i+1]";
	}
	if(defined(${$species2size{$key}}[$i+2])){
	    $count += "${$species2size{$key}}[$i+2]";
	}
	print "$count,";
    }
    print "\n";
}
exit(0);

