#!/usr/bin/perl -w

use strict;

# three hashes for each of the three ontologies
my %molecular;
my %cellular;
my %biological;

# read file

while(<>){
    chomp;
    my @line = split(/\t/);
    # don't want gene names from first column
    if(defined($line[1]) && defined($line[2])){
	# increment counters in three hashes
	$molecular{$line[1]}++  if($line[2] =~ m/Molecular_function/);
	$cellular{$line[1]}++   if($line[2] =~ m/Cellular_component/);
	$biological{$line[1]}++ if($line[2] =~ m/Biological_process/);   
    }
}

# Now print out each count (sort hash numerically by value)

foreach my $key (reverse sort {$molecular{$a} <=> $molecular{$b}} keys %molecular){
    print "Molecular_function,$key,$molecular{$key}\n";
}

foreach my $key (reverse sort {$cellular{$a} <=> $cellular{$b}} keys %cellular){
    print "Cellular_component,$key,$cellular{$key}\n";
}

foreach my $key (reverse sort {$biological{$a} <=> $biological{$b}} keys %biological){
    print "Biological_process,$key,$biological{$key}\n";
}
