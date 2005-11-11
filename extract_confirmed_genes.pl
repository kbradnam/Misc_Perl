#!/usr/local/bin/perl5.6.0 -w

use Ace;

my $db = Ace->connect(-path  => "/wormsrv2/autoace/");

my $query = "Find Sequence; Confirmed_by";
my @confirmed_genes   = $db->fetch(-query=>$query);

open(OUT,">confirmed_genes.dna") || die "Couldn't write to output file\n";

foreach my $seq (@confirmed_genes){
	my $dna = $seq->asDNA();

	my (@type) = $seq->get('Confirmed_by');

	if(defined($type[1])){
		$dna =~ s/(>\w+\.\w+)/$1 Confirmed_by_EST_and_cDNA/;	
	}	
	else{
		$dna =~ s/(>\w+\.\w+)/$1 Confirmed_by_$type[0]/;	
	}	
	print OUT "$dna";
}

close(OUT);

exit(0);


