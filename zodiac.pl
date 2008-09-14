#!/usr/bin/perl
#
# zodiac.pl 
#
# a script to generate molluskan horoscopes from prewritten sentences
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;


# first read data file of pre-written sentences
# will need to store species specific text, plus text specific to any of the ten species
# also another hash will store the dates for each sign

my %species2text;
my %species2dates;

my ($species,$date);

while(my $line = <>){
	chomp($line);
	if ($line =~ m/^#/){
		($species,$date) = split(/,/,$line);
		$species =~ s/# //;
		($species2dates{$species} = $date) unless ($line =~ m/^# all/);
	}
	else{
		push(@{$species2text{$species}},$line);
	}
}


# now loop through each of the 10 species extracting two random sentences per species

my $url = 'http://molluskanzodiac.blogspot.com/2008/09/';

SPECIES: foreach my $species qw(Barnacle Snail Limpet Clam Squid Slug Oyster Scallop Octopus Mussel){

	my $url_name = lc($species);

	print "<a href=\"${url}${url_name}.html\">The $species</a>\n";
	print "<b><i>$species2dates{$species}</b></i>\n";

	my ($rand1,$rand2,$sentence1,$sentence2);

	# if we are dealing with the slug, then we do things differently
	if($species eq "Slug"){
		$rand1 = int(rand(1) * @{$species2text{"$species"}});
		$sentence1 = ${$species2text{"$species"}}[$rand1];
		splice(@{$species2text{"$species"}},$rand1,1);
		
		$rand2 = int(rand(1) * @{$species2text{"$species"}});
		$sentence2 = ${$species2text{"$species"}}[$rand2];
		splice(@{$species2text{"$species"}},$rand2,1);

		print "<P>$sentence1 $sentence2</P>\n";

		next SPECIES;
	}
	
	$rand1 = int(rand(1) * @{$species2text{"all"}});
	$sentence1 = ${$species2text{"all"}}[$rand1];
	splice(@{$species2text{"all"}},$rand1,1);
	
	$rand2 = int(rand(1) * @{$species2text{"all"}});
	$sentence2 = ${$species2text{"all"}}[$rand2];
	splice(@{$species2text{"all"}},$rand2,1);

	print "<P>$sentence1 $sentence2</P>\n";
	
}


exit;


__END__