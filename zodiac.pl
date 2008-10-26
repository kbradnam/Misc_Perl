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

# will want a array of days to substitute any @day@ pattern
my @days = qw(Monday Tuesday Wednesday Thursday Friday Saturday Sunday);


my ($species,$date);

while(my $line = <>){
	chomp($line);
	
	# skip blank lines
	next if ($line =~ m/^$/);
	
	if ($line =~ m/^#/){
		($species,$date) = split(/,/,$line);
		$species =~ s/# //;
		($species2dates{$species} = $date) unless ($line =~ m/^# all/);
	}
	else{
		# if we have any suitable text, then first check to see if it includes @day@
		if($line =~ m/\@day\@/){
			my $rand = int(rand(7));
			$line =~ s/\@day\@/$days[$rand]/;
		}
		# also add in random numbers for each occasion of @number@
		if($line =~ m/\@number\@/){
			my $number = int(rand(100));
			$line =~ s/\@number\@/$number/;
		}
		
		push(@{$species2text{$species}},$line);
	}
}

# now loop through each of the 10 species extracting two random sentences per species

my $url = 'http://molluskanzodiac.blogspot.com/2008/09/';

SPECIES: foreach my $species qw(Barnacle Snail Limpet Clam Squid Slug Oyster Scallop Octopus Mussel){

	my $url_name = lc($species);

	print "<a href=\"${url}${url_name}.html\">The $species</a>\n";
	print "<b><i>$species2dates{$species}</b></i>\n";

	my ($rand1,$rand2,$rand3,$sentence1,$sentence2,$sentence3);

	# if we are dealing with the slug, then we do things differently
	if($species eq "Slug"){
		$rand1 = int(rand(1) * @{$species2text{"$species"}});
		$sentence1 = ${$species2text{"$species"}}[$rand1];
		splice(@{$species2text{"$species"}},$rand1,1);
		
		$rand2 = int(rand(1) * @{$species2text{"$species"}});
		$sentence2 = ${$species2text{"$species"}}[$rand2];
		splice(@{$species2text{"$species"}},$rand2,1);

		$rand3 = int(rand(1) * @{$species2text{"$species"}});
		$sentence3 = ${$species2text{"$species"}}[$rand3];
		splice(@{$species2text{"$species"}},$rand3,1);

		print "<P>$sentence1 $sentence2 $sentence3</P>\n";

		next SPECIES;
	}
	
	$rand1 = int(rand(1) * @{$species2text{"all"}});
	$sentence1 = ${$species2text{"all"}}[$rand1];
	splice(@{$species2text{"all"}},$rand1,1);
	
	$rand2 = int(rand(1) * @{$species2text{"all"}});
	$sentence2 = ${$species2text{"all"}}[$rand2];
	splice(@{$species2text{"all"}},$rand2,1);

	# 3rd sentence is species specific
	# only one sentence taken per species, so need to splice array afterwards
	$rand3 = int(rand(1) * @{$species2text{"$species"}});
	$sentence3 = ${$species2text{"$species"}}[$rand3];

	print "<P>$sentence1 $sentence2 $sentence3</P>\n";
}


exit;


__END__