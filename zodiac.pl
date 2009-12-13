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

# will want a array of colors to substitute any @color@ pattern
my @colors = qw(red orange yellow green blue purple white black pink);

# load up an array of famous people
my @people = load_people();

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
		#if we have any suitable text, then first check to see if it includes @color@
		if($line =~ m/\@color\@/){
			my $rand = int(rand(9));
			$line =~ s/\@color\@/$colors[$rand]/;
		}
		# also add in random numbers for each occasion of @number@
		if($line =~ m/\@number\@/){
			my $number = int(rand(100));
			$line =~ s/\@number\@/$number/;
		}
		# also add in random numbers for each occasion of @person@
		if($line =~ m/\@person\@/){
			my $person = $people[rand(@people)];
			$line =~ s/\@person\@/$person/;
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

	# add hyperlinks to any species names in resulting paragraph
	my $text = check_names("$sentence1 $sentence2 $sentence3");
	print "<P>$text</P>\n";
}

exit;

sub check_names{
	my ($text) = @_;
	foreach my $species qw(Barnacle Snail Limpet Clam Squid Slug Oyster Scallop Octopus Mussel){

		if($text =~ m/$species/){
			my $url_name = lc($species);
			my $new_text = "<a href=\"${url}${url_name}.html\">$species</a>";
			$text =~ s/$species/$new_text/g;
		}
	}
	return($text);
}
sub load_people{
	my @list = (
"David Beckham",
"Tiger Woods",
"Miley Cyrus",
"Bill Cosby",
"Bill Gates",
"George Clooney",
"Margaret Thatcher",
"Thomas Edison",
"Mother Teresa",
"Helen Keller",
"Madonna",
"Jacqueline Kennedy Onasis",
"Tom Brokaw",
"James Taylor",
"Mr. Rogers",
"Isaac Newton",
"Lewis Carrol",
"Andy Rooney",
"General Norman Schwarzkopf",
"Norman Rockwell",
"Pablo Piccaso",
"Paul McCartney",
"Plato",
"Edgar Allen Poe",
"Mae West",
"Ernest Hemingway",
"Vincent Van Gogh",
"W.C.Fields",
"Robin Williams",
"Walt Disney",
"Walter Cronkite",
"William Shakespeare",
"Frank Lloyd Wright",
"Julia Roberts",
"John F. Kennedy, Jr.",
"Terry Bradshaw",
"Gloria Steinem",
"Charles Dickens",
"Thomas Edison",
"Whoopi Goldberg",
"Sigourney Weaver",
"Bill Clinton",
"Dave Letterman",
"Newt Gingrich",
"Jim Carrey",
"Mary Tyler Moore",
"Danny Glover",
"Carol Burnett",
"Paul Harvey",
"Alicia Silverstone",
"Neil Diamond",
"Julia Child",
"George Carlin",
"Valerie Harper",
"John Candy",
"Weird Al Yankovick",
"Marilyn Vos Savant",
"Tom Hanks",
"C. G. Jung",
"William James",
"Henri Mancini",
"Bob Newhart",
"Meryl Streep",
"Benny Goodman",
"Harrison Ford",
"Steve Martin",
"Ronald Regan",
"Dan Aykroyd",
"Susan B. Anthony",
"Arthur Ashe",
"Augustus Caesar",
"Jane Austen",
"William F. Buckley, Jr.",
"Chevy Chase",
"Phil Donahue",
"Peter Jennings",
"Charles Everett Koop",
"C. S. Lewis",
"Roy Rogers",
"Chuck Yeager",
"Jack Nicholson",
"Charlie Brown",
"Oprah Winfrey",
"Paul Newman",
"Pelé",
"Fred Astaire",
"Eddie Murphy",
"Jimmy Conners",
"Michael J. Fox",
"Ross Perot",
"Sean Connery",
"Elizabeth Dole",
"Dick Van Dyke",
"Andy Griffith",
"Peyton Manning",
"Nathaniel Hawthorne",
"Shirley MacLaine",
"Michael Landon",
"John Katz",
"Billy Crystal",
"Carrie Fisher",
"Darth Vader",
"Bob Dylan",
"Carl Sagan",
"Charles Yeager",
"Colin L. Powell",
"Henry A. Kissinger",
"Elvis Presley",
"Madonna",
"Mahatma Gandhi",
"Michael J. Jordan",
"Michele Pfeiffer",
"Doris Day",
"Liberace",
"Elizabeth Taylor",
"Yogi Berra",
"Dan Rather",
"Magic Johnson",
"Michael Jackson",
"John Travolta",
"Tom Cruise",
"Spider Man",
"James Dean",
"Clint Eastwood",
"Ray Charles",
"Jesse Jackson",
"Thomas Jefferson",
"Hank Aaron",
"Mohammad Ali",
"Aristotle",
"Neil Armstrong",
"Lucille Ball",
"Hank Aaron",
"Beethoven",
"Alexander Graham Bell",
"Napoleon",
"George Washington",
"Cleopatra",
"Columbus",
"Dr. Seuss",
"Albert Einstein",
"Eisenhower",
"F Lee Bailey",
"Ben Franklin",
"Sigmund Freud",
"Gandhi",
"Alfred Hitchcock",
"Bob Hope",
"Harry Houdini",
"Martin Luther King",
"John Lennon",
"Leonardo Da Vinci",
"Abraham Lincoln",
"Louis Pasteur",
"Marilyn Monroe",
"Mark Twain",
"Willey Mays",
"Michelangelo",
"Miles Davis",
"Mozart"
);
	return(@list);
}
__END__