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

while(my $line = <DATA>){

	chomp($line);
	
	# skip blank lines
	next if ($line =~ m/^$/);
	
	if ($line =~ m/^#/){
		($species,$date) = split(/,/,$line);
		$species =~ s/# //;
		($species2dates{$species} = $date) unless ($line =~ m/^# all/);
		next;
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

SPECIES: foreach my $species (qw(Barnacle Snail Limpet Clam Squid Slug Oyster Scallop Octopus Mussel)){

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

		print "<P>$sentence1 $sentence2 $sentence3</P><br>\n";

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
	print "<P>$text</P><br>\n";
}

exit;

sub check_names{
	my ($text) = @_;
	foreach my $species (qw(Barnacle Snail Limpet Clam Squid Slug Oyster Scallop Octopus Mussel)){

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
"Mozart",
"Herman Cain",
"Rick Perry",
"Michelle Bachman",
"Rick Santorum",
"Newt Gingrich",
"Mitt Romney"
);
	return(@list);
}

__DATA__
# all,
Dietary choices may be important this week. Consider avoiding foods that are green or yellow in color. 
Now is a good time to live life in the fast lane and be daring and bold. Try wearing one less item of clothing than usual. 
Making more room for music this week will ease current frustrations. The more avant-garde the music the better, and 50's Jazz will particularly prove helpful.
The current problems with your job are partially due to the monotony that surrounds you. Kick start your career by doing something to shock your colleagues and show them your inner beast.
You need to get away from the stress that is currently crushing your spirit. A weekend in an isolation chamber will help you focus.
The numbers 6 or 9 (or possibly 28) hold the key to happiness this week.
Avoid even numbers this week if possible as they will only bring you trouble.
A special number looms large in your life this week, and that number is @number@.
Numbers will prove problematic for you this week. It could be a birthday or other important date, or it could be the lottery. Step wisely when choosing any number.
A man connected with the number @number@ will bring you joy and a woman connected with the number -4 will bring you great sorrow.
Travel this week will broaden the mind, but flatten the wallet.
Ignite your enthusiasm this week by planning a trip, but avoid buses, trains, and planes as these might bring you danger.
A forthcoming trip is causing you much stress, Speak to your doctor for advice.
Don't leave your house on @day@ afternoon, disaster lurks if you step outside.
Take a loved one on a surprise vacation this week and they will be very appreciative, especially if you take them overseas.
Your weight, or the weight of someone important around you, may become a discussion point this week.
You will be attracted to someone in a position of power this week. Do not give in to temptation, make sure they get the cold shoulder.
This week will see you face many important questions. It is important that you answer those questions.
Doubt and uncertainty will cloud your thoughts this week. Try listening to a tall person for advice.
You want what you cannot have. You have what you no longer want. Such is life.
A pet or other animal that is close to you will cause problems this week. Avoid all animals if possible.
If there was ever one week in your life where you should eat cheese, this week is it.
Spend more time not talking to people and your silence will be rewarded.
Think about what you want, and what you need. Are they the same things?
You will breathe more deeply this week when things that you want appear on the horizon, in a shimmering cloud.
Be wise on Thursday, but timid on Friday.
@day@ holds a great surprise for you, unless you already know about it.
You will be troubled by two legs on Tuesday and four legs on Friday.
This week you will be tired. Tired of work. Tired of family and friends. Tired of life. However, you will sleep soundly.
Your friends are being spoons, when all you need is a fork. But being ironic doesn't befit you.
Why do you do what you think you should do when you don't do what you don't think that you should do?
Don't do things that you don't want to do, unless you do want to do the things that you think that you don't want to do.
Are you a lion or a mouse? Now is the time to nail your colors to the flag and decide whether you can squeak or roar.
The rabbit that runs twice as fast, eats twice as slow. Don't be the rabbit that wins a race but ends up hungry.
Clouds are on the horizon. Storm clouds. Storm clouds that will bring rain, hail, thunder, and lightning. Run away.
Even a tiny fly can stop a bullet, if its wings are made of steel. However, your wings are more likely to be made of damp cardboard, which might not be so effective.
There is no difference between what you can do and what you think you can do. The only difference is in your mind, or what you think is in your mind.
Walk briskly this week, because the winds of opposition will try to slow you down. Break through the winds to achieve success, only by breaking wind will you find happiness.
Is there a door opening up in your life? If so then close it, or at most leave it only slightly ajar.
There may be times this week when you will wonder if you will ever make it to Friday unscathed. If you survive until Thursday evening, then you will be fine.
This will be a good week for fun, festivities, and especially fashion. Let your clothes do the talking and don't skimp on the accessories.
Food may be important for you this week, so it might be preferable if you try eating it.
Eat well this week, but don't eat too little, and don't eat too much. Also avoid the wrong types of food and focus on the right types of food.
People will tell you that life can have its ups and downs, but they never tell you to watch out for the sideways.
Avoid photographs this week if you think that your illicit affair may be caught on camera.
This is a great week for trying something completely new such as listening to jazz, ballroom dancing, or invading a neighboring country.
Try to look forward to the future this week, but still keep one eye looking over your shoulder as the past may catch you up and spit in your eye.
Do you want to feel like crap every morning? If the answer is no, then try eating walnuts before bedtime.
It's no use gazing at the stars if your feet are stuck in the mud. Clean your boots and get your life moving forward again.
Up, up, up, up, up, up! That's the direction your life will be heading in this week (terms and conditions may apply).
Think of all the great things that might happen to you this week. They may never happen, but at least you're thinking about them.
If you were a vegetable, you'd probably be a tomato. Watch that you don't get squashed this week.
If there was ever a week to shut the curtains, stay in bed, and hide under the duvet like a frightened kitten, then this is the week...unless you need to go out.
Happy. Happy. Happy. Happy. Sad. Sad. Sad. Sad. Which one of these will you be this week?
You can dance this week if you think that kicking up your heels will make you happier. By the way, it won't.
Love is all around you this week. You will feel it in your fingers. You will feel it in your toes.
You've set your sights high this week, but as the Chinese proverb warns us 'rain always dampens an egg buried in the ground'.
There are many things that you would like to try this week. But remember 'do or do not, there is no try'.
Confront your inner demons this week and arrange for an internal exorcism. 
Are you happy? Are you sad? Are you content? Are you restless? The answers to 3 of these questions will not be revealed this week.
You think that someone is out to get you, you think that they want to see you squashed like a bug. You are wrong. They are wrong. Everyone is wrong.
Smile like a bumblebee in June, and you will be rewarded for your happiness.
A hairy man (or woman) will provide you with a bristly problem this week.
Put some distance between you and a rival. At least @number@ feet, but no more than a mile.
Your enemies are plotting against you. Ignore them, what's the worse that could happen?
Eat well, sleep well, and make sure you put the cat out because you will need a lot of energy to get through this week.
Allergies might prove bothersome this week, especially if you work on a farm or are allergic to milk.
This is the week where you will wish that you could be as slippery as an eel in a mud-wrestling contest.
@day@ afternoon (about 3:15) is the time for making a big decision about your life.
If you keep putting it off (and you know what I mean by 'it'), it will never get done. Sort it out this week once and for all.
A trip to the dry-cleaners could provide the impetus you have been looking for to kick-start your business plans.
Why do people infuriate so much? Could it be because they are all idiots? Probably.
Nobody seems to recognize your genius. You are a jumbo shrimp in a sea of Clams.
Others will spend this week trying to think outside the box. Show them your true genius by turning the box inside out and then thinking inside it.
You work hard but seem to get no reward for your effort. Perhaps this week, you will get effort for your reward.
Be careful not to overexert yourself in the kitchen this week. Remember, too many broths spoil the cook.
If there was ever a week in which you should enroll in a foreign language class, then this is the week.
Tiredness will knock on your door this week, so be prepared to consume vast amounts of energy drinks.
Is there a ray of light at the end of the tunnel? A chance meeting on @day@ with a gynecologist might provide some answers.
You are being driven mad by driving. Don't get mad, get even.
If you ever wanted to place a bet on a big race, then this is the week. A horse whose name begins with the letter G will win big.
Take a second look at what you are wearing. Your friends think that it is time that you burn your wardrobe. Maybe they are right?
You are a genius, only no-one knows it. Maybe you should try telling people.
A foreign fish will play an exciting role in your life this week.
Big developments will occur in the bedroom this week. Make sure your sheets are clean.
Eggs, or products containing eggs, are best avoided this week. Unless you are certain that they are what you want.
Laughter will fill the air this week. But will it be yours? There is only one way to be sure. Rent a good comedy on DVD and watch with a friend.
Train yourself to be mentally stronger and reap the rewards. Especially on @day@ when a chance meeting with a handsome stranger will allow you to think outside the box.
Don't take no for an answer, especially when rancid dairy products are involved.
You have been thinking recently, 'is this the best I can do?'. The answer, sadly, is 'yes'.
Why do you spend so much time waiting for other people to tell you how great you are. Cut out the middle man and start singing your own praises while looking in the mirror.
Life is good at the moment, so be careful not to ruin it all becoming addicted to gambling.
Fish are a big thing in your life at the moment. Catching fish and eating fish are what you are all about.
You are a fighter, not a quitter. Don't let the bastards grind you down.
Why is everyone so keen on cheese these days? You know that steering clear of the yellow stuff is the right thing to do.
A friend in need is a friend indeed...except when they cheat on you behind your back. Keep a careful eye out on those that call themselves your 'friends'.
Ever had to take over the controls of a plane due to an injury to the pilot? This week might provide an occasion to do just that.
Running away from things will not help problems this week. Neither will staying where you are.
Happiness. Happiness. Happiness. Happiness. Happiness. It's a happy week!
Would you accept a taxi ride if the driver was a monkey? Probably not. So be careful of simian chauffeurs this week.
You've always wanted to try drinking a pint of raw eggs...now is the time to try.
Walk faster than the person in front of you if you want to get ahead this week...unless that person is carrying a knife.
Time to remove the 'us' in fuss and put the 'me' in 'medicate'.
If you stayed in bed <b>all</b> week...maybe nobody would know that you were missing.
Take heed of the old sailors warning 'If you drown, you die'.
Ancient mariners used to say that spotting a whale traveling westwards on a Wednesday, meant that you would suffer a bodily discharge on Thursday. Heed these words.
A sailor that can't sail is not a sailor. Likewise a thinker that can't think is not a thinker. Are you a sailor or a thinker?
It could be a good time this week to heed the warning 'Clams, fireworks, and alcohol do not mix well'.
A famous sailor once said "You can kiss a mermaid, but you might still die of scurvy"...these words will have special significance for you this week.
Are olives really 'the Devil's grape'? This is the week where you will find out.
As the old saying goes 'You can hide a shrimp under a shell, but it's still a shrimp, just a shrimp under a shell'. Heed these words this week, especially if you have any run-ins with the police.
If you smoke, then this is a good week to give up. If you don't smoke, then maybe this is a good week to try.
A famous fisherman once noted that while five fish will always feed a family of four, four fish might not feed a family of five. These words will have special meaning for you this week.
Your friends will tell you that you have to make up your mind regarding your big problem. They will tell you that you must sink or swim. Remember though,that there is a third option. Try to achieve a state of neutral buoyancy.
Belief is the key to your problems this week. Belief in the power of a burning flame. Belief in the strength that can only come from catching three green lights in a row. Belief in the proverb that 'Even a lost penguin will find its way home'. It's time to believe.
Paperwork, paperwork, paperwork. The more you finish, the more just keeps piling up on your desk. The solution to your office stress is to buy a box of matches...the rest will become obvious.
A famous sailor once remarked that 'A beached whale is like a boy urinating in a church at a wedding. It doesn't look good, it doesn't smell good, and everyone pretends not to notice, even though they are secretly annoyed. Don't be that beached whale.
A religious fanatic with a speech impediment will cause you much grief this week.
An important financial decision could prove disastrous if you fail to properly understand the intricacies of global macro-economics. Enhance your chances of success by relying on the time-tested tradition of flipping a coin. Heads means 'Buy' and tails means 'Sell'.
Sexual tensions will be further inflamed this week by an inappropriate use of office stationary.
This might be the sort of crazy week where you should try to do the exact <b>opposite</b> of what everyone tries to tell you to do. One exception to this would be if anyone tells you to do the opposite of what you would normally do.
Make sure that you heed the old maritime warning this week: "When whales swim in threes, flatulence comes for thee".
One of the following objects will potentially cause you to have a life-altering event this week: a red car, an unripe avocado, @day@'s edition of your local newspaper, or a vibrating electronic device.
Mishearing the phrase 'Mass perturbation' will prove your undoing this week. Have your excuses at the ready.
It may not make much sense now, but carry an opened umbrella with you on Wednesday (whatever the weather) and you will be thankful that you did.
A man connected with the number @number@ will potentially have a healing effect on your 'little problem' that's been bothering you.
The number @number@ will be a powerful omen for you this @day@, but only if you are in the possession of some dried fruit.
Should you find yourself in a casino this week, then the number @number@ might be the key to a little financial surprise (the surprise might involve the words 'your credit card is no longer valid' so be careful).
A pretty young woman connected to the number @number@ will be involved in a bizarre gardening accident this week. You may or may not know this woman.
Avoid the number @number@ if possible on @day@. The reason for this is unclear, but as a wise sailor once said "I don't mind being swallowed by a whale...as long as I pass out the other end".
Try to spend one day this week in silence. Communicate only with gestures or bodily odors.
Eat, drink, and be merry. But not if you are driving or are a recovering alcoholic. In which case you should just eat.
The number @number@ will have special significance on @day@, but sadly you will never realize just what that significance is, and so it will all be a bit wasted on you.
If you have a pet llama, then you should try to avoid wearing red on @day@. This may seem an unreasonable request, but you really want the violent and bloody death of an innocent llama on your conscience?
Take the time to make some sense of what you want to say. And cast your words away upon the waves. 
If you wear too much make-up on @day@, you could be in for trouble when someone close to you mistakes you for someone even closer to them.
Your watch is making you a prisoner to time. Destroy it. Break it. Smash it up. Be free from the restrictions of a time-delimited schedule...unless you have an important meeting this week.
They say that you should never comment on a woman's age. Maybe you should try to do it this week to see if that saying still holds true.
Avoid cheesecake at all costs this week, except on @day@ where a small slice of cheesecake will be tolerable (but not if it contains unripened fruit).
A CD will be released this week, a CD that you have been waiting a long time to see. You must <b>never</b> buy this CD. If you buy it, you will become more unpopular than you can possibly imagine.
See a penny, pick it up, and all day long you'll have good luck...or will you???
It is written that 'a drunken sailor is a happy sailor', yet it is also written that 'drinking leads to death'. Which one of these sayings do you most believe in?
The ship of your dreams is sailing down the river of despair. It's time to take hold of the tiller of fate, and steer yourself to the calm waters of your future.
If you spot a dead whale (or other cetacean) this week, then beware! This is an omen, an omen of death...or possibly a big sale at your local fish market.
A ship needs a rudder, a ship needs a captain, and a ship needs appropriate health and safety information. Who is the captain of your ship, and who has their hand on the rudder? And most importantly, do you have a life-jacket?
When you walk this week, take only tiny steps. It may take you longer to get where you going, but Rome wasn't built in a day.
This is a good time to give up something, particularly if you have an addiction to any illegal narcotics.
This is a good time to reflect on all the things that you are not. For example, you are not an elephant, nor are you an electric toaster.
Don't look behind you, instead concentrate on what lies ahead. The road that takes you on the longest path is the road that will not take you on the shortest path.
Tell a loved one that you love them this week. Also tell someone you hate that you hate them. Life is all about balance.
Forget what you have learned and instead remember only that which you have yet to learn. If you have never learned anything then you will have that much more to remember and will therefore will become a very wise mollusk indeed.
Hold a dinner party on @day@, but don't invite anyone...that will show them!
Try experiencing the quirkier side of life when you next read a book by only reading the odd-numbered pages.
Remember the saying: you can squeeze the life out of a kitten, but a kitten can't squeeze the life out of you.
The next time that you play poker, you should bet everything you have whenever you see a two and fold whenever you see an ace. This might not actually help you that much but it will keep everybody else on their toes.
A woman bearing gifts might not be the present-carrying-female that she seems. Be wary if she (if it is a she) tries making you any toast.
Wedding bells might be ringing this week, but alas these are very, very quiet wedding bells which have had their clappers lined with velvet. You will have to listen very carefully if you want to hear them.
Take an umbrella with you this week when you go to your 'special' appointment. It won't rain but there will be waterworks.
An accident involving tofu will cause you to dial the emergency services this week. Make sure that you have plenty of warmed milk to hand, and don't worry about the resulting stains.
You will not meet any world leaders this week. Try to deal with this unsatisfying news by remaining calm and not attacking anyone with a sporting accessory.
If you make an appointment on @day@ then it will be cancelled, delayed or postponed. The trick will therefore be to make the appointment for a day that you can't make.
If you start reading a new book this week, but skip over every seventh page, it will lead to an unsatisfying conclusion but you will get the book read that much faster.
They say that 'you are what you eat'. But what if you are a cannibal and ate someone famous...would you become that person?
If you can get away with it, try to eat everything with a spoon on @day@. It will impress a secret admirer.
The color @color@ will be very important to you this week. Especially on @day@, and when connected to the number @number@, and if tomato juice is involved, then let's just say that it will be a day to remember.
Should you wear @color@ on @day@? No, but you'll do it anyway because you have no sense of fashion.
A @color@ car will loom large in your life (or maybe just your rear-view mirror) on @day@. Remember to wear your seat belt and try not to have any small animals in your car on that day.
Why will the color @color@ be important to you this week? The answer to that question may only be revealed when you end up in a police station or supermarket on @day@.
You may say to others that you like cats, but this will be the week where you will be tested on how much you <i>love</i> cats. Particularly when a certain cat could unlock the secret to the whereabouts of a long-lost family member.
Sailors would sometimes avoid wearing the color @color@. They would rather a dolphin spit at them in the eye then wear that color. Heed this advice, particularly on @day@.
A cucumber, a pneumatic drill, and a skateboard. Two of these three items will not give you a major headache this week.
The old sailors motto of 'Kick it. Beat it. Cook it. Eat it.' may have special relevance to you this week when you will be faced with an animal that is in your way.
Self-sufficiency is the name of the game for you this week. If you can avoid buying any food, then so much the better.
A few things to avoid this week if you know what's good for you: cold tea, hot milk, three-legged animals, North Dakota, and books with the word 'fun' in their title.
You will have an important meeting with your boss this week. Be careful, the wrong choice of shoes will prove disastrous to your career.
This is a week that is much less about who you are, but much more about who you could be. You are a kitten but you want to be a tiger. Become the tiger!
Do you go for the unhealthy burger or the healthy salad? This is the type of question that will plague you this week. The solution is to go for neither, and instead choose the poached quails eggs. If they don't have quails eggs then I guess you will go hungry.
Computer problems might cause you headaches this week. Best stick to using a pen and paper.
Something will be hot this week. It could be you, it could be the weather, or it could be some mustard. The heat will be good, just remember to stay cool.
Sometimes it is good to try something new. @day@ will offer you the best chance that you will ever have to try something new that involves cheese.
Embarrassment will loom rather large in your life this week as you are very likely to catch your boss in a somewhat compromising situation involving a small animal and some rubber tubing.
There is a very old tradition that Sailors used to follow when leaving home before embarking on a long voyage. Urinate on three things that you love, and spit on three things that you hate. Only this will ensure a safe trip. Heed these words before undertaking any business travel this week.
You may have heard of the saying 'if you can't beat them, join them', but this is a poor choice in comparison to the original nautical version of this phrase. 'If you can't beat them, then shave their beards off while they sleep'.
Indecision will be your undoing this week. You will say yes, only to then say no. You say will 'large iced latte' only to change your mind to a 'small Americano'. You will say 'I do' only to then have second thoughts and run out the church.
In a year's time you might consider running a marathon or eating a Snickers bar. Either way, this week is when you should start your preparations.
On @day@ your week will take a turn for the worse when you attract the (unwanted) attentions of a born-again reincarnationist. They will try to claim you as their soul-mate. You should run away.
If your boss offers to take you out for a drink this week you should gently decline...unless you want to contract a 'downstairs' disease and be involved in an unpleasant (and protracted) divorce settlement.
Wear a smile on this week because you cannot fail* and everything you do will turn out to be magical and rewarding (* = terms and conditions may apply).
Something involving the color @color@ will be on your mind this week and you are not sure if you need a second opinion about what to do. The solution involves getting a second opinion from a friend as to whether you need to get a second opinion.
@day@ will be a very bad day for you. A <b>very</b> bad day indeed. You might find true love, you may win a large cash sum, and you may even get a promotion. These minor successes will in no way compensate for the badness of the bad thing that will happen to you though.
You might be familiar with the saying that 'you cannot buy success', well this might be a good week to try anyway.
A little rodent problem will cause you a major headache this week. Who knew that rats liked ice-cream?
A famous sailor once remarked 'Life is like jumping overboard without knowing how to swim. You will drown. We all drown. Such is life'. Apply this philosophy to your sales presentation on @day@ this week.
Something about the number @number@ will drive you crazy this week. Luckily, the impending failure of your recent investments on @day@ will keep your mind occupied.
You will see someone this week who looks suspiciously like @person@. This will have no bearing on your life whatsoever.
If you see anybody this week who looks like @person@, you should immediately ask them for the time, but only if their watch is on their right wrist. 
This week you might find yourself inconvenienced in an enclosed space with someone who looks remarkably like @person@. Does this matter? Only time will tell.
Your week will become focused around @day, when the the number @number@, the color @color@ and someone who has a connection to @person@ will potentially change your life, or maybe just your bank balance.
In a parallel universe you were born as @person@. Don't get too excited, because you are still living in <b>this</b> universe.
If you are driving and you see someone who looks like @person@ driving a @color@ car, then it is time to leave town immediately otherwise you will be associated with a very bad smell for many months.
The question everyone will be asking this week is 'are you @person@ in disguise?'. No, I don't know what this means either.
Even when everything is going wrong, and it will go wrong this week, just be thankful that you're not @person@,
If you should happen to bump into anyone who looks like @person@, then this is a good omen. You should immediately go out and rent 'Pretty in Pink' to watch. It will change your life.
If someone should happen to comment that you look just a little bit like @person@, then maybe it's time to consider some heavy duty cosmetic surgery.
There is a 32% probability that someone who looks like @person@ will shower you with unusual gifts on @day@. This will only happen though if you are wearing @color@
Take a deep breath and think to yourself 'Is this really who I am?'. If the answer is 'no', then be afraid, very afraid.
Up for a challenge? Then remove all of the labels from any tins in your house. Meal times will then have an element of surprise and danger about them.
If you have a cat, then consider also getting a dog. If you have a dog, then consider getting a cat. If you already have a cat and a dog, then have you ever thought about owning a moose?
Your challenge for this week is to clear your head of all thoughts concerning sex and mustard.
This is a good week to set sail on a new voyage of discovery and adventure...unless you are feeling tired, in which case you should stay at home.
How can something as simple and harmless as a tube of toothpaste cause so much misery? You will find out this week.
On @day@ you will learn the important difference between a large ukulele and a small guitar.
If you believe in the old addage 'you are what you eat', then you should bear in mind that you eat an awful lot of complete garbage.
Don't be surprised when an accidental slip on a calculator this week could lead to a diplomatic incident involving the French Navy.



# Barnacle,December 2nd - February 19th
You may want to keep a fellow Clam close to your side on @day@.
Watch out for someone saying 'no' to you this week.
Time to get it on with someone this week. It only really matters if they have a pulse.
No, no, no, no, no, no, no! Don't give in to the idiots who are wrong.
Beware an advance from a Limpet this week. They will cling to you like an alcoholic clings to a bottle of cheap whiskey.
Your enemies might tell you that you are not a proper mollusk this week, i.e. that you don't belong in society. Ignore them and you will be more of a mollusk than they could ever be.
A old Snail associate will cross your path this week...very slowly. Be patient, this Snail will provide you with much needed culinary relief.
A collision with a Limpet will literally knock you off your feet this week. Don't spend much time arguing whose fault it was but instead try to reach a consensus that it was due to a stupid Slug that you know.
There is a Clam on the war-path, and that Clam is heading your way. Make like a dead whale, and play dead.
Beware the old saying: 'a Squid in need is a Squid indeed'. It might make no sense, but then again neither does the weather.
This week, you should heed the old nautical expression 'See a Slug, hear a Slug, smell a Slug, hit a Slug'.
This is the week where you will need an Oyster by your side, but there will be none to be found. If you get desperate then try searching at either a bar, brothel, or baptism ceremony. These are all natural haunts for the Oyster.
This is the time to leap to the aid of a Scallop that you work with. They will not thank you for your actions, they may well come to despise you for what you do, but it still needs to be done.
On @day@ you will meet an Octopus who will want to punch your lights out. Did you sleep with their partner behind their back? Only you - or the police - will be able to answer that.
When a Barnacle and a Mussel get together, it's a bit like adding treacle to a slow burning fire. You have been warned.
On a good day, a Barnacle and a Clam can be as an effective a double act as @person@ and Tiger Woods.
Your romantic advances towards a business colleague will suffer a setback on @day@. You will soon get over rejection from this idiot. Especially, as you are still in possession of certain compromising pictures of them using a vacuum cleaner in an 'unnatural' manner.
This might be the week where some vegetarian friends taunt you for not being a 'true' Mollusk. The best way of dealing with these people is to slip some goat blood into their coffee, and then taunt them for not being 'true' vegetarians.

# Snail,February 20th - March 9th
If you see a Barnacle this week, you should probably hit them. They are always trouble.
Limpets are losers so avoid them this week.
You know a Clam who deserves to be punched...twice!
Avoid roller-coasters at <b>all</b> costs on @day@.
This is certainly a week where if you see a Slug, then you should give them a punch on the chin.
If you were a shrimp then you would be an outcast among your Mollusk friends. But you are no shrimp, you are a Snail, and don't you ever forget it!
A Limpet you know well, a new carpet, and a weak bladder will combine with tragic consequences this week.
There is a Clam that is going to do something to you this week which will annoy you greatly. But be prepared by buying a good quality stain remover ahead of time.
A Squid will get in your way this week. If you are in a car, then it is fine to run them down.
Stupid is as stupid does, and as stupidity goes, a run in on @day@ with a Slug will have you reaching for your gun (metaphorically). Shoot down the Slug (metaphorically speaking) before they shoot down your dreams.
When a Snail and an Oyster meet, it is a bit like finding a dead animal in your washing machine. However much you try, the smell just won't go away.
On @day@, walk into the nearest bar after you have finished work and find a friendly Scallop to talk to. If you do not know anyone there, then so much the better.
Would you ever be so stupid to get drunk with a Scallop on a work night, and <b>then</b> go to one of those clubs that your mother warned you about? The answer to this question will be revealed on @day@.
When a Snail and an Octopus get together the results can be hard to predict. So take extra special care on @day@ when you will meet an Octopus in an uncomfortable situation (an industrial-strength stain remover might be required).
Think of a beautiful day where you are happy and carefree. Now think of a fat and sweaty Mussel that you know. They will ruin said beautiful day and an unpleasant bout of flatulence will almost certainly be the cause.
Watch out for a Barnacle in a hurry on @day@. If you time it correctly, you will only end up with a small stain to show for their clumsiness. If you get your timing wrong however, you might be facing a stay in the hospital and you won't be eating solids for a long time.
You generally don't get on with Barnacles, but if you meet one on @day@ who looks at all like @person@, then you should kiss them without hesitation. 
It is imperative that you find a roller coaster to ride on @day@. Your life needs some excitement, and if you want an added kick, don't wear the safety harness.
This is a good time in your life to focus on the things that you really, really, want. Especially if those things involve eggs, cheese, or other dairy products.

# Limpet,March 10th - May 1st
Meet up with an Oyster for a fun time on @day@.
Get your friends to form a circle around you, then they can clap and cheer at your brilliance.
Get your creative juices flowing and write a poem about your favorite cheese.
If you are not in the spotlight this week, then you bloody well should be.
Meet up with a Squid this week for some fun and frolics...beware that alcohol and silicon-based lubricants may be involved.
Take a trip to your local art gallery and prepare to be moved by an unusual pasta-based sculpture.
A Clam in your immediate family will cause trouble by revealing all about your dark secret involving the hamster.
Why do Squids have that annoying habit of saying something at the most inopportune times. If you are speaking at any event this week where there is an opportunity to ask questions, then avoid fielding any such questions from a Squid.
When a Slug comes calling at your door, asking for a little financial favor, tell them in no uncertain terms: "You are a poor excuse for a mollusk, and I would rather force-feed myself to a shark than lend you any money".
You are a good Mollusk, you are a trustworthy Mollusk. So why when an Oyster comes calling at your door on @day@ will you be doubting yourself? I don't know. Do you?
Ever been on holiday with a Scallop before? They'll buy you a lot of drinks but they'll expect certain favors in return. You might not like the sound - or the smell - of those favors.
Sometimes you will try hard to avoid them, you will try your best to pass them by in the street or workplace. But on @day@ there is no escape. You will have to go toe-to-toe with an Octopus. Make sure you have an adequate supply of breath mints.
You will bump into a Mussel on @day@. They will not know you, and you will not know them. You will not talk to them, and they will not talk to you. But it is a meeting of profound importance to your life and career.
Have you ever slept with a Barnacle and regretted it? If not, then this might be the week to try.
'Slow-but-steady' may be the motto of your so-called Snail 'friend'. But what if they are speeding around with you partner behind your back? Don't be heartbroken, just think of how much money they have and then think about that good old word 'blackmail'.
You once knew a Limpet who was vile, repugnant, and had a tendency to sweat heavily. Well bad news for you because that very same family member will be knocking on your door this week.
You've always had two secret role models, but up till now they were so secret you didn't know who they were. Let the truth be revealed, for you secretly covet @person@ and Tiger Woods.
A lively discussion with an old friend will end in one of two ways this week. Either you will resort to bare-knuckle fighting, or you will end up reciting poetry to them. Either way, onlookers will be greatly enthralled.
This is a good week to remember that old nautical expression 'You can make me walk the plank, but I'll drown on my own terms'.



# Clam,May 2nd - June 2nd
You know a Squid who is in trouble this week. Time for a bit of Clam-support.
Whatever anyone says to you, it's not worth telling your boss about...except if you hear a rumor involving mushrooms or anti-wrinkle cream.
Remember, your lips are sealed. If you happened to disclose a certain secret to a certain someone this week, then a certain career (i.e. yours) might be ruined.
You would rather stick a knife in your eye than disclose a less than important secret to your boss.
Hook up with an Octopus on @day@ if you want to see a good time that doesn't involve ambulances.
A Squid in need is a Squid indeed. This is the week to hang out with your Squid buddy and see what pops out of the toaster.
Hook up with a Squid this week in order to relieve those bedroom tensions. Try to avoid using tinned fruit though.
You have a few personal problems at the moment and you might feel that you should turn to a colleague for advice. But asking a Slug for advice is like stepping into a bath full of kerosene and then lighting a firework.
There is an Oyster in town who is looking for a good time. You can join in that good time, but be prepared to run up an substantial credit-card bill, and don't expect to see your shoes again anytime soon.
Get out and enjoy life on @day@, and if you happen to spot a little Scallop who is in need of a good time then so much the better. But remember: Clams and Scallops, good. Clams and Scallops and alcohol, bad.
If you put a Clam and an Octopus together, it is a bit like Laurel and Hardy. There will be much stupidity and much clumsiness. There will also be a lot of pain.
Get in a tussle with a Mussel and they will feel the slam of a Clam.
There is an old saying that goes something like this: 'A Clam, a Clam, a Clam! All I need is a Clam...but a Barnacle might be ok as well'. Heed this warning on @day@.
You may have heard the old sailor's expression 'you can never fail with a Snail'...but you do know that there is an exception to every rule right? Walk very carefully on @day@ when said Snail will try to take you somewhere that a Clam should never go.
You may have heard of the question 'How many Limpets does it take to change a light bulb?', but have you heard of the question 'How many Limpets does it take to change a pacemaker?'. You will this week.
You are not @person@, so don't try to act like them...unless you have a good lawyer of course.
You have an Octopus pal who will need of a shoulder to cry on this week. They are in the wrong, they did the wrong thing, and it will turn out all wrong, but you probably won't want to mention any of that.
Life will be a little bit tough for you this week. Just a little bit though, sort of squidgy-tough rather than hard-tough.


# Squid,June 3rd - July 25th
If a Barnacle, Oyster, or Mussel says anything to you at all this week, don't believe them.
If you have to lie about your age, height, weight, or gender this week, then it's probably for the best.
Tell someone that they look great this week...even if they are pig ugly.
Try applying for a passport using a false identity. It might not work, you might be arrested, but it might be fun trying.
You might know of a Slug who is in trouble this week. But as they are a Slug, you probably won't want to help them.
You may be asked your age this week by a close business colleague...they may be trying to get you into trouble so you should probably lie.
You will see a Slug in trouble this week. You will not care. You are the better Mollusk.
When you and that lovable Oyster colleague of yours get together, then sparks will fly. Unfortunately, that might lead to a charge of arson this week, so best cancel that @day@ night get-together.
Get together for a Scallop this week if you want to have a fun time that involves an activity that is not yet illegal in all countries.
Throw yourself into the (many) arms of an Octopus this week and you will find out whether what they say about an Octopus in an elevator is true.
If you have the time, try to track down a trustworthy Mussel that you know on @day@. Tell them a big secret and see how trustworthy they really are.
Bless your Barnacles, for a Barnacle will come to save the day for you on @day@. You would have never guessed that peanut butter would prove so useful.
On @day@, the color @color@, the number @number@ and a certain little Snail that you know will all combine to create a lot of trouble for you and your pet Yak. You don't have a pet Yak yet, but that's just part of the trouble that you'll be getting into.
Limpets, Limpets everywhere, but not a drop of love to spare. Well this might be the case for you on @day@ when a rabid, potentially-drunk Limpet will cause plenty of trouble for you.
Invite a Clam to dinner this week on @day@. This will be the one day that they can't make, so easy brownie points for you!
On @day@ night you will dream of being @person@. You won't know why until the following @day@ when a chance meeting with an international patent attorney will shed much light on this mystery.
Make like a fox this week and be cunning. Especially when someone is out to deceive you into buying a beef-based product that you really don't need.
This might be a very important week for your career, especially if you carry a jar of mustard with you everywhere you go.

# Slug,July 26th
Your miserable existence will take a further turn for the worse this week, so be prepared to sink to new lows.
Why do you try so hard, when everything you do fails?
Now is not the time to shed a tear. Now is the time to weep uncontrollably.
You have tried so hard, and accomplished so little. Now is the time to give up.
Your friends talk about you behind your back. Are they pathetic...or are you?
Everything you try to achieve ends up being surrounded in failure, perhaps you should consider early retirement?
You have nothing to offer anyone this week, so it's business as usual.
Too many cooks spoil the broth, but if you are making the broth, then you will spoil it all by yourself.
There is an elephant in the room. You are the elephant.
In a week where everything that can go wrong, will go wrong, you just have to accept that this is largely your fault.
Have you looked outside recently? If you have you will have noticed that it has been dull and gloomy for some time. A bit like you.
You can cry, you can weep, you can rant and rage, you can demand attention. You can do all of these things and more, but the bottom line is that maybe you deserve it.
There will be good news on Wednesday this week. However, it will turn out to be very bad news by Friday.
Be careful what you choose to eat this week...there is a lot of food poisoning about.
Romance looms large this week. But not for you unfortunately.
This week, you should be wary of the hapless idiot...especially when the idiot in question is you.
Have you ever truly been happy? Probably not.
Your friends will gather closely around you this week, so please take steps to lessen your foul odor.
You need to talk to people to tell them how you really feel about things. They desperately want to know how you feel. Well, maybe not desperately. Actually, they don't really want to know how you feel...or even if you are still drawing breath.
The sound of thunder will hang over you until you can put a smile on your face. As you are one of the most miserable people around, this may not be easy.
The person that you have a secret crush on does not feel the same way about you. If they knew how you truly felt, then they would probably be violently ill.
You might want to take a second look in the mirror at some point this week...just to confirm that you really are that ugly.
Remember, things can go only get better...actually for you they can probably still get quite a bit worse.
You need to go on a low-sodium diet to improve your health...pity this won't improve your looks though.
There is a chance that things will go well for you this week...remember though, there is also a chance that pigs might fly.
You may have heard of the saying "Don't worry, be happy"...well, that doesn't apply to miserable idiots like you.
A friend will come to you seek your advice on a sensitive subject this week. They will also come to deeply regret asking you about anything because your advice sucks.
You have dandruff, do something about it!
You may be feeling down. You may be feeling that nothing good ever happens to you. But don't worry. Just remember, that 99.9% of the rest of the population are much happier than you. So at least it all balances out!
Happiness. Joy. Financial success. Just another three things that you will not experience this week.
In your hour of need, an Oyster that you know will have all the answers to your problems this week. However, they are not going to tell you any of the answers.
You will receive a call this week with fantastic news about a possible love interest. Unfortunately, it will be a wrong number.
Hanging out with a Scallop on @day@ might gain you some attention as you bask in the aura of Mr/Mrs Popular. However, they will hate you for this unwanted association and your evening may well end with the threat of extreme physical violence.
With such a tragic life, with an existence full of misery, you may think you are a suitable candidate for the Guinness World Record of 'Most miserable life'. Don't think about applying for this record however...you will be rejected.
Life is looking good for you this week...actually, that's not exactly true. More likely, life is looking very bad for you.
Is it possible for everyone you know to violently dislike you? Yes. It is.
The number @number@ will be important for you this week. This will possibly be an amount of money that you will lose, or the number of days you might be held for questioning by the police.
A long lost family member will appear in your life once again this week. You will be overcome with emotion at meeting up with this person. That is until you find out that they have only tracked you down to ask you for money.
Days to avoid this week include Thursday, Saturday and Sunday. Also Monday might be bad and Friday has an outside chance of being a miserable day. Wednesday is not looking too good either. But Tuesday will be ok...except if you have to talk to anyone in which case it will be a very bad day indeed.
If everything goes to plan this week then you will be a very happy Mollusk indeed. Chance are though, that it will fall to pieces...again!
Improve your popularity this week by a) not saying anything to anyone and b) wearing a bag over your head.
It's ok, your complete failure to achieve anything of significance in life is not entirely your fault...oh wait a minute, yes it is.
You will go to an auction on @day@. You will pay too much for something that you won't be able to sell and which you will take an instant disliking too the moment after you buy it. You are an idiot.
You will be followed about by a bad smell everywhere that you go this week. This is not much of a mystery, the smell is you.
This week your colleagues will be trying to heed the words of the old nautical expression 'If you see a Slug, run for your lives'.
Things will be mostly crap for you this week, but on the plus side of things, you will already know exactly what this feels like.
It's a tough life being a Slug. Nobody likes you, nobody wants to be around you, and nobody can stand your personal hygiene problems. Are you just misunderstood? Actually, no.
This is going to be a very good week...not for you personally, but you can't have everything.
In Roman times, ancient mariners had a special word for people who are Slugs. That word translates from the original latin to 'eternal failure'.
Just give up making any sort of plans this week. They will all fail so best stay in bed.
You will be very popular this week and will receive lots of mail. Oh, actually they are all overdue bills as you have forgot to pay off your utility bills...you idiot.
You've been thinking about having some minor cosmetic surgery done, but here's a word of warning...if you polish a turd, it's still a turd.
You know the old saying 'Don't worry, be happy'? Well you <i>will</i> worry, and you <i>won't</i> be happy. Such is the life of a Slug.
Do you remember that when you were young, that your parents said 'When you grow up, you have the potential to do anything you want to in life'? They were lying. You only have the potential to be a failure.
Did you know that 'Slug' is very nearly an anagram of 'ugly'. This is quite fitting as your grim features are enough to put a dying dog off its food.
Want some advice? Trying to be popular is never going to work. An alternative solution would be to crawl under a large rock and stay there.
On @day@, the number @number@ will signify bad news. Really. Bad. News. 
Your Slug-like nature will mean that you will suffer twice as much as normal this week when a rival colleague will attempt to <i>literally</i> rub salt into old wounds.
One more week on the planet, means another week of learning and discovery; it also means that you're one week closer to your death.

# Oyster,July 27th - August 19th
If the level of your confidence was a country, it would be Australia.
Get the guys or girls around your place on @day@ for a lurve fest.
Get some attention this week by wearing 7 items of clothing on Monday, and then remove an item each day
Does it really count as adultery if you don't tell anyone?
Invite a Squid over this week for some mollusk-on-mollusk action.
Is is really vanity if you pay to put an advert in a national newspaper to point out to everyone how beautiful you are?
Take note of the old saying 'An Oyster and a Scallop is like quarter pounder and cheese...only without the cheese'.
You will get romantically entangled with an Octopus this week. They will regret it, but the quantities of alcohol involved mean that you won't remember anything so don't worry too much about it.
How many times do you get a Mussel trying to chat you up over a drink and a hot dog? Well this is the week where a Mussel with a point to prove will try to ply you with hot dogs and beer. Just go easy on the mustard!
Ever get stuck in an elevator with a Barnacle? Well be prepared for that eventuality on @day@. Also be prepared for a very bad body odor problem.
Ever hear the joke about the Oyster and the Snail who lived next door to each other. They drove each other to drink. Then they drove each other to hard drugs. Then they became the best of friends and started playing Scrabble together on a regular basis. Let that be a lesson to you (if you live next door to any Snails).
Do you know a Limpet? Do you <b>want</b> to know a Limpet? If the answer is yes, then on @day@ night make your way to where the cool people go. And take lots of loose change with you.
Take extra special care on @day@ because your life might be changed forever by a chance encounter with a Clam. The Clam will demand one of the following: money, sex, or citrus fruit. If you can meet their demands, then things will work out well for you. If you can't, then you will spend the rest of your life regretting it.
You will fall in love with a Squid on 6:45 am on @day@. By 7:15 you will realize that actually they are quite repulsive.
On @day@, your day will be swiftly ruined by an odious Slug that you know. You can't prevent what they are going to do, the only thing you can do is feel a small degree of satisfaction when you sue them for every penny they've got.
This week, if you meet a Squid that looks at all like @person@ then you might be in for some fun times. If however, you meet a Squid that resembles Tiger Woods, then you will almost certainly become violently ill before the end of the week.
Other Mollusks would say that 'one-on-one is fun', but you are an Oyster, in which case you should adhere to the 'eight-on-eight is great' school of bedroom philosophy.
Look yourself in the mirror on @day@ and say to yourself "I'm an Oyster, an Oyster, an Oyster!". If you don't say this, no-one else will.



# Scallop,August 20th - October 1st
Beware, Clams are plotting against you! And even if they are not actually plotting, they are probably thinking about plotting. And even if they are not thinking, they will be.
If someone offers you any food this week, then beware! It might be spiked with pepper. You should no longer trust this person, even if you are married to them.
Eat anything you want this week, but avoid the kung po chicken at all costs.
Your sex-life could be greatly improved by judicious use of peanut butter this week. Naturally, 'Crunchy' would be better than 'Smooth'.
This week you may take any life-threatening actions that come your way. But whatever you do, go easy on the chili sauce.
A Clam you know will offer to cook for you this week. Be careful, they might have ulterior motives, and they will certainly try spiking your food with Tabasco sauce.
A distantly-related Octopus will offer an interesting opportunity to you this week. Whether to accept that offer will depend heavily on a) whether you trust your wife and b) how quickly you are prepared to learn Korean.
A female Mussel friend will give you something very precious this week, try not to blow the moment by commenting on her oversized rear.
You will come to the defense of a Barnacle this week when a common friend insults them for "not being a true Mollusk".
Take a Scallop and a Snail. Two very similar Mollusks who are also so entirely different. On @day@ you will find out just how similar or different you are when you will be inadvertently stuck in a toilet cubicle with said Snail.
On @day@ you might want to try playing Limpet limbo, but only if you know any sexually-charged Limpets. Otherwise stay at home with a good book.
A Clam that you know will try to kill you this week. Well maybe they are just plotting the act at this stage. Actually, they might only be thinking about it. On second thoughts, it's more of a vague intention. So don't worry about it too much. Just be careful around them if they are holding any sharp objects.
Squids are all around you this week. Try not to get smothered in their tentacles. One Squid in particular will try to make romantic advances towards you. If you can't smell their hideous body odor, then you are a perfect match. If you can, then you are not.
The best thing you can do to help a Slug in trouble this week is remind them what a failure they are and that you would help, only they will probably be in trouble again next week so why bother?
Hot fudge sauce will be your downfall this week, and the reason for this is that you will believe the foolish advice of a Scallop that you know. Believe me, hot fudge sauce is never the solution to problems in the bedroom.
On @day@ just remind yourself that you are lucky to not have been born a Slug.
That stranger who you keep seeing in your neighborhood, the one who looks a bit like Tiger Woods, well you can rest easy because they're not @person@ at all.
A stupid Clam friend that you know will prove very bothersome on @day@. Just ignore them. Unless they start removing clothes in which case you should just run away.
You might find it useful to spend some this week in the company of seagulls. Just make sure you wear appropriate headwear.


# Octopus,October 1st - October 29th
Someone will swear at you this week. You will not be happy, in fact you will be livid. In these scenarios, physical retribution is only fair.
You know which way is north and that ain't no lie. Use this information to your advantage on @day@ when a navigationally-challenged colleague will seek your guidance.
If you hear so much as one mention of the F-word from a friend or colleague, then forcefully wash their mouth out with soap (or battery acid).
Tell a loved one that you are going to take up base jumping. You're not going to do this of course, but it's good to keep people on their toes.
Your map reading skills might just help you save a stranded puppy this week.
If you hear just one more person swear within a 20-foot radius of you, then it is time to tear up the map and get out of this town.
A portly Mussel that you work with will literally get in your way this week. You might want to tactfully suggest that the fat lump of lard should go on a diet.
You will find yourself in one of those situations where time is of the essence this week. However, a portly Barnacle involved in a roller-skating experiment is going to ensure that your scheduling goes out the window.
If a Snail tries buying your affections by spending vast amounts on money on you, then don't fall for it. It may make you happy, but happiness is not everything...at least not when the Snail in question has spent time inside for attempted manslaughter.
When you and a Limpet get together on @day@, sparks will literally fly. That's what you get when a chance encounter with a welder goes horribly wrong.
Turn up on time for a meeting with a Clam on @day@ and experience the 'Clocktopus Effect' - a beneficial outcome that will have arisen because you were on time.
When a Squid and an Octopus meet it's full-on tentacle action. So if you are out and about on @day@, then make sure you take enough moisturizer.
This will be a week full of stress and angst for you. Try releasing that angst by finding a Slug that lives in your street. Wait for them to leave their home and then paint the words 'I am better than you' on their doors and windows. You will feel much better.
When an Oyster that you know comes around to visit you on @day@ and asks if you can help them out with a little financial problem, be very careful. Offer them drugs. Offer them sex. But do not offer them money!
The letters F, Y, and K will all be very important to you this week, especially in conjunction with a Scallop wearing @color@. Be especially cautious if they offer you a hot-dog, but don't offer you any mustard.
Given the choice, you might think that you would have preferred to be born as @person@, but the reality is that you would end up spending a lot more money on lubrication products.
You are starting to tire of a colleague's constant profanity in the workplace. It would be great if they were to 'accidentally' be punched in the throat. Well one can dream.
Wake up at 3:00 AM on @day@ to remind yourself why it's such a bad idea to get up at 3:00 AM.


# Mussel,October 30th - December 1st
Tick tock, someone will be running late for a meeting with you. They are lazy fools.
To be on the safe side, arrive 5 hours early for your special work meeting this week.
You will kill yourself if you arrive late for work this week so purchase 7 alarm clocks to be on the safe side.
Make some sweet love in the afternoon...about 3:43 pm.
Remember, it is always better to arrive early. Arriving late is a sign of a drunken loser.
There is a time and a place for everything. This week, that time will more often that not be 8:22 am.
You will see a Barnacle in considerable distress this week. If they are left-handed, you should step in to help, otherwise keep walking.
A casual comment by a Snail acquaintance of yours might make you think twice before making that important purchase this week. Don't worry. As long as they have it in red, things will turn out just fine.
You will be asked to look after a Limpet this week. That may be a good thing but it may be a bad thing. Be especially careful on @day@ when said Limpet might ask you to do something which could be considered illegal in many countries.
'Wham, bam, thank you Clam'...that might be a motto for you to learn this week as Clam-antics in the bedroom will get you all worked up.
There is a Squid that you really like. There is a Squid that really likes you. Unfortunately they are half your age and live on the other side of the world.
You might be feeling low this week, things might not be going so well for you. There is a silver lining to your cloud of depression though. It could be worse, you could be a Slug.
Try relaxing on @day@ evening in the company of an Oyster. Just make sure you don't let them consume too much alcohol else they might leave you with an embarrassing stain to clear up.
You might get some advice about this week from a Scallop about which orifice is most suitable for a particular pursuit that you might try on @day@. Please get a second opinion from someone else before you embark on said pursuit.
You know an Octopus who is almost the perfect person. Polite, charming, attractive, and financially independent. Sadly, you chose to marry their poor, ugly, and alcoholic cousin instead.
You're a Mussel, so that's good. But you know a Barnacle who resembles @person@ a little too much for your liking, so that's not so good. Well that's life I guess, it's all about balance.
When you say meet me at 3:47 pm, you of course mean 'meet me at 3:47 pm'. So when a stupid Slug that you know turns up at 3:49 pm, you are entitled to walk away and never talk to them again.
What you lack in wisdom, you make up for in strength. So maybe this is a good week to settle an argument with a fist fight.


