#!/usr/bin/perl -w 

# blackjack

use strict;

# set up initial pack of 52 cards
my @pack =qw(AH 2H 3H 4H 5H 6H 7H 8H 9H 10H JH QH KH 
	     AC 2C 3C 4C 5C 6C 7C 8C 9C 10C JC QC KC 
	     AD 2D 3D 4D 5D 6D 7D 8D 9D 10D JD QD KD 
	     AS 2S 3S 4S 5S 6S 7S 8S 9S 10S JS QS KS);

my %cards2points;
$cards2points = ("AH" => "A",
		 "2H" => "2",
		 "3H" => "3",
		 );

# how many decks in shoe?
my $decks = 6;

# set up inititial shoe
my @shoe;

for (my $i=1; $i <= $decks; $i++){
    @shoe = (@shoe,@pack);
}


# shuffle routine
&shuffle_shoe;

print "@shoe\n";






###########################################################
#
# Shuffle routine, will re-order cards in shoe
#
###########################################################


sub shuffle_shoe{
    
    # how many cards in total
    my $size = $decks * 52;

    # need to fill a new shoe when shuffling
    my @new_shoe;

    my $number;
    
    while (@shoe){
	
	# choose random position in old shoe and move that card to new shoe
	$number = int(rand($size));
	push(@new_shoe,$shoe[$number]);    

	# now remove that card from the old shoe (delete array position)
	splice(@shoe,$number,1,@shoe[$number..-1]);

	# decrement number of cards left in old shoe
	$size--;
    }
    
    # finally set new shoe to replace shoe
    @shoe = @new_shoe;
    @new_shoe = ();
}
