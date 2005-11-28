#!/usr/bin/perl -w 

# blackjack

use strict;

my @pack;
my %ace_lookup;
my $decks = 6;
my @shoe;

# how many cards to be dealt from shoe before reshuffle
my $stop;

# set up inititial shoe based on how many decks there are
# shuffle shoe and set penetration
&setup;


# deal
my $card_count = 0;

# two arrays to store dealers and players cards
my @dealer;
my @player;

# two variables to keep track of total score of cards for
# dealer and player
my $d_score;
my $p_score;

my $card;

while($card_count < $stop){

    #deal the first four cards
    &initial_deal;

    # reset counters and arrays
    $d_score=$p_score=0;
    @dealer=@player=();
}
print "Should stop at $stop\n";


exit(0);



##############################################################
#
# setup routine, build shoe, shuffle it and set penetration
#
##############################################################


sub setup{
    # set up initial pack of 52 cards
    @pack =qw(AH 2H 3H 4H 5H 6H 7H 8H 9H 10H JH QH KH 
	     AC 2C 3C 4C 5C 6C 7C 8C 9C 10C JC QC KC 
	     AD 2D 3D 4D 5D 6D 7D 8D 9D 10D JD QD KD 
	     AS 2S 3S 4S 5S 6S 7S 8S 9S 10S JS QS KS);

    # make shoe from packs of cards
    for (my $i=1; $i <= $decks; $i++){
	@shoe = (@shoe,@pack);
    }

    # shuffle routine
    &shuffle_shoe;

    # Calculate the degree of penetration, need to know at what point to stop
    # dealing from shoe. Penetration will be between 56% (worst case scenario,
    # remove 45% of shoe) and 85% (best case scenario, remove 15% of shoe)
    my $penetration = 56+int(rand(30));
    $stop=length(@shoe)*$penetration;
    
    #
    %ace_lookup = (1  => [1,11],
		   2  => [2,12], 
		   3  => [3,13],	   
		   4  => [4,14],
		   5  => [5,15],
		   6  => [6,16],
		   7  => [7,17],
		   8  => [8,18],	   
		   9  => [9,19],
		   10 => [10,20],
		   11 => [11,21],
		   12 => [12],
		   13 => [13],
		   14 => [14],
		   15 => [15],
		   16 => [16],
		   17 => [17],
		   18 => [18],
		   19 => [19],
		   20 => [20],
		   21 => [21],
		   22 => [22]);       	   
}


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


sub initial_deal{

    #deal two cards for dealer
    push(@dealer,shift(@shoe),shift(@shoe));
    $card_count += 2;
    my @dealer_score = &score("D");
       
    print "Dealer) $dealer[0]  ?\t @dealer_score \t\t ($dealer[0],$dealer[1])\n";

    
    #deal two cards for player
    push(@player,shift(@shoe),shift(@shoe));
    $card_count += 2;
    my @player_score = &score("P");
			 
    print "Player) $player[0],$player[1] \t @player_score\n\n";
}


sub score{
    # score for player or dealer?
    my $type = shift;

    my @hand;
    (@hand=@dealer) if ($type eq "D"); 
    (@hand=@player) if ($type eq "P");
   
    my $score=0;
    my @scores;
    my $ace_count = 0;

    foreach my $card (@hand){
      
	$card =~ s/[CDHS]//;
	$card =~ s/[JQK]/10/;
    
	if($card eq "A"){	
	    $ace_count++;
	}
	else{
	    $score += $card;
	}
    }

    # if there are aces, we need to return an array of all possible scores
    if($ace_count>0){
	# look up possible scores in ace_lookup and for each possible score,
	# add to @scores array
	my @ace_scores = @{$ace_lookup{$ace_count}};
	foreach my $ace_score (@ace_scores){
	    push(@scores,$score+$ace_score);
	}	
	return(@scores);
    }
    # no aces, just return straight-forward score
    else{
	return($score);
    }    
}


