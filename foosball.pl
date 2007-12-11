#!/usr/bin/perl -w
#
# foosball.pl
#
# A script to manage the GC Foosball league
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;

# two main hashes for storing data, the stats hash will contain details of games played [0], 
# games won [1] and games drawn [2]. The ratings hash will just contain the floating point rating of
# each player
my %player2ratings;
my %player2stats;

# another simple hash to allow quicker player input, enter their rank in the table to look
# up their name
my %rank2player;

# track maximum all-time score and player
my $max = 0;
my $max_player;

my $game_counter = 0;


# in command-line mode?
if(defined($ARGV[3])){
	$player2ratings{"A"} = $ARGV[0];
	$player2ratings{"B"} = $ARGV[1];
	my $a_score = $ARGV[2];
	my $b_score = $ARGV[3];
	$player2stats{"A"}[0] = 1;
	$player2stats{"B"}[0] = 1;
	&calculate_ratings("A","$a_score","B","$b_score","n");
	exit;
}

# read results from a csv file
&read_results;

#&calculate_ratings("Ian Korf","1","Victoria Whitworth","10","y");

# print two tables, one for people who have played more than 5 games, one for the rest
my ($high,$low);

foreach my $key (reverse sort {$player2ratings{$a} <=> $player2ratings{$b}} keys %player2ratings){
	my $formatted =  sprintf("%.2f",$player2ratings{$key});

	if($player2stats{$key}[0] >= 10){
		$high .= "\"$key\",$player2stats{$key}[0],$player2stats{$key}[1],$player2stats{$key}[2],$formatted\n";	
	}
	else{
		$low .= "\"$key\",$player2stats{$key}[0],$player2stats{$key}[1],$player2stats{$key}[2],$formatted\n";	
	}
}
print "$high\n\n$low\n\n";

print "Max score so far achieved is $max by $max_player\n";
exit(0);


sub read_results{
	open(IN,"$ARGV[0]") || die "Can't open input file\n";
	my @details;

	print "\n\n";
	while(<IN>){
		$game_counter++;
		chomp;
		@details = split(/,/);
		
		# if competitive match field not specified, set to 'n'
		($details[4] = "n") if (!defined($details[4]));			
		
		&calculate_ratings("$details[0]","$details[1]","$details[2]","$details[3]","$details[4]");
	}	
}
close(IN);



sub calculate_ratings{
	# take details of players names, scores, and whether match was a competitive game or now (y/n)
	my $p1 = shift;
	my $score1 = shift;
	my $p2 = shift;
	my $score2 = shift;
	my $competitive = shift;

	# add players to hash if they are new and set rating and player stats to zero
	if(!defined($player2ratings{$p1})){
		$player2ratings{$p1} = 0;
		$player2stats{$p1}[0] = 0;
		$player2stats{$p1}[1] = 0;
		$player2stats{$p1}[2] = 0;
	}
	if(!defined($player2ratings{$p2})){
		$player2ratings{$p2} = 0;
		$player2stats{$p2}[0] = 0;
		$player2stats{$p2}[1] = 0;
		$player2stats{$p2}[2] = 0;
	}

	# need to assign score of higher ranking player to $p1_score and lower ranking 
	# to $p2_score as this affects the results. If rankings are equal then it doesn't
	# make any difference
	my ($p1_score, $p2_score, $p1_name, $p2_name);
	
	if($player2ratings{$p1} >= $player2ratings{$p2}){
		$p1_score = $score1;
		$p2_score = $score2;
		$p1_name = $p1;
		$p2_name = $p2;
	}
	else{
		$p1_score = $score2;
		$p2_score = $score1;
		$p1_name = $p2;
		$p2_name = $p1;
	}
	
	# calculate difference in ratings and score
	my $ratings_diff = $player2ratings{$p1_name} - $player2ratings{$p2_name};
	my $score_diff = abs($p1_score - $p2_score);

	# stats for ELO rating
	my $k;   # weighting constant, depending on competitive or friendly and margin of victory
	my $g;   # goal coefficent based on margin of victory
	my $w1;  # result of game (1 for winner, 0.5 each for a draw, 0 for loser)
	my $w2;  # result of game (1 for winner, 0.5 each for a draw, 0 for loser)
	my $we1; # expected result, based on a formula which uses $dr
	my $we2; # expected result, based on a formula which uses $dr
	my $dr1; # difference in ranking (can be positive or negative depending on player)
	my $dr2; # need one score for each player

	# player 1 scores and statistics should always be the player that had the higher rating to
	# begin with (and hence is expected to win). If ratings are the same it doesn't matter
	# who is assigned to player 1.

	# calculate $we1 (higher ranking player) and $we2 (lower ranking player)
	# These two scores should sum to 1
	$we1 = sprintf("%.3f",1/(10**(-$ratings_diff/400)+1));
	$we2 = sprintf("%.3f",1/(10**( $ratings_diff/400)+1));

	print "====================== Game $game_counter ==========================\n";
	print "$p1_name: played $player2stats{$p1_name}[0], rated $player2ratings{$p1_name}, Win expectation = $we1\n";
	print "$p2_name: played $player2stats{$p2_name}[0], rated $player2ratings{$p2_name}, Win expectation = $we2\n\n";


	# calculate $w (i.e. who won?) and add wins and draws to player2stats
	if ($p1_score == $p2_score){
		($w1 = $w2 = 0.5); 
		print "$p1_name ties with $p2_name ($p1_score - $p2_score)\n";
		$player2stats{$p1_name}[2]++;
		$player2stats{$p2_name}[2]++;
	}
	if ($p1_score > $p2_score){
		$w1 = 1;
		$w2 = 0;
		print "$p1_name beats $p2_name ($p1_score - $p2_score)\n";
		$player2stats{$p1_name}[1]++;
	}
	elsif($p1_score < $p2_score){
		$w1 = 0;
		$w2 = 1;
		print "$p2_name beats $p1_name ($p2_score - $p1_score)\n";
		$player2stats{$p2_name}[1]++;
	}

	# add game count to player2stats
	$player2stats{$p1_name}[0]++;
	$player2stats{$p2_name}[0]++;

	
	# calculate $k and $g
	$g = 1   if ($score_diff <= 1);
	$g = 1.5 if ($score_diff == 2);
	$g = ((11+$score_diff)/8) if ($score_diff > 2);
	$k = 20 if ($competitive eq "n");
	$k = 30 if ($competitive eq "y");

#	print "\n";	
#	print "Match status weighting (k) = $k, goal difference weighting (g) = $g\n";

	# now can calculate final ratings
	my $rating1 = sprintf("%.2f",$k*$g*($w1-$we1));
	my $rating2 = sprintf("%.2f",$k*$g*($w2-$we2));

	# add these points to their current ratings
	$player2ratings{$p1_name}+= $rating1;
	$player2ratings{$p2_name}+= $rating2;
	print "\n";
	print "$p1_name gains $rating1 points, new total $player2ratings{$p1_name}\n";
	print "$p2_name gains $rating2 points, new total $player2ratings{$p2_name}\n";	
	
	print "========================================================\n\n\n";

	# update max rating details if necessary
	if($player2ratings{$p1_name} > $max){
		$max = $player2ratings{$p1_name};
		$max_player = $p1_name;
	}
	if($player2ratings{$p2_name} > $max){
		$max = $player2ratings{$p2_name};
		$max_player = $p2_name;
	}
}
