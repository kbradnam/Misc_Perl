#!/usr/bin/perl
#
# foosball_league.pl
#
# A script to manage the GC Foosball singles league
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;

# two main hashes for storing data, the stats hash will contain details of games played [0], 
# games won [1] and games drawn [2], goals scored [3] and goals conceded [4]
my %player2stats;

# hash for final scores for each player (3 points for a win, 1 point for a draw)
my %player2scores;

# read results from a csv file
my $game_counter = 0;
&read_results;

my $table;
foreach my $key (reverse sort {$player2scores{$a} <=> $player2scores{$b}} keys %player2scores){

	# set win and draw count to zero if not defined
	($player2stats{$key}[1] = 0) if (!defined($player2stats{$key}[1]));
	($player2stats{$key}[2] = 0) if (!defined($player2stats{$key}[2]));
	
	my $score =  $player2scores{$key};
	$table .= "\"$key\",$player2stats{$key}[0],$player2stats{$key}[1],$player2stats{$key}[2],";
	$table .= "$player2stats{$key}[3],$player2stats{$key}[4],$score\n";	

}
print "$table\n";

exit(0);


sub read_results{
	open(IN,"$ARGV[0]") || die "Can't open input file\n";
	while(<IN>){
		$game_counter++;
		chomp;
		my($p1,$s1,$p2,$s2)= split(/,/);			
		
		# increment games played counter
		$player2stats{$p1}[0]++;
		$player2stats{$p2}[0]++;
		
		# increment wins or draw counters
		$player2stats{$p1}[1]++ if ($s1 > $s2);
		$player2stats{$p2}[1]++ if ($s2 > $s1);
		$player2stats{$p1}[2]++ if ($s1 == $s2);
		$player2stats{$p1}[2]++ if ($s1 == $s2);
				
		# increment goals scored and conceded
		$player2stats{$p1}[3]+= $s1;
		$player2stats{$p2}[3]+= $s2;
		$player2stats{$p1}[4]+= $s2;
		$player2stats{$p2}[4]+= $s1;

		# set win and draw count to zero if not defined
		($player2stats{$p1}[1] = 0) if (!defined($player2stats{$p1}[1]));
		($player2stats{$p1}[2] = 0) if (!defined($player2stats{$p1}[2]));
		($player2stats{$p2}[1] = 0) if (!defined($player2stats{$p2}[1]));
		($player2stats{$p2}[2] = 0) if (!defined($player2stats{$p2}[2]));
		
		# can now calculate total score (3 points for a win, 1 point for a draw
		$player2scores{$p1} = (3 * $player2stats{$p1}[1]) + (1 * $player2stats{$p1}[2]);
		$player2scores{$p2} = (3 * $player2stats{$p2}[1]) + (1 * $player2stats{$p2}[2]);
	}	
	close(IN);

}


