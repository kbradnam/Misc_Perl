#!/usr/bin/perl
#
# tw.pl
#
# A lightweight twitter client that uses cURL
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;

my $user;      # get info about a specified user
my $recent;    # show recent tweets by friends
my $following; # show who you are following
my $friends;   # show who is following you

GetOptions ("user=s"    => \$user,
		    "recent"    => \$recent,
		 	"following" => \$following,
			"friends"   => \$friends
);


# resource file containing user name and password 
my $twitter_rc = "$ENV{HOME}/.twitterrc";

# read config file
open(TWIT,"<$twitter_rc") or die("Failed to open $twitter_rc: $!");
my $tw_user = <TWIT>;
my $tw_pass = <TWIT>;
chomp($tw_user,$tw_pass);
close(TWIT);

my $tweet;
if($ARGV[0]){
	$tweet = "$ARGV[0]";
	die "Tweet is longer than 140 characters\n" if (length($tweet) > 140);
	`curl -s -u $tw_user:$tw_pass -d status=\"$tweet\" http://twitter.com/statuses/update.json`;
	exit;
}
elsif($user){
	my $output = `curl -s http://twitter.com/users/show.xml?screen_name=$user`;
	my ($name) = $output =~ m/<name>(.*)<\/name>/;
	my ($bio) = $output =~ m/<description>(.*)<\/description>/;
	my ($friend_count) = $output =~ m/<friends_count>(.*)<\/friends_count>/;
	my ($follow_count) = $output =~ m/<followers_count>(.*)<\/followers_count>/;
	my ($last_tweet) = $output =~ m/<text>(.*)<\/text>/;

	print "$name\n";
	print "Follows: $friend_count\n";
	print "Followers: $follow_count\n";
	print "Bio: $bio\n";
	print "Last tweet: $last_tweet\n";
	exit;
}
elsif($recent){
	my @output = `curl -s -u $tw_user:$tw_pass http://twitter.com/statuses/friends_timeline.xml`;
	my ($name,$last_tweet);
	foreach my $line (@output){
		print "$line\n";
		if ($line =~ m/<name>.*<\/name>/){
			($name) = $line =~ m/<name>(.*)<\/name>/ ;
		}
		if ($line =~ m/<text>.*<\/text>/){
			($last_tweet) = $line =~ m/<text>(.*)<\/text>/ 
			
		}
		if($line =~ m/<\/status>/){
			print "$name:\t$last_tweet\n";
		}
	}
	exit;
}
