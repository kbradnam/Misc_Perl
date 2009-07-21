#!/usr/bin/perl
#
# This program is designed to implement twittering.  Isn't
# that awesome?  It's intended for monitoring purposes, 
# and similar activities.
#
# This program doesn't require XML parsers cuz I am cheating
# and not using one.  I know what to look for, so I'm 
# simply parsing out the response.
#
# Written by Gabriel Cain, <gabriel@dreamingcrow.com>
# Licensed under the GNU GPL v2.

use strict;
use lib '/home/gabriel/lib';
use Data::Dumper;
use Getopt::Std;
use Twitter;

# optionally support curses
eval { use Term::ANSIColor qw(:constants); $Term::ANSIColor::AUTORESET = 1; };

my $colorEnabled = ($@ ) ? 0 : 1;

if( $colorEnabled ) {
	print BLUE, "\t\t\t *** ANSI COLOR ENABLED ***\n", RESET;
}

# Set our options
# -f = alternate config file
# -d = Remove last twitter
# -D <name> = d <Name> message
# -F <name> = follow <Name>
# -r = get recent stuff
# -x = no color

# Process command-line args
my %o;
getopts('xhrdD:F:f:', \%o);

if( $o{x} ) {
	$colorEnabled = 0;
}


# Read our configuration from $HOME/.twitterrc
# first line = username, 2nd line = password

my $tw_rc = $o{f} || "$ENV{HOME}/.twitterrc";

if( ! -f $tw_rc ) {
	doCreateTwitterRC();
}

die("$tw_rc: $!") unless -f $tw_rc;	
open(T,"<$tw_rc") or die("Failed to open $tw_rc: $!");
	
my $tw_user = <T>;
my $tw_pass = <T>;
close(T);

# Eat useless newlines.
chomp( $tw_user, $tw_pass );

# Test for -d (delete twit)
if( $o{d}) {
	print "Sorry, delete isn't implemented yet.\n";
	exit 0;
}

# Direct & empty is contraindicated
if(  length $o{D} == 0 && defined $o{D} ) {
	print "Empty direct recipient not ok.\n";
	exit;
}

# fail if we're not fetching & no message supplied
if( ($#ARGV < 0 && ! $o{r}) ||  $o{h} ) {
	print << "END";
-=-=-=-=-=-=< Twitter Command Line Client >-=-=-=-=-=-=-
Written by Gabriel Cain, <gabriel\@dreamingcrow.com>

Usage: 
	twitter.pl <message> 
	twitter.pl -r
	twitter.pl [options]

Options:
	-f			Read alternate config file, 
				default is \$HOME/.twitterrc
	-d			Remove last twitter (not implemented yet)
	-D <name> <message>	Send a direct message to a user
	-F <name> 		Follow a user
	-r			Get friends timeline (twenty most recent)
	-x			Do not use colored output (forces)

This program is licensed under the GNU GPL v2.  Enjoy it.  

END
	exit;
}

# Build twitter message.
my $twit;
$twit = "d $o{D} " if $o{D}; # direct message if -D <name> 
$twit .= join ' ', @ARGV;

# Log messages make the baby twitter cry
if( length $twit > 140 ) { 
	warn("Warning: Message longer than 140 characters\n");
}

# Follow
if( defined $o{F} && length $o{F} > 1 ) {
	$twit = "follow $o{F}"
}

my $res;
my $T = new Twitter( $tw_user, $tw_pass );
die("Object couldn't be made\n") unless $T;

if( $o{r}) {
	$res = $T->whatsup();
}
else {
	$res = $T->say( $twit );
}

# Check response, and if we failed, complain about it.
if( ! $res->is_success ) {
	print "Request failed: ".$res->status_line."\n";
	exit;
}

# If we're reading our "friends" page, then we need to spit
# back a list of what our friends say.

if( $o{r}) {
	# Parse content for the items and dates.
	my @items;
	
	my @content = split /<item>/, $res->content;
	foreach my $c (@content) {
		# examine each <item> .. </item> block for what we
		# need.
		my $i = {};
		$c =~ /<description>(.*)<\/description>/mi;
		$i->{status} = $1;
		$c =~ /<pubDate>(.*)<\/pubDate>/mi;
		$i->{date} = $1;
		next if $i->{date} =~ /Twitter updates from/;
		push @items, $i;
	}

	# Spit them out in reverse order -- newest at bottom
	foreach ( reverse @items ) {
		my $a = sprintf "%s", $_->{date};
		my ($u,$s) = split /:/, $_->{status}, 2;

		if( $colorEnabled ) {
			print GREEN, $a, RESET, ":", RED, $u, 
				RESET, ":", CYAN, $s, "\n", RESET;
		}
		else {
			print $a, ":", $u, ":", $s, "\n";
		}
	}
	
}

	

sub doCreateTwitterRC {
	# create $HOME/.twitterrc
	my ($l,$p);
	open( R, ">$ENV{HOME}/.twitterrc") or die("Couldn't create .twitterrc.\n");
	$| = 1;
	print "What's your twitter login? ";

	chomp( $l = <> );
	print "What's your twitter password? ";
	chomp( $p = <> );

	# should have login/passwd
	print R "$l\n$p\n";
	close(R);
}

