# Written by Gabriel Cain, <gabriel@dreamingcrow.com>
# Licensed under the GNU GPL v2.

package Twitter;

use strict;
use LWP::UserAgent;
use HTTP::Cookies;

# This package implements a twitter object that is
# used to connect to twitter, and "do stuff"

my $twitter_post = q{http://twitter.com/statuses/update.json};
my $twitter_friends_timeline = q{http://twitter.com/statuses/friends_timeline.rss};

# Constructor.  
# new ( <user>, <pass>, [cookie jar] )
#
# returns undef or new Twitter obj.

sub new {
	my $class = shift;
	my $self = {};
	bless $self, $class;
	
	# get user/pass
	my ($u, $p) = (shift,shift);
	
	# and optional cookie file.
	my $cookiefile = shift || "$ENV{HOME}/.twitter_cookies";
	
	# Create object for LWP
	$self->{ua} = LWP::UserAgent->new;
	
	# If we fail, eat shit and die
	return undef unless $self->{ua};
		
	# set some parameters up.
	$self->{ua}->timeout(10);

	# Basic auth is for the win.
	$self->{ua}->credentials( "www.twitter.com:80", "Twitter API", $u, $p );
	$self->{ua}->credentials( "twitter.com:80", "Twitter API", $u, $p );
		
	# It's fun to store our cookies.
	$self->{ua}->cookie_jar( HTTP::Cookies->new( 
		file => $cookiefile, autosave => 1 ) 
		);
	
	# Return created object.
	return $self;
}

sub whatsup {
	my $self = shift;
	
	return $self->{ua}->get( $twitter_friends_timeline );
}

# say ( something )
# returns:
#	0 - no message to say
#	1 - said it
# 	some other value: an error

sub say {
	my $self = shift;
	my $msg = shift;
	
	# bail on empty messages
	return undef unless $msg;
	
	# Post it
	$self->{ua}->post( 
		$twitter_post, { status => qq{$msg} }
	);
	
}


1; # eof #####
