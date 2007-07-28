#!/usr/bin/perl
# Quick dictionary lookup using regular expressions
use strict;
use warnings;

# get words from dictionary
my @words = `cat /usr/share/dict/words`;

my $answer = 1;
while(1){
		# ask for pattern unless one specified on command line
		if($ARGV[0]){
			$answer = $ARGV[0];
			$ARGV[0] = "";
		}
		else{
	        print "Enter pattern for regular expression:";			
	        chomp($answer = <STDIN>);
		}
        last unless(length $answer && defined $answer);
        print  grep eval {m/$answer/}, @words;       
        print "Continuing after error: $@\n\n" if ($@);
}
exit(0);