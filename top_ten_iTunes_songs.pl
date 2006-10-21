#!/usr/bin/perl -w
#
# A script to check the iTunes music library for problems
# (missing music files, duplicate artwork, no year etc.)
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;

my $library = "/Users/Nod/Music/iTunes/iTunes Music Library.xml";
print "Using library: $library\n";


# The basic variables, it's pretty obvious what they will store
my ($artist,$album,$rating);
my $label; 				# stores album name + artist joined by '@@@' (to prevent confusion when different artists have the same album name)
my %artists;  			# hash structure for all ratings by same artist
my %albums;				# hash structure for all ratings for each album
my %artists2average;	# hash which will be used for sorting artists by rating
my %albums2average;		# hash which will be used for sorting albums by rating
my $chart =10;			# how many songs to print in the final 'charts', defaults to 10 if not specified on command line
$chart = $ARGV[0] if ($ARGV[0]);
my %ratings2counts; # rating is key, value is count

###################################################################################
#
# Main loop parses iTunes XML file and extracts album, artist, and rating details
#
###################################################################################

open(IN,"<$library") || die "Couldn't open library file\n";

while(<IN>){
    last if (m/<key>Playlists<\/key>/);    

    if (m/<key>Artist<\/key><string>(.*)<\/string>$/){
	$artist = $1;
    }
    if (m/<key>Album<\/key><string>(.*)<\/string>$/){
		$album = $1;
	    $label = $album."@@@".$artist;
		
		# keep track of overal song count for artists and albums
		$artists{$artist}[0]++;
		$albums{$label}[0]++;

	}
    if (m/<key>Rating<\/key><integer>(.*)<\/integer>$/){
		$rating = $1;
		$ratings2counts{$rating}++;
		
		# store total of all ratings so far ([1]), number of songs rated ([2])
		$artists{$artist}[1] += $rating;	
		$artists{$artist}[2] ++;
		$artists2average{$artist} = $artists{$artist}[1] / 	$artists{$artist}[2];
	
	
		$albums{$label}[1] += $rating;
		$albums{$label}[2] ++;
		$albums2average{$label} = $albums{$label}[1] / $albums{$label}[2];
    }
}
close(IN);


# now sort albums by average rating

my $counter = 0;
my $final_rating;

print "\nALBUMS\n";

my $title;
foreach my $key (sort {$albums2average{$b} <=> $albums2average{$a}} (keys(%albums2average))){
	my $percent = $albums{$key}[2] / $albums{$key}[0];
	if(($percent >=0.75) && ($albums{$key}[2] >5)){
		$counter++;	
    	# need to convert to a 1-5 star rating
    	$final_rating = sprintf("%.2f",$albums2average{$key}/20);
	
		# need to separate out artist and album info
		$title = $key;
		$title =~ s/@@@/, /;
	    print "$counter) $final_rating - $title ($albums{$key}[2] songs out of $albums{$key}[0] rated)\n";
	}
    last if ($counter == $chart);
}



# now do same thing for artists

$counter=0;
print "\n\nARTISTS\n";

foreach my $key (sort {$artists2average{$b} <=> $artists2average{$a}} (keys(%artists2average))){

	my $percent = $artists{$key}[1] / $artists{$key}[0];
	if(($percent >=0.50) && ($artists{$key}[2] >10)){
		$counter++;	
    	# need to convert to a 1-5 star rating
    	$final_rating = sprintf("%.2f",$artists2average{$key}/20);
	    print "$counter) $final_rating - $key ($artists{$key}[2] songs out of $artists{$key}[0] rated)\n";
	}
    last if ($counter == $chart);
}


# Finally, print counts of each rating
print "\n\nRATING COUNTS\n";

foreach my $key (sort {$a <=> $b} (keys(%ratings2counts))){
    my $rating = sprintf("%.0d",$key/20);
    print "$rating - $ratings2counts{$key}\n";
}
exit(0);
