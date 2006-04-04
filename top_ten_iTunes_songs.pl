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

open(IN,"<$library") || die "Couldn't open library file\n";

my %album2ratings;
my %album2rating_count;
my %album2artist;
my %album2song_count;

my ($artist,$album,$rating);

while(<IN>){
    last if (m/<key>Playlists<\/key>/);    

    if (m/<key>Artist<\/key><string>(.*)<\/string>$/){
	$artist = $1;
    }
    if (m/<key>Album<\/key><string>(.*)<\/string>$/){
	$album = $1;
	($album2artist{$album} = $artist) if ($artist);
	$album2song_count{$album}++;
    }
    if (m/<key>Rating<\/key><integer>(.*)<\/integer>$/){
	$rating = $1;
	$album2ratings{$album} += $rating;
	$album2rating_count{$album}++;
    }
}
close(IN);


# now need a hash to store average rating for each album 
my %album2average; 

foreach my $key (keys %album2ratings){

    # want at least 5 songs on album and at least 75% of album rated
    my $percent = $album2rating_count{$key}/$album2song_count{$key};
    if(($album2song_count{$key} >= 5) && ($percent >= 0.75)){
	
	my $average_rating = sprintf("%.2f",$album2ratings{$key} / $album2rating_count{$key});

	$album2average{$key} = $average_rating;	
    }
}


# now sort albums by average rating

# how many songs to show
my $chart = 10;
my $counter = 0;

foreach my $key (sort sort_by_rating (keys(%album2average))){
    $counter++;
    print "$counter) $album2average{$key} - $key, $album2artist{$key} ($album2rating_count{$key} out of $album2song_count{$key} songs rated)\n";

    last if ($counter ==$chart);
}



exit(0);


# sort subroutine to sort by hash values
sub sort_by_rating{
    $album2average{$b} <=> $album2average{$a};
}
