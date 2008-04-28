#!/usr/bin/perl 
                               
use strict;
use warnings;
use Getopt::Long;

########################
# Command line options
########################

my $local;      # option to *NOT* ignore hits from local network (128.120.136.* unless -local specified)
my $images;     # option to *NOT* ignore image files unless -images specified
my $forums;     # ignore to *NOT* hits to the discussion forums

GetOptions ("local"  => \$local,
			"images" => \$images,
			"forums" => \$forums);

open(LOG,"/var/log/apache2/access_log") || die "could not open log file\n";

# hash to keep scores
my %pages2hits;
my %pages2sites;

while(my $temp=<LOG>){
    chomp($temp);
	
	# only want GET requests
	next unless ($temp =~ m/GET/);
	
	# only want successful (200) requests
	next unless ($temp =~ m/ 200 /);
	
	# ignore hits from server console
	next if ($temp =~ m/^::1/);

	# ignore favicons
	next if ($temp =~ m/favicon\.ico /);
	
	# ignore style sheets
	next if ($temp =~ m/\.css /);

	# ignore hits from genome center network
	unless($local){
		next if ($temp =~/128\.120\.136/);
	}
	
	unless($images){
		next if ($temp =~ m/(\.gif|\.jpg|\.jpeg|\.png)/i);
	}

	unless($forums){
		next if ($temp =~ m/\/Forums\//i);
	}

	# tidy up some URLs
    $temp =~ s/.*%7E/\/~/;
    $temp =~ s/.*%7e/\/~/;           
    $temp =~ s/index\.html//;
	$temp =~ s/(\/~\w*) /$1\/ /;  
	$temp =~ s/(.*) - -.*GET* (.*) HTTP.*/$2 $1/;
	$temp =~ s/-//g;
    $temp =~ s/\`//g;
    $temp =~ s/\|//g;
    
	my ($url,$site) = split(/ /,$temp);

    $pages2hits{$url}++;
	push(@{$pages2sites{$url}},$site);
       
}
close(LOG);

# sort array
foreach my $key (sort by_number keys %pages2hits) {
	print "$pages2hits{$key}\t$key\n";
}

sub by_number { $pages2hits{$b} <=> $pages2hits{$a}; }
exit(0);
