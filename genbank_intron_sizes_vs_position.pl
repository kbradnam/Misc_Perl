#!/usr/bin/perl -w
#
# genbank_intron_sizes_vs_position.pl
#
# a script to extract intron sizes for all genes in all organisms and average the size according to
# whether it is the 1st, 2nd, 3rd, nth intron. Do this for all species combined, and also for 
# separate species
#
# by Keith Bradnam
#
# created 11th November 2005
# drastic rewrite 4th April 2007 (to take account of genbank_info_dump.pl)
#
# Last updated by: $Author$
# Last updated on: $Date$    
#
#######################################################################################


use strict;
use Getopt::Long;
use List::Util qw(sum);

########################
# Command line options
########################
my $limit;     # set a cut-off for how many introns are needed per species
my $min;       # Only match genes with a minimum number of introns
my $max;       # only match genes with a maximum number of introns
my $input;	   # the name of the intron dump file from genbank_info_dump.pl

GetOptions("limit=i"  => \$limit,
	       "min=i"    => \$min,
	       "max=i"    => \$max,
		   "input=s" => \$input);

die "Must use -input option to specify a file containing intron info from genbank\n" if (!defined($input));

#set default for $limit,$min, and $max if not specified
$limit = 100  if (!defined($limit));
$min   = 1     if (!defined($min));
$max   = 10000 if (!defined($max));


########################
# misc variables
########################

# simple hash for keeping track of all species: key is species name, values are just set to '1'
my %all_species;

# hash for counting how many introns occur in any given species for any given number of introns
# e.g. how many C. elegans genes have four introns
my %species2count;

# Now store everything else in one complex hash structure where key is species name, values are an array of arrays of arrays!
# first array covers number of introns in a gene 
# e.g {Caenorhabditis elegans}[4] is for worm genes with exactly 4 introns 
# 
# second array covers the intron position within the former array:
# e.g. {Homo sapiens}[3][2] would refer to the 2nd intron of all human genes that have exactly 3 introns
#
# the final array is a list of all intron sizes in the category of the former array
# e.g. {Bos taurus}[7][1] = (34,52,244,121)  - four cow genes have 7 introns exactly and th lengths of the first introns are 34,52 etc
my %species2introns;



############################################################
#
#  M A I N   L O O P
#
#  process genbank dump file to put info into hashes
#
############################################################

open (IN,"<$input") || die "Can't open $input file\n";

while (<IN>) {
	# split line and extract relevant information to separate variables
	# need to know how long each line is to work out how many introns there are
	chomp;
	
	# skip the very few occurrences of negative (or zero) sized introns and move to next line in file
	next if (m/,-/);
	next if (m/,0/);
	
	my @details = split(/,/);
	my $length = @details;
	my ($species,$accession,$cds) = (@details[0..2]);
	my @introns = (@details[3..$length-1]);
	
	# might have occasional problem with some species names in genbank
	# names should be alphanumerical characters, spaces, and maybe a period and dash
	next if ($species !~ m/[\w\s\.\-]+/);
	
	# should always have at least four entries per line 
	next if ($length < 4);
	
	# add to species2count hash
	$species2count{$species}{$length-3}++;
	
	# put all introns into species specific hash of arrays of arrays of arrays!
	for(my $i = 3; $i<@details;$i++){		
		push(@{$species2introns{$species}[$length-3][$i-2]},$details[$i]); 
	}
	$all_species{$species} = 1;
}

close(IN) || die "Couldn't close $input\n";

# now loop through big data structure to get stats. Species first...

foreach my $species (sort(keys(%all_species))){

	# now loop through number of introns in a gene 1..x
	for (my $i=1; $i < @{$species2introns{$species}};$i++){
	
		# skip if $min and/or $max has been set
		next if ($i < $min);
		next if ($i > $max);
		
		# do we actually have data for the current number of introns, and do we have
		# more genes than defined by $limit
		if(defined($species2count{$species}{$i}) && ($species2count{$species}{$i}>$limit)){
			#print "$species - $species2count{$species}{$i} genes with $i introns\n";
		
			# now loop through each intron position to calculate average length, stdev etc. at each position
			my (@means,@std_errs);
			for (my $j=1; $j < @{$species2introns{$species}[$i]};$j++){
				
				# calculate statistics and load to temporary arrays
				my $n = @{$species2introns{$species}[$i][$j]};
				my $mean = sum(@{$species2introns{$species}[$i][$j]})/$n;
				my $std_dev = sqrt(sum(map {($_ - $mean) ** 2} @{$species2introns{$species}[$i][$j]}) / ($n-1));
				my $std_err = $std_dev/sqrt($mean);
				push(@means,sprintf("%.2f",$mean));
				push(@std_errs,sprintf("%.2f",$std_err));
				#print "Intron position $j: @{$species2introns{$species}[$i][$j]}\n";
				#print "Mean = $mean Stdev = $std_dev Std Err =  $std_err\n";
				}
			# now print results
			print "$species,$i,$species2count{$species}{$i},";
			for(my $k=0;$k<@means;$k++){
				print "$means[$k],";
			}
			print "\n";
			print "$species,$i,$species2count{$species}{$i},";
			for(my $k=0;$k<@std_errs;$k++){
				print "$std_errs[$k],";
			}
			print "\n\n";
			
		}
	}
}

	
exit(0);
