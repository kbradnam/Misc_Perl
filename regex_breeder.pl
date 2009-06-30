#!/usr/bin/perl
#
# regex_breeder.pl
#
# A script to 
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use List::Util qw(sum);

my $file;               # input file
my $n;                  # population size
my $generations;        # number of generations
my $max_motifs;         # maximum number of motifs in each motif set
my $min_r;              # minimum correlation value needed from motifsets in starting population
my $max_pattern_length; # maximum pattern length in regular expression
my $starting_motif;     # optionally seed one motif in starting generation 

GetOptions ("n=i"              => \$n,
			"max_motifs=i"     => \$max_motifs,
			"file=s"           => \$file,
			"min_r=f"          => \$min_r,
			"max_pattern=i"    => \$max_pattern_length,
			"starting_motif=s" => \$starting_motif
			);


# set defaults
$n = 10000              if (!$n);
$generations = 1000     if (!$generations);
$max_motifs = 2         if (!$max_motifs);
$min_r = 0.25           if (!$min_r);
$max_pattern_length = 3 if (!$max_pattern_length);
die "Specify intron file with -file option\n" if (!$file);





# turn on autoflush
$| = 1;

# big multi-dimensional array to hold all of the data
# plus a second array to hold copy of data for mutating
my @motifset;
my @copy;

# want to keep track of sequences from input file, plus their expression level
my @seqs;
my @expression;
my %seq_to_id;

# need list of all IUPAC bases plus a translation to turn those into regular expressions
# overload the list of bases with most common ones (these will be selected at random)
my @bases = qw(A A A A A A A A C C C C C C C C G G G G G G G G T T T T T T T R R R Y Y Y S S S W W W K M B D H V N);
my %bases_to_regex = (
	A => 'A',
	C => 'C',
	G => 'G',
	T => 'T',
	R => '(A|G)',
	Y => '(C|T)',
	S => '(C|G)',
	W => '(A|T)',
	K => '(T|G)',
	M => '(C|A)',
	B => '(C|G|T)',
	D => '(A|G|T)',
	H => '(A|C|T)',
	V => '(A|C|G)',
	N => '(A|C|G|T)'
);

# read sequence files and also extract expression values
read_file();

# one time set-up
print "Creating starting population\n";
create_starting_population($n);

# mutate half of the population 
mutate_copy();

# add mutated copies back to original set
@motifset = (@motifset,@copy);


#############################
#
#     M A I N   L O O P
#
#############################

my %motifsets_to_scores;

for (my $i = 1; $i <= $generations; $i++){
	print "\nGeneration $i\n";
	
	# wipe relevant hashes and arrays
	%motifsets_to_scores = ();
	@copy = ();

	# an array to keep the winning n/2 survivors of selection
	my @survivors;

	# grab scores from each motifset
	score_motifsets($i);
	
	# now need to remove the worst $n/2 motifsets
	my @sorted = (reverse sort {$motifsets_to_scores{$a} <=> $motifsets_to_scores{$b}} (keys(%motifsets_to_scores)));
 
	# simple counter
	my $c = 1;

	my $print_limit = 10;
	
	# loop through scores and make a new motifset array of the top n/2 survivors
	foreach my $motifset (@sorted){
		
		# print details of top few motifs
		my $number_of_motifs = @{$motifset[$motifset]};
		print "$c) $number_of_motifs motifs r = $motifsets_to_scores{$motifset}\n" if ($c <= $print_limit);
		for my $j (0 .. $number_of_motifs-1){

			# get motif out of data structure
			my $number_of_bases = scalar @{$motifset[$motifset][$j]};

			# assemble motif as a string out of @motifset data, also keep a regex version
			my $motif;
			my $motif_regex;
			for my $k (0 .. $number_of_bases-1){
				$motif .= $motifset[$motifset][$j][$k]{base};
				$motif_regex .= $bases_to_regex{$motifset[$motifset][$j][$k]{base}};

				# only print out min and max if max is greater than 1 (otherwise it's just one base)
				if($motifset[$motifset][$j][$k]{max} > 1){
					$motif .= "{" . $motifset[$motifset][$j][$k]{min} . ",";
					$motif .= $motifset[$motifset][$j][$k]{max} . "}";

					$motif_regex .= "{" . $motifset[$motifset][$j][$k]{min} . ",";
					$motif_regex .= $motifset[$motifset][$j][$k]{max} . "}";
				}

				# add details of this motifset to two new arrays, survivors + copy
				$survivors[$c-1][$j][$k]{base} = $motifset[$motifset][$j][$k]{base};
				$survivors[$c-1][$j][$k]{min} = $motifset[$motifset][$j][$k]{min};
				$survivors[$c-1][$j][$k]{max} = $motifset[$motifset][$j][$k]{max};

				$copy[$c-1][$j][$k]{base} = $motifset[$motifset][$j][$k]{base};
				$copy[$c-1][$j][$k]{min} = $motifset[$motifset][$j][$k]{min};
				$copy[$c-1][$j][$k]{max} = $motifset[$motifset][$j][$k]{max};

			}
			print "\tmotif $j ($number_of_bases nt): $motif\t$motif_regex\n" if ($c <= $print_limit);
		}
				
		# stop at the halfway point
		last if ($c == ($n/2));
		$c++;
	}	
	
	mutate_copy();
	
	# add mutated copies back to original along with the best survivors from above
	@motifset = (@survivors,@copy);
}



exit;

sub read_file{
	open(IN,"<$file") || die "Couldn't open $file file\n";

    my $fasta = new FAlite(\*IN);

    # loop through each sequence in target file
    while(my $entry = $fasta->nextEntry) {

 		# get and trim header
    	my $header = $entry->def;
    	$header =~ m/.* x(\d{1,2}\.\d)/;
    	my $expression = $1;

		die "$header\n" if (!$expression);
    	$header =~ s/ .*//;
		$header =~ s/>//;

        my $seq = uc($entry->seq);
		$seq_to_id{$seq} = $header;
		
		# add sequence and expression values to arrays
		push(@seqs,$seq);
		push(@expression,$expression);
	}
	close(IN);
}

sub create_starting_population{
	my $n = shift;
	$n /= 2;
	
	MOTIF: for my $i (0 .. $n-1){		
		my $number_of_motifs = int(rand($max_motifs));
		for my $j (0 .. $number_of_motifs){

			my $motif_length = int(rand(6))+1;
			
			for my $k (0 .. $motif_length){
				
				# choose a random base (using IUPAC codes)
				my $base = $bases[int(rand(@bases))];
				
				# start off with minimum pattern length of 1 
				my $min = 1;
				
				# and choose maximum number, bias this to be mostly 1, 5% chance of it being 2
				my $increment = int(rand(1.05));
				my $max = $min + $increment;

				$motifset[$i][$j][$k]{base} = $base;
				$motifset[$i][$j][$k]{min}  = $min;
				$motifset[$i][$j][$k]{max}  = $max;
				
				# now test score of this motif
				# same information is repeated into a different set which will be mutated
				$copy[$i][$j][$k]{base} = $base;
				$copy[$i][$j][$k]{min}  = $min;
				$copy[$i][$j][$k]{max}  = $max;
			}
		}
		# get score for this motifset
		my $r = score_motifset($i);
#		print "motifset $i, r = $r\n";

		# reject and try again if below some threshold (to ensure a good starting population)
		if ($r < $min_r){
			splice(@motifset,$i,1);
			splice(@copy,$i,1);
			redo MOTIF;
		}
		else{
			foreach my $trigger qw(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1){
				print (($trigger*100) . "% ") if ((($i+1)/$n) == $trigger);
			}
		}
	}
	print "\n";
	
	# if -starting_motif has been specified, we overwrite first_motif of motif set
	if($starting_motif){
		# erase existing motif at that position
		@{$motifset[0][0]} = ();
		@{$copy[0][0]} = ();

		my $counter = 0;
		my @motif = split(//,$starting_motif);
		
		while(@motif){
			my ($min,$max);			
			my $base = shift(@motif);
			my $next_base;
			if($motif[0]){
				$next_base = $motif[0];
			}
			else{
				$next_base = "END";
			}
			
			# are both positions regular characters with no length pattern?
			if ($base =~ m/[ACGTRYSWKMBDHVN]/ && $next_base =~ m/[ACGTRYSWKMBDHVN]/){		

				# if so then just set min and max length to be equal to 1
				$motifset[0][0][$counter]{base} = $base;
				$motifset[0][0][$counter]{min}  = 1;
				$motifset[0][0][$counter]{max}  = 1;
				$copy[0][0][$counter]{base} = $base;
				$copy[0][0][$counter]{min}  = 1;
				$copy[0][0][$counter]{max}  = 1;
		
				$counter++;
			}
			else{
				shift(@motif);			# get rid of '{'
				$min = shift(@motif);
				shift(@motif);			# get rid of ','
				$max = shift(@motif);
				shift(@motif);			# get rid of '}'	
				
				$motifset[0][0][$counter]{base} = $base;
				$motifset[0][0][$counter]{min}  = $min;
				$motifset[0][0][$counter]{max}  = $max;
				$copy[0][0][$counter]{base} = $base;
				$copy[0][0][$counter]{min}  = $min;
				$copy[0][0][$counter]{max}  = $max;
				
				$counter++;
			}
		}
	}
}

sub mutate_copy{
	
	for (my $i=0; $i < @copy; $i++){
		
		###################################
		# 1) Add 1 new motifs to a motifset
		###################################

		if(rand(1) < 0.25){
			
			my $number_of_motifs = scalar @{$copy[$i]};

			# can only add a motif if we have less than max number of $max_motifs
			if($number_of_motifs < $max_motifs){
				# next position in array will be equal to size of array
	 			my $next_position = scalar @{$copy[$i]};

				# create random details for new motif
				my $motif_length = int(rand(6));

				for my $k (0 .. $motif_length){

					# choose a random base (using IUPAC codes)
					my $base = $bases[int(rand(@bases))];

					# choose mininum number of copies of that base 
					my $min = 1;

					# and choose maximum number
					my $increment = int(rand(1.05));
					my $max = $min + $increment;

					$copy[$i][$next_position][$k]{base} = $base;
					$copy[$i][$next_position][$k]{min}  = $min;
					$copy[$i][$next_position][$k]{max}  = $max;
				}				
			}
		}

		###################################
		# 2) Remove 1 motif from motifset
		###################################

		if(rand(1) < 0.25){
			my $number_of_motifs = scalar @{$copy[$i]};

			# can't remove motif if there is only one
			last if ($number_of_motifs == 1);
			
			# choose which motif to remove
			my $rand = int(rand($number_of_motifs));
			splice(@{$copy[$i]},$rand,1);
		}
		
		#########################################
		# 3) Add 1 base to one motif in motifset
		#########################################
		
		if(rand(1) < 0.25){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base added
			my $motif = int(rand($number_of_motifs));
			
			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));
			
			# choose a random base (using IUPAC codes)
			my $base = $bases[int(rand(@bases))];

			# choose mininum number of copies of that base 
			my $min = int(rand(2))+1;

			# and choose maximum number
			my $increment = int(rand(2));
			my $max = $min + $increment;
	
			# to simplify things we will just add base at end of motif (for now)
			$copy[$i][$motif][$number_of_bases]{base} = $base;
			$copy[$i][$motif][$number_of_bases]{min}  = $min;
			$copy[$i][$motif][$number_of_bases]{max}  = $max;
		}

		################################################
		# 4) Remove 1 base from one motif in motifset
		################################################
		
		if(rand(1) < 0.3){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base added
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to remove new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};

			# can only do this if you have at least 2 bases
			
			if($number_of_bases >1){
				my $position = int(rand($number_of_bases));

				# add details of new base to @copy
				splice(@{$copy[$i][$motif]},$position,1);
			}
		}
		
		################################################
		# 5) Substitute 1 base from one motif in motifset
		################################################

		if(rand(1) < 0.75){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base substituted
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));

			# choose a random base (using IUPAC codes)
			my $base = $bases[int(rand(@bases))];

			# add details of new base to @copy
			$copy[$i][$motif][$position]{base} = $base;			
		}

		###################################################
		# 6) Increment minimum pattern length for one base 
		###################################################

		if(rand(1) < 0.1){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base substituted
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));

			my $min = $copy[$i][$motif][$position]{min};
			my $max = $copy[$i][$motif][$position]{max};

			# increment value
			# can only do this if incremented value of min is less than max
			$copy[$i][$motif][$position]{min}++ if ($copy[$i][$motif][$position]{min}+1 < $max);
		}

		################################################
		# 7) Decrease minimum pattern length for one base
		################################################

		# deliberately make it more likely to reduce pattern length than increase it
		if(rand(1) < 0.6){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base substituted
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));
			
			my $min = $copy[$i][$motif][$position]{min};
			my $max = $copy[$i][$motif][$position]{max};

			# decrease value
			# can only do this if min is greater than 1
			$copy[$i][$motif][$position]{min}-- if ($min > 1);
		}


		################################################
		# 8) Increment maximum pattern length for one base
		################################################

		if(rand(1) < 0.5){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base substituted
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));
			
			my $min = $copy[$i][$motif][$position]{min};
			my $max = $copy[$i][$motif][$position]{max};

			# increment value (as long as this doesn't exceed $max_pattern_length )
			if($copy[$i][$motif][$position]{max} < $max_pattern_length){
				$copy[$i][$motif][$position]{max}++;				
			}
			die "$i $motif $position $min $max *\n" if (!defined($max) || !defined($min));
		}
		
		################################################
		# 9) Decrease maximum pattern length for one base
		################################################

		if(rand(1) < 0.5){
			my $number_of_motifs = scalar @{$copy[$i]};

			# choose a motif to have a base substituted
			my $motif = int(rand($number_of_motifs));

			# choose a position in that motif to add new base
			my $number_of_bases = scalar @{$copy[$i][$motif]};
			my $position = int(rand($number_of_bases));

			my $min = $copy[$i][$motif][$position]{min};
			my $max = $copy[$i][$motif][$position]{max};

			# decrease value
			# can only do this if max will still be greater than $min
			$copy[$i][$motif][$position]{max}-- if (($copy[$i][$motif][$position]{max}-1 >= $min));	
		}
	}
}


sub score_motifsets{
	my $generation = shift;
	for my $i (0 .. $n-1){

		my $number_of_motifs = scalar @{$motifset[$i]};

		# total count will be the count of all motifs in each motifset
		my $total_count = 0;

		# need to store set of total counts for each sequence (to compare against expression values)
		my @counts;
		
		# may need to remove motifs from a motif set if they do not occur in any sequence
		# otherwise there is no selection pressure to remove them as they won't make any contribution 
		# to the correlation score
		my @motifs_to_remove = ();

		for my $j (0 .. $number_of_motifs-1){
			
			# separate counter just for each motif
			my $motif_count = 0;

			# get motif out of data structure
			my $number_of_bases = scalar @{$motifset[$i][$j]};

			#die "No bases: i = $i, j = $j, num_motifs $number_of_motifs\n" if ($number_of_bases == 0);
			
			# assemble motif as a string out of @motifset data
			my $motif;
			for my $k (0 .. $number_of_bases-1){
				$motif .= $bases_to_regex{$motifset[$i][$j][$k]{base}};
				$motif .= "{" . $motifset[$i][$j][$k]{min} . ",";
				$motif .= $motifset[$i][$j][$k]{max} . "}";
			}

			# check to see if motif starts or ends with an N, if so can remove this base (but only if motif has more than one base)
			if($motifset[$i][$j][0]{base} eq "N"){
				$number_of_bases = scalar @{$motifset[$i][$j]};
				splice(@{$motifset[$i][$j]},0,1) if ($number_of_bases > 1);								
			}
			$number_of_bases = scalar @{$motifset[$i][$j]};
			if($motifset[$i][$j][$number_of_bases-1]{base} eq "N" && ($number_of_bases > 1)){
				$number_of_bases = scalar @{$motifset[$i][$j]};
				splice(@{$motifset[$i][$j]},$number_of_bases-1,1);
			}

			# now loop through each sequence
			for (my $l=0; $l < @seqs; $l++){
				my $seq = $seqs[$l];
				
				my $id = $seq_to_id{$seq};
				
				# count motif in sequence
				my $count = 0;
				$count = $seq =~ s/($motif)/$1/g;
				($count = 0) if (!$count);
	
				# add to total_count for that sequence
				$total_count += $count;
				$motif_count += $count;
				$counts[$l] += $count;			
			}
			
			# should remove any motif from motifset if it doesn't occur in any sequence at all
			if ($motif_count == 0){
				# add to an array to remove after processing all motifs in motifset
				push(@motifs_to_remove,$j)
				
			}
	
		}
#		print "Motifset $i contains $total_count total motifs\n";
		
		# can now calculate correlation (r) for that motifset
		my $r = r2(\@expression,\@counts);
		
		# add this to the hash so that we can later sort all motifsets by there score (and get rid of the lowest ones)
		$motifsets_to_scores{$i} = $r;
		
		# clean up
		# are there any motifs to remove from this motifset?
		# if this is the only motif in the motifset we can leave it as it will be selected against
		# otherwise just remove one motif with zero counts. Ideally should remove them all but this is getting
		# tricky to code. Remaining zero-count motifs will be cleaned in subsequent generations
		$number_of_motifs = scalar @{$motifset[$i]};
		if(@motifs_to_remove && ($number_of_motifs > 1)){
			splice(@{$motifset[$i]},$motifs_to_remove[0],1);
		}		
	}	
}



# similar routine but just for scoring a single motifset
sub score_motifset{
	my $i = shift;
	
	my $number_of_motifs = scalar @{$motifset[$i]};

	# total count will be the count of all motifs in each motifset
	my $total_count = 0;

	# need to store set of total counts for each sequence (to compare against expression values)
	my @counts;
	
	# may need to remove motifs from a motif set if they do not occur in any sequence
	# otherwise there is no selection pressure to remove them as they won't make any contribution 
	# to the correlation score
	my @motifs_to_remove = ();

	for my $j (0 .. $number_of_motifs-1){
		
		# separate counter just for each motif
		my $motif_count = 0;

		# get motif out of data structure
		my $number_of_bases = scalar @{$motifset[$i][$j]};

		#die "No bases: i = $i, j = $j, num_motifs $number_of_motifs\n" if ($number_of_bases == 0);
		
		# assemble motif as a string out of @motifset data
		my $motif;
		for my $k (0 .. $number_of_bases-1){
			$motif .= $bases_to_regex{$motifset[$i][$j][$k]{base}};
			$motif .= "{" . $motifset[$i][$j][$k]{min} . ",";
			$motif .= $motifset[$i][$j][$k]{max} . "}";
		}

		# check to see if motif starts or ends with an N, if so can remove this base (but only if motif has more than one base)
		if($motifset[$i][$j][0]{base} eq "N"){
			$number_of_bases = scalar @{$motifset[$i][$j]};
			splice(@{$motifset[$i][$j]},0,1) if ($number_of_bases > 1);								
		}
		$number_of_bases = scalar @{$motifset[$i][$j]};
		if($motifset[$i][$j][$number_of_bases-1]{base} eq "N" && ($number_of_bases > 1)){
			$number_of_bases = scalar @{$motifset[$i][$j]};
			splice(@{$motifset[$i][$j]},$number_of_bases-1,1);
		}

		# now loop through each sequence
		for (my $l=0; $l < @seqs; $l++){
			my $seq = $seqs[$l];
			
			my $id = $seq_to_id{$seq};
			
			# count motif in sequence
			my $count = 0;
			$count = $seq =~ s/($motif)/$1/g;
			($count = 0) if (!$count);

			# add to total_count for that sequence
			$total_count += $count;
			$motif_count += $count;
			$counts[$l] += $count;			
		}
	}	
	# can now calculate correlation (r) for that motifset
	my $r = r2(\@expression,\@counts);
	return($r);

}

#######################################################
#
#          C a l c u la t e   R  - s q u a r e d
#
#######################################################

sub r2{
	# two arrays of data to regress against each other
    my $x = shift;
    my $y = shift;

    # need a bunch of (fairly standard) of statistics from both arrays
    my $n       = @{$x};
    my $sum_x   = sum(@{$x});
    my $sum_y   = sum(@{$y});
    my $sum_xy  = sum(map {${$x}[$_] * ${$y}[$_]} (0..$n-1));
    my $sum_x2  = sum(map {${$x}[$_]**2} (0..$n-1));
    my $sum_y2  = sum(map {${$y}[$_]**2} (0..$n-1));

    # can't calculate r2 if sum_x or sum_y = 0; so return 0 instead
    if($sum_x == 0 || $sum_y == 0){
    	return("0");
    }
    # calculate r and return r2
    else{
        my $r = (($n * $sum_xy) - ($sum_x * $sum_y))/sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2));        
    	return(sprintf("%.8f",$r));                  
    }
}