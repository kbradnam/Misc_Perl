#!/usr/bin/perl
#
# regex_breeder.pl
#
# A script to breed sets of motifs (using regular expressions) and see which best explains
# known expression values for various introns
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use List::Util qw(sum);
use DataBrowser;

my $file;               # input file
my $n;                  # population size
my $generations;        # number of generations
my $max_motifs;         # maximum number of motifs in each motif set
my $min_r;              # minimum correlation value needed from motifsets in starting population
my $max_pattern_length; # maximum pattern length in regular expression
my $starting_motif;     # optionally seed one motif in starting generation 
my $exclude_motif;      # optionally specify a motif pattern that can never occur 
my $mortality;          # what fraction of motifsets die in 1st generation (mortality will increase)
my $untouchables;       # what fraction of motifsets in 1st generation should be saved from mutation (the untouchables)
my $verbose; 			# turn on extra output

GetOptions ("n=i"                   => \$n,
			"max_motifs=i"          => \$max_motifs,
			"file=s"                => \$file,
			"min_r=f"               => \$min_r,
			"max_pattern=i"         => \$max_pattern_length,
			"starting_motif=s"      => \$starting_motif,
			"exclude_motif=s"       => \$exclude_motif,
			"generations=i"         => \$generations,
			"mortality=f"           => \$mortality,
			"untouchables=f"        => \$untouchables,
			"verbose"               => \$verbose
			);


# set defaults
$n = 300                  if (!$n);
$generations = 1000       if (!$generations);
$max_motifs = 5           if (!$max_motifs);
$min_r = 0.25             if (!$min_r);
$max_pattern_length = 9   if (!$max_pattern_length);
$mortality = 0.33         if (!$mortality);
$untouchables = 0.25      if (!$untouchables);

die "Specify intron file with -file option\n" if (!$file);
die "Values for -mortality ($mortality) and -untouchables ($untouchables) options are incompatible. Must sum to < 1\n" if ((1 - $mortality) <= $untouchables);

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
my @bases = qw(A C G G T R Y S W K M B D H V N);
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
	N => '(A|C|G|T|N)'
);

# read sequence files and also extract expression values
read_file();

# one time set-up
print "Creating starting population\n";
create_starting_population($n);


#############################
#
#     M A I N   L O O P
#
#############################

my %motifsets_to_scores;

# want to keep track of how many generations in a row that the top and bottom motifsets
# have the same value of r (in which case we can stop the breeder)
my $r_watch = 0;

# will work with a copy of $mortality that will increase
my $current_mortality = $mortality;

for (my $i = 1; $i <= $generations; $i++){

	# decrease number of survivors each generation (but stop when you get to 5)
	# use a simple formula to decrease number of survivors in a linear fashion
	my $percent_completed = ($i / $generations);
	my $survivors = int(($n * $untouchables) * (1 - $percent_completed)) + 1;
	$survivors = 5 if ($survivors < 5);

	# mortality will increase with each generation
	$current_mortality = $percent_completed;

	# but set a minimum and maximum amount of mortality 
	$current_mortality = $mortality if ($current_mortality < $mortality);
	$current_mortality = 0.9 if ($current_mortality > 0.9);
	my $dead = int($current_mortality * $n);
	print "\n===================================================\n";
	print "GENERATION $i\n";
	print "N = $n, Untouchables = $survivors, Dead = $dead\n";
	print "===================================================\n";
	
	# mutate all motifsets (except any untouchables)
	mutate_copy($i);
	
	# wipe relevant hashes and arrays
	%motifsets_to_scores = ();

	# grab scores from each motifset
	score_motifsets();

	# sort scores
	my @sorted = (reverse sort {$motifsets_to_scores{$a} <=> $motifsets_to_scores{$b}} (keys(%motifsets_to_scores)));


	print "\nTop motifs\n----------\n";
	# print details of the top $survivors motifsets
	for (my $j = 0; $j < @sorted; $j++){
		if ($j < $survivors){
			# make these top scoring motifs untouchable
			$motifset[$sorted[$j]]{untouchable} = 1;			

			# only ever print details for top 5 scoring motifsets
			if($j < 5){
				print_motifset($sorted[$j]);
				print "\n";				
			}
		} else {
			last;
		}
	}

	# now just out motifset that is the lowest scoring but surviving
	my $offset = int((1 - $current_mortality) * $n);
	print "Last surviving motif (position $offset of $n)\n-------------------------------------------\n";
	print_motifset($sorted[$offset-1]);
	
	# is the highest scoring motifset producing the same r value as the lowest?
	# if so we will increment a counter, otherwise reset it
	if($motifset[$sorted[0]]{r} == $motifset[$sorted[$offset-1]]{r}){
		$r_watch++;
	} else{
		$r_watch = 0;
	}
	
	# If we have too many generations in a row where the top and bottom motifsets are giving
	# the same score of r, then we should stop the simulation
	if($r_watch == 5){
		die "\nStopping simulation after $i generations because all motifsets have had the same r value for 5 generations in a row\n";
	}
	
	# death - remove the $mortality percentage of the lowest scoring motifsets
	kill_motifsets(\@sorted, $i);
	print "\n";

	# reproduction, add new motifsets from survivors
	reproduction();
}

exit;



##########################################################################
#
# 
#              T  H  E     S  U  B  R  O  U  T  I  N  E  S
#  
#
###########################################################################

sub read_file{
	open(IN,"<$file") || die "Couldn't open $file file\n";

    my $fasta = new FAlite(\*IN);

    # loop through each sequence in target file
    while(my $entry = $fasta->nextEntry) {
 		# get and trim header
    	my $header = $entry->def;
    	$header =~ m/.* x(\d{1,2}\.\d)/;
    	my $expression = $1;

		die "Can't find expression value in: $header\n" if (!$expression);
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


###########################################################################
#
# CREATE_STARTING_POPULATION
#
# Create $n motifsets
#
###########################################################################

sub create_starting_population{
	my ($n) = @_;

	# will create $n motifset objects
	
	MOTIF: foreach my $i (0 .. $n - 1){

		# first decide how many motifs can be in each motifset
		my @motifs;
		my $number_of_motifs = int(rand($max_motifs)) + 1;

		foreach my $j (1 .. $number_of_motifs){
			my $motif = create_single_motif();
			push(@motifs, $motif);
		}
		
		# choose starting distance from TSS from which motifs will only be considered
		# initially make this choice quite simple
		my @distances = qw(100 250 500 1000 5000);
		my $distance = $distances[int(rand(@distances))];
		
		
		# choose whether motifs can be placed on the leading strand or on both strands
		my @strands = qw(1 2);
		my $strand = $strands[int(rand(@strands))];
		
		push @motifset, {
			motifs   => \@motifs,
			strand   => $strand,
			distance => $distance,
			untouchable => 0
		};
		
		
		# print motifs
#		print "\nmotifset $i\n";
#		print "Strand = $motifset[$i]{strand}\n";
#		print "Distance threshold = $motifset[$i]{distance}\n";

#		print "Number of motifs = $number_of_motifs\n";
		my $motif_count = 1;
		foreach my $motif (@{$motifset[$i]{motifs}}){
#			print "\tMotif $motif_count\n";
			# how to loop through motifs at this point?
			foreach my $pos (@{$motif}){
#				print "\t\t${$pos}{base}\{${$pos}{min},${$pos}{max}\}\n";
			}
			$motif_count++;
		}
	}
#	print "$motif->[0]{base}\n";

	# if -starting_motif has been specified, we overwrite first motif of motifset with a single motif
	if($starting_motif){

		# erase first motifset but add default values for strand (1) and distance (5000)
        $motifset[0] = ();
		$motifset[0]{distance} = 500;
		$motifset[0]{r} = 0;
		$motifset[0]{untouchable} = 1;
		$motifset[0]{strand} = 1;
		${$motifset[0]}{motifcounts}[0] = 0;
		
        my $counter = 0;
        my @motif = split(//,$starting_motif);
         
        while(@motif){
	     	my ($min, $max);                 
	        my $base = shift(@motif);
	        my $next_base;
	        if($motif[0]){
	        	$next_base = $motif[0];
	        } else{
	        	$next_base = "END";
	        }
		
	        # are both positions regular characters with no length pattern?
	        if ($base =~ m/[ACGTRYSWKMBDHVN]/ && $next_base =~ m/[ACGTRYSWKMBDHVN]/){               
            	# if so then just set min and max length to be equal to 1
                $motifset[0]{motifs}[0][$counter]{base} = $base;
                $motifset[0]{motifs}[0][$counter]{min}  = 1;
                $motifset[0]{motifs}[0][$counter]{max}  = 1;

                $counter++;
	        } else{
		
                shift(@motif);                  # get rid of '{'
                $min = shift(@motif);

				# now we might have a single value for min and max (e.g. Y{3} rather than Y{3,3})
				# treat accordingly
				
				if($motif[0] eq '}'){
					$max = $min;
				} else{
	                shift(@motif);                  # get rid of ','
	                $max = shift(@motif);
					
				}
                shift(@motif);                  # get rid of '}'        
        
                $motifset[0]{motifs}[0][$counter]{base} = $base;
                $motifset[0]{motifs}[0][$counter]{min}  = $min;
                $motifset[0]{motifs}[0][$counter]{max}  = $max;
               
                $counter++;
	         }
		}
	}
}

sub create_single_motif{
	my $motif_length = int(rand(6))+1;
	
	my @motif;
	foreach my $i (0 .. $motif_length){
		
		# choose a random base (using IUPAC codes)
		my $base = $bases[int(rand(@bases))];
		
		# decide minimum and maximum number of repetitions for that base
		
		# choose maximum number: between 1 and 3
		my $max = int(rand(3)) + 1;
		
		# min is anynumber up to the max
		my $min = int(rand($max)) + 1;
		
		push @motif, {
			base => $base,
			min  => $min,
			max  => $max,
		};
	}
	return(\@motif);
}

###########################################################################
#
# KILL_MOTIFSETS
#
# Empty a fraction of motifsets with the lowest r2 scores
#
###########################################################################

sub kill_motifsets{
	my ($sorted_motifsets, $i) = @_;
		
	print "\nKilling motifsets\n=================\n" if ($verbose);
	my $to_die = $n * $current_mortality;
	print "\n$to_die motifsets must die!\n\n" if ($verbose);

	my @doomed = @{$sorted_motifsets}[$n-$to_die..$n-1];
	for (my $s = 0; $s < @doomed; $s++){
		print "$s) Killing motifset $doomed[$s] r = $motifsets_to_scores{$doomed[$s]}\n" if ($verbose);
		$motifset[$doomed[$s]] = undef;
	}
}



###########################################################################
#
# REPRODUCTION
#
# Fill missing motifsets by randomly choosing from surviving motifsets
#
###########################################################################

sub reproduction{
	print "\nReproduction\n============\n" if ($verbose);
	
	# loop through all motifsets
	for my $i (0 .. $n-1){
		
		# is the current motifset dead?
		if(!defined($motifset[$i])){
			print "Motifset $i is dead\n" if ($verbose);
			print "\tfinding suitable replacement\n" if ($verbose);

			# stay in while loop until we find a suitable replacement motif
			while(1){
				my $rand = int(rand($n));
				if(defined($motifset[$rand])){
					print "\treplacing motifset $i with motifset $rand\n" if ($verbose);
					
#					browse($motifset[$rand]);

					# now copy over from an existing motifset					
					copy_motifset($rand,$i);
				
					# you might have copied an untouchable motifset, in which case we should make the copy
					# 'touchable', i.e. up for possibility of mutation
					$motifset[$i]{untouchable} = 0;					

					last;
				}
			}
		} else {
			if($motifset[$i]{untouchable} == 1){
				print "Motifset $i is alive (and untouchable)\n" if ($verbose);
				
			} else {
				print "Motifset $i is alive\n" if ($verbose);
				
			}
		}
	}
}

###########################################################################
#
# COPY_MOTIFSET
#
# To deal with references need to copy motifsets by traversing data structure
#
###########################################################################

sub copy_motifset{
	my ($old, $new) = @_;
	# will copy $old motifset to $new

	# copy over any hash key-value pairs at top level of motifset object
	${$motifset[$new]}{distance}    = ${$motifset[$old]}{distance};
	${$motifset[$new]}{strand}      = ${$motifset[$old]}{strand};
	${$motifset[$new]}{untouchable} = ${$motifset[$old]}{untouchable};
	${$motifset[$new]}{r}           = ${$motifset[$old]}{r};
	
	# now loop through all motifs in $old motifset
	for (my $m = 0; $m < @{$motifset[$old]{motifs}}; $m++){

		# and now loop through all bases in each motif of $old motifset and copy to $new motifset
		for (my $n = 0; $n < @{$motifset[$old]{motifs}[$m]}; $n++){
			$motifset[$new]{motifs}[$m][$n]{base} = $motifset[$old]{motifs}[$m][$n]{base};
			$motifset[$new]{motifs}[$m][$n]{min}  = $motifset[$old]{motifs}[$m][$n]{min};
			$motifset[$new]{motifs}[$m][$n]{max}  = $motifset[$old]{motifs}[$m][$n]{max};
		}		
	}
}



###########################################################################
#
# MUTATE_COPY
#
# Iterate through all untouchable motifsets and perform various mutations
#
###########################################################################

sub mutate_copy{
	my ($generation) = @_;
	
	print "\nMutating motifsets\n==================\n" if ($verbose);

	for (my $i = 0; $i < @motifset; $i++){
		print "$i)" if ($verbose);

		# skip the untouchables and the dead
		if (defined($motifset[$i]) && ($motifset[$i]{untouchable} == 1)){
			print " - this motifset is untouchable\n" if ($verbose);
			next;
		} else{
			print " - this motifset is up for mutation\n" if ($verbose);
		}
		
		# choose how many mutations to make (potentially all possible types of mutation could be applied to a single motifset)
		my $number_of_mutations = int(rand(12)+1);
		print "Motifset $i will undergo $number_of_mutations mutations\n" if ($verbose);
		
		# now loop through each mutation
		for(my $p = 0; $p < $number_of_mutations; $p++){

			# now choose which mutation will happen this time
			my $mutation_type = int(rand(12)+1);

			if($mutation_type == 1){			
				# 1) Add 1 new motifs to a motifset
				my ($number_of_motifs) = get_mutation_info($i);

				# can only add a motif if we have less than max number of $max_motifs
				if($number_of_motifs < $max_motifs){
					print "\tADDING MOTIF!\n" if ($verbose);				
					my $motif = create_single_motif();
					push @{$motifset[$i]{motifs}}, $motif;
				}
			} 
			
			elsif ($mutation_type == 2){			
				# 2) Remove 1 motif from motifset

				my ($number_of_motifs) = get_mutation_info($i);

				# can't remove motif if there is only one
				next if ($number_of_motifs == 1);
			
				print "\tREMOVING MOTIF!\n" if ($verbose);
			
				# choose which motif to remove
				my $rand = int(rand($number_of_motifs));
				splice(@{$motifset[$i]{motifs}},$rand,1);
			} 
			
			elsif($mutation_type == 3){			
				# 3) Add 1 base to one motif in motifset

				my ($number_of_motifs, $motif, $number_of_bases, $position) = get_mutation_info($i);

				print "\tADDING BASE TO MOTIF $motif AT POSITION $position\n" if ($verbose);
						
				# choose a random base (using IUPAC codes)
				my $base = $bases[int(rand(@bases))];

				# choose maximum number: between 1 and 3
				my $max = int(rand(3)) + 1;

				# min is any number up to the max
				my $min = int(rand($max)) + 1;

				my $new_pos = {
					base => $base,
					min  => $min,
					max  => $max
				};

				# add new base at randomly chosen position
				splice(@{$motifset[$i]{motifs}[$motif]},$position,0,$new_pos);
			}
		
			elsif($mutation_type == 4){			
				# 4) Remove 1 base from one motif in motifset
				my ($number_of_motifs, $motif, $number_of_bases, $position) = get_mutation_info($i);

				# can only do this if you have at least 2 bases
				if($number_of_bases >1){
					print "\tREMOVING BASE $position FROM MOTIF $motif\n" if ($verbose);

					# add details of new base to @copy
					splice(@{$motifset[$i]{motifs}[$motif]},$position,1);				
				}
			}
		
			elsif($mutation_type == 5){			
				# 5) Substitute 1 base from one motif in motifset
				my ($number_of_motifs, $motif, $number_of_bases, $position) = get_mutation_info($i);

				print "\tSUBSTITUTING BASE $position FROM MOTIF $motif\n" if ($verbose);

				# choose a random base (using IUPAC codes)
				my $base = $bases[int(rand(@bases))];

				# add details of new base to @copy
				$motifset[$i]{motifs}[$motif][$position]{base} = $base;				
			}

			elsif($mutation_type == 6){			
				# 6) Increment minimum pattern length for one base 
				my ($number_of_motifs, $motif, $number_of_bases, $position, $min, $max) = get_mutation_info($i);

				# increment value
				# can only do this if incremented value of min is less than max
				if ($min + 1 < $max){
					print "\tINCREASING MIN QUANTIFER FOR BASE $position IN MOTIF $motif\n" if ($verbose);
					$motifset[$i]{motifs}[$motif][$position]{min}++;
				}
			}

			elsif($mutation_type == 7){			
				# 7) Decrease minimum pattern length for one base
				my ($number_of_motifs, $motif, $number_of_bases, $position, $min, $max) = get_mutation_info($i);

				# decrease value, can only do this if value of min is greater than one
				if ($min > 1){
					print "\tDECREASING MIN QUANTIFER FOR BASE $position IN MOTIF $motif\n" if ($verbose);
					$motifset[$i]{motifs}[$motif][$position]{min}--;
				}
			}

			elsif($mutation_type == 8){			
				# 8) Increment maximum pattern length for one base
				my ($number_of_motifs, $motif, $number_of_bases, $position, $min, $max) = get_mutation_info($i);

				# increment value (as long as this doesn't exceed $max_pattern_length )
				if($max < $max_pattern_length){
					print "\tINCREASING MAX QUANTIFER FOR BASE $position IN MOTIF $motif\n" if ($verbose);
					$motifset[$i]{motifs}[$motif][$position]{max}++;				
				}
				die "$i $motif $position $min $max *\n" if (!defined($max) || !defined($min));
			}
		
			elsif($mutation_type == 9){			
				# 9) Decrease maximum pattern length for one base

				my ($number_of_motifs, $motif, $number_of_bases, $position, $min, $max) = get_mutation_info($i);

				# decrease value
				# can only do this if max will still be greater than $min
				if ($max - 1 > $min){
					print "\tDECREASING MAX QUANTIFER FOR BASE $position IN MOTIF $motif\n" if ($verbose);
					$motifset[$i]{motifs}[$motif][$position]{max}--;
				}
			}
			
			elsif($mutation_type == 10){			
				# 10) Increase threshold distance
			
				# incremement by varying amounts
				my @distances = (1, 5, 10, 100, 250, 500);
				my $increment = $distances[rand(@distances)];
				print "\tINCREASING DISTANCE THRESHOLD by $increment\n" if ($verbose);
				$motifset[$i]{distance} += $increment;
			}
		
			elsif($mutation_type == 11){			
				# 11) Decrease threshold distance
			
				# decrease by varying amounts
				my @distances = (1, 5, 10, 100, 250, 500);
				my $decrease = $distances[rand(@distances)];

				# can't decrease below a minimum (arbitraily chosen for now at 25 nt)
				if (($motifset[$i]{distance} - $decrease) >= 25){
					print "\tDECREASING DISTANCE THRESHOLD by $decrease\n" if ($verbose);
					$motifset[$i]{distance} -= $decrease;
				}
			}
		
			elsif($mutation_type == 12){			
				# 12) Change strand count
			
				# simply switch it from whatever it currently is
				if ($motifset[$i]{strand} == 1){
					print "\tINCREASING STRAND TO 2\n" if ($verbose);
					$motifset[$i]{strand} = 2;
				} else {
					print "\tDECREASING STRAND TO 1\n" if ($verbose);
					$motifset[$i]{strand} = 1;
				}
			}	
			else{
				# nothing happening
				# shouldn't come here I think
			}
		}
	}
}

###########################################################################
#
# GET_MUTATION_INFO
#
# Short subroutine to return some basic info that many of the mutation
# steps will repeatedly ask for
#
###########################################################################

sub get_mutation_info{
	# grab the motifset number
	my ($i) = @_;
	my $number_of_motifs = @{$motifset[$i]{motifs}}; 
	my $rand_motif = int(rand($number_of_motifs));
	my $number_of_bases = @{$motifset[$i]{motifs}[$rand_motif]};
	my $rand_position = int(rand($number_of_bases));
	my $min = $motifset[$i]{motifs}[$rand_motif][$rand_position]{min};
	my $max = $motifset[$i]{motifs}[$rand_motif][$rand_position]{max};

	# return details (which will only be for one randomly chosen base of one randomly chosen motif)
	return($number_of_motifs, $rand_motif, $number_of_bases, $rand_position, $min, $max);
}


###########################################################################
#
# SCORE_MOTIFSETS
#
# Calculate an r2 value for how well each motifset explains expression values
#
###########################################################################


sub score_motifsets{
	print "\nScoring motifsets\n=================\n" if ($verbose);
	for my $i (0 .. $n-1){
		
#		print "\nScoring motifset $i\n==================\n" ;
		my $number_of_motifs = @{$motifset[$i]{motifs}};
		
		# total count will be the count of all motifs in each motifset
		my $total_count = 0;
		
		# first want to trim the target sequences that we might be matching with all motifs. This will happen 
		# if we are only scoring motifs in the first X nt of each sequence
		# we will create a copy of @seqs to do this
		my @seqs_copy;
		my @rev_seqs_copy;
		
		# reset details of the survivors (i.e. the untouchables will be up for mutation if they are no longer top-scoring)
		$motifset[$i]{untouchable} = 0;
		
		# reset counts of all motifs from last generation
		$motifset[$i]{motifcounts} = ();
		
		my $distance_threshold = $motifset[$i]{distance};
		my $strand             = $motifset[$i]{strand};
#		print "Distance threshold = $distance_threshold nt\n";
#		print "Strand = $strand\n";
		
		# load up arrays with copies of sequences that will be used to match motifs
		# we have to make copies as we may be modifying them
		foreach my $seq (@seqs){
			my $copy = substr($seq, 0, $distance_threshold);
			push(@seqs_copy, $copy);
			if ($strand == 2){
				my $revcomp = Keith::revcomp($seq);
				push(@rev_seqs_copy, $revcomp);
			}
		}

		# need to store set of total counts for each sequence (to compare against expression values)
		my @counts;

		# now loop through the trimmed sequences (which we can modify if necessary)
		for (my $j = 0; $j < @seqs_copy; $j++){
#			print "\tSequence: $j\n";
			# keep track of all motifs found in this sequence
			my $all_motifs_in_seq = 0;

			my $motif_counter = 0;
			foreach my $motif (@{$motifset[$i]{motifs}}){

				# convert motif to a regex string
				my $motif_regex = motif_to_regex($motif);

				# count occurrence of regex in trimmed sequence and replace sequence with 6 Ns (to stop subsequent motifs matching in same spot)
				my ($count, $revcount) = (0, 0);
				$count = $seqs_copy[$j] =~ s/($motif_regex)/NNNNNN/g;
				($count = 0) if (!$count);
				# now consider reverse strand matches if necessary
				if ($strand == 2){
					$revcount = $rev_seqs_copy[$j] =~ s/($motif_regex)/NNNNNN/g;
					($revcount = 0) if (!$revcount);
				}

				# add counts from both strands to current seq total
				$all_motifs_in_seq += ($count + $revcount);
				
				# add current motif count to array in motifset object
				$motifset[$i]{motifcounts}[$motif_counter] += ($count + $revcount);
				
#				print "\t\tMotif $motif_counter ($motif_regex) occurs $count times (+ strand) and $revcount times (- strand)\n";
#				print "\t\tMotif $motif_counter occurs $count times (+ strand) and $revcount times (- strand) (running total = $motifset[$i]{motifcounts}[$motif_counter])\n";

				$motif_counter++;
			}
#			print "\tSeq $j: found $all_motifs_in_seq motifs in total in this sequence\n\n";

			# add counts to array 
			$counts[$j] = $all_motifs_in_seq;
			
		}		
		
		# can now calculate correlation (r) for that motifset
		my $r = r2(\@expression,\@counts, $i);
		
#		print "r is $r\n\n";
		
		# add this to a separate hash so that we can later sort all motifsets by there score (and get rid of the lowest ones)
		# also add it to the motifset hash (not using this for now, but might later on)
		$motifsets_to_scores{$i} = $r;
		$motifset[$i]{r} = $r;
		
		# As a final check, we want to see if any motifs don't occur at all in any of the sequences
		# If this is the case, they will not be contributing to the r2 value and we can remove the motif from the motifset		
		# only worth doing this if there is more than one motif, if there is just one motif, it will end up with a zero r2
		# score and will be selected out later on
		if ($number_of_motifs > 1){

			for(my $m = 0; $m < $number_of_motifs; $m++){

				# was the current motif not seen in any sequence?
				if ($motifset[$i]{motifcounts}[$m] == 0){
					print "REMOVING MOTIF $m as it does not occur in any sequence\n" if ($verbose);
					splice(@{$motifset[$i]{motifs}},$m,1);
					$number_of_motifs = @{$motifset[$i]{motifs}}; 
				}
			}
			
		}
	}	
}

###########################################################################
#
# PRINT_MOTIFSET
#
# subroutine to print out all of the details from a single motifset
#
###########################################################################

sub print_motifset{
	my ($i) = @_;

	my $distance_threshold = ${$motifset[$i]}{distance};
	my $strand             = ${$motifset[$i]}{strand};
	my $r                  = ${$motifset[$i]}{r};
	
	my $formatted_r = sprintf("%.4f", $r);

	my $number_of_motifs   =  @{$motifset[$i]{motifs}};

	print "id=$i: r = $formatted_r, motifs = $number_of_motifs, distance threshold = $distance_threshold nt, strands = $strand\n";

	for (my $m = 0; $m < $number_of_motifs; $m++){
		my $number_of_bases = @{$motifset[$i]{motifs}[$m]};
		my $motif_count = ${$motifset[$i]}{motifcounts}[$m];
		print "Motif $m:\t";
		printf 'n =%4s ', "$motif_count";
		print "\t";
		for (my $n = 0; $n < $number_of_bases; $n++){
			my $base = $motifset[$i]{motifs}[$m][$n]{base};
			my $min = $motifset[$i]{motifs}[$m][$n]{min};
			my $max = $motifset[$i]{motifs}[$m][$n]{max};
			
			# if min and max are both equal to 1 then no point printing the min and max
			# equally if min is equal to max, then can just print one value
			if($min == 1 && $max == 1){
				print "$base";				
			} elsif ($min == $max){
				print "$base\{$max\}";	
				
			} else {
				print "$base\{$min,$max\}";	
			}			
		}		
		print "\n";
	}
}

sub motif_to_regex{
	my ($motif) = @_;
	my $regex;
	foreach my $pos (@{$motif}){
#		print "\t\t${$pos}{base}\{${$pos}{min},${$pos}{max}\}\n";
		# convert each IUPAC code to a regex and add to a new string (along with repetition quantifiers)
		$regex .= ($bases_to_regex{${$pos}{base}} . "\{${$pos}{min},${$pos}{max}\}");
	}
	return($regex);
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
	my $i = shift;

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
		if ($sum_x2 == 0){
			print_motifset($i);
			die "sum_x2 is 0\n";
		}
		if ($sum_y2 == 0){
			print_motifset($i);
			die "sum_y2 is 0\n";
		}
		my $numerator = (($n * $sum_xy) - ($sum_x * $sum_y));
		my $denominator = sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2));
		
		# want to avoid divide by zero errors, so cheat a bit
		$denominator = 0.01 if ($denominator == 0);
		
        my $r = $numerator/$denominator;        
    	return(sprintf("%.8f",$r));                  
    }
}