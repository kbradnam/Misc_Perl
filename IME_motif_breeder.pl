#!/usr/bin/perl
#
# IME_r_squared.pl
#
# A script to take a set of motifs found by Nested MICA, and test them against introns that have
# known expression values in Arabidopsis
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;
use List::Util qw(sum);

########################
# Command line options
########################

my $verbose;    # show details for motif count/density in each sequence
my $cycles;		# how many cycles to run for
my $n;			# how many motifs
my $mortality;  # what fraction of motifs will die each cycle
my $limit;      # what threshold do we start outputting successful motifs (in NestedMica format)
my $ime_data;   # file with experimental intron data in
my $threshold;  # choose a different threshold to count motifs
my $min;		# minimum size to generate motifs
my $max;		# maximum size to generate motifs
my $moredeath;  # mortality rate autoincrements
my $print_count;# can seed the counter of motifs which will appear in the names of motifs in XMS files

GetOptions ("verbose"       => \$verbose,
			"cycles=i"      => \$cycles,
			"n=i"           => \$n,
			"mortality=f"   => \$mortality,
			"limit=f"       => \$limit,
			"ime_data=s"    => \$ime_data,
			"threshold=f"   => \$threshold,
			"min=i"		    => \$min,
			"max=i"		    => \$max,
			"moredeath"     => \$moredeath,
			"print_count=i" => \$print_count
);

# set some defaults
$cycles = 100 if (!$cycles);
$n = 100 if (!$n);
$mortality = 0.9 if (!$mortality);
$limit = 0.8 if (!$limit);
$threshold = 3 if (!$threshold);
$min = 4 if (!$min);
$max = 15 if (!$max);
$print_count = 0 if (!$print_count);

# in -moredeath mode, mortality starts v. low but increases by 1% per generation
if($moredeath){
	print "Using -moredeath, mortality rate will autoincrement each generation\n";
	$cycles = 100;
	$mortality = 0.25;
}

die "Please specify a file with IME intron sequences in with -ime_data option\n" if (!$ime_data);


#############################
# global variables
#############################

# motifs will be stored in an array of array of hashes
# 1st array is motif number, 2nd array is position within motif
# final hash has keys which are A,C,G, or T, and values which are nucleotide frequencies
my @motifs;

# how many motifs to generate
my $num_motifs = $n;

# keep track of best score
my $r2_max = 0;

# Need sets of expected nucleotide frequencies to compute log likelihood scores
my %expected = ("A" => "0.2713","C" => "0.1534", "G" => "0.1701","T" => "0.4051");

# key is motif number, value is r2 score of motif
my %motif2score;

# keep track of motifs that survive each generation, these will be spared from mutation
# motif index is key, value is 1 (true)
my %survivors;

##############################################################
##############################################################
#
#    S T A R T 
#
##############################################################
##############################################################

my @seqs;
my @expression_scores;

&read_sequences;

print "-- Generating $num_motifs Motifs --\n" if ($verbose);
&generate_motifs;

for(my $i = 0; $i < $cycles; $i++){
	print "\n=== Cycle $i ===\n";
	print "-- Scoring Motifs --\n" if ($verbose);
	&score_motifs($i);
	print "-- Killing Motifs --\n" if ($verbose);
	&death;
	print "-- Mutating Motifs --\n" if ($verbose);
	&mutate;
	print "-- Growing & Shrinking Motifs --\n" if ($verbose);
	&grow_and_shrink;
	print "-- Recombining Motifs --\n" if ($verbose);
	&recombine;	
}

# Final scoring cycle
print "-- Scoring Motifs --\n" if ($verbose);
&score_motifs($cycles);
exit(0);


#############################################
#
#
#          S u b r o u t i n e s 
#
#
##############################################



# simple routine to generate a small number of motifs
# initial base frequencies chosen randomly

sub generate_motifs{

	for(my $i = 0; $i < $num_motifs; $i++){

		# motifs will be random length between $min and $max nt
		my $rand = int(rand($max-$min))+$min;

		for (my $j = 0; $j < $rand; $j++){

			# set random frequencies for each base (they must sum to 1.0)
			# all frequencies calculated to 5 dp
			my $a = sprintf("%.4f",rand(1));
			my $c = sprintf("%.4f",rand(1) * (1 - $a));
			my $g = sprintf("%.4f",rand(1) * (1 - $a - $c));
			my $t = sprintf("%.4f",1 - $a - $c -$g);

			# Remove minus signs should they appear (-0.0000 will sometimes appear)
			($a =~ s/-//) if ($a =~ m/-/);
			($c =~ s/-//) if ($c =~ m/-/);
			($g =~ s/-//) if ($g =~ m/-/);
			($t =~ s/-//) if ($t =~ m/-/);
			
			$motifs[$i][$j]{"A"} = $a;
			$motifs[$i][$j]{"C"} = $c;
			$motifs[$i][$j]{"G"} = $g;
			$motifs[$i][$j]{"T"} = $t;
		}
	}	
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

	# calculate r and return r2
	# can't calculate r2 if sum_x or sum_y = 0; so return a token low value instead
	if($sum_x == 0 || $sum_y == 0){
		return("0.0666");
	}
	else{
		my $r = (($n * $sum_xy) - ($sum_x * $sum_y))/sqrt( (($n * $sum_x2) - $sum_x**2) * (($n * $sum_y2) - $sum_y**2)); 	
		return(sprintf("%.5f",$r**2));
	}
}


########################################################################
#
#          L o a d   I n t r o n   d a t a   f r o m    f i l e 
#
#########################################################################
sub read_sequences{
	
	open(DATA,"<$ime_data") || die "Couldn't open $ime_data file\n";

	my $fasta = new FAlite(\*DATA);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {

		# get and trim header
		my $header = $entry->def;
		$header =~ m/.* ([\d\.]+) /;
		my $expression = $1;
		$header =~ s/ .*//;
	    my $seq = uc($entry->seq);
	
		# add intron sequences expression values to arrays
		push(@seqs,$seq);
		push(@expression_scores,$expression);
	}
	close(DATA);	
}

##############################################################################################
#
# calculate r2 score for all motifs
#
##############################################################################################

sub score_motifs{

	# will need to know what cycle we are in if a high scoring motif is found
	my $cycles = shift;
	
	# first loop over all potential motifs
	for(my $i = 0; $i < $num_motifs; $i++){


		# for each motif we evaluate, we will track of motif count in each of the 12 introns
		my @counts;
		
		# get motif length
		my $motif_length = @{$motifs[$i]};
		
		# will track final r2 score for each motif
		my $motif_r2;
						
		# loop through each of the 12 experimental introns, scoring them for motif		
		for (my $j = 0; $j < @seqs; $j++){
			# tidy up intron sequence to remove spaces
			my $intron = $seqs[$j];
			$intron =~ s/\s+//g;
			my $intron_length = length($intron);
			
			# will count motifs above threshold
			my $motif_count = 0;


		    # loop through intron sequence in windows equal to motif size
		    for (my $k = 0; $k < $intron_length-$motif_length + 1; $k++) {

				# extract a window of sequence, split it, and place in array	
				my @window = split(//,substr($intron, $k, $motif_length));

				# Calculate motif score for current window of sequence by looping over each base in window
				my $score = 0;
				for(my $l = 0; $l < @window; $l++){
					
					my $base = $window[$l];
					
					my $motif_base_frequency = $motifs[$i][$l]{$base};
					my $expected_base_frequency = $expected{$base};
					
					# should observed frequencies go below 0 or above 1 then we cheat
					($motif_base_frequency = 0.0001)   if ($motif_base_frequency <= 0);
					($motif_base_frequency = 1)        if ($motif_base_frequency >1);

					my $log = log($motif_base_frequency/$expected_base_frequency);
					$score += $log;

				}
				# add to motif count if we've found some sequence that is above the threshold...
				if ($score > $threshold){
					$motif_count++;
				}

			}

			# add these values to array for later stats analysis
			push(@counts,$motif_count);
		}
		my $count_r2 = &r2(\@expression_scores,\@counts);
		if($count_r2 > $r2_max){
			$r2_max = $count_r2;
			print "\n* New r2 max: motif $print_count = $r2_max, length $motif_length\n";
			
			# print out motif if score is above $limit
			if($r2_max > $limit){
				&print_motif($cycles,$i,$r2_max,$motif_length);
			}

		}		
		# add r2 to %motif2score, this should replace existing value if one is there
		$motif2score{$i} = $count_r2;
	}
}


##############################################################################################
#
# kill weak motifs, clone good motifs to fill remaining spots
#
##############################################################################################

sub death{
	
	# first clear suvivors hash from last generation
	%survivors = ();
	
	# mortality rate must increase if using -moredeath mode
	# but make sure that it stops at 0.95
	$mortality += 0.01 if ($moredeath && ($mortality != 0.95));
	my $formatted = sprintf("%.2f",$mortality);
	print "\nMortality rate = $formatted\n\n" if ($moredeath);
	
	# sort the r2 values and take the best proportion (as defined by using the $mortality factor)
	# first calculate the max number of motifs to keep.
	my $dead = $num_motifs * $mortality;
	my $alive = $num_motifs - $dead;
	print "Mortality factor = $mortality, so $dead motifs must die and $alive will live!\n" if ($verbose);

	my @sorted = (reverse sort {$motif2score{$a} <=> $motif2score{$b}} (keys(%motif2score)));
		
	# now work through list of motifs replacing those that will die because of low scores with 
	# replacement motifs with high scores 
	my $count = 0;
	foreach my $i (@sorted){
		if ($count < $alive){
			my $motif_length = @{$motifs[$sorted[$count]]};
			print "$count) $i - $motif2score{$i} - length $motif_length\n" if ($count<10);
			print "...\n" if ($count == 10);
			print "$count) $i - $motif2score{$i} - length $motif_length\n" if ($count>$alive-10);

			# add to survivors hash
			$survivors{$sorted[$count]} = 1;
		}
		else{
			#print "$count) $i - $motif2score{$i} - DIES\n" if ($verbose);
			# now need to swap this replace this array element with one that lives			
			#first randomly choose one of the surviving motif numbers 
			# the index number of @sorted for these will be between 0 and $alive
			my $rand = int(rand($alive));
			
			# now replace @motifs element $i with $sorted[$rand]
			&clone($i,$sorted[$rand])
		}
		$count++;
	}
	print "\n" if ($verbose);

}


##############################################################################################
#
# mutate motif base frequencies
#
##############################################################################################

sub mutate{
	
	for(my $i = 0; $i < $num_motifs; $i++){
		next if ($survivors{$i});
		
		my $motif_length = @{$motifs[$i]};
 		my $rand = rand(1);
		
		# 20% motifs (not including survivors) acquire very small mutations (0.0001% change in nt frequencies)
		# loop through each position in motif
		
		if($rand <=0.20){
				for (my $j = 0; $j < $motif_length; $j++){
				# send base frequencies to change_bases subroutine. The first value will specify what the extent
				# of any mutation will be
				my @new_base_freqs = &change_bases('0.0001',$motifs[$i][$j]{"A"},$motifs[$i][$j]{"C"},$motifs[$i][$j]{"G"},$motifs[$i][$j]{"T"});
				$motifs[$i][$j]{"A"} = $new_base_freqs[0];
				$motifs[$i][$j]{"C"} = $new_base_freqs[1];
				$motifs[$i][$j]{"G"} = $new_base_freqs[2];
				$motifs[$i][$j]{"T"} = $new_base_freqs[3];
			}
		}		
		# Another 20% of motifs acquire additional medium mutations (0.001% change in nt frequencies)
		elsif($rand <=0.4){
			# loop through each position in motif
			for (my $j = 0; $j < $motif_length; $j++){
				# send base frequencies to change_bases subroutine. The first value will specify what the extent
				# of any mutation will be
				my @new_base_freqs = &change_bases('0.001',$motifs[$i][$j]{"A"},$motifs[$i][$j]{"C"},$motifs[$i][$j]{"G"},$motifs[$i][$j]{"T"});
				$motifs[$i][$j]{"A"} = $new_base_freqs[0];
				$motifs[$i][$j]{"C"} = $new_base_freqs[1];
				$motifs[$i][$j]{"G"} = $new_base_freqs[2];
				$motifs[$i][$j]{"T"} = $new_base_freqs[3];
			}
		}
		# Another 20% of motifs acquire major mutations (0.01% change in nt frequencies)
		elsif($rand <=0.6){			
			# loop through each position in motif
			for (my $j = 0; $j < $motif_length; $j++){
				# send base frequencies to change_bases subroutine. The first value will specify what the extent
				# of any mutation will be

				my @new_base_freqs = &change_bases('0.01',$motifs[$i][$j]{"A"},$motifs[$i][$j]{"C"},$motifs[$i][$j]{"G"},$motifs[$i][$j]{"T"});		
				$motifs[$i][$j]{"A"} = $new_base_freqs[0];	
				$motifs[$i][$j]{"C"} = $new_base_freqs[1];
				$motifs[$i][$j]{"G"} = $new_base_freqs[2];
				$motifs[$i][$j]{"T"} = $new_base_freqs[3];
			}
		}
		# Another 20% of motifs acquire very major mutations (0.1% change in nt frequencies)
		elsif($rand <=0.8){			
			# loop through each position in motif
			for (my $j = 0; $j < $motif_length; $j++){
				# send base frequencies to change_bases subroutine. The first value will specify what the extent
				# of any mutation will be

				my @new_base_freqs = &change_bases('0.1',$motifs[$i][$j]{"A"},$motifs[$i][$j]{"C"},$motifs[$i][$j]{"G"},$motifs[$i][$j]{"T"});		
				$motifs[$i][$j]{"A"} = $new_base_freqs[0];	
				$motifs[$i][$j]{"C"} = $new_base_freqs[1];
				$motifs[$i][$j]{"G"} = $new_base_freqs[2];
				$motifs[$i][$j]{"T"} = $new_base_freqs[3];
			}
		}
		# This means 20% of motifs shouldn't have mutations, but these still have a chance to grow, shrink or recombine
	}	
}
##############################################################################################
#
# change base frequencies up or down
#
##############################################################################################

sub change_bases{
	my $mutagenicity = shift;
	
	my %bases;
	$bases{'A'} = shift;
	$bases{'C'} = shift;
	$bases{'G'} = shift;
	$bases{'T'} = shift;
	
	my $sum = $bases{'A'} + $bases{'C'} +  $bases{'G'} + $bases{'T'};
	
	my $base_count = 0;
	my $running_total = 0;
	
	# shuffle order of bases
	my @bases = qw (A C G T);
	for (my $i = @bases; --$i; ) {
	 	my $j = int rand ($i+1);
	    next if $i == $j;
	    @bases[$i,$j] = @bases[$j,$i];
	 }
	
	# keep track of which way the mutations go (useful for debugging)
	my %mutations;
	
	# loop through each base
	LOOP: foreach my $key (@bases){
		$base_count++;
	
#		$sum = $bases{'A'} + $bases{'C'} +  $bases{'G'} + $bases{'T'};	
#		print "$base_count)     A = $bases{'A'}\tC = $bases{'C'}\tG = $bases{'G'}\tT = $bases{'T'}\tSUM = $sum\n";

		# if this is 4th base then frequency is already decided (all bases must sum to 1)
		if($base_count == 4){
			($bases{$key} = sprintf("%.4f",1 - $bases{'C'} - $bases{'G'} - $bases{'T'})) if ($key eq 'A');
			($bases{$key} = sprintf("%.4f",1 - $bases{'A'} - $bases{'G'} - $bases{'T'})) if ($key eq 'C');
			($bases{$key} = sprintf("%.4f",1 - $bases{'A'} - $bases{'C'} - $bases{'T'})) if ($key eq 'G');
			($bases{$key} = sprintf("%.4f",1 - $bases{'A'} - $bases{'C'} - $bases{'G'})) if ($key eq 'T');

			# occasionally need to remove the negative sign for -0.00000 values
			($bases{$key} = "0.0000") if ($bases{$key} eq "-0.0000");
		#	$sum = $bases{'A'} + $bases{'C'} +  $bases{'G'} + $bases{'T'};
		#	print "$base_count) $key   A = $bases{'A'}\tC = $bases{'C'}\tG = $bases{'G'}\tT = $bases{'T'}\tSUM = $sum\n\n";

			# now exit the loop, can't risk mutating the last base now
			last LOOP;
		}
		
		# randomly decide whether base will mutate upwards or downwards in frequency or stay the same
		my $rand = rand(1);
		my $direction = "plus";
		$direction = "same"  if ($rand <0.666);
		$direction = "minus" if ($rand <0.333);


		# can't mutate to above 1 or below 0
		if($direction eq "plus"){
			$mutations{$key} = "+";
			if (($bases{$key} + $mutagenicity + $running_total) >1){
				$bases{$key} = sprintf("%.4f",1 - $running_total);
			}
			else{
				$bases{$key} = sprintf("%.4f",$bases{$key} + $mutagenicity)	;				
			}
			$running_total += $bases{$key};
		}
		elsif($direction eq "minus"){
			$mutations{$key} = "-";
			if (($bases{$key} - $mutagenicity) <0){
				$bases{$key} = 0;
			}
			else{
				$bases{$key} = sprintf("%.4f",$bases{$key} - $mutagenicity);				
				# have to scale back if we have exceeded 1
				($bases{$key} = sprintf("%.4f",1 - $running_total)) if ($bases{$key} + $running_total > 1);
			}
			$running_total += $bases{$key};
			
		}	
		else{
			$mutations{$key} = "=";
			# have to scale back if we have exceeded 1
			($bases{$key} = sprintf("%.4f",1 - $running_total)) if ($bases{$key} + $running_total > 1);
			$running_total += $bases{$key};
		}
		

	#	$sum = $bases{'A'} + $bases{'C'} +  $bases{'G'} + $bases{'T'};
	#	print "$base_count) ${key}$mutations{$key}  A = $bases{'A'}\tC = $bases{'C'}\tG = $bases{'G'}\tT = $bases{'T'}\tSUM = $sum\t$running_total\n\n";
		
	}

	die "A $bases{'A'}\n" if ($bases{'A'} =~ m/-/);
	die "C $bases{'C'}\n" if ($bases{'C'} =~ m/-/);
	die "G $bases{'G'}\n" if ($bases{'G'} =~ m/-/);
	die "T $bases{'T'}\n" if ($bases{'T'} =~ m/-/);

	return($bases{'A'},$bases{'C'},$bases{'G'},$bases{'T'});
}


# simple subroutine to replace one element of @motifs with another
sub clone{
	my $replace = shift;
	my $keep = shift;
	
	my $motif_length = @{$motifs[$keep]};
	
	# wipe the array element that we are going to replace
	$motifs[$replace] = ();
	
	for(my $j=0;$j<$motif_length;$j++){
		foreach my $base qw(A C G T){
			$motifs[$replace][$j]{$base} = $motifs[$keep][$j]{$base};
		}
	}	
}


##############################################################################################
#
# grow and shrink motifs
#
##############################################################################################

sub grow_and_shrink{
	
	for(my $i = 0; $i < $num_motifs; $i++){
		next if ($survivors{$i});
		
		my $motif_length = @{$motifs[$i]};
 		my $rand = rand(1);

		# 15% of motifs gain a base
		if($rand <=0.15){
			# choose random position to add base at
			$rand = int(rand($motif_length));
			# copy bases after this position along one
			for(my $j = $motif_length;$j>$rand;$j--){
				$motifs[$i][$j]{"A"} = $motifs[$i][$j-1]{"A"};
				$motifs[$i][$j]{"C"} = $motifs[$i][$j-1]{"C"};
				$motifs[$i][$j]{"G"} = $motifs[$i][$j-1]{"G"};
				$motifs[$i][$j]{"T"} = $motifs[$i][$j-1]{"T"};
			}
			
			# Now set random frequencies for new base (they must sum to 1.0) at position $rand
			my $a = sprintf("%.4f",rand(1));
			my $c = sprintf("%.4f",rand(1) * (1 - $a));
			my $g = sprintf("%.4f",rand(1) * (1 - $a - $c));
			my $t = sprintf("%.4f",1 - $a - $c -$g);

			$motifs[$i][$rand]{"A"} = $a;
			$motifs[$i][$rand]{"C"} = $c;
			$motifs[$i][$rand]{"G"} = $g;
			$motifs[$i][$rand]{"T"} = $t;					
		}
		# A different 15% of motifs will lose a base
		elsif($rand <=0.3 && $motif_length >1){
			
			# choose random position to remove
			$rand = int(rand($motif_length));
			splice(@{$motifs[$i]},$rand,1);
			$motif_length = @{$motifs[$i]};		
		}
	}
}


##############################################################################################
#
# recombine motifs
#
##############################################################################################

sub recombine{
	
	for(my $i = 0; $i < $num_motifs; $i++){
		next if ($survivors{$i});
		
		my $motif_length = @{$motifs[$i]};
 		my $rand = rand(1);

		# 10% of motifs swap a position with another motif
		if($rand <=0.1){
			
			# choose motif to swap with, can't be self or survivor
			my $partner = 0;
			while($partner == 0){
				$partner = int(rand($num_motifs));
				$partner = 0 if ($partner == $i);
				$partner = 0 if ($survivors{$partner});
			}
			
			# choose random positions to swap base in donor and partner motif
			my $motif1_pos = int(rand($motif_length));
			my $partner_motif_length = @{$motifs[$partner]};
			my $motif2_pos = int(rand($partner_motif_length));
			
			# need to store partner base frequencies in tmp variables before overwriting
			my ($a,$c,$g,$t);
			$a = $motifs[$partner][$motif2_pos]{'A'};
			$c = $motifs[$partner][$motif2_pos]{'C'};
			$g = $motifs[$partner][$motif2_pos]{'G'};			
			$t = $motifs[$partner][$motif2_pos]{'T'};

			# now overwrite partner motif with position from current motif (motif $i)
			$motifs[$partner][$motif2_pos]{'A'} = $motifs[$i][$motif1_pos]{'A'};
			$motifs[$partner][$motif2_pos]{'C'} = $motifs[$i][$motif1_pos]{'C'};
			$motifs[$partner][$motif2_pos]{'G'} = $motifs[$i][$motif1_pos]{'G'};
			$motifs[$partner][$motif2_pos]{'T'} = $motifs[$i][$motif1_pos]{'T'};
			
			# now replace current motif from tmp frequencies
			$motifs[$i][$motif1_pos]{'A'} = $a;
			$motifs[$i][$motif1_pos]{'C'} = $c;
			$motifs[$i][$motif1_pos]{'G'} = $g;
			$motifs[$i][$motif1_pos]{'T'} = $t;
		}
	}
}

##############################################################################################
#
# print motif - write to a NestedMica XMS file
#
##############################################################################################

sub print_motif{
	my $cycle = shift;
	my $i = shift;
	my $r2 = shift;
	my $length = shift;
	my $date = `date`;
	
	my $seqs = @seqs;
		
	open(OUT,">${print_count}_c${cycle}_l${length}_r2_${r2}_${seqs}_seqs.xms") || die "Couldn't create output file\n";
	print OUT "<motifset xmlns=\"http://biotiffin.org/XMS/\">\n";
	print OUT "  <prop>\n";
	print OUT "    <key>creator.name</key>\n";
	print OUT "    <value>nminfer</value>\n";
	print OUT "  </prop>\n";
	print OUT "  <prop>\n";
	print OUT "    <key>creator.version</key>\n";
	print OUT "    <value>0.8.0</value>\n";
	print OUT "  </prop>\n";
	print OUT "  <prop>\n";
	print OUT "    <key>input</key>\n";
	print OUT "    <value>IME_motif_breeder.pl</value>\n";
	print OUT "  </prop>\n";
	print OUT "  <prop>\n";
	print OUT "    <key>date</key>\n";
	print OUT "    <value>$date</value>\n";
	print OUT "  </prop>\n";
	print OUT "  <motif>\n";
	print OUT "    <name>IME Breeder $print_count</name>\n";
	print OUT "    <weightmatrix alphabet=\"DNA\" columns=\"$length\">\n";
	
	for (my $k =0; $k<$length;$k++){
		my $a = sprintf("%.5f",$motifs[$i][$k]{'A'});
		my $c = sprintf("%.5f",$motifs[$i][$k]{'C'});
		my $g = sprintf("%.5f",$motifs[$i][$k]{'G'});
		my $t = sprintf("%.5f",$motifs[$i][$k]{'T'});
	
		print OUT "      <column pos=\"$k\">\n";
		print OUT "        <weight symbol=\"adenine\">$a</weight>\n";
		print OUT "        <weight symbol=\"cytosine\">$c</weight>\n";
		print OUT "        <weight symbol=\"guanine\">$g</weight>\n";
		print OUT "        <weight symbol=\"thymine\">$t</weight>\n";
		print OUT "      </column>\n";
	}
	print OUT "    </weightmatrix>\n";
	print OUT "    <threshold>0.0</threshold>\n";
	print OUT "    <prop>\n";
	print OUT "      <key>nmica.history_thread</key>\n";
	print OUT "      <value>412</value>\n";
	print OUT "    </prop>\n";
	print OUT "  </motif>\n";
	print OUT "</motifset>\n";
	close(OUT);
	
	# increment count of printed motif
	$print_count++;
	
}
