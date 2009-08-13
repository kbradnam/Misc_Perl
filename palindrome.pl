#!/usr/bin/perl
#
# palindrome.pl
#
# A script to find palindromes in a DNA sequence
# modified from an original script written by Jules J. Berman 2002
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;

########################
# Command line options
########################

my $filename;
my $min_loop;
my $max_loop;
my $min_regex; # minimum size of regular expression to search for
my $max_regex; 
my $allow_dimers;   # allow palindromes which are all just a poly-dimer (e.g. ATAT*atata*TATA), ignored by default
my $gc_pair_min; # specify a minimum number of G-C pairs in the palindrome unit, useful if you are expecting these to form stem loop
my $summary; 	 # print out summary statistics for all sequences in file
my $individual;  # print out statistics for each sequence in file
my $show_match;  # show matching palindromes
my $print_pre_match; # print FASTA output of sequences that occur before the regex
my $print_match;     # print FASTA output of sequences of the matching palindrome (two matches + loop)
my $print_post_match; # print FASTA output of sequences that occur after the regex

GetOptions ("filename=s"    => \$filename,
            "min_loop=i"    => \$min_loop,
            "max_loop=i"    => \$max_loop,
            "min_regex=i"   => \$min_regex,
            "max_regex=i"   => \$max_regex,             
			"allow_dimers"  => \$allow_dimers,
			"gc_pair_min=i" => \$gc_pair_min,
			"summary"       => \$summary,
			"individual"    => \$individual,
			"show_match"    => \$show_match,
			"print_pre_match" => \$print_pre_match,
			"print_match" => \$print_match,
			"print_post_match" => \$print_post_match
);

# set defaults 
$min_loop  = 2   if (!$min_loop);
$max_loop  = 25  if (!$max_loop);
$min_regex = 6   if (!$min_regex);
$max_regex = 20  if (!$max_regex);
$gc_pair_min = 0 if (!$gc_pair_min);

die "Usage: palindrome.pl -filename <filename> <options>\n" if (!$filename);

# generate list of regexes
my @regexes;

for my $n ($min_regex..$max_regex){
	my $re = qr /[cagt]{$n}/;
	$regexes[$n-$min_regex] = $re;
}
# reverse array so that we will start with the largest regex first
@regexes = reverse(@regexes);


# keep track of number of sequences in file & total number of palindromes, plus number of redundant (removed) matches
my $seq_count = 0;
my $total_palindrome_count = 0;
my $total_seq_length = 0;
my $redundant_counter = 0;

open (FASTA, "$filename") or die"Cannot open $filename\n";
my $fasta = new FAlite(\*FASTA);

while(my $entry = $fasta->nextEntry) {

	# need to also count palindromes in each sequence
	my $palindrome_count = 0;
	
	$seq_count++;
 	my $header = $entry->def;
	$header =~ s/^>//;
	
    my $seq = lc($entry->seq);
	$total_seq_length += length($seq);
	
	my ($distance_to_TSS) = $header =~ m/.*_i\d+_(\d+)_\w+/;

	# keep track of all regexes that match each sequence, might want to eliminate duplicates (e.g. ACAGT vs ACAGTA) later on
	my @regexes_in_seq = ();

	# loop through each regex
    foreach my $regex (@regexes){
		# loop through each match to the regex
		REGEX: while ($seq =~ m/$regex/g){
         	my $match = $&; # the match to the regex

			# count how man G or C nucleotides there are in the match
			my $gc_pair_count = $match =~ tr/[cg]/[cg]/;

			# skip to next REGEX if there are not enough
			next REGEX if ($gc_pair_count < $gc_pair_min);

			# store some more details of the match
			my $endline = $'; # the match of everything after regex
			my $match_position = length($`) + 1; # store the position of the match within the sequence
			my $distance_to_match_from_TSS = $distance_to_TSS + $match_position;
			
			# make opposite pattern to match (reverse & complement)
         	my $revmatch = reverse($match);
         	$revmatch =~ tr/cagt/gtca/;
			my $loop;
			
			# now look for match to reverse complement sequence within remainder of $seq
			# allowing for 0-15 bases in the middle
         	if ($endline =~ /^([cagt]{$min_loop,$max_loop})($revmatch)/){
				my ($loop, $match_2) = ($1, $2);

            	my $palindrome = uc($match) . "*" . $loop . "*" . uc($match_2);

				# format a version to use if -print_match is specified
				my $print_palindrome = lc($palindrome);
				$print_palindrome =~ s/\*//g;
				
				# want to be to ignore any palindrome which is just a poly-dimer (unless -allow_dimers specified)
				unless($allow_dimers){
					# use an all lower case version  of palindrome for testing
					my $tmp = $match . $loop . $match_2;
					foreach my $dimer qw(ac ag at ca cg ct ga gc gt ta tc tg){
						# jump to next REGEX if the palindrome is all dimers
						next REGEX if ($tmp =~ m/^($dimer){2,1000}$/);
					}
				}				
				push(@regexes_in_seq, $match);
				# if we have more than one matching regex in this sequence, check that they are not redundant, if so, just keep the longest 
				if (@regexes_in_seq > 1){
					my $redundant = 0;
					($redundant, @regexes_in_seq) = remove_redundant_regexes(@regexes_in_seq);
					# skip to next regex if there was a redundant one
					next REGEX if ($redundant > 0);
				}
				
				# if we are here, we can count palindrome details and print if necessary
				$palindrome_count++;
				$total_palindrome_count++; 
				print "$total_palindrome_count) $header $distance_to_match_from_TSS $match_position ",length($loop),"\t$palindrome\n" if ($show_match);
				print ">$header\n$print_palindrome\n" if ($print_match);
				# get pre and post match for the palindrome
				$seq =~ m/$print_palindrome/;
				my $pre_match = $`;
				my $post_match = $';
				print ">$header\n$pre_match\n" if ($print_pre_match);
				print ">$header\n$post_match\n" if ($print_post_match);
				
            }
		}
	}

	my $number_of_palindromes = $palindrome_count;
	my $palindromes_per_1000nt = sprintf("%.2f",($palindrome_count/(length($seq)/1000)));
	print "$header Contains $number_of_palindromes palindromes, $palindromes_per_1000nt palindromes per 1000 nt\n" if ($individual && $number_of_palindromes >0);
}
close FASTA;

# print out overall stats?
if($summary){
	my $palindromes_per_1000nt = sprintf("%.2f",($total_palindrome_count/($total_seq_length/1000)));

	print "$filename contains $total_palindrome_count palindromes from $seq_count sequences ($total_seq_length nt total)";
	print ": $palindromes_per_1000nt palindromes per 1000 nt of sequence. ";
	print "$redundant_counter redundant palindromes were removed";
	if ($allow_dimers){
		print "all-dimer palindromes allowed\n";
	}
	else{
		print "\n";
	}

}
exit;

# subroutine to remove redundant (overlapping) regular expressions
# i.e. if ACGTTGA .. TCAACGT is a palindrome, then we can ignore the shorter CGTTGA .. TCAACG palindrome
sub remove_redundant_regexes{
	my @regexes = @_;
	# flag variable to say whether we encountered a redundant regex
	my $redundant = 0;
	OUTER: for (my $i = 0; $i < @regexes; $i++){
		INNER: for (my $j = $i+1; $j < @regexes; $j++){

			# remove redundant (shorter) match if two regexes match
			if ($regexes[$i] =~ m/$regexes[$j]/){
				$redundant++;
				$redundant_counter++;
				splice(@regexes,$j,1);
				redo(OUTER);
			}
		}
	}
	return($redundant,@regexes);
}
