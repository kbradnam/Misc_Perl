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
my $report;     # only show matching regexes that occur $limit times
my $no_dimers;   # ignore palindromes which are all just a poly-dimer (e.g. ATAT*atata*TATA)

GetOptions ("filename=s"  => \$filename,
            "min_loop=i"  => \$min_loop,
            "max_loop=i"  => \$max_loop,
            "min_regex=i" => \$min_regex,
            "max_regex=i" => \$max_regex,             
            "report=i"    => \$report,
			"no_dimers"    => \$no_dimers
);

# set defaults 
$min_loop  = 0  if (!$min_loop);
$max_loop  = 15 if (!$max_loop);
$min_regex = 5  if (!$min_regex);
$max_regex = 20 if (!$max_regex);
$report    = 0  if (!$report);

my %palindromes;

# generate list of regexes
my @regexes;

for my $n ($min_regex..$max_regex){
	my $re = qr /[cagt]{$n}/;
	$regexes[$n-$min_regex] = $re;
}

# keep track of number of sequences in file & total number of palindromes
my $seq_count = 0;
my $palindrome_count = 0;
 
open (FASTA, "$filename") or die"Cannot open $filename\n";
my $fasta = new FAlite(\*FASTA);

while(my $entry = $fasta->nextEntry) {
	$seq_count++;
 	my $header = $entry->def;
    my $seq = lc($entry->seq);

	# loop through each regex
    foreach my $regex (@regexes){
		# loop through each match to the regex
		REGEX: while ($seq =~ m/$regex/g){
         	my $match = $&; # the match to the regex
			my $endline = $'; # the match of everything after regex
			
			# make opposite pattern to match (reverse & complement)
         	my $revmatch = reverse($match);
         	$revmatch =~ tr/cagt/gtca/;
			
			# now look for match to reverse complement sequence within remainder of $seq
			# allowing for 0-15 bases in the middle
         	if ($endline =~ /^([cagt]{$min_loop,$max_loop})($revmatch)/){
            	my $palindrome = uc($match) . "*" . $1 . "*" . uc($2);

				# want to be able to potentially filter out any palindrome
				# which is just a poly-dimer
				my $dimer_check = 0;
				
				if ($no_dimers){
					# use an all lower case version  of palindrome for testing
					my $tmp = $match . $1 . $2;
					foreach my $dimer qw(ac ag at ca cg ct ga gc gt ta tc tg){
						if ($tmp =~ m/^($dimer){2,50}$/){
							#print "poly-dimer: $palindrome\n"; 
							$dimer_check = 1;						
						}
					}
				}
	
            	$palindromes{$palindrome}++ unless ($dimer_check);
				$palindrome_count++ unless ($dimer_check); 
            }
		}
	}
}
close FASTA;


# now loop through hash printing out final counts of all palindromes
foreach my $key (sort {$palindromes{$a} <=> $palindromes{$b}} keys (%palindromes)) {
	print "$key => $palindromes{$key}\n" if ($palindromes{$key} > $report);
}

my $number_of_palindromes = scalar(keys(%palindromes));
my $palindromes_per_sequence = sprintf("%.2f",($palindrome_count/$seq_count));

print "$filename contains $number_of_palindromes different palindromes ($palindrome_count total) from $seq_count sequences";
print ": $palindromes_per_sequence palindromes per sequence. ";
if ($no_dimers){
	print "poly-dimer check ON\n";
}
else{
	print "poly-dimer check OFF\n";
}

exit;
