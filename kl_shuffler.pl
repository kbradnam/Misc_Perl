#!/usr/bin/perl
#
# kl_shuffler.pl 
#
# a script to calculate the KL distance between two files (A & B) and then:
# 1) calculate the *intra* sequence KL distance for file A (randomly split sequences in file A into two new files)
# 2) calculate the *intra* sequence KL distance for file B
# 3) calculate the distance between a shuffled version of A, and a shuffled version of B
# 4) repeat steps 1-3 N times
# 5) see how many times the intra KL distances exceed the known distance

# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
kl_distance.pl - a program for calculating K-L distance between two FASTA files
and then comparing this distance with K-L distances derived from intra-file comparisons
of sequences in both files.

Usage: <fasta file 1> <fasta file 2> <word size> <number of shuffles>\n" unless @ARGV == 4;


my ($FASTA1,$FASTA2,$WORD,$SHUFFLES) = @ARGV;
die "Please specify a word size of 10 or less\n" if ($WORD >10);



# First get the KL distance *BETWEEN* the two files
my $comparison_distance = kl($FASTA1,$FASTA2);
print "\nInter-sequence K-L distance between file 1 and file 2 using word size $WORD = $comparison_distance\n\n";



# will store sequences from both files in array to save time later on
my (@seqs_from_1st_file,@seqs_from_2nd_file);
my $file_counter=0;


# Now need to calculate the the KL distance *WITHIN* each of the two files.

foreach my $file ($FASTA1,$FASTA2){
	$file_counter++;
	
	# now read the first input file and store all the sequences in an array
	# each element of array will be a sequence from input file
	my @seqs_from_file;
	
	open(FILE,"$file") || die "Can't open $file\n";

	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file and add to array
	while(my $entry = $fasta->nextEntry) {  
	        push(@seqs_from_file,uc($entry->seq));
	}
	close(FILE);

	# need count how many times distance is exceeded
	my $c=0;

	# copy sequence arrays for later on
	if($file_counter ==1){
		@seqs_from_1st_file = @seqs_from_file;
	}
	else{
		@seqs_from_2nd_file = @seqs_from_file;
	}

	# now make lots of pairs of files by randomly selecting sequence from @seqs
	for (my $i=1;$i<=$SHUFFLES;$i++){
		# always want to keep original array unspoilt
		my @seqs = @seqs_from_file;

		#now want two output files
		# need two tmp file names which won't conflict with anything else
		my $tmp1 = "tmpseq1.$i";
		my $tmp2 = "tmpseq2.$i";
		
		open(OUT1,">$tmp1") || die "can't create $tmp1\n";
		open(OUT2,">$tmp2") || die "can't create $tmp2\n";

		while(@seqs){
			# choose a random position
			my $rand1 = int(rand(1)*@seqs);

			# write output files, removing sequence from @seqs array 
			print OUT1 ">$i\n$seqs[$rand1]\n";
			splice(@seqs,$rand1,1);

			# repeat for 2nd randomly chosen sequence (if there is still something in @seqs)
			if(@seqs){
				my $rand2 = int(rand(1)*@seqs);
				print OUT2 ">$i\n$seqs[$rand2]\n";
				splice(@seqs,$rand2,1);
			}
		}
		close OUT1;
		close OUT2;

		# now calculate kl distance from pair of new files
		my $distance = kl("$tmp1","$tmp2");
		#print "$i $distance";

		# has this exceeded first distance? If so, increment counter
		$c++ if ($distance > $comparison_distance);

		# tidy up and remove tmp files
		unlink($tmp1);
		unlink($tmp2);
	}

	print "Inter-sequence distance exceeded $c times out of $SHUFFLES intra-sequence comparisons for file $file\n";
}


# need new counter
my $c = 0;

for (my $i=1;$i<=$SHUFFLES;$i++){
	
	# need local copies of arrays for each shuffling run
	my @seqs1 = @seqs_from_1st_file;
	my @seqs2 = @seqs_from_2nd_file;
	
	# need two tmp file names which won't conflict with anything else
	my $tmp1 = "tmpseq1.fasta";
	my $tmp2 = "tmpseq2.fasta";

	# now want two output files
	open(OUT1,">$tmp1") || die "can't create $tmp1\n";
	open(OUT2,">$tmp2") || die "can't create $tmp2\n";

	# loop through each sequence in both input files
	# these have been stored in @seqs_from_1st_file and @seqs_from_2nd_file arrays

	while(my $seq = shift(@seqs1)) {
		print OUT1 ">$i\n";
        # randomize sequence
        my $random = randomize($seq);
        print OUT1 "$random\n";
	}
	while(my $seq = shift(@seqs2)) {
		print OUT2 ">$i\n";
        # randomize sequence
        my $random = randomize($seq);
        print OUT2 "$random\n";
	}

	close OUT1;
	close OUT2;
	
	# now calculate kl distance from pair of new files
	my $distance = kl("$tmp1","$tmp2");

	# has this exceeded first distance? If so, increment counter
	$c++ if ($distance > $comparison_distance);

	# tidy up and remove tmp files
	unlink($tmp1);
	unlink($tmp2);
}

print "Inter-sequence distance exceeded $c times out of $SHUFFLES when comparing shuffled versions of both files\n\n";


exit(0);




################################################

sub kl{
	my $file1 = shift;  
	my $file2 = shift;

	# get frequencies of words in each sequence file
	my $freq1     = &frequency_table($file1,$WORD);
	my $freq2     = &frequency_table($file2,$WORD);

	# calculate reciprocal K-L distances
	my $distance1 = &kl_distance($freq1,$freq2);
	my $distance2 = &kl_distance($freq2,$freq1);

	# take average of both distances
	my $distance = ($distance1 + $distance2) /2;

	# convert score to bits
	$distance /= log(2);
	return($distance);
}





sub frequency_table{
	my ($file,$w) = @_;

	# keep track of all words
	my $total_words = 0;

	# open file to read sequence and break into words
	# then count work frequencies to return values

	open(FILE,"$file") || die "Can't open $file\n";

	my $fasta = new FAlite(\*FILE);

	# keep track of each word, initially set to 1 to avoid zero errors
	my %count = &word_table($w,1);
	
	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {  
	        my $seq = uc($entry->seq);
			
			# loop through sequence in windows equal to motif size
		    for (my $i = 0; $i < length($seq)-$w+1; $i++) {
				
				# extract a window of sequence, split it, and place in array    
		    	my $word = substr($seq, $i, $w);
		
				# check that word is a valid sequence word, count that words and increment total counter
				if (exists $count{$word}){
					$count{$word}++;
					$total_words++;
				}
	        }
	}
	close(FILE);
	
	# check (and warn) if too many words do not exist in the input sequence files
	# it is useful to know if many of the different possible words only exist as pseudocounts
	# Also need to convert counts to frequencies
	my %freq;
	
	# counter for how many words only exist as pseudocounts
	my $only_pseudocounts = 0;
	
	foreach my $word (keys %count){
		$freq{$word} = $count{$word}/$total_words;
		$only_pseudocounts++ if ($count{$word} == 1);
	}

	my $total_keys = keys(%count);
	my $percentage = sprintf("%.0f",$only_pseudocounts/$total_keys *100);
	
	# Only want to warn if a certain proportion of words only exist with counts of 1 (only pseudocounts)
	# The K-L distance may be less meaningful if it based on many comparisons of 1 vs 1 counts
	# The threshold level (10%) is arbitrary but should help to give a clue as to whether using a smaller word size would be more approrpriate

	my $threshold = 10;
	
	if ($only_pseudocounts && ($percentage >= $threshold)){
		print "\nWARNING: in file $file, $only_pseudocounts out of a possible $total_keys words ($percentage%) do not exist at all\n";
		print "Maybe consider using a smaller word size in order to calculate a more reliable K-L distance\n"; 		
	}

	return \%freq;
}

#####################################################
# pseudocount all words
#####################################################

sub word_table{
	my ($w,$pseudocount) = @_;
	my @alphabet = qw(A C G T);
	
	my %table;
	
	# populate tables with pseudocounts, do this for each possible word size
	if($w == 1){
		foreach my $c1 (@alphabet){
			$table{"$c1"} = $pseudocount;
		}
	}
	
	elsif($w == 2){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				$table{"$c1$c2"} = $pseudocount;
			}
		}
	}
	
	elsif($w == 3){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					$table{"$c1$c2$c3"} = $pseudocount;
				}
			}
		}
	}
	
	elsif($w == 4){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						$table{"$c1$c2$c3$c4"} = $pseudocount;
					}
				}
			}
		}
	}
	
	elsif($w == 5){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							$table{"$c1$c2$c3$c4$c5"} = $pseudocount;
						}
					}
				}
			}
		}
	}
	
	elsif($w == 6){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								$table{"$c1$c2$c3$c4$c5$c6"} = $pseudocount;
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 7){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									$table{"$c1$c2$c3$c4$c5$c6$c7"} = $pseudocount;
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 8){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										$table{"$c1$c2$c3$c4$c5$c6$c7$c8"} = $pseudocount;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 9){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										foreach my $c9 (@alphabet){
											$table{"$c1$c2$c3$c4$c5$c6$c7$c8$c9"} = $pseudocount;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 10){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										foreach my $c9 (@alphabet){
											foreach my $c10 (@alphabet){
												$table{"$c1$c2$c3$c4$c5$c6$c7$c8$c9$c10"} = $pseudocount;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return %table;	
}

###############################################################


sub kl_distance{
	my ($p,$q) = @_;
	
	my $distance = 0;
	
	foreach my $w (keys %{$p}){
		$distance += $p->{$w} * log ($p->{$w} / $q->{$w});
	}
	return $distance;
}



sub randomize{
    my $seq = shift;
    my $length = length($seq);
    
    # split sequence into an array
    my @seq = split(//,$seq);
    
    # want a random sequence array to store new sequence
    my @random;

    # loop through input sequence, randomly moving k-mers into a new randomized sequence
    while(@seq){
    
		 # choose random coordinate in sequence
		 my $rand = int(rand(1) * @seq);
        
		# add chosen random nt to second array
		push(@random,$seq[$rand]);    
		
	    # and remove from @seq
	    splice(@seq,$rand,1);
    }
    my $random = join('',@random);
    
    return $random;
}

__END__

