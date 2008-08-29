#!/usr/bin/perl
#
# kl_shuffler.pl 
#
# a script to calculate the KL distance between two files (A & B) and then see  
# how many times this distance exceeds the distance of half of the sequences in A vs other half of sequences in A
# can run multiple iterations of shuffling sequences in A
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
kl_distance.pl - a program for calculating K-L distance between two FASTA files
and then comparing this distance with K-L distances derived from intra-file comparisons
of sequences in file 1.

Usage: <fasta file 1> <fasta file 2> <word size> <number of shuffles>\n" unless @ARGV == 4;


my ($FASTA1,$FASTA2,$WORD,$SHUFFLES) = @ARGV;

die "Please specify a word size of 10 or less\n" if ($WORD >10);



my $comparison_distance = kl($FASTA1,$FASTA2);
print "K-L distance between file 1 and file 2 = $comparison_distance\n";



# now read the first input file and store all the sequences in an array
# each element of array will be a sequence from input file
my @file1_seqs;
open(FILE,"$FASTA1") || die "Can't open $FASTA1\n";

my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file and add to array
while(my $entry = $fasta->nextEntry) {  
        push(@file1_seqs,uc($entry->seq));
}
close(FILE);


# count how many times distance is exceeded
my $c=0;

# need two tmp file names which won't conflict with anything else
my $tmp1 = "tmpseq1.fasta";
my $tmp2 = "tmpseq2.fasta";

# now make lots of pairs of files by randomly selecting sequence from @seqs
for (my $i=1;$i<=$SHUFFLES;$i++){
	# always want to keep original array unspoilt
	my @seqs = @file1_seqs;

	#now want two output files
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
	print "$i $distance";
	
	# has this exceeded first distance?
	if($distance > $comparison_distance){
		print "\tDISTANCE EXCEEDED\n";
		$c++;
	}
	else{
		print "\n";
	}

	# tidy up and remove tmp files
	unlink($tmp1);
	unlink($tmp2);
}

print "\ncomparison distance ($comparison_distance) exceeded $c times out of $SHUFFLES\n";



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
	
	# convert counts to frequencies
	my %freq;
	
	foreach my $word (keys %count){
		$freq{$word} = $count{$word}/$total_words;
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



__END__

