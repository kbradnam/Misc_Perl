#!/usr/bin/perl
#
# kmer_frequencies.pl 
#
# A script to just count and print frequencies of any specified Kmer in a FASTA file
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
kmer_frequencies.pl - a script to count the frequency of any kmer in a FASTA file
Usage: <fasta file 1> <word size>\n" unless @ARGV == 2;


my ($FASTA,$WORD) = @ARGV;

die "Please specify a word size of 10 or less\n" if ($WORD >10);


# get frequencies of words in each sequence file
frequency_table($FASTA,$WORD);


exit(0);

################################################
#
#
#   T H E   S U B R O U T I N E S
#
#
################################################

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
	
	foreach my $word (sort(keys %count)){
		$freq{$word} = $count{$word}/$total_words;
		my $frequency =  sprintf("%.4f",$freq{$word});
		print lc($word),"\t$frequency\n";
		$only_pseudocounts++ if ($count{$word} == 1);
	}

	my $total_keys = keys(%count);
	my $percentage = sprintf("%.0f",$only_pseudocounts/$total_keys *100);
	
	# Only want to warn if a certain proportion of words only exist with counts of 1 (only pseudocounts)
	# The threshold level (5%) is arbitrary but should help to give a clue as to whether using a smaller word size would be more approrpriate

	my $threshold = 5;
	
	if ($only_pseudocounts && ($percentage >= $threshold)){
		print "\nWARNING: in file $file, $only_pseudocounts out of a possible $total_keys words ($percentage%) do not exist at all\n";
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




__END__

