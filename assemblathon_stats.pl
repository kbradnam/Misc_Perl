#!/usr/bin/perl
#
# assemblathon_stats.pl
#
# A script to calculate a basic set of metrics from a genome assembly
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use List::Util qw(sum max min);

###############################################
# 
#  C o m m a n d   l i n e   o p t i o n s
#
###############################################

my $limit; # limit processing of data to first $limit sequences (for quick testing)
my $graph; # produce some output ready for Excel or R

GetOptions ("limit=i" => \$limit,
			 "graph"  => \$graph);

# set defaults
$limit = 1000000000 if (!$limit);



# check we have a suitable input file
die "Usage: assemblathon_stats.pl <assembly_file>\n" unless (@ARGV == 1);
my ($file) = @ARGV;



###############################################
# 
#  S o m e   G l o b a l   v a r i a b l e s
#
###############################################

my $n_limit = 25;						# how many Ns are we using to split scaffolds into contigs?
my $contig_file = "tmp_contigs$$.fa";	# might need to create a temp output file for contigs
my $known_genome_size = 112500000;  	# the approximate genome size of species A
my $scaffolded_contigs = 0;				# how many contigs that are part of scaffolds (sequences must have $n_limit consecutive Ns)
my $scaffolded_contig_length = 0;		# total length of all scaffolded contigs
my $unscaffolded_contigs = 0;			# how many 'orphan' contigs, not part of a scaffold
my $unscaffolded_contig_length =0;		# total length of all contigs not part of scaffold
my $w = 50;								# formatting width for output
my %data;								# data structure to hold all sequence info key is either 'scaffold', 'contig' or intermediate', values are seqs & length arrays

	
	
	
# make first loop through file, capture some basic info and add sequences to arrays
process_FASTA($file);

print "\n<-- Information for assembly \'$file\' -->\n\n";

# produce scaffold statistics
sequence_statistics('scaffold');

# produce a couple of intermediate statistics based on scaffolded contigs vs unscaffolded contigs
sequence_statistics('intermediate');	

# finish with contig stats
sequence_statistics('contig');

exit(0);



##########################################
#
#
#    S  U  B  R  O  U  T  I  N  E  S
#
#
##########################################


##########################################
#    M A I N  loop through FASTA file
##########################################

sub process_FASTA{
	
	my ($seqs) = @_;
	
	open(my $input, "<", "$seqs") or die "Can't open $seqs\n";
	my $fasta = new FAlite(\*$input);

	# want to keep track of various contig + scaffold counts
	my $seq_count = 0;

	while(my $entry = $fasta->nextEntry){
	    my $seq = uc($entry->seq);
		my $length = length($seq);
		$seq_count++;

		# everything gets pushed to scaffolds array
		push(@{$data{scaffold}{seqs}},$seq);	
		push(@{$data{scaffold}{lengths}},$length);	
		
		# if there are not at least 25 consecutive Ns in the sequence we need to split it into contigs
		# otherwise the sequence must be a contig itself and it still needs to be put in @contigs array
		if ($seq =~ m/N{$n_limit}/){
			
			# add length to $scaffolded_contig_length
			$scaffolded_contig_length += $length;
			
			# loop through all contigs that comprise the scaffold
			foreach my $contig (split(/N{25,}/, $seq)){
				$scaffolded_contigs++;
				my $length = length($contig);				
				push(@{$data{contig}{seqs}},$contig);	
				push(@{$data{contig}{lengths}},$length);	
			}
		} else {
			# must be here if the scaffold is actually just a contig (or is a scaffold with < 25 Ns)
			$unscaffolded_contigs++;
			$unscaffolded_contig_length += $length;
			push(@{$data{contig}{seqs}},$seq);	
			push(@{$data{contig}{lengths}},$length);	
		}	
		# for testing, just use a few sequences
		last if ($seq_count >= $limit);
		
	}
	close($input);
}


##########################################
#    Calculate basic assembly metrics
##########################################

sub sequence_statistics{
	my ($type) = @_;
	
	print "\n";
	
	# there are just a couple of intermediate level statistics to print
	if($type eq 'intermediate'){
		my $total_size = sum(@{$data{scaffold}{lengths}});
		
		# now calculate percentage of assembly that is accounted for by scaffolded contigs
		my $percent = sprintf("%.1f",($scaffolded_contig_length / $total_size) * 100);
		printf "%${w}s %10s\n","Percentage of assembly in scaffolded contigs", "$percent%";		

		# now calculate percentage of assembly that is accounted for by scaffolded contigs
		$percent = sprintf("%.1f",($unscaffolded_contig_length / $total_size) * 100);
		printf "%${w}s %10s\n","Percentage of assembly in unscaffolded contigs", "$percent%";
		return();
	}
	
	
	# n
	my $count = scalar(@{$data{$type}{lengths}});
	printf "%${w}s %10d\n","Number of ${type}s", $count;
	
	
	# more contig details (only for contigs)
	if ($type eq 'contig'){
		printf "%${w}s %10d\n","Number of contigs in scaffolds",$scaffolded_contigs;
		printf "%${w}s %10d\n","Number of unscaffolded contigs",$unscaffolded_contigs;
	}


	# total size of sequences
	my $total_size = sum(@{$data{$type}{lengths}});
	printf "%${w}s %10d\n","Total size of ${type}s", $total_size;


	# For scaffold data only, can caluclate the percentage of known genome size 
	if ($type eq 'scaffold'){		
		my $percent = sprintf("%.1f",($total_size / $known_genome_size) * 100);
		printf "%${w}s %10s\n","Percentage of known genome size", "$percent%";		
	}
		
	
	# longest and shortest sequences
	my $max = max(@{$data{$type}{lengths}});
	my $min = min(@{$data{$type}{lengths}});	
	printf "%${w}s %10d\n","Longest $type", $max;
	printf "%${w}s %10d\n","Shortest $type", $min;
	
	
	# find number of sequences above certain sizes
	my %sizes_to_shorthand = (500	  => '500',
							  1000    => '1K',
							  10000   => '10K',
							  100000  => '100K',
							  1000000 => '1M');

	foreach my $size qw(500 1000 10000 100000 1000000){
		my $matches = grep { $_ > $size } @{$data{$type}{lengths}};
		my $percent = sprintf("%.1f", ($matches / $count) * 100);
		printf "%${w}s %10d %5s%%\n","Number of ${type}s > $sizes_to_shorthand{$size} nt", $matches, $percent;
	}
	
	
	# mean sequence size
	my $mean = $total_size / $count;
	printf "%${w}s %10d\n","Mean $type size", $mean;


	# median sequence size
	my $median = ${$data{$type}{lengths}}[$count/2];
	printf "%${w}s %10d\n","Median $type size", $median;
	
	
	
	##################################################################################
 	# 
	# N50 values
	#
	# Includes N(x) values, NG(x) (using known genome size)
	# and L(x) values (number of sequences larger than or equal to N50 sequence size)
	##################################################################################

	# keep track of cumulative assembly size (starting from smallest seq)
	my $running_total = 0;
	
	# want to store all N50-style values from N1..N100. First target size to pass is N1
	my $n_index = 1;
	my $ng_index = 1;
	my @n_values;
	my @ng_values;
	
	my $i = 0;
	foreach my $length (reverse sort{$a <=> $b} @{$data{$type}{lengths}}){
		$i++;
		$running_total += $length;

		# check the current sequence and all sequences shorter than current one
		# to see if they exceed the current NX value
		while($running_total > int (($n_index / 100) * $total_size)){	
			if ($n_index == 50){
				printf "%${w}s %10d\n","N50 $type length", $length;
				printf "%${w}s %10d\n","L50 $type count", $i;
			}
			$n_values[$n_index] = $length;
			$n_index++;
		}		

		# now do the same for NG values, using known genome size
		while($running_total > int (($ng_index / 100) * $known_genome_size)){	
			if ($ng_index == 50){
				printf "%${w}s %10d\n","NG50 $type length", $length;
				printf "%${w}s %10d\n","LG50 $type count", $i;				
			}
			$ng_values[$ng_index] = $length;
			$ng_index++;
		}		

	}
	# add final value to @n_values and @ng_values which will just be the shortest sequence
#	$n_values[100] = $min;
#	$ng_values[100] = $min;
	


	# base frequencies
	my %bases;

    my $seq = join('',@{$data{$type}{seqs}});  
	my $length = length($seq);

    # count mononucleotide frequencies
    $bases{A} = ($seq =~ tr/A/A/); 
    $bases{C} = ($seq =~ tr/C/C/); 
    $bases{G} = ($seq =~ tr/G/G/); 
    $bases{T} = ($seq =~ tr/T/T/); 
    $bases{N} = ($seq =~ tr/N/N/);
	
	my $base_count = 0;
	foreach my $base qw (A C G T N){
		my $percent = sprintf("%.2f", ($bases{$base} / $length) * 100);
		printf "%${w}s %10s\n","%$base", $percent;
		$base_count += $bases{$base};
	}
    # calculate remainder ('other) in case there are other characters present
	my $other = $length - $base_count;
	my $percent = sprintf("%.2f", ($other / $length) * 100);
	printf "%${w}s %10s\n","%other", $percent;


	# anything to dump for graphing?
	if($graph){
		my $file_name = "${type}_graph.csv";
		my $out;
		
		# append to file if it exists, else write new one
		if (-e $file_name){
			open($out, ">>", "$file_name") or die "Can't append to $file_name\n";			
		} else{
			open($out, ">", "$file_name") or die "Can't create $file_name\n";
			print $out join (',',"Assembly",1..99), "\n";
		}
		
		# CSV file, with filename in first column
		print $out "$file,";
		
		for (my $i = 1; $i < 100; $i++){
			# higher NG values might not be present if assembly is poor
			if (defined $ng_values[$i]){
				print $out "$ng_values[$i],";	
			} else{
				print $out "0,";
			}
		}
		print $out "\n";
		close($out);
	}
}

