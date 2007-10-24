#!/usr/bin/perl -w
#
# genbank_intron_analysis.pl
#
# a script to extract intron sizes for all genes in all organisms and average the size according to
# whether it is the 1st, 2nd, 3rd, nth intron as well as non-first and last categories. Then compare
# 1st intron sizes to other size categories
#
# by Keith Bradnam
#
# created 23rd October 2007 (but largely based on genbank_intron_size_vs_position.pl)
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
my $limit;     # set a cut-off for how many introns are needed per species (at first intron position)
my $input;	   # the name of the intron dump file from genbank_info_dump.pl
my $min;       # set min length for introns, any CDS with any intron below that length is ignored
my $summary;   # only print summary of size differences between first and non-first introns

GetOptions("limit=i"  => \$limit,
		   "input=s"  => \$input,
		   "min=i"    => \$min,
			"summary" => \$summary);

die "Must use -input option to specify a file containing intron info from genbank\n" if (!defined($input));

#set default for $limit if not specified
$limit  = 100   if (!defined($limit));

########################
# misc variables
########################

# simple hash for keeping track of all species: key is species name, values are just set to '1'
my %all_species;

# keep track of CDSs with introns shorter than value set by $min
my $short_counter = 0;

my %species2first_intron;
my %species2second_intron;
my %species2third_intron;
my %species2fourth_intron;
my %species2fifth_intron;
my %species2first_intron_multi;
my %species2last_intron;
my %species2last_intron_multi;
my %species2non_first_intron;
my %species2single_intron;




############################################################
#
#  M A I N   L O O P
#
#  process genbank dump file to put info into hashes
#
############################################################

open (IN,"<$input") || die "Can't open $input file\n";

LINE: while (<IN>) {
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
	
	# check all intron sizes if -ignore option specified, skip CDSs with short intron
	if($min){
		foreach my $i (@introns){
			if ($i < $min){
				$short_counter++;
				next LINE;
			}
		}
	}
	
	# might have occasional problem with some species names in genbank
	# names should be alphanumerical characters, spaces, and maybe a period and dash
	next if ($species !~ m/[\w\s\.\-]+/);
	
	# should always have at least four entries per line 
	next if ($length < 4);
	
	# keep track of species
	$all_species{$species} = 1;
	
	# now add data to relevant hashes
	my $intron = shift(@introns);
	push(@{$species2first_intron{$species}},$intron);
	
	# any more introns to come? If not add to last intron hash and single intron hash
	if(!@introns){
		push(@{$species2last_intron{$species}},$intron);	 
		push(@{$species2single_intron{$species}},$intron);	 
		next LINE;
	}
	
	# if we are here then there must be more introns
	# so can also add first intron to hash of first introns in multi-intron genes
	push(@{$species2first_intron_multi{$species}},$intron);

	# get 2nd intron
	$intron = shift(@introns);
	push(@{$species2second_intron{$species}},$intron);
	push(@{$species2non_first_intron{$species}},$intron);
	
	# any more introns to come? If not add to last intron hash (and last intron of multi-intron genes hash)
	if(!@introns){
		push(@{$species2last_intron{$species}},$intron);
		push(@{$species2last_intron_multi{$species}},$intron);	
		next LINE;
	}

	# get 3rd intron
	$intron = shift(@introns);
	push(@{$species2third_intron{$species}},$intron);
	push(@{$species2non_first_intron{$species}},$intron);
	
	# any more introns to come? If not add to last intron hash (and last intron of multi-intron genes hash)
	if(!@introns){
		push(@{$species2last_intron{$species}},$intron);
		push(@{$species2last_intron_multi{$species}},$intron);	
		next LINE;
	}
	
	# get 4th intron
	$intron = shift(@introns);
	push(@{$species2fourth_intron{$species}},$intron);
	push(@{$species2non_first_intron{$species}},$intron);
	
	# any more introns to come? If not add to last intron hash (and last intron of multi-intron genes hash)
	if(!@introns){
		push(@{$species2last_intron{$species}},$intron);
		push(@{$species2last_intron_multi{$species}},$intron);	
		next LINE;
	}

	# get 5th intron
	$intron = shift(@introns);
	push(@{$species2fifth_intron{$species}},$intron);
	push(@{$species2non_first_intron{$species}},$intron);
	
	# any more introns to come? If not add to last intron hash (and last intron of multi-intron genes hash)
	if(!@introns){
		push(@{$species2last_intron{$species}},$intron);
		push(@{$species2last_intron_multi{$species}},$intron);	
		next LINE;
	}
	
	# not counting individual introns past 5th so can treat remainder a little differently
	while(@introns){
		$intron = shift(@introns);
		push(@{$species2non_first_intron{$species}},$intron);
		if(!@introns){
			push(@{$species2last_intron{$species}},$intron);
			push(@{$species2last_intron_multi{$species}},$intron);
		}
	}

}

close(IN) || die "Couldn't close $input\n";

if($summary){
	print "Species,N,mean first intron size,mean non-first intron size,%diff,Significance\n";
}

# now loop through big data structure to get stats. Species first...
foreach my $species (sort(keys(%all_species))){
	
	# only want to look at species where there is some non-first intron data to compare with 1st intron
	# i.e. skip species if there are not enough second introns (>= $limit)
	next if (!defined(@{$species2second_intron{$species}}) || (@{$species2second_intron{$species}} < $limit));
	
	my ($n,$mean,$sd);
	my ($n1a,$mean1a,$sd1a);
	my ($n1m,$mean1m,$sd1m);
	my ($n2,$mean2,$sd2);
	my ($n3,$mean3,$sd3);
	my ($n4,$mean4,$sd4);
	my ($n5,$mean5,$sd5);
	my ($n_nf,$mean_nf,$sd_nf);
	my ($n_last,$mean_last,$sd_last);
	my ($n_single,$mean_single,$sd_single);
	
	# just print summary data (for excel) if required
	if($summary){
		if(defined(@{$species2first_intron_multi{$species}}) && (@{$species2first_intron_multi{$species}} >= $limit)){
			$n1m = @{$species2first_intron_multi{$species}};
			$mean1m = sprintf("%.2f",sum(@{$species2first_intron_multi{$species}})/$n1m);
			$sd1m = sqrt(sum(map {($_ - $mean1m) ** 2} @{$species2first_intron{$species}}) / ($n1m-1));

			$n_nf = @{$species2non_first_intron{$species}};
			$mean_nf = sprintf("%.2f",sum(@{$species2non_first_intron{$species}})/$n_nf);
			$sd_nf = sqrt(sum(map {($_ - $mean_nf) ** 2} @{$species2first_intron{$species}}) / ($n_nf-1));

			# what percentage of average first intron length are non-first introns
			my $percentage_diff = sprintf("%.2f",$mean_nf/$mean1m);

			my $z = sprintf("%.2f",($mean1m - $mean_nf)/sqrt(($sd1m**2/$n1m)+($sd_nf**2/$n_nf)));
			my $sig;
			$sig = "NS";
			$sig = "*"     if ($z >= 1.96);
			$sig = "**"    if ($z >= 2.58);
			$sig = "***"   if ($z >= 3.29);
			$sig = "****"  if ($z >= 3.89);
			$sig = "*****" if ($z >= 4.41);
			print "$species,$n1m,$mean1m,$mean_nf,$percentage_diff,$sig\n";
		}	
		next;
	}
	
	# proceed if not in summary mode
	print "\n=== $species ===\n";	

	if(defined(@{$species2first_intron{$species}}) && (@{$species2first_intron{$species}} >= $limit)){
		$n1a = @{$species2first_intron{$species}};
		$mean1a = sprintf("%.2f",sum(@{$species2first_intron{$species}})/$n1a);
		$sd1a = sqrt(sum(map {($_ - $mean1a) ** 2} @{$species2first_intron{$species}}) / ($n1a-1));
		print "1st introns (all): n = $n1a, mean = $mean1a nt\n";	
	}
	
	if(defined(@{$species2first_intron_multi{$species}}) && (@{$species2first_intron_multi{$species}} >= $limit)){
		$n1m = @{$species2first_intron_multi{$species}};
		$mean1m = sprintf("%.2f",sum(@{$species2first_intron_multi{$species}})/$n1m);
		$sd1m = sqrt(sum(map {($_ - $mean1m) ** 2} @{$species2first_intron{$species}}) / ($n1m-1));
		print "1st introns (multi-intron genes): n = $n1m, mean = $mean1m nt\n";		
	}

	if(defined(@{$species2second_intron{$species}}) && (@{$species2second_intron{$species}} >= $limit)){
		$n2 = @{$species2second_intron{$species}};
		$mean2 = sprintf("%.2f",sum(@{$species2second_intron{$species}})/$n2);
		$sd2 = sqrt(sum(map {($_ - $mean2) ** 2} @{$species2first_intron{$species}}) / ($n2-1));
		print "2nd introns: n = $n2, mean = $mean2 nt\n";		
	}

	if(defined(@{$species2third_intron{$species}}) && (@{$species2third_intron{$species}} >= $limit)){
		$n3 = @{$species2third_intron{$species}};
		$mean3 = sprintf("%.2f",sum(@{$species2third_intron{$species}})/$n3);
		$sd3 = sqrt(sum(map {($_ - $mean3) ** 2} @{$species2first_intron{$species}}) / ($n3-1));
		print "3rd introns: n = $n3, mean = $mean3 nt\n";		
	}
	
	if(defined(@{$species2fourth_intron{$species}}) && (@{$species2fourth_intron{$species}} >= $limit)){
		$n4 = @{$species2fourth_intron{$species}};
		$mean4 = sprintf("%.2f",sum(@{$species2fourth_intron{$species}})/$n4);
		$sd4 = sqrt(sum(map {($_ - $mean4) ** 2} @{$species2first_intron{$species}}) / ($n4-1));
		print "4th introns: n = $n4, mean = $mean4 nt\n";		
	}
	
	if(defined(@{$species2fifth_intron{$species}}) && (@{$species2fifth_intron{$species}} >= $limit)){
		$n5 = @{$species2fifth_intron{$species}};
		$mean5 = sprintf("%.2f",sum(@{$species2fifth_intron{$species}})/$n5);
		$sd5 = sqrt(sum(map {($_ - $mean5) ** 2} @{$species2first_intron{$species}}) / ($n5-1));
		print "5th introns: n = $n5, mean = $mean5 nt\n";		
	}
	
	if(defined(@{$species2non_first_intron{$species}}) && (@{$species2non_first_intron{$species}} >= $limit)){
		$n_nf = @{$species2non_first_intron{$species}};
		$mean_nf = sprintf("%.2f",sum(@{$species2non_first_intron{$species}})/$n_nf);
		$sd_nf = sqrt(sum(map {($_ - $mean_nf) ** 2} @{$species2first_intron{$species}}) / ($n_nf-1));
		# will often get data for this category but where there is not enough genes for 1st intron data,
		# so only want to print this data if there is something that we can compare it to (1st intron from
		# multi-intron genes)
		print "non-1st introns: n = $n_nf, mean = $mean_nf nt\n" if (@{$species2first_intron_multi{$species}} >= $limit);		
	}

	if(defined(@{$species2last_intron{$species}}) && (@{$species2last_intron{$species}} >= $limit)){
		$n = @{$species2last_intron{$species}};
		$mean = sprintf("%.2f",sum(@{$species2last_intron{$species}})/$n);
		print "Last introns: n = $n, mean = $mean nt\n";		
	}

	if(defined(@{$species2last_intron_multi{$species}}) && (@{$species2last_intron_multi{$species}} >= $limit)){
		$n_last = @{$species2last_intron_multi{$species}};
		$mean_last = sprintf("%.2f",sum(@{$species2last_intron_multi{$species}})/$n_last);
		$sd_last = sqrt(sum(map {($_ - $mean_last) ** 2} @{$species2first_intron{$species}}) / ($n_last-1));
		print "Last introns (multi-intron genes): n = $n_last, mean = $mean_last nt\n";		
	}

	if(defined(@{$species2single_intron{$species}}) && (@{$species2single_intron{$species}} >= $limit)){
		$n_single = @{$species2single_intron{$species}};
		$mean_single = sprintf("%.2f",sum(@{$species2single_intron{$species}})/$n_single);
		$sd_single = sqrt(sum(map {($_ - $mean_single) ** 2} @{$species2first_intron{$species}}) / ($n_single-1));
		print "Single intron genes: n = $n_single, mean = $mean_single nt\n";		
	}

	# now make z-tests
	my $z;
	print "\n";
	
	# 1st (multi) vs 2nd
	if(($n1m >= $limit) && ($n2 >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean2)/sqrt(($sd1m**2/$n1m)+($sd2**2/$n2)));
		print "Z-test 1st(multi) vs 2nd = $z\n";
	}
	# 1st (multi) vs 3rd
	if(($n1m >= $limit) && $n3 && ($n3 >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean3)/sqrt(($sd1m**2/$n1m)+($sd3**2/$n3)));
		print "Z-test 1st(multi) vs 3rd = $z\n";
	}
	# 1st (multi) vs 4th
	if(($n1m >= $limit) && $n4 && ($n4 >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean4)/sqrt(($sd1m**2/$n1m)+($sd4**2/$n4)));
		print "Z-test 1st(multi) vs 4th = $z\n";
	}
	# 1st (multi) vs 5th
	if(($n1m >= $limit) && $n5 && ($n5 >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean5)/sqrt(($sd1m**2/$n1m)+($sd5**2/$n5)));
		print "Z-test 1st(multi) vs 5th = $z\n";
	}
	# 1st (multi) vs non-first
	if(($n1m >= $limit) && $n_nf && ($n_nf >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean_nf)/sqrt(($sd1m**2/$n1m)+($sd_nf**2/$n_nf)));
#		print "Z-test 1st(multi) vs non-first = $z\n";
		if($z >= 2.58){
			print "Z-test 1st(multi) vs non-first = $z *\n";
		}
		else{
			print "Z-test 1st(multi) vs non-first = $z\n";
		}
	}
	# 1st (multi) vs last (multi introns)
	if(($n1m >= $limit) && $n_last && ($n_last >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean_last)/sqrt(($sd1m**2/$n1m)+($sd_last**2/$n_last)));	
		print "Z-test 1st(multi) vs last(multi) = $z\n";
	}
	# 1st (multi) vs single introns
	if(($n1m >= $limit) && $n_single && ($n_single >= $limit)){
		$z = sprintf("%.2f",($mean1m - $mean_single)/sqrt(($sd1m**2/$n1m)+($sd_single**2/$n_single)));
		print "Z-test 1st(multi) vs single introns = $z\n";
	}
}


print "\n$short_counter CDSs with at least one introns shorter than $min bp length were ignored\n" if ($min);
	
exit(0);
