#!/usr/bin/perl
#
# get_trace_reads.pl
#
# A script to download and process reads from the NCBI trace archive
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

$SIG{'INT'} = 'INT_handler';


my @taxa = ("Arabidopsis thaliana","Xenopus laevis","Drosophila melanogaster","Homo sapiens");
#my @taxa = ("Arabidopsis thaliana","Drosophila melanogaster");

# script name that does the actual fetching of data (supplied by NCBI)
my $prog = "query_tracedb.pl";

# minimum number of bases that you need in a sequence after clipping to keep it
my $min_bases = 25;

# keep track of how many bases are clipped and how many traces are rejected for too many low quality bases
# do this for all species (grand total), and keep separate species totals
my $total_bases = 0;
my $total_clipped_bases = 0;
my $rejected_traces = 0;
my $species_total_bases;
my $species_clipped_bases;
my $species_rejected_traces;

# need an output file name for each species that we work with
my $file;


my $date = `date`; 
chomp($date);

SPECIES: foreach my $species (@taxa){

	print "\n------------------------------------------------------\n$date\n\n";
	
	# reset species-level counters
	$species_total_bases = 0;
	$species_clipped_bases = 0;
	$species_rejected_traces = 0;
	
	# store some of the query criteria in a separate string as we will want to reuse it later
	my $query = "species_code = '$species' and (trace_type_code = 'WGS' or trace_type_code = 'WCS' or trace_type_code = 'SHOTGUN' or trace_type_code = '454' or trace_type_code = 'CLONEEND') and source_type = 'GENOMIC'";

	# use 'query count' commands to get count of how many records match $query
	my $count = `$prog query count \"$query\"`;
	$count =~ s/\s+//g;
		

	# want to just check (and report) how many records are being filtered out at this first stage
	my $count_all = `$prog query count \"species_code = '$species'"`;
	$count_all =~ s/\s+//g;
	

	# open a FASTA output file for all reads for that species
	my $species_file_name = $species;
	$species_file_name =~ s/\s+/_/g;


	# can only download 40,000 records at a time
	my $pages = int($count/40000)+1;

	# there might be no records after above filtering and so we need to skip to next species
	if($count == 0){
		print "${species_file_name}: 0 reads remain after filtering (from $count_all possible reads)\n";
		next SPECIES;
	}
	else{
		my $percent_retained = sprintf("%.1f",($count/$count_all)*100);
		print "${species_file_name}: Attempting to fetch $count reads out of $count_all possible (${percent_retained}%), $pages pages\n";		
	}

	# if we have any sequences to fetch then need to make an output file
	open(OUT,">${species_file_name}_trace_reads.fa") or die "Can't open output file for $species\n";

	
	$pages = 4;
	
	for (my $i=0;$i<$pages;$i++){
		$file = "${species_file_name}_trace_read_data${i}";
		my $command = "(echo -n \"retrieve_gz fasta xml 0b\"\; $prog \"query page_size 40000 page_number $i binary $query\") | $prog > ${file}.gz"; 
#		my $command = "(echo -n \"retrieve fasta xml 0b\"\; $prog \"query page_size 40000 page_number $i binary $query\") | $prog > ${file}"; 
	
		print "${species_file_name}: Processing page ",$i+1,"/$pages\n";
		#print  "$command\n";

		# run the command to grab the sequences and xml files
		system("$command") && die "Can't execute $command\n";

		# now want to unpoack the file
		system("gunzip ${file}.gz") && die "Can't unzip archive\n";

		# process sequences based on info in xml file
		process_file($file);
	}
	close(OUT);
	
	#########################
	#
	#
	# FTP code to go here
	#
	#
	#########################
	
	my $percent_clipped = sprintf("%.1f",($species_clipped_bases/$species_total_bases)*100);
	print "$species_file_name: Processed $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "$species_file_name: $species_rejected_traces traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n";
	
	$date = `date`; 
	chomp($date);
	print "\n$date\n------------------------------------------------------\n";
}

# Final print out of stat
my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);
print "\n\n======================================================\n\n";
print "TOTAL: Processed $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
print "TOTAL: $rejected_traces traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n\n";
print "======================================================\n\n";

exit;

#########################
#
#
#  S U B R O U T I N E S
#
#
#########################




# subroutine to load file and clip sequences based on the following fields
# clip_quality_left
# clip_quality_right
# clip_vector_left
# clip_vector_right

sub process_file{
	# first want to split the file into two to make it easier to deal with the FASTA and XML separately
	# use UNIX split, with -p option to specify a regular expression
	
	my $split_file_prefix = "trace_read_split_file";
	
	system("split -p \"^<\\?xml\" $file $split_file_prefix") && die "Couldn't split file\n";
	
	# now remove the original file
	unlink($file) || die "Can't remove $file\n";
	
	# want to store fasta headers, fasta sequences, and seq lengths in hashes all tied to TI number of each sequence
	my %ti_to_seq;
	my %ti_to_header;
	my %ti_to_seqlength;


	# parse FASTA file first
	open(FASTA, "${split_file_prefix}aa") or die "Can't open fasta file";
	my $FA = new FAlite (\*FASTA);
	while (my $entry = $FA->nextEntry) {
		my $header = $entry->def;
		$header =~ m/ti\|(\d+) /;
		my $ti = $1;
		die "No ti in $header\n" if (!$ti);
		$ti_to_header{$ti} = $header;
		$ti_to_seq{$ti} = $entry->seq;
		$ti_to_seqlength{$ti} = length($entry->seq);
		
		# update global and species specific stats
		$total_bases += length($entry->seq);
		$species_total_bases += length($entry->seq);
		
	}
	close(FASTA);
	unlink("${split_file_prefix}aa") || die "Can't remove ${split_file_prefix}aa\n";
		 
	# now parse xml	
	# looking to find up to four different fields, and we will only consider the most extreme value for either left/right field
	# set artificially high value for right-values as we want to see if there are values that are lower than this
	my $clip_left=0;
	my $clip_right=1000000;
	my $vector_left=0;
	my $vector_right=1000000;
	my $left;
	my $right;
	

	# use a UNIX command-line combo to get all relevant information on one line for each trace
	my $command = "grep -E \"<ti>|<clip_quality_(left|right)>|<clip_vector_(left|right)>\" ${split_file_prefix}ab | sed 's/ti>\$/ti>\@/'  | tr '\n' '\#' | tr '\@' '\n'";
	open(XML,"$command |") or die "Can't open pipe\n";

	# can now parse line-by-line in a straightforward while loop
TRACE: while (my $line = <XML>) {
		next if ($line !~ /<ti>/);	
		
		# first want to get ti number and look up length
		(my $ti = $1)        if($line =~ m/<ti>(\d+)</);
		my $length = $ti_to_seqlength{$ti};
	
		# looking at each line to find up to four different fields, and we will only consider the most extreme value for either left/right field
		# set clip right-values to length of sequence as we want to see if there are values that are lower than this
		my $clip_left=0;
		my $clip_right=$length;
		my $vector_left=0;
		my $vector_right=$length;
		my $left=0;
		my $right=$length;


		($clip_left = $1)    if($line =~ m/<clip_quality_left>(\d+)</);
		($clip_right = $1)   if($line =~ m/<clip_quality_right>(\d+)</);
		($vector_left = $1)  if($line =~ m/<clip_vector_left>(\d+)</);
		($vector_right = $1) if($line =~ m/<clip_vector_right>(\d+)</);
		#print "\nTI:*$ti* LENGTH: $length CLIP_LEFT:$clip_left CLIP_RIGHT:$clip_right VECTOR_LEFT:$vector_left VECTOR_RIGHT:$vector_right\n";
		
		
		# now assign $left and $right to the most extreme values discovered
		($left = $clip_left)     if ($clip_left > $left);
		($left = $vector_left)   if ($vector_left > $left);
		($right = $clip_right)   if ($clip_right < $right);
		($right = $vector_right) if ($vector_right < $right);	
		#print "Length = $length Left = $left, right = $right\n";
		
		# now clip sequence if necessary
		my $seq = $ti_to_seq{$ti};

		if($left > 0 || $right < $length){
			
			# want to keep track total number of clipped bases and species specific clipped bases
			my $clipped_bases = $left + ($length - $right - 1);
			$total_clipped_bases += $clipped_bases;
			$species_clipped_bases += $clipped_bases;

			#print "Clipping $clipped_bases bases\n";		
			my $remaining_bases = $length - $clipped_bases;

			# use substr to do the actual clipping
			$seq = substr($seq,$left-1,$remaining_bases);
			
			# die if we don't have any sequence for some reason
			die "No seq\n$ti\tlength=$length\tleft=$left\tright=$right\n" if (!$seq);
		
			#######################			
			#
			# Add check in case remaining sequence is below some useful limit?
			#
			#######################
			
			if($remaining_bases < $min_bases){
				$rejected_traces++;
				$species_rejected_traces++;
				print "$ti has $remaining_bases bases after clipping\n";
				next(TRACE);
			}
			
			
			# modify FASTA header
			$ti_to_header{$ti} .= " CLIPPED: $clipped_bases nt";
		}
		my $tidied_seq = tidy_seq($seq);
		
		# now send to output file
		print OUT "$ti_to_header{$ti}\n$tidied_seq\n";
	}
	
	# clean up
	close(XML);
	unlink("${split_file_prefix}ab") || die "Can't remove ${split_file_prefix}ab\n";
	
}

# signal event handler in case of interrupts (Ctrl+C)
sub INT_handler {
	
	# print final statistic of how many bases were clipped
	$date = `date`; 
	chomp($date);
	
	my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);
	print "\n\n======================================================\n\n";
	print "SCRIPT INTERRUPTED at $date\n";
	print "TOTAL: Processed $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "TOTAL: $rejected_traces traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n\n";
	print "======================================================\n\n";
    exit(0);
}


sub tidy_seq{
#adds a new line character every 60 bases  
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;
    
    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
}