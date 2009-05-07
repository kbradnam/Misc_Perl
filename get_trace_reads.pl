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
use Getopt::Long;
use LWP::UserAgent;
use HTTP::Request::Common 'POST';

$ENV{'LANG'}='C'; # copied from query_tracedb.pl - not sure what this is doing
$ENV{'LC_ALL'}='C'; # copied from query_tracedb.pl - not sure what this is doing
$SIG{'INT'} = 'INT_handler'; # for handling signal interrupts

my $min_bases; # minimum number of bases that you need in a sequence after clipping to keep it
my $max_pages; # maximum number of pages to download
my $page_size; # how many records to try to download in each network query (max is 40,000)
my $max_n;	   # what is the maximum percentage of N's allowed in the clipped sequence
my $verbose;   # turn on extra output - e.g. errors for traces that were rejected
my $program;   # can specify path to query_tracedb.pl program

GetOptions ("min_bases:i" => \$min_bases,
			"max_pages:i" => \$max_pages,
            "page_size:i" => \$page_size,
            "max_n:f"     => \$max_n,
			"verbose"     => \$verbose,
			"program:s"   => \$program);

# set defaults if not specified on command line
$min_bases = 20    if (!$min_bases);
$max_pages = 2     if (!$max_pages);
$page_size = 40000 if (!$page_size);
$max_n = 5         if (!$max_n);
$program = "query_tracedb.pl" if (!$program);


$SIG{'INT'} = 'INT_handler';



my @taxa = ("Plasmodium falciparum","Arabidopsis thaliana","Xenopus laevis","Drosophila melanogaster","Homo sapiens","Caenorhabditis japonica", "Procavia capensis");
#@taxa = ("Plasmodium falciparum");
#my @taxa = ("Arabidopsis thaliana");

# script name that does the actual fetching of data (supplied by NCBI)
# using glob incase '~' was used to specify location
my $prog = glob($program);

# need to keep track of:
# 1) total traces parsed
# 2) how many bases are clipped
# 3) how many traces are rejected for being too short after clipping low quality bases, 
# 4) how many traces are rejected for containing too many Ns
# 5) how many other errors there were (e.g. when clip left coordinate is greater than total length of sequence)
# Do this for all species (grand total), and keep separate species totals
my $total_traces = 0;
my $total_bases = 0;
my $total_clipped_bases = 0;
my $total_rejected_high_n = 0;
my $total_rejected_too_short = 0;
my $total_errors = 0;

my $species_total_bases;
my $species_clipped_bases;
my $species_rejected_high_n;
my $species_rejected_too_short;
my $species_errors;


# want to store fasta headers, fasta sequences, and seq lengths in hashes all tied to TI number of each sequence	
my %ti_to_seq;
my %ti_to_header;
my %ti_to_seqlength;



# first thing to do is test whether connection to trace server is down, may as well exit now if it is
check_network();

SPECIES: foreach my $species (@taxa){
	my $date = `date`; 
	chomp($date);
	print "\n------------------------------------------------------\n$date\n\n";
	
	# reset species-level counters
	$species_total_bases = 0;
	$species_clipped_bases = 0;
	$species_rejected_high_n = 0;
	$species_rejected_too_short = 0;
	$species_errors = 0;
	
	# store some of the query criteria in a separate string as we will want to reuse it later
	my $query = "species_code = '$species' and (trace_type_code = 'WGS' or trace_type_code = 'WCS' or trace_type_code = 'SHOTGUN' or trace_type_code = '454' or trace_type_code = 'CLONEEND') and source_type = 'GENOMIC'";

	#print "$prog query count \"$query\"\n";

	# Connection check
	check_network();

	# use 'query count' commands to get count of how many records match $query
	my $count = `$prog query count \"$query\"`;
	$count =~ s/\s+//g;
		

	# want to just check (and report) how many records are being filtered out at this first stage
	my $count_all = `$prog query count \"species_code = '$species'"`;
	$count_all =~ s/\s+//g;
	
	# make a string of the species name but without spaces
	my $species_file_name = $species;
	$species_file_name =~ s/\s+/_/g;

	# can only download 40,000 records at a time
	my $pages = int($count/40000)+1;

	# there might be no records after above filtering and so we need to skip to next species
	if($count == 0){
		print "${species_file_name}: 0 reads remain after filtering (from $count_all possible reads)\n";
		print "------------------------------------------------------\n";
		next SPECIES;
	}
	else{
		my $percent_retained = sprintf("%.1f",($count/$count_all)*100);
		print "${species_file_name}: Attempting to fetch $count reads out of $count_all possible (${percent_retained}%), $pages pages\n";		
	}

	# if we have any sequences to fetch then need to make an output file
	open(OUT,">${species_file_name}_trace_reads.fa") or die "Can't open output file for $species\n";
	
	
	PAGES: for (my $i=0;$i<$pages;$i++){

		if ($i >= $max_pages){
			print "${species_file_name}: Stopping because max_pages limit of $max_pages has been reached\n";
			last PAGES;
		}
		
		my $file = "${species_file_name}_trace_read_data${i}";
	
		# clear hashes
		%ti_to_seq = ();
		%ti_to_header = ();
		%ti_to_seqlength = ();
	
		print "${species_file_name}: Processing page ",$i+1,"/$pages\n";

		
		##########################################
		#
		# G R A B   T R A C E   S E Q U E N C E S 
		#
		##########################################

		# Connection check
		check_network();
				
		# Form command that will grab FASTA sequences of specified reads
		my $command = "(/bin/echo -n \"retrieve_gz fasta 0b\"\; $prog \"query page_size $page_size page_number $i binary $query\") | $prog"; 

		# send command through a pipe to gunzip and then output will be sent to FAlite module
		open(FASTA,"$command | /usr/bin/gunzip -c | ") or die "Can't open pipe: $? $!\n";
		my $FA = new FAlite (\*FASTA);
		while (my $entry = $FA->nextEntry) {
			die "$entry\n" if ($entry =~ m/Couldn't connect to TRACE server/i);			
			my $header = $entry->def;
			$header =~ m/ti\|(\d+) /;
			my $ti = $1;
			die "No ti in $header\n" if (!$ti);
			$ti_to_header{$ti} = $header;
			$ti_to_seq{$ti} = $entry->seq;
			
			die "Can't get sequence for $ti\n" if (!defined($ti_to_seq{$ti}));
			
			$ti_to_seqlength{$ti} = length($entry->seq);

			# update global and species specific stats
			$total_bases += length($entry->seq);
			$species_total_bases += length($entry->seq);
		}	
		close(FASTA);
	
	
		################################
		#
		# G R A B   T R A C E   I N F O
		#
		################################
		
		
		# now make new command for getting trace info for each read
		$command = "(/bin/echo -n \"retrieve_gz info 0b\"\; $prog \"query page_size $page_size binary page_number $i $query\") | $prog"; 
		
		# Connection check
		check_network();
	
	
		# send it through a pipe to send it to perl, use gunzip -c to uncompress data in place
		open(INFO,"$command | /usr/bin/gunzip -c | ") or die "Can't open pipe\n";

		# change input record separator when processing this pipe
		# a 'state' field should be the last field in each entry
		$/ = "state:\t";
		
		# can now parse trace read information straightforward while loop
		
 		TRACE: while (my $record = <INFO>) {
		
			next if ($record !~ /ti:/);	
			$total_traces++;
		
			# first want to get ti number and look up length
			(my $ti = $1) if ($record =~ m/ti:\t(\d+)/);
			
			my $length = $ti_to_seqlength{$ti};
		
			# looking at each line to find up to four different fields, and we will only consider the most extreme value for either left/right field
			# set clip right-values to length of sequence as we want to see if there are values that are lower than this
			my $clip_left=0;
			my $clip_right=$length;
			my $vector_left=0;
			my $vector_right=$length;
			my $left=0;
			my $right=$length;

			($clip_left = $1)    if($record =~ m/clip_quality_left:\t(\d+)/);
			($clip_right = $1)   if($record =~ m/clip_quality_right:\t(\d+)/);
			($vector_left = $1)  if($record =~ m/clip_vector_left:\t(\d+)/);
			($vector_right = $1) if($record =~ m/clip_vector_right:\t(\d+)/);
		 	#print "\nTI:*$ti* LENGTH: $length CLIP_LEFT:$clip_left CLIP_RIGHT:$clip_right VECTOR_LEFT:$vector_left VECTOR_RIGHT:$vector_right\n";

			# now assign $left and $right to the most extreme values discovered
			($left = $clip_left)     if ($clip_left > $left);
			($left = $vector_left)   if ($vector_left > $left);

			# some clip_right values are greater than the length of the sequence, in which case we can change to
			# set them to the sequence length, i.e. no right clipping
			($clip_right   = $length) if ($clip_right   > $length);
			($vector_right = $length) if ($vector_right > $length);
		
			# for the right-fields, need to remember that a zero value just means no clip information is present
			($right = $clip_right)   if (($clip_right < $right)   && ($clip_right   != 0));
			($right = $vector_right) if (($vector_right < $right) && ($vector_right != 0));	
			#print "Length = $length Left = $left, right = $right\n";
		
		
			# have some basic sanity checks to catch errors		
			my $error_test = 0;
		
			$error_test = 1 if ($clip_left   > $length);
			$error_test = 1 if ($vector_left > $length);
			$error_test = 1 if ($left        > $right);
		
			if($error_test){
				$species_errors++;
				$total_errors++;
				print "ERROR: TI:*$ti* LENGTH: $length CLIP_LEFT:$clip_left CLIP_RIGHT:$clip_right VECTOR_LEFT:$vector_left VECTOR_RIGHT:$vector_right\n" if ($verbose);
			
				# no point going any further
				next TRACE;
			}
		
			# now clip sequence if necessary
			my $seq = $ti_to_seq{$ti};
			if($left > 0 || $right < $length){
			
				# want to keep track total number of clipped bases and species specific clipped bases
				my $clipped_bases = $left + ($length - $right - 1);
				$total_clipped_bases += $clipped_bases;
				$species_clipped_bases += $clipped_bases;

				# how many bases will be left after clipping?
				my $remaining_bases = $length - $clipped_bases;
				
				# Add check in case remaining sequence is below some useful limit?		
				if($remaining_bases < $min_bases){
					$total_rejected_too_short++;
					$species_rejected_too_short++;
					print "ERROR: $ti has $remaining_bases bases after clipping\n" if ($verbose);
					next(TRACE);
				}
				
				# use substr to do the actual clipping and modify corresponding FASTA header with clipping info
				$seq = substr($seq,$left-1,$remaining_bases);
				$ti_to_header{$ti} .= " CLIPPED: $clipped_bases nt";
			}
			
			# die if we don't have any sequence for some reason
			die "No seq\n$ti\tlength=$length\tleft=$left\tright=$right\n" if (!$seq);
		
			# now check whether remaining sequence contains more than $max_n Ns in which case we should ignore it
			my $n = ($seq =~ tr/N/N/);
			my $percent_n =($n / length($seq)*100);
			if($percent_n > $max_n){
				$total_rejected_high_n++;
				$species_rejected_high_n++;
				#print "ERROR: $ti contains $percent_n Ns, more than ${max_n}% Ns in its sequence\n";
				
				print "ERROR: $ti contains more than ${max_n}% Ns in its sequence\n" if ($verbose);
				next(TRACE);
			}
			
			# If we have survived to this point, then we have a valid sequence which we can tidy and then print to output
			my $tidied_seq = tidy_seq($seq);
			print OUT "$ti_to_header{$ti}\n$tidied_seq\n";
		}
	
		# clean up
		close(INFO);
		
		# reset input record separator
		$/ = "\n";
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
	print "$species_file_name: Processed $species_total_bases nt of which $species_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "$species_file_name: WARNING - $species_rejected_too_short traces were rejected for being too short after vector/quality clipping\n" if ($species_rejected_too_short > 0);
	print "$species_file_name: WARNING - $species_rejected_high_n traces were rejected for containing too many unknown bases (Ns)\n" if ($species_rejected_high_n > 0);
	print "$species_file_name: WARNING - $species_errors traces were rejected for containing errors (inconsistant information)\n" if ($species_errors > 0);
	
	print "------------------------------------------------------\n";
}

my $date = `date`; 
chomp($date);

print "\nFINISHED RUN AT $date\n";

# Final print out of stat
my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);

print "\n\n======================================================\n\n";
print "TOTAL: Processed $total_traces traces containing $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
print "TOTAL: $total_rejected_too_short traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n";
print "TOTAL: $total_rejected_high_n traces were rejected for containing too many unknown bases (>%$max_n)\n";
print "TOTAL: $total_errors traces were rejected for containing errors (inconsistant information)\n\n";
print "======================================================\n\n";

exit;

#########################
#
#
#  S U B R O U T I N E S
#
#
#########################


# signal event handler in case of interrupts (Ctrl+C)
sub INT_handler {
	
	# print final statistic of how many bases were clipped
	$date = `date`; 
	chomp($date);
	
	my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);
	print "\n\n======================================================\n\n";
	print "SCRIPT INTERRUPTED at $date\n";
	print "TOTAL: Processed $total_traces traces containing $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "TOTAL: $total_rejected_too_short traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n";
	print "TOTAL: $total_rejected_high_n traces were rejected for containing too many unknown bases (>%$max_n)\n";
	print "TOTAL: $total_errors traces were rejected for containing errors (inconsistant information)\n\n";
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


sub check_network{
                
    # try to retrieve FASTA of TI number 1 as a test of network
    my $request = POST 'http://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=raw', [query=>'retrieve fasta 1'];
    my $response =  LWP::UserAgent->new->request($request);
	
	# first check that connection actually worked
    if ($response->is_error){
   		die "Connection attempt to trace server resulted in an error\n";
    }
	
	# now check that there wasn't an 001 or 002 error, 
	# i.e. connection worked but trace server returned an error code rather than sequence or info
	if($response->content =~ m/^00\d:/i) {
		die "Connection to trace server returned an error code: ",$response->content,"\n";
	}
}

