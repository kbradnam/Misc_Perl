#!/usr/bin/perl
#
# ftp_trace_reads.pl
#
# A script to download trace reads of selected species from NCBI trace archive
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use Net::FTP;


# Need this or else things won't work!
$ENV{FTP_PASSIVE} = 1;

my $max_files;    # maximum number of files to download for any species
my $debug;        # whether to turn on debugging in FTP module
my $timeout;      # set timeout value for Net::FTP module
my $sleep;        # how long to sleep for before retrying ftp download
my $max_attempts; # how many attempts to download one file before giving up
my $min_traces;   # what is the minimum number of traces that each species needs to have to continue processing it
my $species_list; # optionally specify a file which contains species (quicker than looking up via separate script)
my $prog;         # specify path to a program that will produce list of (eukaryotic) species
my $verbose;      # turn on extra output (usefulf for debugging)

GetOptions ("max_files:i"    => \$max_files,
			"debug"          => \$debug,
			"timeout:i"      => \$timeout,
			"sleep:i"        => \$sleep,
			"max_attempts:i" => \$max_attempts,
			"min_traces:i"   => \$min_traces,
			"species_list:s" => \$species_list,
			"prog:s"         => \$prog,
			"verbose"        => \$verbose);

# set defaults if not specified on command line
$max_files = 2     if (!$max_files);
$timeout = 180     if (!$timeout);
$sleep = 10        if (!$sleep);
$max_attempts = 5  if (!$max_attempts);
$min_traces = 1000 if (!$min_traces);
$prog = glob("~/Work/bin/find_eukaryotes_in_trace_archive.pl -min_traces $min_traces") if (!$prog);

if (!$debug){
	$debug = 0;
}
else{
	$debug = 1;
}

#########################################
# GET LIST OF EUKARYOTES IN TRACE ARCHIVE
#########################################

my @taxa = get_species_names();






########################
# BASIC FTP settings
########################
my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "krbradnam\@ucdavis.edu";
my $root = "pub/TraceDB";

#my @taxa = ("Plasmodium falciparum","Arabidopsis thaliana","Xenopus laevis","Drosophila melanogaster","Homo sapiens","Caenorhabditis japonica", "Procavia capensis");
#my @taxa = ("Drosophila melanogaster","Homo sapiens","Caenorhabditis japonica", "Procavia capensis");


# keep track of how many species did not have an exact file name match on FTP site
my $missing_counter = 0;
my $species_counter = 0;


my $ftp;


SPECIES: foreach my $species (@taxa){
	$species_counter++;

	# get species name in correct format (should already be lower case)
	chomp($species);	
	$species = lc($species);
	$species =~ s/ /_/g;
	
	print "\nAttempting to fetch sequences for $species\n";
	print "======================================================\n";
	my $dir = "$root/$species";

	$ftp = Net::FTP->new($host, Debug => $debug, Timeout => $timeout)  or die "Cannot connect to $host: $@\n",$ftp->message;
	$ftp->login($user,$password) or die "Cannot login ", $ftp->message, "\n";
	$ftp->binary;


	# need to find out how many files are in directory, grab all FASTA files and use array index of last file
	# also keep count of how many species we don't get an exact name match for
	my @fasta;
	unless(@fasta = $ftp->ls("$dir/fasta.$species.[0-9]*.gz")){
		print "MISSING: $species\n" if ($verbose);
		$missing_counter++;
		next SPECIES;	
	}		
	
	my $last_file = $fasta[-1];
	$last_file =~ m/fasta.$species.*(\d+).gz/;
	my $number_of_files = $1;

	FILE: for (my $i=1;$i<=$number_of_files;$i++){
		
		# break out of loop if we have exceeded max number of pages
		if ($i > $max_files){
			print "Maximum number of files ($max_files) has been exceeded, skipping to next species\n";
			last FILE;
		}	
		
		# grab files in pairs, fasta + clip (they should pair up) 
		my $fasta_return = get_files($dir,$species,$i,$number_of_files,"fasta",1);
		print "get_files (FASTA) failed for $species $i\n" if (!$fasta_return);
		
		my $clip_return = get_files($dir,$species,$i,$number_of_files,"clip",1);
		print "get_files (CLIP) failed for $species $i\n" if (!$clip_return);
		
	}
	# tidy up
	$ftp->quit or die "Can't quit FTP",$ftp->message;
}



print "\n$missing_counter species (out of $species_counter) could not be found on FTP site, might be due to slight variations in species names\n\n";

exit;



sub get_files{
	my ($dir,$species,$i,$number_of_files,$type,$attempt) = @_;

	my $number = sprintf("%03d", $i);

	print "Getting $type number $i of $number_of_files, attempt number $attempt\n";

	# potentially have to fetch two differently formatted file names depending on whether there are no quality scores or not
	my $file1 = "$type.$species.$number.gz";
	my $file2 = "$type.$species.qualityless.$number.gz";
	
	
	if(defined $ftp->size("$dir/$file1")){

		# check whether file exists locally (and is same size)
		my $size = $ftp->size("$dir/$file1");
		
		if(-e $file1 && (-s $file1 == $size)){
			print "$file1 exists locally - skipping\n\n";
			return(1);
		}
		elsif(-e $file1 && (-s $file1 != $size)){
			print "$file1 exists locally but is a different size, will download again\n";
		}

		# use eval statement to catch any timeouts
		eval{$ftp->get("$dir/$file1") or die "Can't grab $file1\n",$ftp->message};	
		if ($@ =~ /Timeout/){
		   # catching code goes here
			print "$@\n";
			return(0) if ($attempt > $max_attempts);
			print "Sleeping for $sleep seconds, and then retrying\n";
			sleep($sleep);
			get_files($dir,$species,$number,$i,$type,++$attempt);
		}		
	}
	
	# now try qualityless file name 
	elsif(defined $ftp->size("$dir/$file2")){
		
		# check whether file exists locally (and is same size)
		my $size = $ftp->size("$dir/$file2");
		
		if(-e $file2 && (-s $file2 == $size)){
			print "$file2 exists locally - skipping\n\n";
			return(1);
		}
		elsif(-e $file2 && (-s $file2 != $size)){
			print "$file2 exists locally but is a different size, will download again\n";
		}
		
		eval{$ftp->get("$file2") or die "Can't grab $file2\n",$ftp->message};
		if ($@ =~ /Timeout/){
		   # catching code goes here
			print "$@\n";
			return(0) if ($attempt > $max_attempts);
			print "Sleeping for $sleep seconds, and then retrying\n";
			sleep($sleep);
			get_files($dir,$species,$number,$i,$type,++$attempt);
		}		
	}
	
	# or give up
	else{
		print "ERROR: $type file $i can not be found\n";
		return(0);
	}
	

	# if we get here, then we should have downloaded a file sucessfully
	return(1);
}



sub get_species_names{
	print "Fetching list of eukaryotes by using $prog\n";

	my @taxa;
	if($species_list){
		open (IN,"<$species_list") or die "-Can't find file specified by -species_file: $species_list\n";
		while(<IN>){
			push(@taxa,$_);
		}
		close(IN);
	}
	else{
		@taxa = `$prog`or die "Can't run $prog\n";	
	}
	return(@taxa);
}