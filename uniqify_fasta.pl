#!/usr/bin/perl
#
# uniqify_fasta.pl
#
# A script to process FASTA files to add a unique numerical prefix to the header
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use Keith;

# check we have at least one thing specified 
die "Usage: uniqify_fasta.pl <fasta file(s)>\n" unless (@ARGV);

my (@input_files) = @ARGV;

# will add a unique integer to each FASTA header
my $n = 0;

foreach my $file (@input_files){

	# only process *.fa or *.fasta files
	next unless ($file =~ m/\.fa$/ or $file =~ m/\.fasta$/);

	print STDERR "Processing $file\n";

	# open files
	my $output_file = "/tmp/$$.tmp";
	open(my $input, "<", "$file") or die "Can't open $file\n";
	open(my $out, ">", "$output_file") or die "Can't create $output_file\n";
	my $fasta = new FAlite(\*$input);
	while(my $entry = $fasta->nextEntry){
		$n++;
	    my $seq = Keith::tidy_seq($entry->seq);
	    my $def = $entry->def;
		$def =~ s/>(.*)/>$n $1/;
		print $out "$def\n$seq\n";
	}
	close($input);
	close($out);
	# now rename tmp file to original file
	rename($output_file, $file) or die "Can't rename $output_file to $file\n";
}

exit(0);


