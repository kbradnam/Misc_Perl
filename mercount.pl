#!/usr/bin/perl -w
#
# mercount.pl
#
# A script to count mers in FASTA files
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use Getopt::Long;

########################
# Command line options
########################

my $mer;      # what size oligo to count

GetOptions ("mer=i" => \$mer);

# sanity checks
if(!$mer){
	die "Please specify the size of oligo that you wish to count with the -mer option\n";
}

if($mer>10){
	die "Are you out of your mind? Do you know how long this script will take to run!\n";
}

# simple hash, key is n-mer sequence, value is count
my %mers;

# Two arrays to hold current line and end of line if not an exact
# multiple of $mer
my (@seq,@end);

# string to hold current oligo 
my $seq;


# input file should contain fasta sequences.  Multiple sequences will be
# concatenated together and treated as one long sequence
open (IN, "<$ARGV[0]") || die "Failed to open input file\n\n"; 
while(<IN>){
    next if (m/^>/);
	chomp;
	@seq = split(//);

	# if we had some left from last line, add it now
	if(@end){
		@seq = (@end,@seq); 
	}
	
	# loop through each line counting mers	
	while(scalar(@seq)>=$mer){
		# grab oligo and make sure it is upper case
		$seq = uc(join ('',@seq[0..$mer-1]));
		# add to count of hash and reduce @seq by 1 bp
		$mers{$seq}++;
		shift(@seq);
	}	

	# need to know if there was sequence (< $mer bp) left on the end of the line
	# if so, put it in @end
	if(scalar(@seq)<$mer){
		@end = @seq;
	}
	else{
		@end = ();
	}
}
close(IN) || die "Couldn't close input file\n";  
 
# now print counts to output file
open(OUT,">$ARGV[0].${mer}.out") || die "This output file is being stubborn\n";
foreach my $key ((sort {$mers{$b} <=> $mers{$a,}} (keys %mers))){
	print OUT "$key $mers{$key}\n";
}
close(OUT) || die "Couldn't close damned output file\n";


exit(0);