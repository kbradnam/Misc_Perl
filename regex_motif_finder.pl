#!/usr/bin/perl
#
# regex_motif_finder.pl
#
# A script that uses specified regular expressions to find motifs in FASTA files
# Behaves a lot like where_is_motif.pl. Calculates motif density for any specified
# motif and can find individual motifs specified on the command line or multiple
# motifs listed in a file
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;

########################
# Command line options
########################

my $motif;      # nested MICA *.xms file with (single) motif
my $target;     # file with sequences in which to find motif
my $scores;     # Show scores
my $seqs;       # Show motif sequences in output (one sequence per motif in each intron)
my $msummary;	# show motif count and percentage for all sequences combined
my $file;		# switch to using input file of motifs (overides $motif)
my $flanking;   # inspect nucleotide frequencies at flanking position
my $distance;   # specify a distance threshold, and only consider motifs up to that distance

GetOptions ("motif=s"    => \$motif,
	       "target=s"    => \$target,
	       "scores"      => \$scores,
	       "seqs"        => \$seqs,	    	    
		   "msummary"    => \$msummary,
		   "file=s"	     => \$file,
		   "flanking"    => \$flanking,
		   "distance=i"  => \$distance
) or die "\n";

# are we using correct command-line options?
&pre_flight_checks;

# open and read motif file is $file is specified
my @motifs;
my %motifs_to_counts;

# hash to store counts of nucleotides before and after motif position
my %before;
my %after;

# will extract 1 nt either side of motif if -flanking specified - this might change in future
my $flank = 1;

# are we working with many motif regular expressions in a file?
if($file){
	open(FILE,"<$file") || die "Couldn't open $file file\n";
	while (my $entry = <FILE>){
		next if ($entry =~ m/^$/);
		chomp($entry);
		push(@motifs,$entry)
	}
	close(FILE);	
}
# otherwise just add single motif to @motifs
else{
	$motifs[0] = $motif;
}

##############################################################


open(TARGET,"<$target") || die "Couldn't open $target file\n";

my $fasta = new FAlite(\*TARGET);

# keep track of:
# 1) total length of all sequences in input file 
# 2) total length of all motifs in input file
# 3) number of sequences in input file
# 4) number of motifs in all sequences

my ($total_seq_length,$total_motif_length,$seq_count,$total_motif_count);


# loop through each sequence in target file
SEQ: while(my $entry = $fasta->nextEntry) {
	$seq_count++;
	# get header
	my $header = $entry->def;

    my $seq = lc($entry->seq);

	#trim sequences as necessary depending if $distance is being used
	($seq = substr $seq, 0, $distance) if ($distance);

    my $length = length($seq);
	$total_seq_length += $length;
		
	# will want to store counts from multiple motifs in one variable for later	
	my $output_text;

	
	foreach my $motif (@motifs){

		# convert motif to regex (it might already be in regex format)
		# e.g. RYN is motif, (A|G)(C|T)(A|C|G|T) would be the regex
		$motif = convert_to_regex($motif);
		$motif = lc($motif);
		
		# skip to next sequence if motif is longer than sequence
		my $motif_length = length($motif);
		next SEQ if ($motif_length > $length);
		
		# extract flanking positions and add to two hashes
		my @flanks = $seq =~ m/(\w{$flank}${motif}\w{$flank})/g if ($flanking);
		foreach my $seq (@flanks){
			$before{substr($seq,0,$flank)}++;
			$after{substr($seq,-$flank,$flank)}++;
		}

	 	# use regexes to find motif
		my $count = $seq =~ s/($motif)/ \U$1 /g;
		($count = 0) if (!$count);

		# count motifs
		$output_text .= "$count ";
		$total_motif_count += $count;
		$motifs_to_counts{$motif} += $count;
	}

	if ($scores){
		print "\n$header D=$length nt MOTIFS: $output_text\n" 
		
	}
	print "$seq\n" if ($seqs);
	

	# count how many bases are in motif, add to running total
	my $motif_sequence = ($seq =~ tr /A-Z/A-Z/);
	$total_motif_length += $motif_sequence;
	my $percent_motif = sprintf("%.3f",($total_motif_count / $length) * 100);	
}

close(TARGET) || die "Couldn't close $target\n";

# print motif summary if requested
if($msummary){
	my $percent_motif = sprintf("%.3f",($total_motif_length/$total_seq_length) * 100);
	print "\nSUMMARY:\n"; 
	print "number_of_sequences: $seq_count total_sequence_length: $total_seq_length\n";
	print "number_of_motifs: $total_motif_count total_motif_length: $total_motif_length motif_density: $percent_motif%\n\n";

	foreach my $motif (@motifs){
		print "$motifs_to_counts{$motif} ",uc($motif),"\n";
	}
	print "\n\n";
}


# do we need to print details of flanking nucleotides?
if($flanking){
	print "BEFORE motif:\n";

	foreach my $key (sort keys (%before)){
	
		print "$key ",sprintf("%.3f",$before{$key}/$total_motif_count), " $before{$key}\n";
	}
	print "\nAFTER motif:\n";

	foreach my $key (sort keys (%before)){
	
		print "$key ",sprintf("%.3f",$after{$key}/$total_motif_count), " $after{$key}\n";
	}
}

exit(0);


#############################################
#
#
#          S u b r o u t i n e s 
#
#
##############################################

sub pre_flight_checks{
	# check that both command line options are specified
	die "Need to specify both -motif and -target options\n" if(!$motif && !$target);

	# check that only -file or -motif mode is used
	die "Can only specify -motif *OR* -file options\n" if($motif && $file);


	# check files exist
	die "$target does not seem to exist\n" if (! -e $target);
	die "$file does not seem to exist\n" if ($file && ! -e $file);
}

# if input is in IUPAC format (RYSWN etc), need to convert to regexes
sub convert_to_regex{

	my ($motif) = @_;

	# first convert any IUPAC characters to regex equivalents
	my %bases_to_regex = (
		R => '(A|G)',
		Y => '(C|T)',
		S => '(C|G)',
		W => '(A|T)',
		K => '(T|G)',
		M => '(C|A)',
		B => '(C|G|T)',
		D => '(A|G|T)',
		H => '(A|C|T)',
		V => '(A|C|G)',
		N => '(A|C|G|T)'
	);

	foreach my $code (keys %bases_to_regex){
		$motif =~ s/$code/$bases_to_regex{$code}/g;
	}
	return($motif);
}
