#!/usr/bin/perl -w
#
# seqzip.pl
#
# A script to create compressed FASTA files which are still usable
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use IO::Handle;
# turn on autoflush
# $|=1;

# input file should contain fasta sequences.  Multiple sequences will be
# concatenated together and treated as one long sequence
open (IN, "<$ARGV[0]")     || die "Failed to open input file\n\n"; 
open (OUT,">$ARGV[0].cna") || die "This output file is being stubborn\n";

# need to know when we move to a new sequence
my $newseq = 0;

# variable to store sequence to be compressed
my $seq;

# need to store end of each line to potentially prepend
# to start of next line. We always have three possible endlines so we need a few variables
my ($endline,$endline1,$endline2,$endline3);

while(<IN>){
	chomp;
	
    if (m/^>/){
		# if we reach a new sequence first print out $endline from
		# previous sequence (if it exists)
		
		# Also need to reduce any multiple characters if possible
		if($endline){
			$endline = &compress($endline);
			print OUT "$endline\n";
		}
		# then print header line and set $newseq to remember that we are in new sequence
		print OUT "$_\n";
		$newseq = 1;
	}
	elsif(m/^(\s)*$/){
		# skip if blank
		next;
	}
	else{
		
		# will prepend end of previous line unless we are in new sequence
		unless($newseq == 1){	
			$_ = $endline.$_;
		}
		$newseq = 0;
		# make lower case just to make sure that we can more easily
		# recognise polyN runs later on (which will be in uppercase)
		tr/A-Z/a-z/;

		# need to take whatever repeat ends the line (even if it is just a single nucletotide). 
		# This could also be the whole line look for poly-N runs first, but then look for poly-XY runs, 
		# and then poly-XYZ runs
		m/(a+|c+|g+|t+|n+)$/i;
		$endline1 = $1;

		m/((ac)+|(ag)+|(at)+|(ca)+|(cg)+|(ct)+|(ga)+|(gc)+|(gt)+|(ta)+|(tc)+|(tg)+)$/i;
		$endline2 = $1;

		m/((aac)+|(aag)+|(aat)+|(aca)+|(acc)+|(acg)+|(act)+|(aga)+|(agc)+|(agg)+|(agt)+|(ata)+|(atc)+|(atg)+|(att)+)$/i;
		m/((caa)+|(cac)+|(cag)+|(cat)+|(cca)+|(ccg)+|(cct)+|(cga)+|(cgc)+|(cgg)+|(cgt)+|(cta)+|(ctc)+|(ctg)+|(ctt)+)$/i;
		m/((gaa)+|(gac)+|(gag)+|(gat)+|(gca)+|(gcc)+|(gcg)+|(gct)+|(gga)+|(ggc)+|(ggt)+|(gta)+|(gtc)+|(gtg)+|(gtt)+)$/i;
		m/((taa)+|(tac)+|(tag)+|(tat)+|(tca)+|(tcc)+|(tcg)+|(tct)+|(tga)+|(tgc)+|(tgg)+|(tgt)+|(tta)+|(ttc)+|(ttg)+)$/i;		
		$endline3 = $1;
		
		# now see which endline was the longest
		$endline = $endline1 if (length($endline1) >= length($endline2));
		$endline = $endline2 if (length($endline2) >= length($endline1));
		$endline = $endline3 if (length($endline3) >= length($endline));
		
		# remove $endline from $_ (this might remove whole line)
		s/($endline)$//;
			
		# skip line now if $endline has swallowed up everything and $_ is blank
		# but process $endline if $endline is getting too big
#		print length($endline)."\n";

		if(m/^(\s)*$/){
			if(length($endline)>5000){
					$_ = $endline;
					$endline = "";
			}
			else{next;}
		}
			
		# compress sequence and print
		$seq = &compress($_);
		print OUT "$seq\n";	
		OUT->autoflush(1);	
	}
}
# print whatever is left from the end of the last sequence in the file
$endline = &compress($endline);	
print OUT "$endline\n";

close(IN)  || die "Couldn't close input file\n";  
close(OUT) || die "Couldn't close damned output file\n";

# The main subroutine that replaced poly-N nucleotide or poly-XY dinucleotide runs
sub compress{
	my $seq = shift;
	
	# replace poly-N mononucleotide runs (only if there are at least 4)
	$seq =~ s/a{4,}|c{4,}|t{4,}|g{4,}|n{4,}/length($&)."{".uc(substr($&,0,1))."}"/ge;

	# replace 3 or more identical dinucleotides with a number
	$seq =~ s/(ac){3,}|(ag){3,}|(at){3,}/(length($&)\/2)."{".uc(substr($&,0,2))."}"/ge;
	$seq =~ s/(ca){3,}|(cg){3,}|(ct){3,}/(length($&)\/2)."{".uc(substr($&,0,2))."}"/ge;
	$seq =~ s/(ga){3,}|(gc){3,}|(gt){3,}/(length($&)\/2)."{".uc(substr($&,0,2))."}"/ge;
	$seq =~ s/(ta){3,}|(tc){3,}|(tg){3,}/(length($&)\/2)."{".uc(substr($&,0,2))."}"/ge;

	# replace 3 or more triplet repeats
	$seq =~ s/(aac){3,}|(aag){3,}|(aat){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(aca){3,}|(acc){3,}|(acg){3,}|(act){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(aga){3,}|(agc){3,}|(agg){3,}|(agt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(ata){3,}|(atc){3,}|(atg){3,}|(att){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;

	$seq =~ s/(caa){3,}|(cac){3,}|(cag){3,}|(cat){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(cca){3,}|(ccg){3,}|(cct){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(cga){3,}|(cgc){3,}|(cgg){3,}|(cgt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(cta){3,}|(ctc){3,}|(ctg){3,}|(ctt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;

	$seq =~ s/(gaa){3,}|(gac){3,}|(gag){3,}|(gat){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(gca){3,}|(gcc){3,}|(gcg){3,}|(gct){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(gga){3,}|(ggc){3,}|(ggt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(gta){3,}|(gtc){3,}|(gtg){3,}|(gtt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	
	$seq =~ s/(taa){3,}|(tac){3,}|(tag){3,}|(tat){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(tca){3,}|(tcc){3,}|(tcg){3,}|(tct){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(tga){3,}|(tgc){3,}|(tgg){3,}|(tgt){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;
	$seq =~ s/(tta){3,}|(ttc){3,}|(ttg){3,}/(length($&)\/3)."{".uc(substr($&,0,3))."}"/ge;

	return($seq);
}
exit(0);