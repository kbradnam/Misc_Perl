#!/usr/bin/perl -w
#
# genbank_info_dump.pl
#
# a script to process all genbank sequence files and extract information on intron and exon sizes
# this script mostly dumps info into two csv files (for exons and introns) to be processed by other scripts. 
#
# by Keith Bradnam
#
# created 11th November 2005 as another script, huge revision 29th March 2007
#
# Last updated by: $Author$
# Last updated on: $Date$    
#
#######################################################################################


use strict;
use Cwd;
use Getopt::Long;
use List::Util qw(sum);

########################
# Command line options
########################
my $release;   	# which release of genbank to query against?  Specify an integer
my $test; 		# simple test mode which just uses files in current working directory
GetOptions("release=i"=> \$release,
			"test"    => \$test);

die "Must use -release option to specify a GenBank release number to query against\n" if (!defined($release) && !$test);

################
# Paths etc.
################

# specify path (which changes if in test mode)
my $path;
$path = "/Volumes/genbank/${release}" if ($release);
$path = getcwd if ($test);

die "$path directory does not exist.\n" if (! -e "$path");

# have list of all GenBank files that we will be interested in (invertebrates, mammals, plants,
# primates, rodents, and vertebrates
my @inv = glob("$path/gbinv*.seq");
my @mam = glob("$path/gbmam*.seq");
my @pln = glob("$path/gbpln*.seq");
my @pri = glob("$path/gbpri*.seq");
my @rod = glob("$path/gbrod*.seq");
my @vrt = glob("$path/gbvrt*.seq");
my @htg = glob("$path/gbhtg*.seq");

# get list of wgs files not kept in standard division files
my @wgs = glob("${path}/wgs/wgs.*.gbff");

# combine everything together
my @files = (@inv,@mam,@pln,@pri,@rod,@vrt,@htg,@wgs);
	 
die "No *.seq files found in $path\n" if (@files == 0);    
	
########################
# misc variables
########################

# Use array to keep track of 6 stats:
# 1) total number of genbank entries processed
# 2) number of entries with at least one CDS feature
# 3) total number of CDS features
# 4) number of unusable partial CDSs
# 5) number of single bp exons at end of CDS (presumably part of UTR)
# 6) numer of usuable CDSs (can be just one exon)
# 7) number of usable CDSs with at least one intron
my @entries;

# three main things to extract per entry 
my ($mol,$accession,$species);

# true/false flags to know if we are on reverse strand or not 
my $reverse_strand;

# need to keep track of which CDS within a variable you are at
my $cds_counter;

# two variables to store exon and intron lengths
my ($exons,$introns);


############################################################
#
#  M A I N   L O O P
#
#  loop through GenBank files extracting data
#
############################################################


# two output streams for exons and introns
open(EXON,">genbank_r${release}_exons.csv") || die "Can't open exon output file\n";
open(INTRON,">genbank_r${release}_introns.csv") || die "Can't open intron output file\n";

# Change record delimiter to split on //
# newlines needed to avoid // in URLs for example
$/ = "\n//\n";

while (my $file = shift(@files)){

    print "Processing $file\n";
    open (FILE,"<$file") || die "Can't open $file\n"; 
    
    ENTRY: while (<FILE>) {
	
		# count entry & reset CDS counter
		$entries[0]++;
		$cds_counter = 0;
	
		# skip to start of first record if at the beginning of a file
		if(m/^GB\S+\.SEQ/){
	    	m/LOCUS\s+\S+\s+\d+ bp\s+(\S+)\s+/ || die "1) No molecule field found for:\n\n$_\n"; 
	    	$mol  = $1;
		}
		else{
	    	m/^LOCUS\s+\S+\s+\d+ bp\s+(\S+)\s+/ || die "2) No molecule field found for:\n\n$_\n"; 
	    	$mol  = $1;
		}       
		
		# check molecule type, ignore mRNAs as there won't be any introns
		next ENTRY unless ($mol eq "DNA");
		
		# skip if there is no CDS entry anywhere
		# this will no doubt match other things but will generally help speed things up 		
		next ENTRY unless (m/CDS/);

		
		# count entries with CDS
		$entries[1]++;

		# insert dividers for splitting remainder of entry
		s/\n(\S+)/\nSPLITHERE\n$1/g;
	
		# now loop through rest of single GenBank entry
		for (split(/SPLITHERE\n/)){

	    	# get accession
	    	if (m/ACCESSION\s+(\S+)/){
				$accession = $1;
	    	}

	    	# get species name
	    	if (/^SOURCE/) {
				m/ORGANISM\s+(\S+\s+\S+)/;
				$species = $1;
				# occasionally certain entries are structured in a way that gives problems, simplest to ignore these
				# these entries usually only contain the genus and no species suffix, such that $species always ends
				# up with 'Eukaryota' from the next line in the genbank entry
				next ENTRY if ($species =~ m/Eukaryota/);
	    	}
	    	    
	   		# split up each feature object
	    	if (/^FEATURES/) {
				s/\n\s{5}(\S+)/ZZZZ$1/g;
		
				# loop through each feature
	      		FEATURE: for (split(/ZZZZ/)) {
		    
		  			# skip if not a CDS, remove the CDS feature name as well
		  			next unless (s/^CDS\s+//);
		  			
					# now we are at a CDS, so count it (whether it is good or bad) and note it's number
					# relative to whole entry
					$entries[2]++;
					$cds_counter++;			
					
					# are we on reverse strand?
					if(m/complement/){
						$reverse_strand = 1;
					}
					else{
						$reverse_strand = 0;
					}
					
		  			# insert dividers for splitting into qualifiers and substitute excess space
		  			chomp ;
		  			s/\n\s{21}\//ZZZZ/g;
		  			s/\n\s{21}//g;
		  
		  			# now only want first part of CDS feature which will be just the coordinates
		  			my @qualifiers = split (/ZZZZ/); # split into qualifiers
		  			my $cds = shift (@qualifiers); # pull out location               				
							 
		  			$cds =~ s/\(//g;
		  			$cds =~ s/\)//g;
		  			$cds =~ s/join//g;
		  			$cds =~ s/complement//g;
		  	
		  			# ignore any CDS entry which contains references to other accessions, e.g.
		  			# AY437144.1:26..80
		  			next if ($cds =~ m/[A-Z]/);
		  		  	
					# ignore partial CDSs (with < or > in the CDS description) but keep count of them
					if ($cds =~ m/[<>]/){
						$entries[3]++;
						next;
					}
					
					# sometimes there are single base exons at the start or end of a CDS 
		  			# these only get one coordinate rather than the same coordinate twice (joined by ..)
		  			# E.g. complement(join(7218853,7218909..7219111))
					#
		  			# these are likely where the exon is a part of a 5' or 3' UTR. Post-processing
					# scripts can decide whether they want to use these or not, but we need to change
					# this into a format where we still calculate the length
					if (($cds =~ m/^\d+,\d+/) || ($cds =~ m/\d+,\d+$/)){
						$entries[4]++; 
						# just duplicate the single coordinate
						$cds =~ s/^(\d+),(\d+)/$1\.\.$1,$2/;
						$cds =~ s/(\d+),(\d+)$/$1,$2\.\.$2/;
					}
					
		  			# to find internal occurances, look for three consecutive digits (no ..)
		  			next if ($cds =~ m/\d+,\d+,\d+/);		   		 
			
					# now want to find out both exon and intron sizes
					$exons = $introns = $cds;	
					$exons =~ s/(\d+)\.\.(\d+)/$2-$1+1/ge;
					
					# set introns to be nothing if we have a single exon gene
					$introns = "" if ($cds !~ m/,/);


					# replace intron coords with intron lengths
					$introns =~ s/(\d+),(\d+)/$2-$1-1/ge;

					# skip any entry with zero or negative intron/exon lengths
					next if ($introns =~ m/\.\.0/);
					next if ($introns =~ m/\.\.-/);
					next if ($exons =~ m/\.\.0/);
					next if ($exons =~ m/\.\.-/);	
						
					# now we have removed most of the problem CDSs we can treat what is left
					# as valid entries for further processing
					$entries[5]++;

					# increment intron counter if CDS had introns
					$entries[6]++ if (defined($introns));
	
					# remove first and last exon coordinates and change '..' to commas
					$introns =~ s/^\d+\.\.([\d,\.]+)\.\.\d+$/$1/;
					$introns =~ s/\.\./,/g;

					# final challenge is to reverse order of coordinates if on the reverse strand
					if($reverse_strand == 1){
						my @int = split(/,/,$introns);
						$introns = join(",",reverse(@int));
						my @ex = split(/,/,$exons);
						$exons = join(",",reverse(@ex));						
					}

					# print out details unless we are dealing with single exon genes (so no introns)
					print EXON "$species,$accession,$cds_counter,$exons\n";
					print INTRON "$species,$accession,$cds_counter,$introns\n" if ($introns);
	      		}
	    	}
		}
    }
    close(FILE) || die "Can not close file\n";
}

print "\n\nProcessed $entries[0] GenBank entries, $entries[1] of which contains CDSs:\n\n";
print "$entries[2] CDS features in total\n";
print "$entries[3] partial CDS features which can't be used\n";
print "$entries[5] usable full-length CDS features\n";
print "$entries[6] usable CDS features with at least one intron\n";

close(EXON) || die "Can't close exon file\n";
close(INTRON) || die "Can't close intron file\n";

exit(0);
