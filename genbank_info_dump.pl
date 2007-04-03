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
use Getopt::Long;
use List::Util qw(sum);

########################
# Command line options
########################
my $release;   # which release of genbank to query against?  Specify an integer

GetOptions("release=i"=> \$release);

die "Must use -release option to specify a GenBank release number to query against\n" if (!defined($release));

################
# Paths etc.
################

# specify path to GenBank release
my $path = "/Volumes/genbank/${release}";
#$path = glob("~/Desktop");

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
	     
#@files = glob("$path/gbinv1.seq");


########################
# misc variables
########################

# Use array to keep track of 4 basic stats:
# 1) total number of genbank entries processed
# 2) number of entries with at least one valid CDS feature
# 3) number of entries with at least one valid exon in a CDS feature
# 4) numer of entries with undefined start/end coordinates of CDS
# 5) number of entries with at least one intron in a CDS feature
my @entries;

# three main things to extract per entry 
my ($mol,$accession,$species);

# true/false flags to know if we are on reverse strand or not and
# whether CDS has undetermined start/end coordinates of an exon
my ($reverse_strand,$undetermined);

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
open(EXON,">genbank_dump_r${release}_exon.csv") || die "Can't open exon output file\n";
open(INTRON,">genbank_dump_r${release}_intron.csv") || die "Can't open intron output file\n";

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
	    	}
	    	    
	   		# split up each feature object
	    	if (/^FEATURES/) {
				s/\n\s{5}(\S+)/ZZZZ$1/g;
		
				# loop through each feature
	      		FEATURE: for (split(/ZZZZ/)) {
		    
		  			# skip if not a CDS, remove the CDS feature name as well
		  			next unless (s/CDS\s+//);
		  		
				
					# are we on reverse strand?
					if(m/complement/){
						$reverse_strand = 1;
					}
					else{
						$reverse_strand = 0;
					}
					
					# set flag for undetermined start/end of exons
					$undetermined = 0;
		
		  			# insert dividers for splitting into qualifiers and substitute excess space
		  			chomp ;
		  			s/\n\s{21}\//ZZZZ/g;
		  			s/\n\s{21}//g;
		  
		  			# now only want first part of CDS feature which will be just the coordinates
		  			my @qualifiers = split (/ZZZZ/); # split into qualifiers
		  			my $cds = shift (@qualifiers); # pull out location               
					
					# count all CDSs
					$entries[2]++;
							 
		  			$cds =~ s/\(//g;
		  			$cds =~ s/\)//g;
		  			$cds =~ s/join//g;
		  			$cds =~ s/complement//g;
		  	
		  			# ignore any CDS entry which contains references to other accessions, e.g.
		  			# AY437144.1:26..80
		  			next if ($cds =~ m/[A-Z]/);
		  		  		
		  			# check that there are always pairs of exon coordinates, e.g. avoid scenarios like
		  			# this in accession AE003538 (has a final single exon start coordinate but with no end
		  			# CDS             complement(join(290022..290246,290389..290873,
		  			#                 290960..291199,291264))
		  			# bit of a fudge, only looking for pair of coordinates at start and end of CDS					
		  			next unless ($cds =~ m/^[\d<]+\.\.\d+/);
		  			next unless ($cds =~ m/\d+\.\.[>\d]+$/);
					
		  			# to find internal occurances, look for three consecutive digits (no ..)
		  			next if ($cds =~ m/\d+,\d+,\d+/);		   		 

					# now we have removed most of the problem CDSs we can treat what is left
					# as valid entries for further processing
					# increment total CDS counter, and current CDS couunter for this entry
					$entries[3]++;
					$cds_counter++;					

					# print "CDS: $cds\n";				
			
					# some CDSs don't know the exact ends of exons (denoted by use of < and > characters)
					# so have to exclude these entries when calculating exon lengths
			 		if ($cds =~ m/[<>]/){
						$entries[4]++;
						$undetermined = "1";
						# remove these characters so that we can still calculate intron sizes
						$cds =~ s/[<>]//g;
			  		}
					
					# now want to find out both exon and intron sizes
					$exons = $introns = $cds;	
					$exons =~ s/(\d+)\.\.(\d+)/$2-$1+1/ge;
					
					# set introns to be nothing if we have a single exon gene
					# else advance intron counter
					if ($cds !~ m/,/){
						$introns = "";
					}
					else{
						$entries[5]++;
					}

					# replace intron coords with intron lengths
					$introns =~ s/(\d+),(\d+)/$2-$1-1/ge;

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

					# print out details unless we are dealing with exons with undetermined ends
					# or single exon genes
					print EXON "$species,$accession,$cds_counter,$exons\n" unless ($undetermined == 1);
					print INTRON "$species,$accession,$cds_counter,$introns\n" if ($introns);
	      		}
	    	}
		}
    }
    close(FILE) || die "Can not close file\n";
}

print "\n\nProcessed $entries[0] GenBank entries:\n\n";
print "$entries[1] entries contained $entries[2] CDS features\n";
print "$entries[3] usable CDS features in total ($entries[4] of which have undetermined exon start/end coordinates)\n";
print "$entries[5] usable CDS features with at least one intron\n";

close(EXON) || die "Can't close exon file\n";
close(INTRON) || die "Can't close intron file\n";

exit(0);
