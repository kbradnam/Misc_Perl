#!/usr/bin/perl -w
#
# genbank_intron_counter.pl
#
# a script to count the number of CDSs with and without introns in each species in GenBank
# and also count average number of introns per gene per species. Only does this for entries
# with molecule type of DNA
#
# by Keith Bradnam
#
# created 11th November 2005
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

die "Must specify a GenBank release number to query against\n" if (!defined($release));


########################
# misc variables
########################

# Hash for keeping track of all species and related data
# Key is species name, value is an array which contains:
# 0 - Total number of CDSs
# 1 - Total number of introns
# 2 - Number of CDSs with no introns (single exon genes)
# 3 - CDSs with multiple exons but where the start or end may not be known

my %all_species;


# specify path to GenBank release
my $path = "/Volumes/GenBank/genbank${release}";
#$path = glob("~keith");


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


############################################################
#
#  M A I N   L O O P
#
#  loop through GenBank files extracting data
#
############################################################

# count number of genbank entries processed
my $entries;

# Change record delimiter to split on //
# newlines needed to avoid // in URLs for example
$/ = "\n//\n";

while (my $file = shift(@files)){

    print "Processing $file\n";
    open (FILE,"<$file") || die "Can't open $file\n"; 
    
    ENTRY: while (<FILE>) {
	
		my ($locus,$mol,$accession,$species);
		$entries++;
	
		# skip to start of first record if at the beginning of a file
		if(m/^GB\S+\.SEQ/){
	    	m/LOCUS\s+(\S+)\s+\d+ bp\s+(\S+)\s+/ || die "1) No LOCUS field found for:\n\n$_\n"; 
	    	$locus = $1;
	    	$mol  = $2;
		}
		else{
	    	m/^LOCUS\s+(\S+)\s+\d+ bp\s+(\S+)\s+/ || die "2) No LOCUS field found for:\n\n$_\n"; 
	    	$locus = $1;
	    	$mol  = $2;
		}       
		
		# check molecule type, ignore mRNAs as there won't be any introns
		next ENTRY unless ($mol eq "DNA");
		
		# skip if there is no CDS entry anywhere
		# this will no doubt match other things but will generally help speed things up 		
		next ENTRY unless (m/CDS/);

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
				# skip if not human
				#next ENTRY unless ($species eq "Homo sapiens");
				
	    	}
	    	    
	   		# split up each feature object
	    	if (/^FEATURES/) {
				s/\n\s{5}(\S+)/ZZZZ$1/g;
		
				# loop through each feature
	      		FEATURE: for (split(/ZZZZ/)) {
		    
		  			# skip if not a CDS, remove the CDS feature name as well
		  			next unless (s/CDS\s+//);
		  
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
		  		  
		  			# check that there are always pairs of exon coordinates, e.g. avoid scenarios like
		  			# this in accession AE003538 (has a final single exon start coordinate but with no end
		  			# CDS             complement(join(290022..290246,290389..290873,
		  			#                 290960..291199,291264))
		  			# bit of a fudge, only looking for pair of coordinates at start and end of CDS
		  			next unless ($cds =~ m/^[\d<]+\.\.\d+/);
		  			next unless ($cds =~ m/\d+\.\.[>\d]+$/);

		  			# to find internal occurances, look for three consecutive digits (no ..)
		  			next if($cds =~ m/\d+,\d+,\d+/);		   		 
	
					# store number of introns
					my $introns;

		  			# some CDSs don't know the exact ends of exons (denoted by use of < and > characters)
		  			# so will treat these entries separately
		  			if ($cds =~ m/[<>]/){
			  			$introns = $cds =~ tr/,/,/;
			  			# keep count of these CDSs which have at least 1 intron
						if($introns > 0){
							$all_species{$species}[0] += 0;
							$all_species{$species}[1] += 0;
							$all_species{$species}[2] += 0;
							$all_species{$species}[3]++;
						}
		  			}
			 		else{
		  				# can now trust that we have a complete, correct CDS entry so can increment
		  				# cds counter
		  				$all_species{$species}[0]++;
		  
		  				# number of commas will in $cds will be the number of introns		  
		  				$introns = $cds =~ tr/,/,/;

		  				# now increment total number of introns, and number of CDSs with no introns (if there are any)
		  				$all_species{$species}[1]+= $introns;
		  
		  				if($introns == 0){
		  					$all_species{$species}[2]++;
		  				}    
					}
	      		}
	    	}
		}
    }
    close(FILE) || die "Can not close file\n";
}

print "\n\nProcessed $entries GenBank entries\n\n";

open(OUT,">species_CDS_intron_count.csv") || die "Can't open output file\n";
print OUT "Species,CDS_count,Intron_count,Single_exon_CDS_count,CDSs_with_introns_but_missing_ends\n";

foreach my $species (sort(keys(%all_species))){
	# Only print entries with at least 100 CDSs?
	if($all_species{$species}[0] >= 100){
		
		# reset third counter to zero if there were no data for this category
		($all_species{$species}[3] = 0) if (!defined($all_species{$species}[3]));

		print OUT "$species,$all_species{$species}[0],$all_species{$species}[1],$all_species{$species}[2],$all_species{$species}[3]\n";     
	}    
}

close(OUT) || die "Can't close file\n";

exit(0);
