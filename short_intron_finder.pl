#!/usr/bin/perl -w
#
# short_intron_finder.pl
#
# Simple script to report on which entries in GenBank have suspiciously short introns
#
# by Keith Bradnam
#
# 26th October 2005
#
#######################################################################################


use strict;
use Getopt::Long;

########################
# Command line options
########################
my $release;   # which release of genbank to query against?  Specify an integer

GetOptions("release=i"=> \$release);

die "Must specify a GenBank release number to query against\n" if (!defined($release));


########################
# misc variables
########################

# counters to keep track of short exons/introns (key is size, value is count)
my %intron_warnings;
my %species2intron_warnings;
my $negative_introns = 0;
my $intron_counter;

# specify path to GenBank release
my $path = "/GenBank/genbank${release}";

die "$path directory does not exist.\n" if (! -e "$path");


# have list of all GenBank files that we will be interested in (invertebrates, mammals, plants,
# primates, rodents, and vertebrates
my @files = ("inv1", "inv2", "inv3", "inv4", "inv5", "inv6", "inv7","inv8",
	     "mam1", "mam2",
	     "pln1", "pln2", "pln3", "pln4", "pln5", "pln6", "pln7", "pln8", "pln9", "pln10", "pln11",
	     "pln12","pln13", "pln14", "pln15", "pln16",
	     "pri1", "pri2", "pri3", "pri4", "pri5", "pri6", "pri7", "pri8", "pri9", "pri10", "pri11",
	     "pri12","pri13", "pri14", "pri15", "pri16", "pri17", "pri18", "pri19", "pri20", "pri21",
	     "pri22","pri23", "pri24", "pri25", "pri26", "pri27", "pri28","pri29",
	     "rod1", "rod2", "rod3", "rod4", "rod5", "rod6", "rod7", "rod8", "rod9", "rod10", "rod11",
	     "rod12","rod13", "rod14", "rod15", "rod16", "rod17", "rod18","rod19","rod20","rod21",
	     "vrt1", "vrt2", "vrt3", "vrt4", "vrt5", "vrt6", "vrt7", "vrt8", "vrt9");

# get list of wgs files
my @more_files = glob("${path}/wgs/wgs.*.gbff");

# combine into one array
@files = (@files,@more_files);

#@files = ("inv1", "inv2", "inv3", "inv4", "inv5", "inv6", "inv7","inv8");
#@files = ("inv1");



############################################################
#
#  M A I N   L O O P
#
#  loop through GenBank files extracting data
#
############################################################

# Change record delimiter to split on //
# newlines needed to avoid // in URLs for example
$/ = "\n//\n";


while (my $file = shift(@files)){

    print "$file\n";
    
    # different file opening routines depending on name of file
    if($file =~ m/gbff$/){
	open (FILE,"<$file") || die "Can't open gb{$file}.seq\n";
    }
    else{
	open (FILE,"<${path}/gb${file}.seq") || die "Can't open ${path}/gb${file}.seq\n";
    }
    
    ENTRY: while (<FILE>) {
	
	my ($locus,$mol,$accession,$species);
	
	# skip to start of first record if at the beginning of a file
	if(m/^GB\S+\.SEQ/){
	    m/LOCUS\s+(\S+)\s+\d+ bp\s+(\S+)\s+/ || die "No LOCUS field found for:\n\n$_\n"; 
	    $locus = $1;
	    $mol  = $2;
	}
	else{
	    m/^LOCUS\s+(\S+)\s+\d+ bp\s+(\S+)\s+/ || die "No LOCUS field found for:\n\n$_\n"; 
	    $locus = $1;
	    $mol  = $2;
	}       
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
		    

		    # match at least two coding exons to find entries with at least one intron
		    next unless ($cds =~ m /\.\.\d+,\d+\.\.\d+/);
		    
		    # to find internal occurances, look for three consecutive digits (no ..)
		    next if($cds =~ m/\d+,\d+,\d+/);		   
		    
		    # some entries use > & < internally of the first and last exon, ignore these entries altogether
		    next if ($cds =~ m/</);
		    next if ($cds =~ m/>/);

		    
                    # now get intron sizes
		    my @coords;
		    
		    # split on ..  character to get sets of introns
		    @coords = split(/\.\./,$cds);

		    # now loop through exon coordinates to get sizes of adjacent introns and exons
		  COORDS: for(my $i=0; $i<@coords; $i++){
		      
		      #skip single exons
		      next COORDS unless ($coords[$i] =~ m/,/);

		      
		      # calculate intron size
		      $coords[$i] =~ m/^(\d+),(\d+)$/;
		      my $exon_end = $1;
		      my $next_exon_start = $2;
		      my $intron_size = $next_exon_start-$exon_end-1;
#		      print " $intron_size\n" if ($intron_size < 10);		      
		      $intron_warnings{$intron_size}++ if ($intron_size < 11);
		      $negative_introns++ if($intron_size <0);
		      $species2intron_warnings{$species}++ if ($intron_size <11);
		      $intron_counter++;
		  }		    
		}
	    }
	}
    }
    close(FILE) || die "Can not close file\n";
}



# print some stats
print "\n";

print "$negative_introns introns at < 0 bp\n";
foreach my $i (0..10){
    if(defined($intron_warnings{$i})){
	print "$intron_warnings{$i} introns at $i bp\n";
    }
    else{
	print "0 introns at $i bp\n";
    }
}

print "\n";
foreach my $species (sort(keys(%species2intron_warnings))){
    print "$species $species2intron_warnings{$species}\n" if ((defined($species2intron_warnings{$species})) && ($species2intron_warnings{$species} > 10));
}

print "\nTotal number of introns processed: $intron_counter\n";

exit(0);
