#!/usr/bin/perl -w
#
# genbank_intron_sizes_vs_position.pl
#
# a script to extract intron sizes for all genes in all organisms and average the size according to
# whether it is the 1st, 2nd, 3rd, nth intron.
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
my $limit;     # set a cut-off for how many introns/exons are needed per species

GetOptions("release=i"=> \$release,
	   "limit=i"  => \$limit);

die "Must specify a GenBank release number to query against\n" if (!defined($release));

#set default for $limit if not specified
$limit = 1000 if (!defined($limit));


########################
# misc variables
########################

# simple hash for keeping track of all species (regardless of whether you are dealing with exons or introns)
# values are just set to '1'
my %all_species;

# For each species will also want to keep an array (which will be sorted) of exons, introns, flanking exons
# and flanking introns.  Key is species, value is an array of sizes
my %species2introns;

# another couple of hashes for storing only first introns and other introns
my %species2first_introns;
my %species2other_introns;


# specify path to GenBank release
my $path = "/GenBank/genbank${release}";
#$path = glob("~keith");


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

# want to write some output to a file as we go (for scatterplots of exon/intron vs flanking intron/exon sizes)
#open(CSV,">introns_vs_exons_r${release}.csv") || die "Couldn't open introns vs exons output file\n";

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
		  
		  
		  # match CDSs with at least two coding exons (i.e one intron)
		  next unless ($cds =~ m /\.\.\d+,\d+/);
		  
		  # or even look for two intron genes at least
		  next unless ($cds =~ m /\.\.\d+,\d+\.\.\d+,\d+/);

		  # some entries use > & < internally of the first and last exon, remove these symbols if they are terminal exons
		  next if ($cds =~ s/^<//);
		  next if ($cds =~ s/>$//);

		  
		  # now check that there are always pairs of exon coordinates, e.g. avoid scenarios like
		  # this in accession AE003538 (has a final single exon start coordinate but with no end
		  # CDS             complement(join(290022..290246,290389..290873,
		  #                 290960..291199,291264))
		  # bit of a fudge, only looking for pair of coordinates at start and end of CDS
		  next unless ($cds =~ m/^\d+\.\.\d+/);
		  next unless ($cds =~ m/\d+\.\.\d+$/);
		  # to find internal occurances, look for three consecutive digits (no ..)
		  next if($cds =~ m/\d+,\d+,\d+/);		   
		  
		  # now get intron and exon sizes depending on which command line options were used
		  my @coords;
		  
		  # split on ..  character to get sets of exons		    
		  @coords = split(/\.\./,$cds);

		  # need to know how many introns in gene
		  my $intron_counter=0;
		  # now loop through exon coordinates to get sizes of adjacent introns and exons
		COORDS: for(my $i=0; $i<@coords; $i++){

		    # double check that we haven't got a single coordinate intron
		    next if ($coords[$i] !~ m/^\d+,\d+$/);
				  		 
		    $intron_counter++;
		    # maybe we don't really care if genes have more than 10 or so introns
		    # really just want to see size difference in first half a dozen
		    next FEATURE if ($intron_counter > 10);

		    # calculate intron size
		    $coords[$i] =~ m/^(\d+),(\d+)$/;
		    my $exon_end = $1;
		    my $next_exon_start = $2;
		    my $intron_size = $next_exon_start-$exon_end-1;
		    		    
		    # sanity check 1
		    # some exons appear butt-ended and some produce negative intron sizes, so
		    # always wise to check and just ignore spurious exon/introns
		    next FEATURE if ($intron_size <1);
		    
		    # add intron size to hash
		    push(@{$species2introns{$species}[$intron_counter]},$intron_size);
			
		    $all_species{$species} = 1;
		    
		    # If this is the first intron, add sizes to a separate key in hash
		    # else add all other positions to another separate key in hash
		    if($intron_counter == 1){
			push(@{$species2first_introns{$species}[$intron_counter]},$intron_size);			
		    }
		    else{
			push(@{$species2other_introns{$species}[$intron_counter]},$intron_size);
		    }
		}		    
	      }
	    }
	}
    }
    close(FILE) || die "Can not close file\n";
}


# hash
# species name is key, value is an array of average intron sizes at position N
my %species2means;


# another hash to keep track of count of introns in each species
my %species2count;

open(OUT,">species2intron_size.csv") || die "Can't open output file\n";

foreach my $species (sort(keys(%all_species))){
       
    # what's the max number of introns for any gene in that species?
    my $intron_positions = scalar(@{$species2introns{$species}});


    for(my $i=1; $i<$intron_positions;$i++){

	# calculate mean intron size at position N
	my $sum = sum(@{@{$species2introns{$species}}[$i]});
	my $introns = scalar(@{@{$species2introns{$species}}[$i]});
	# add number of introns to count
	$species2count{$species}+= $introns;
	my $mean = $sum/$introns;

	# for simplicity store the means in a new array (hash of arrays)
	$species2means{$species}[$i] = $mean;
    }
}




# Now work through new array to print average intron size at each intron position in each species
foreach my $species (sort(keys(%all_species))){

    # need to count all introns in each species and only print those above a threshold (set by $limit)
    if($species2count{$species} > $limit){
	print OUT "$species,$species2count{$species},$species,";
	for(my $i =1;$i<(@{$species2means{$species}});$i++){
	    my $mean = sprintf "%.1f",${$species2means{$species}}[$i];
	    print OUT "$mean,";
	}
	print OUT "\n";
    }
}

print OUT "\n";



# now print out a table of how many introns at each size in each species (we want to be careful if we
# are calculating an average size based on only a few introns)
foreach my $species (sort(keys(%all_species))){

    # need to count all introns in each species and only print those above a threshold (set by $limit)
    if($species2count{$species} > $limit){

	print OUT "$species,,$species,";
	# what's the max number of introns for any gene in that species?
	my $intron_positions = scalar(@{$species2introns{$species}});
	
	for(my $i=1; $i<$intron_positions;$i++){	    
	    # calculate mean intron size at position N
	    my $introns = scalar(@{@{$species2introns{$species}}[$i]});
	    print OUT "$introns,";
	}
	print OUT "\n";
    }
}
print OUT "\n";

close(OUT) || die "Can't close file\n";


exit(0);
