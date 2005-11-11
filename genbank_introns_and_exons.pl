#!/usr/bin/perl -w
#
# genbank_introns_and_exons.pl
#
# a script to extract intron (and exon) sizes from the CDS features of a GenBank entry
#
# by Keith Bradnam
#
# created 26th July 2005
#
# Last updated by: $Author$
# Last updated on: $Date$    
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


#######################################
# misc variables and data structures
#######################################

# species name is key, 1st value is either 'INTRON' or 'EXON' which is a hash key in turn
# second hash values are an array of counts of introns/exons at different sizes
my %species2count;

# species name is key, 1st value is either 'INTRON' or 'EXON' which is a hash key in turn
# secoond hash value is the cumulative number of bases for introns/exons in that species
# can be later used to calculate average intron/exon size per species in conjunction with
# %species2count
my %species2size;

# species name is key, 1st value is either 'INTRON' or 'EXON' which is a hash key in turn
# second hash value is count of total number of introns/exons in that species
my %species2sumsize;

# simple tracker for all species names (there may be some species with exon but not intron data
# and vice versa)
my %all_species;


# Change record delimiter to split on //
# newlines needed to avoid // in URLs for example
$/ = "\n//\n";


# specify main part of path, maybe this should be a command line option?
my $path = "/GenBank/genbank${release}";

die "/GenBank/genbank${release} directory does not exist.\n" if (! -e "/GenBank/genbank${release}");


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



while (my $file = shift(@files)){

    print "$file\n";
    
    # different file opening routines depending on name of file
    if($file =~ m/gbff$/){
	open (FILE,"<$file") || die "Can't open $file\n";
    }
    else{
	open (FILE,"<${path}/gb${file}.seq") || die "Can't open gb{$file}.seq\n";
    }

    while (<FILE>) {
	
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
		/ORGANISM\s+(\S+\s+\S+)/;
		$species = $1;
		# want to only look at Eukaryotes
		next unless (m/Eukaryota/);
	    }
	    
	    
	    # split up each feature object
	    if (/^FEATURES/) {
		s/\n\s{5}(\S+)/ZZZZ$1/g;
		
#	    print "\n$locus $mol $species\n";
		
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

		    

		    # Need to avoid CDS features which might be part of a UTR exon.  E.g. the start codon
		    # (and hence position 1 of the CDS) might be at the end of the first exon in a gene.
		    # To be sure of getting real exons sizes, need to only use CDSs features which describe
		    # at least three exons. The outer two exons must then always be excluded.  But can still
		    # use two exon CDSs to calculate intron size.

		    # Skip anything which doesn't have at least two exons (nnn..nnn,nnn..nnn)
		    next unless ($cds =~ m /\.\..*\.\./);

		    # now can ignore first and last exon altogether, only remove one of two coordinates
		    # as will need the other to calculate intron size
		    $cds =~ s/^\d+/@@@/;
		    $cds =~ s/^<\d+/@@@/;
		    $cds =~ s/\d+$/@@@/;
		    $cds =~ s/>\d+$/@@@/;

	    
		    # now check that there are always pairs of exon coordinates, e.g. avoid scenarios like
		    # this in accession AE003538 (has a final single exon start coordinate but with no end
		    # CDS             complement(join(290022..290246,290389..290873,
                    #                 290960..291199,291264))
		    # bit of a fudge, only looking for pair of coordinates at start and end of CDS
		    next unless ($cds =~ m/^\@\@\@\.\.\d+/);
		    next unless ($cds =~ m/\d+\.\.\@\@\@$/);
		    # to find internal occurances, look for three consecutive digits (no ..)
		    next if($cds =~ m/\d+,\d+,\d+/);		   
		    		    
		    # some entries use > & < internally of the first and last exon, ignore these entries altogether
		    next if ($cds =~ m/</);
		    next if ($cds =~ m/>/);
		    
		    # now get intron and exon sizes
		    my @coords;
		    
		    # split on , character to get exon start & stop coordinates, introns will be predicted based on exon coordinates	    
		    @coords = split(/,/,$cds);

		    for(my $i=0; $i<@coords; $i++){

			# skip first exon if we don't know the true start, i.e. might be UTR exon
			next if ($coords[$i] =~ m/@@@\.\./);
			
			
			# calculate exon size
			my ($exon_start,$exon_end,$exon_size);
			
			if($coords[$i] =~ m/^(\d+)\.\.(\d+)$/){			   
			    $exon_start = $1;
			    $exon_end = $2;
			    $exon_size = $exon_end-$exon_start+1;

			    # SANITY CHECK 1: some exons appear butt-ended and some produce negative intron sizes, so
			    # always wise to check and just ignore spurious exon/introns
			    next FEATURE if ($exon_size <1);

			    # increment exon and species counters
			    $all_species{$species} = 1;
			    $species2count{"EXON"}{$species}++;
			    $species2sumsize{"EXON"}{$species} += $exon_size;
			    ${$species2size{"EXON"}{$species}}[$exon_size]++;

			}
			# exception for last exon (which might be a 3' UTR exon).  We only need exon start coordinate
			# in order to calulate preceding intron size.
			elsif($coords[$i] =~ m/^(\d+)\.\.\@\@\@$/){
			    $exon_start = $1;
			}

			# now work backwards to work out preceding intron size
			# (note different way of calculating size compared to exons)
			my $preceding_intron = $coords[$i-1];
			$preceding_intron =~ s/[\d@]+\.\.//;
			my $preceding_intron_size = $exon_start - $preceding_intron -1;

			# SANITY CHECK 2
			next FEATURE if ($preceding_intron_size <1);
			
			# increment intron and species counters
			$all_species{$species} = 1;		       
			$species2count{"INTRON"}{$species}++;
			$species2sumsize{"INTRON"}{$species} += $preceding_intron_size;
			${$species2size{"INTRON"}{$species}}[$preceding_intron_size]++;

		    }
		}		
	    }
	}
    }
    close(FILE);
}

# open output file (which will contain both exon and intron data)
open(OUT,">exon_intron_genbank${release}.csv") || die "Couldn't open output file\n";


# List of things which we will query and filter on for both introns and exons
my ($min,$max,$interval,$average,$mode);

# big wrapping loop to separately process exons and then introns
foreach my $type ("EXON","INTRON"){

    if($type eq "EXON"){
	$min = 1; $max = 750; $interval = 9;
    }

    elsif($type eq "INTRON"){
	$min = 1; $max = 750; $interval = 3;
    }

    # print header information for csv file
    print OUT "$type,Species,N,Average,Mode,Species,";       
    my $start = $min;
    while($start < $max){
	my $end = $start+$interval-1;
	print OUT "$start-$end,";   
	$start+=$interval;
    }
    print OUT "\n";

    
    # now get the actual intron/exon data, cycliung through each species
    foreach my $key (keys(%all_species)){

	# ignore species with less than 500 exons/introns
	next if (!defined($species2count{$type}{$key}));
	next if ($species2count{$type}{$key} < 500);

	print OUT "$type,$key,$species2count{$type}{$key},";

	# now calculate and print average and modal exon/intron size
	$average = $species2sumsize{$type}{$key} / $species2count{$type}{$key};
	$mode = &mode($key,$type);
	print OUT "$average,$mode,$key,";


	# now print counts of introns/exons (need to sum up counts in each interval size)
	for(my $i =$min; $i <$max; $i+=$interval){
	    my $count=0;

	    my $bin_counter = 0;	    
	    # loop through all possible sizes in each interval range
	    while($bin_counter < $interval){
		if(defined(${$species2size{$type}{$key}}[$i+$bin_counter])){
		    $count += "${$species2size{$type}{$key}}[$i+$bin_counter]";
		}
		$bin_counter++;
	    }
	    print OUT "$count,";	    
	}		
	print OUT "\n";
    }
}

close(OUT)   || die "Couldn't close output file\n";

exit(0);

    
# simple subroutine to calculate the modal value(s) from an array of integers

sub mode {
    my $species = shift;
    my $type = shift;
    my $max = 0;
    my $mode = 0;
    my $counter = 0;
    my @list;

    
    # loop through array of intron/exon sizes
    foreach(my $size = 0; $size <  @{$species2size{$type}{$species}}; $size++) {
	
	# are number of introns greater than max?
	if(defined(${$species2size{$type}{$species}}[$size]) &&  ${$species2size{$type}{$species}}[$size] >= $max ){
	    
	    # want to track multiple modes when count of intron/exons at a given size is same
	    # as another size
	    if(${$species2size{$type}{$species}}[$size] == $max ){
		push(@list,$mode,$size);
	    }
	    # but wipe array if new intron/exon size becomes new mode
	    elsif(${$species2size{$type}{$species}}[$size] > $max ){
		@list = ();		
	    }
	    
	    # if so change the mode to reflect this
	    $mode = $size;
	    
	    # now change max value that we need to beat to set a new mode
	    $max = ${$species2size{$type}{$species}}[$size];  
	}
    }


    print "*** $species has mutiple modes for $type : @list\n" if (defined($list[0]));	

    # return single value if one mode, or set of joined values if multiple modes
    if(defined($list[0])){
	my $modes = join("_",@list);
	return($modes);
    }
    else{
	return($mode);
    }
    
}
