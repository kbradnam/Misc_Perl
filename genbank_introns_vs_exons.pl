#!/usr/bin/perl -w
#
# genbank_introns_vs_exons.pl
#
# a script to extract intron (and exon) sizes from the CDS features of a GenBank entry
# only extract when there are at least three coding exons in a gene, in which case discard
# flanking exons (which may be part of UTR exons).  Output is a list of pairs of:
# i) intron size + average flanking exon size
# ii) exon size + average flanking intron size
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

# For each species, will want to keep track of number of introns and exons
# key will be species name, value will be another hash key (either EXON or INTRON)
# and final value will be count of exons (with flanking introns) or introns (with flanking exons)
my %species2count;

# For each species will also want to keep an array (which will be sorted) of exons, introns, flanking exons
# and flanking introns.  Key is species, value is an array of sizes
my %species2exons;
my %species2introns;
my %species2flanking_exons;
my %species2flanking_introns;

# For each possible intron/exon size want to keep a list of all the flanking exons/introns that occur
# E.g. keep list (an array) of all intron sizes that flank exons of size 47 bp.
my %introns2flanking_exons;
my %exons2flanking_introns;

# counters to show total number of exons & introns in all species
my $total_exon_count = 0;
my $total_intron_count = 0;


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
		    

		    # match at least three coding exons as terminal exons might be part of longer UTR exons
		    # hence coordinates are not telling you the full length of the entire exon
		    next unless ($cds =~ m /\.\.\d+,\d+\.\.\d+,\d+\.\./);
		    
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

		    
                    # now get intron and exon sizes depending on which command line options were used
		    my @coords;
		    
		    # split on ,  character to get sets of exons		    
		    @coords = split(/,/,$cds);

		    # now loop through exon coordinates to get sizes of adjacent introns and exons
		  COORDS: for(my $i=0; $i<@coords; $i++){
			
			# skip first exon if we don't know the true start, i.e. might be UTR exon
			next COORDS if ($coords[$i] =~ m/@@@/);

			# calculate exon size
			$coords[$i] =~ m/^(\d+)\.\.(\d+)$/;
			my $exon_start = $1;
			my $exon_end = $2;
			my $exon_size = $exon_end-$exon_start+1;

			# sanity check 1
			# some exons appear butt-ended and some produce negative intron sizes, so
			# always wise to check and just ignore spurious exon/introns
			next FEATURE if ($exon_size <1);

			
			# now work backwards to work out preceding intron size
			my $preceding_intron_end = $coords[$i-1];
			$preceding_intron_end =~ s/[\d@]+\.\.//;
			my $preceding_intron_size = $exon_start - $preceding_intron_end -1;
			
			# sanity check 2
			next FEATURE if ($preceding_intron_size <6);
					       	
			# and now work forward to calculate next intron size (there should always
			# be a next intron because we have selected CDSs features with at least three exons
			# Note different way of calculating intron size (with a -1 rather than +1) compared to exons
			my $next_intron = $coords[$i+1];
			$next_intron =~ s/\.\.[\d@]+//;
			my $next_intron_size = $next_intron - $exon_end -1;

			# sanity check 3
			next FEATURE if ($next_intron_size <6);

				
			# calculate average of flanking introns
#			my $average_intron_size = ($preceding_intron_size + $next_intron_size)/2;
#			print CSV "$species,EXON,$exon_size,$preceding_intron_size,$next_intron_size,$average_intron_size\n";

			# increment exon, flanking intron and species counters
			$total_exon_count++;
			$all_species{$species} = 1;

			push(@{$species2exons{$species}},$exon_size);
			push(@{$species2flanking_introns{$species}},$preceding_intron_size,$next_intron_size);
			push(@{$exons2flanking_introns{$species}{$exon_size}},$preceding_intron_size,$next_intron_size);

			
			# If there is another exon to come, then we can calculate intron size vs
			# flanking exons sizes (the opposite of above)
			if($coords[$i+1]){
			    # again we need to check that next exon isn't a possible UTR exon, if it is we move to the next feature
			    next FEATURE if ($coords[$i+1] =~ /@@@/);

			    # now calcuate next exon size
			    $coords[$i+1] =~ m/^(\d+)\.\.(\d+)$/;
			    my $next_exon_start = $1;
			    my $next_exon_end = $2;
			    my $next_exon_size = $next_exon_end-$next_exon_start+1;

			    # sanity check 4
			    next FEATURE if ($next_exon_size <1);
			    
			    # calculate average of flanking exons
#			    my $average_exon_size = ($exon_size + $next_exon_size)/2;
#			    print CSV "$species,INTRON,$next_intron_size,$exon_size,$next_exon_size,$average_exon_size\n";

			    # increment intron, flanking exon and species counters			    
			    $total_intron_count++;
			    $all_species{$species} = 1;

			    push(@{$species2introns{$species}},$next_intron_size);
			    push(@{$species2flanking_exons{$species}},$exon_size,$next_exon_size);
			    push(@{$introns2flanking_exons{$species}{$next_intron_size}},$exon_size,$next_exon_size);
			}

			# Now want to skip the next exon altogether so that we are not double counting exons and introns
			# i.e. in a five exon gene (four introns) we will always exclude exons 1 & 5 as they are terminal exons
			# so can compare:
			# exon 2 with introns 1+2
			# then skip to exon 4 (as exon 3 is flanked by intron 2 which has already been counted)
			# exon 4 with introns 3+4
			$i++;
		    }		    
		}
	    }
	}
    }
    close(FILE) || die "Can not close file\n";
}

#close(CSV) || die "Couldn't close introns vs exons file\n";


# print some stats
print "\n";
print "Total exons = $total_exon_count\n";
print "Total introns = $total_intron_count\n\n";

	
# Now need to get the 25% of shortest & longest exons/introns (for each species) and then get the associated
# flanking introns/exons and perform a Z-test to compare the means of the flanking introns/exons.
# E.g. is the average length of introns that flank the shortest exons in worm significantly different from
# the length of introns that flank the longest exons

open(OUT,">z_test_r${release}_restricted.csv") || die "Can't open z-test output file\n";

foreach my $species (sort(keys(%all_species))){
  TYPE: foreach my $type ("EXON","INTRON"){
      
      # sort array of exons/introns to get them in ascending size order
      my @sorted;
      my $size;
      
      if($type eq "EXON"){
	  # can only proceed if there are exons for that species
	  next TYPE if (!defined($species2exons{$species}[0]));
	  @sorted = sort{$a <=> $b;} @{$species2exons{$species}};	  
	  $size = @sorted;
	  # only show values when there are more exons/introns than the specified threshold
	  next TYPE unless ($size >$limit);      

      }
      elsif($type eq "INTRON"){
	  # can only proceed if there are introns for that species
	  next TYPE if (!defined($species2introns{$species}[0]));
	  @sorted = sort{$a <=> $b;} @{$species2introns{$species}};
	  $size = @sorted;
	  # only show values when there are more exons/introns than the specified threshold
	  next TYPE unless ($size >$limit);      
      }
      
      # want to find out positions of 25% and 75% quartiles in @sorted array, so just multiply array length by 0.25 and 0.75
      # need to subtract 1 to convert to array coordinates.  Can then extract that (as an array slice) to get two new arrays
      # which just have the first quartile of exon/intron sizes and then the last quartile
      my @lower_range = @sorted[0..sprintf("%.0f",$size * 0.25)-1];
      my @upper_range = @sorted[sprintf("%.0f",$size* 0.75)-1..$size-1];

      # only want unique values in @lower_range and @upper_range as we are looking up associated values in %introns2flanking_exons
      # and %exons2flanking_introns which have *all* associated values at any given exon/intron size.
      # if you didn't remove duplicate values in @sorted then you will double count all associated exon & intron sizes
      my %seen = ();
      my @unique_lower = grep { ! $seen{$_} ++ } @lower_range;
      my @unique_upper = grep { ! $seen{$_} ++ } @upper_range;
      
	  ;
      # For each exon/intron size in the first & last quartiles, need to get the associated flanking intron/exons and
      # for each of those new arrays, calculate N, mean, and standard deviation so as to be able to perform a z-test
      my (@short,$short_count,$short_mean,$short_stdev);
      my (@long,$long_count,$long_mean,$long_stdev);

#      print "@sorted\n";
#      print "@lower_range\n";
#      print "@unique_lower\n";

      if($type eq "INTRON"){	    
	  foreach my $intron_size (@unique_lower){
	      foreach my $flanking_exon (@{$introns2flanking_exons{$species}{$intron_size}}){
		  push(@short,$flanking_exon);
	      }
	  }
	  foreach my $intron (@unique_upper){
	      foreach my $flanking_exon (@{$introns2flanking_exons{$species}{$intron}}){
		  push(@long,$flanking_exon);
	      }
	  }
      }		    
      if($type eq "EXON"){		
	  foreach my $exon (@unique_lower){
	      foreach my $flanking_intron (@{$exons2flanking_introns{$species}{$exon}}){
		  push(@short,$flanking_intron);
	      }
	  }			  
	  foreach my $exon (@unique_upper){
	      foreach my $flanking_intron (@{$exons2flanking_introns{$species}{$exon}}){
		  push(@long,$flanking_intron);
	      }
	  }		
      }
      
      # z-test
      $short_count = @short;
      $short_mean  = sprintf("%.1f",sum(@short)/$short_count);
      $short_stdev = sprintf("%.1f",sqrt(sum(map {($_ - $short_mean) ** 2} @short) / ($short_count-1)));
      
      $long_count = @long;
      $long_mean  = sprintf("%.1f",sum(@long)/$long_count);
      $long_stdev = sprintf("%.1f",sqrt(sum(map {($_ - $long_mean) ** 2} @long) / ($long_count-1)));
      
      my $z = sprintf("%.2f",($short_mean - $long_mean)/sqrt(($short_stdev**2 / $short_count) + ($long_stdev**2 / $long_count)));
      
      my $critical;
      # $critical = 1.96; # P < 0.05
      $critical = 2.58; # P < 0.01
      # $critical = 3.29; # P < 0.001
      # $critical = 3.89; # P < 0.0001
      # $critical = 4.41; # P < 0.00001
      
      
#      if(abs($z) > $critical){
#	  print OUT "$species : Z = $z ($size ${type}s)\n";
#	  print OUT "Lower:  N,mean,Stdev = $short_count,$short_mean,$short_stdev\t";
#	  print OUT "Higher: N,mean,Stdev = $long_count,$long_mean,$long_stdev\n\n";	  	
#      }
      print OUT "$type,$species,$z,$size,$short_count,$long_count,$short_mean,$long_mean,$short_stdev,$long_stdev\n";      
  }
}

close(OUT) || die "Can't close file\n";


exit(0);
