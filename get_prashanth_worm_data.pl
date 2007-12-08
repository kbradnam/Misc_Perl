#!/usr/bin/perl
#
# get_prashanth_worm_data.pl
#
# a script to extract some data from C.elegans GFF files needed by Prashanth
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;
use Getopt::Long;

########################
# Command line options
########################

my $test;      # test mode, only test one chromosome for speed
my $confirmed; # should specify a file of confirmed gene IDs
my $out1; 	   # specifies gff like output of just some features
my $out2;	   # specifies vertical output, 1 line per nt with all features

GetOptions ("test"        => \$test,
			"confirmed=s" => \$confirmed,
			"out1"        => \$out1,
			"out2"        => \$out2);

#############
# Check options 
#############
die "Use -out1 OR -out2\n" if ($out1 && $out2);
die "Must specify -out1 or -out2 output options\n" if (!$out1 && !$out2);
die "Must specify a file of confirmed gene IDs with -confirmed option\n" if (!$confirmed);

#############
# Paths etc 
#############

my $dbdir       = "/Korflab/Data_sources/WormBase/WS180";    # Database path
my $gffdir      = "$dbdir/CHROMOSOMES";                      # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 
@chromosomes = qw (I) if ($test);

###################
# Misc. variables
###################

my $seq;              # will hold chromosome sequence as a string...
my @dna;              # ...which will then get split into an array...
my @dna2;			  #...and copied. This copy will store a hash of the information at each base
my @dna3;			  # ...this copy will store details of genes and which strand they are on
my $chr_length;       # obvious really
my $intergenic;       # will store sequence of an intergenic region
my $intergenic_start; # start coord
my $intergenic_end;   # stop coord
my $intergenic_size;  # length
my $previous;         # strandedness of gene (+ or -) 5' to intergenic region
my $next;             # strandedness of gene (+ or -) 3' to intergenic region



#########################################
# read list of confirmed genes into hash
#########################################

my %genes;

open(GENES,"<$confirmed") || die "Can't open confirmed genes file: $confirmed\n";
 while(<GENES>){
	chomp;
	$genes{$_} = 1;
}
close(GENES);


####################################
# Main loop through each chromosome
####################################

foreach my $chromosome (@chromosomes) {
    
    # get chromosome sequence, load into $seq string
	if($test){
		open (DNA, "<CHROMOSOME_I.dna") || die "Failed to open dna file\n\n";
	}
    else{
		open (DNA, "<$gffdir/CHROMOSOME_${chromosome}.dna") || die "Failed to open dna file\n\n";
    }
    $seq = "";
    while(my $tmp =<DNA>){
        chomp($tmp);
        # skip header line
        next if ($tmp =~ m/^>/);
        $seq .= $tmp;
    }
    close(DNA);
   
    # to make things easier to calculate, will add a character to left end of $seq and @dna, such that position 1 of dna
    # sequence becomes array element 1 (rather than zero).  Add 'S' for 'Start'
    $seq = "S".$seq;

    # split to an array, make a copy and calculate size
    @dna=();
    @dna = split(//,$seq);
	@dna3 = @dna;
    $chr_length = scalar(@dna)-1;

    #################################################
    # GFF step: extract details of genes
    #################################################

    if($test){
		open (GFF_in, "<CHROMOSOME_I.gff") || die "Failed to open gff file\n\n";
	}
    else{
		open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
	}
	
    while (<GFF_in>) {
		chomp;
	
		# skip header info and process each line
		next if m/^\#/;
		my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split /\t/;
		my $length = $stop-$start+1;
		
		# confirmed intron data will match the following
		# GFF_source = curated
		# GFF_feature = exon
		# Comment contains CDS name that matches a confirmed gene list
			
		# process gene entries
		if (($source eq "gene") && ($feature eq "gene")){
			&annotate("G",$start,$stop);
			
			# if this gene is confirmed, annotate that as a separate thing
			$comment =~ s/Gene \"(.*)\"//;
			my $id = $1;
			if ($genes{$id}){
				annotate("+",$start,$stop);
				print "$chromosome Confirmed_gene $start $stop $length\n" if ($out1);
			}
			
			# substitute the appropriate part of @dna3 and with pluses (+) for forward genes
			# and minuses (-) for reverse genes.	
			my @replace = ("$strand") x $length;	
			splice(@dna3,$start,$length, @replace);
		}
		# process intron entries
		elsif(($source eq "Coding_transcript") && ($feature eq "intron")){
			&annotate("I",$start,$stop);
			# print details for confirmed introns
			print "$chromosome Confirmed_intron $start $stop $length\n" if ($comment =~ m/Confirmed_/ && $out1);
		}
		# process exon entries
		elsif(($source eq "Coding_transcript") && ($feature eq "exon")){
			&annotate("E",$start,$stop);
		}
		# process protein-coding exon entries
		elsif(($source eq "Coding_transcript") && ($feature eq "coding_exon")){
			&annotate("P",$start,$stop);
		}
		# process 5' UTR entries
		elsif(($source eq "Coding_transcript") && ($feature eq "five_prime_UTR")){
			&annotate("5",$start,$stop);
		}
		# process 3' UTR entries
		elsif(($source eq "Coding_transcript") && ($feature eq "three_prime_UTR")){
			&annotate("3",$start,$stop);
		}
		# process RepeatMasker repeats
		elsif(($source eq "RepeatMasker") && ($feature eq "repeat_region")){
			&annotate("R",$start,$stop);
		}
		# process tandem and inverted repeats 
		elsif(($source eq "tandem") && ($feature eq "tandem_repeat")){
			&annotate("L",$start,$stop);
		}
		elsif(($source eq "inverted") && ($feature eq "inverted_repeat")){
			&annotate("L",$start,$stop);
		}
		# skip if we get to this point
		else{
			next;
		}
		
    }
    close(GFF_in);
	    
    ########################################
    # Main loop through chromosome sequence
    ########################################
    
    # need a marker to be able to move past first non-coding (telomeric) sequence
    my $flag = 0;
	
	# for the 3' telomere calculation, need to keep track of the end of hte last gene
	my $previous_gene_end;
	
    SEQ: for (my $i = 1; $i<$chr_length;$i++){

		# want to ignore the first stretch of non-coding sequence at the telomere
		# as this is not strictly intergenic sequence
		if($flag == 0){
	    	while(($dna3[$i] ne "+") && ($dna3[$i] ne "-")){
				$i++;
	    	}
	    	# now must be in the first coding region so can change $flag value and print out
			# telomere coordinates
	    	$flag = 1;
			my $j = $i-1;
			print "$chromosome Telomere 1 $j\n" if($out1);
			&annotate("T",1,$j);
		}	

		# keep skipping forwards until you enter an intergenic region, i.e. not a plus or minus
		if (($dna3[$i] eq "+") || ($dna3[$i] eq "-")){
			$previous_gene_end = $i;
			next SEQ;	
		}

	    $intergenic_start = $i;

	    # what direction was previous gene?
	    $previous = $dna3[$i-1];
	    
        # now start counting until you leave intergenic region
	    $intergenic_size = 0;

	    while(($dna3[$i] ne "+") && ($dna3[$i] ne "-")){
			$intergenic_size++;
			$i++;
			# Need to quit if we are at the end of the chromosome so
			# last telomeric region doesnt't get counted as intergenic
			if ($i == $chr_length){
				$previous_gene_end++;
				print "$chromosome Telomere $previous_gene_end $i\n" if($out1);
				&annotate("T",$previous_gene_end,$i);
				last SEQ;
			}
	    }

	    # calculate end coordinate of intergenic region
	    $intergenic_end = $i-1;
	    
	    # what is direction of next gene?
	    $next = $dna3[$i];
		
		# work out orientation of genes
		my $symbol;
		($symbol = "T") if ($previous eq $next);
		($symbol = "C") if ($previous eq "+" && $next eq "-");
		($symbol = "D") if ($previous eq "-" && $next eq "+");
		
		print "$chromosome Intergenic $intergenic_start $intergenic_end $intergenic_size $symbol $previous$next\n" if($out1);
		&annotate("N",$intergenic_start,$intergenic_end);
	}
	    
	# print vertical output (unless -out1 specified)
	exit unless ($out2);
	
	foreach my $i (1..$chr_length){
		my @keys = keys (%{$dna2[$i]});
		print "$i) $dna[$i] @keys \n";
	}
}

# subroutine that will annotate the virtual sequence with characters at each base
sub annotate{
	my $type = shift;
	my $start = shift;
	my $stop = shift;
		
	foreach my $i ($start..$stop){
		$dna2[$i]{$type} = 1;
	}
}

exit(0);

