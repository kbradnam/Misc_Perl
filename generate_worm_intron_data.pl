#!/usr/bin/perl
#
# generate_worm_intron_data.pl
#
# a script to generate intron FASTA files with IME data in header from worm GFF files
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;
use FAlite;


#############
# Paths etc 
#############

my @chromosomes = qw( I II III IV V X );                     
@chromosomes = qw (chr1.gff);     


#####################################################################
# First loop through each chromosome to just get transcript details
# and to find out confirmation status of CDS
#####################################################################

foreach my $chromosome (@chromosomes) {
    
	
	###################
	# Misc. variables
	###################

	my %transcript2start;   # transcript ID is key, transcript start coord is value
	my %transcript2strand;  # transcript ID is key, strand is the value 
	my %transcript2introns; # transcript ID is key, value is an array of intron coords
	my %cds2status;	        # CDS ID is key, confirmation status is value
	my %cds2start;		    # need start and stop coordinates of each CDS
	my %cds2stop;

	
#	open (GFF, "<CHROMOSOME_${chromosome}.gff") || die "Failed to open dna file\n\n";
	open (GFF, "<$chromosome") || die "Failed to open GFF file\n\n";

    while(my $tmp =<GFF>){
        
		# skip header lines
		next if ($tmp =~ m/#/);
		
		# split to get details
		chomp($tmp);
		my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split (/\t/,$tmp);
        
		# only want lines that contain protein_coding_primary_transcript_details
        if ($feature eq "protein_coding_primary_transcript"){					
			
			# get transcript ID from $comment and skip if transcript is an alternative splice form
			# i.e. keep *.a versions and *.n.1
			$comment =~ m/\"(.*)\"/;
			my $id = $1;
			
			# skip alternative CDS isoforms of genes
			next if ($id =~ m/[b-z]$/);
			next if ($id =~ m/[b-z]\.[0-9]*$/);

	
			# skip alternative transcript isoforms of genes
			next if ($id =~ m/\.[0-9]\.[2-9]$/);

			# finally, skip alternative transcript isoforms of a form of CDS variants
			next if ($id =~ m/a\.[2-9]$/);
			
			# now add start coord and strand to hash
			# use stop coordinate for reverse strand genes
			($transcript2start{$id} = $start) if ($strand eq "+");
			($transcript2start{$id} = $stop) if ($strand eq "-");
			
			$transcript2strand{$id} = $strand;
		}
		elsif(($source eq "curated") && ($feature eq "CDS")){
			# grab confirmation status
			$comment =~ m/CDS \"(.*?)\".*Status \"(.*?)\"/;
			$cds2status{$1} = $2;
			
			# add coords to hash but reverse for negative strand
			($cds2start{$1} = $start) if ($strand eq "+");
			($cds2stop{$1} = $stop)   if ($strand eq "+");
			($cds2start{$1} = $stop)  if ($strand eq "-");
			($cds2stop{$1} = $start)  if ($strand eq "-");


		}
		else{
			next;
		}
    }
    close(GFF);
  
	# now do 2nd loop through each chromosome to just get intron details
    	
#	open (GFF, "<CHROMOSOME_${chromosome}.gff") || die "Failed to open dna file\n\n";
	open (GFF, "<$chromosome") || die "Failed to open GFF file\n\n";

    while(my $tmp =<GFF>){
        
		# skip header lines
		next if ($tmp =~ m/#/);
		
		# split to get details
		chomp($tmp);
		my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split (/\t/,$tmp);
        
		# only want lines that contain protein_coding_primary_transcript_details
        next unless (($source eq "Coding_transcript") && ($feature eq "intron"));					
		
		# grab transcript ID
		$comment =~ m/Transcript \"(.*?)\"/;
		my $id = $1;
			
		# skip introns from alternative CDS isoforms of genes
		next if ($id =~ m/[b-z]$/);
		next if ($id =~ m/[b-z]\.[0-9]*$/);

		# skip introns from alternative transcript isoforms of genes
		next if ($id =~ m/\.[0-9]\.[2-9]$/);

		# finally, skip introns from alternative transcript isoforms of a form of CDS variants
		next if ($id =~ m/a\.[2-9]$/);	
			
		# transcript ID will sometimes be longer than the CDS ID. E.g. Transcript=ZK270.2a.2 CDS=ZK270.2a
		# so need to delete end of transcript ID for some lookups of CDS related data
		my $id2 = $id;
		($id2 =~ s/\.[0-9]+$//) if ($id2 =~ m/\..*\./);
		
		# if this is a confirmed gene, then we can store the intron data, reverse coords for - strand genes
		if($cds2status{$id2} eq "Confirmed"){				
			push(@{$transcript2introns{$id}},"$start:$stop")  if ($strand eq "+");				
			push(@{$transcript2introns{$id}},"$stop:$start")  if ($strand eq "-");				
		}
		
	}
    close(GFF);
  

	# get chromosome sequence, load into $seq string
#    open (DNA, "<CHROMOSOME_${chromosome}.dna") || die "Failed to open dna file\n\n";
    open (DNA, "<chr1.dna") || die "Failed to open dna file\n\n";
	
    my $seq;

    while(my $tmp =<DNA>){
        chomp($tmp);
        # skip header line
        next if ($tmp =~ m/^>/);
        $seq .= $tmp;
    }
    close(DNA);


	#####################################################################
	# Now loop through all transcripts, processing each intron in turn
	#####################################################################

	foreach my $key (keys %transcript2introns){

	#	print "${key}_\n";
	#	print "Start: $transcript2start{$key} Strand: $transcript2strand{$key}\n";

		# transcript ID will sometimes be longer than the CDS ID. E.g. Transcript=ZK270.2a.2 CDS=ZK270.2a
		# so need to delete end of transcript ID for some lookups of CDS related data
		my $id2 = $key;
		($id2 =~ s/\.[0-9]+$//) if ($id2 =~ m/\..*\./);
	#	print "CDS start: $cds2start{$id2} CDS stop: $cds2stop{$id2}\n";

		# sort introns into ascending (or descending) coordinates
		my @introns = sort (@{$transcript2introns{$key}});
		(@introns = reverse @introns) if ($transcript2strand{$key} eq "-");

		# count order of introns within each transcript
		my $counter = 1;


		foreach my $intron (@introns){

			my ($start,$stop) = split(/:/,$intron);

			# initially assume that intron is in CDS and only change if it is an UTR
			my $status = "CDS";

			# now flag status if intron is in UTR, calculate distance to TSS and get sequence
			# treat differently depending on strand
			my $distance;			
			my $sequence;
			
			if($transcript2strand{$key} eq "+"){
				$status = "5UTR" if ($start < $cds2start{$id2});
				$status = "3UTR" if ($start > $cds2stop{$id2});
				$distance = ($start - $transcript2start{$key});
				$sequence = substr($seq,$start-1,$stop-$start+1);
			}
			else{
				$status = "5UTR" if ($start > $cds2start{$id2});
				$status = "3UTR" if ($start < $cds2stop{$id2});
				$distance = ($transcript2start{$key} - $start);
				
				# need to reverse complement sequence
				my $tmp_seq = substr($seq,$stop-1,$start-$stop+1);
				$sequence = reverse $tmp_seq;
			    $sequence =~ tr/ACGTacgt/TGCAtgca/;
				
			}
			# format sequence for printing
			my $output_seq = &tidy_seq($sequence);
			
			print ">${key}_i${counter}_${distance}_${status} $start-$stop\n";
			print "$output_seq\n";
			$counter++;
		}	
	}

}


exit(0);


sub tidy_seq{
	# adds a new line character every 60 bases  
    my $seq = shift;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;

	# how many newlines to add?
    my $to_add = int(length($seq)/60);
    
    my $output_seq;	
    my ($start,$end);
    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
	# add remaining sequence (if it exists)
    $output_seq .= substr($seq,$start);

    return ($output_seq);
}

