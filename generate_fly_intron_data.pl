#!/usr/bin/perl
#
# generate_fly_intron_data.pl
#
# a script to generate intron FASTA files with IME data in header from fly GFF files
#
# Last updated by: $Author$
# Last updated on: $Date$
#
#############################################################################################

use strict;
use warnings;
use FAlite;


###################
# Misc. variables
###################
 
# store path to Drosophila chromosome files
my $dir = "/Korflab/Genomes/Sequences/Drosophila_melanogaster/r5.4";

my %transcript2start;   # transcript ID is key, transcript start coord is value
my %transcript2strand;  # transcript ID is key, strand is the value 
my %transcript2name;   # transcript ID is key, gene name is the value 
my %transcript2score;   # transcript ID is key, score is the value 
my %transcript2introns; # transcript ID is key, value is an array of intron coords
my %protein2start;		# need start and stop coordinates of each protein
my %protein2stop;


my @chromosomes = qw (dmel_mitochondrion_genome 2RHet 2L X 3L XHet Uextra 4 YHet U 2LHet 2R 3LHet 3R 3RHet);

foreach my $chr (@chromosomes){

	open (GFF, "<dmel-${chr}-r5.4.gff") || die "Failed to open $chr\n\n";

	while(my $tmp =<GFF>){

		# skip header lines
		next if ($tmp =~ m/#/);
		
		# skip sequence files at end of file	
		last if ($tmp =~ m/^>/);
		
		# exit if you get to FASTA sequence
		last if	($tmp =~ m/##FASTA/);

		# split to get details
		chomp($tmp);
		my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split (/\t/,$tmp);
	
		# only want lines that contain mRNA details
	    if (($source eq "FlyBase") && ($feature eq "mRNA")){					

			# get transcript ID, gene name and score from $comment 
			$comment =~ m/ID=(FBtr.*?);Name=(.*?);/;
			my $id = $1;
			my $name = $2;
			
			# skip alternative CDS isoforms of genes
			next if ($name =~ m/R[B-Z]$/);
						
			# need a score and need a score of at least 8 (the 'Strongly supported' category)
			my $score;
			
			if ($comment !~ m/.*score=(\d+);/){
				next;
			}
			else{
				$score = $1;
				next if ($score < 8);
				$transcript2score{$id} = $score;
			}
			
			# keep track of gene name
			$transcript2name{$id} = $name;

			# now add start coord and strand to hash
			# use stop coordinate for reverse strand genes
			($transcript2start{$id} = $start) if ($strand eq "+");
			($transcript2start{$id} = $stop)  if ($strand eq "-");

			$transcript2strand{$id} = $strand;
		}
		elsif(($source eq "FlyBase") && ($feature eq "protein")){

			# get ID
			$comment =~ m/;Derives_from=(FBtr.*?);/;
			# add coords to hash but reverse for negative strand
			($protein2start{$1} = $start) if ($strand eq "+");
			($protein2stop{$1} = $stop)   if ($strand eq "+");
			($protein2start{$1} = $stop)  if ($strand eq "-");
			($protein2stop{$1} = $start)  if ($strand eq "-");
		}
		else{
			# just skip everything else in this run
			next;
		}
	}
	close(GFF);
	
}


# now do 2nd loop through file to just get intron details

foreach my $chr (@chromosomes){
		
	# first empty hash
	%transcript2introns = (); 
	
	open (GFF, "<dmel-${chr}-r5.4.gff") || die "Failed to open $chr file\n\n";

	while(my $tmp =<GFF>){

		# skip header lines
		next if ($tmp =~ m/#/);

		# exit if you get to FASTA sequence
		last if	($tmp =~ m/##FASTA/);

		# skip sequence files at end of file	
		last if ($tmp =~ m/^>/);
			
		# split to get details
		chomp($tmp);
		my ($location,$source,$feature,$start,$stop,$score,$strand,$phase,$comment) = split (/\t/,$tmp);

		# only want lines that contain protein_coding_primary_transcript_details
	    next unless (($source eq "FlyBase") && ($feature eq "intron"));					

		# get transcript ID, gene name and score from $comment 
		$comment =~ m/Parent=(FBtr.*?);/;
		my $id = $1;


		# if this is intron is from a high scoring transcript (with a cDNA)
		# then we can store the intron data, reverse coords for - strand genes
		if($transcript2start{$id}){				
			push(@{$transcript2introns{$id}},"$start:$stop")  if ($strand eq "+");				
			push(@{$transcript2introns{$id}},"$stop:$start")  if ($strand eq "-");				
		}
	}
	close(GFF);

	# now loop through chromosome sequences 
	# read sequence to one variable
	open (DNA, "<dmel-${chr}-chromosome-r5.4.fasta") || die "Failed to open $chr dna file\n\n";
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

		# transcript ID will sometimes be longer than the CDS ID. E.g. Transcript=ZK270.2a.2 CDS=ZK270.2a
		# so need to delete end of transcript ID for some lookups of CDS related data
		my $id2 = $key;
		($id2 =~ s/\.[0-9]+$//) if ($id2 =~ m/\..*\./);
	
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
				$status = "5UTR" if ($start < $protein2start{$id2});
				$status = "3UTR" if ($start > $protein2stop{$id2});
				$distance = ($start - $transcript2start{$key});
				$sequence = substr($seq,$start-1,$stop-$start+1);
			}
			else{
				$status = "5UTR" if ($start > $protein2start{$id2});
				$status = "3UTR" if ($start < $protein2stop{$id2});
				$distance = ($transcript2start{$key} - $start);
				
				# need to reverse complement sequence
				my $tmp_seq = substr($seq,$stop-1,$start-$stop+1);
				$sequence = reverse $tmp_seq;
			    $sequence =~ tr/ACGTacgt/TGCAtgca/;
				
			}
			# format sequence for printing
			my $output_seq = &tidy_seq($sequence);
			
			print ">${key}_i${counter}_${distance}_${status} $start-$stop ($transcript2strand{$key})\n";
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

