#!/usr/bin/perl
#
# frag_directory_generator.pl
#
# A script to generate symlinks to files in /Chanlab/Data/FRAG_project/EXP/* from /Chanlab/Data/FRAG_project/FRAG/*
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;

###################################################### 
#              Command-line options                  #
###################################################### 

my $usage = "
$0 <csv file containing EXP->FRAG ID information from Google Doc>

optional arguments: 
        --help
";

my ($help);

die $usage unless GetOptions (
            "help"            => \$help
            );

die $usage if ($help);
die $usage unless (@ARGV == 1);


# hard-coded path to where directories will be made...shouldn't need to change
my $path_prefix = "/Chanlab/Data/FRAG_project/FRAG";

# link targets
my $target_dir = "/Chanlab/Data/FRAG_project/EXP";


######################################
# Process CSV file
######################################

# should look like this

# Experiment,FRAG ID,Plant,Wild type parent,Seeds?,Barcode/Index,Comai/BiOO barcode reference,Old Sample ID
# EXP001,FRAG00001A,5,2x Wa-1,n,CTCTG,F4,5
# EXP001,FRAG00002A,6,2x Wa-1,n,TCAAG,H2,6
# EXP001,FRAG00003A,16,2x Wa-1,n,GCGTT,G4,16
# EXP001,FRAG00004A,18,2x Wa-1,n,TGATG,H4,18

while(<>){

	# skip header line
	next if ($_ =~ m/^Experiment/);

	# grab main identifiers, ignore rest (@stuff)
	my ($exp, $frag_id, $plant_id, @stuff) = split(/,/, $_);

	# FRAG dir (e.g. FRAG12345) will contain FRAG IDs (e.g. FRAG12345A)
	my $frag_dir = $frag_id;
	$frag_dir =~ s/[A-Z]$//;
			
	print "EXP = $exp\tFRAG_dir = $frag_dir\tFRAG = $frag_id\n";

	# does FRAG prefix directory already exist? If not create it
	unless (-e "$path_prefix/$frag_dir"){
		mkdir("$path_prefix/$frag_dir") or die "Can't make directory $path_prefix/$frag_dir";
	} 
	
	
	# does the full frag name directory exist?
	unless (-e "$path_prefix/$frag_dir/$frag_id"){
		mkdir("$path_prefix/$frag_dir/$frag_id") or die "Can't make directory $path_prefix/$frag_dir/$frag_id";
	}


	# can now make symbolic links to BAM, BAI, FQ etc files

	chdir("$path_prefix/$frag_dir/$frag_id") or die "Can't changedir to $path_prefix/$frag_dir/$frag_id";

	my $bai_file     = "$target_dir/$exp/BAI_files/$frag_id.bai.gz";
	my $bam_file     = "$target_dir/$exp/BAM_files/$frag_id.bam.gz";
	my $fq_file      = "$target_dir/$exp/FQ_files/$frag_id.fq.gz";
	my $pileup_file  = "$target_dir/$exp/Pileup_files/$frag_id.pileup.gz";
	my $sai_file     = "$target_dir/$exp/SAI_files/$frag_id.sai.gz";
	my $sam_file     = "$target_dir/$exp/SAM_files/$frag_id.sam.gz";

	foreach my $file ($bai_file, $bam_file, $fq_file, $pileup_file, $sai_file, $sam_file){
		unless(-e $file){
			system("ln -s $file .") and die "Can't make sym link to $file";
		}
	}

	# Now process image links (if a PL12345 style ID is present)
	# does Images subdirectory already exist? If not create it
	unless (-e "$path_prefix/$frag_dir/$frag_id/Images"){
		mkdir("$path_prefix/$frag_dir/$frag_id/Images") or die "Can't make directory $path_prefix/$frag_dir/$frag_id/Images";
	} 
	chdir("$path_prefix/$frag_dir/$frag_id/Images") or die "Can't changedir to $path_prefix/$frag_dir/$frag_id/Images";
	if ($plant_id =~ m/PL\d{5}/){
			my $target = "$target_dir/$exp/Images/$plant_id.JPG";
			system("ln -s $target .") and die "Can't make sym link to $target";	
	}
	

}

exit;
