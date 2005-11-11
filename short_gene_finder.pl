#!/usr/local/bin/perl5.8.0 -w
#
# find_short_genes.pl  
# 
# by Keith Bradnam, aged 12 and half
#
# A script to find (and classify) potentially short, spurious genes (<100 aa)
#
# Last updated by: $Author$     
# Last updated on: $Date$     


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;
use Getopt::Long;


##############################
# misc variables             #
##############################

my $cutoff = 100;        # what length genes do we want to look at, cutoff is max length in amino acids
my $maintainers = "All";  # who will receive log file
my $ws_version   = &get_wormbase_version_name;



##############################
# command-line options       #
##############################

my $help;                # Help/Usage page
my $build;               # flag to say whether this is a build script or run outside of build
my $verbose;             # turn on extra output
my $debug;               # For sending output to just one person
my $database;            # which database to use?


GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose,
            "debug=s"    => \$debug,
            "help"       => \$help,
	    "build"      => \$build
            );


# Help pod if needed
exec ('perldoc',$0) if ($help);


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}





#####################################################################
#
# Open log file, output file, and database connection
#
#####################################################################

my $log = Log_files->make_build_log();  
$log->write_to(&runtime,": starting script\n===============\n\n");
$log->write_to("Looking for spurious genes shorter or equal to $cutoff amino acids in length\n");
print "Looking for spurious genes shorter or equal to $cutoff amino acids in length\n" if ($verbose);


# choice of output file location depends on if you are running this as part of build
if($build){
  open(OUT,">/wormsrv2/autoace/CHECKS/short_spurious_genes.$ws_version.csv") || die "Can't open output file\n";
}
else{
  open(OUT,">/wormsrv2/tmp/short_spurious_genes.$ws_version.csv") || die "Can't open output file\n";
}


# open a database connection and grab list of genes
my $db = Ace->connect(-path  =>  "$database")  || die "Cannot open $database",Ace->error;
my @genes = $db->fetch(-class => "Predicted_gene");




################################################
#
# M A I N   L O O P  -  loop through each gene
#
################################################

foreach my $gene (@genes){
  
  # get protein information, and make sure that there is a protein object!
  my ($protein) = $gene->at("Visible.Corresponding_protein");
  next if !($protein);  
  $protein = $db->fetch(Protein => "$protein") || die "Cannot fetch protein\n";
  my $length = $protein->at("Peptide")->right(2);

  # ignore proteins longer than cutoff value
  if ($length > $cutoff){
    $protein->DESTROY();
    next;
  }


  # Get confirmed/partially_confirmed etc. status
  my $status;
  if   ($gene->at("Properties.Coding.Confirmed_by")){$status = "Confirmed";}	
  elsif($gene->at("Visible.Matching_cDNA"))         {$status = "Partially_confirmed";}
  else                                              {$status = "Predicted";}			
  

  # get lab
  my $lab = $gene->From_laboratory;


  # get RNAi info, status defaults to 'N/A' when there are no RNAi results
  # otherwise status is set to WT or non-WT
  my $rnai_result = "N/A";
  my @RNAi = $gene->RNAi_result;
  foreach my $item (@RNAi){
    $rnai_result = "WT";
    my $RNAi = $db->fetch(RNAi => "$item");
    my $phenotype = $RNAi->Phenotype;
    if ($phenotype ne "WT"){
      $rnai_result = "non-WT";
      $RNAi->DESTROY();
      last;
    }
    $RNAi->DESTROY();
  }

  
  # look for associated PCR products that do or do not amplify
  # status is 'N/A' if there are no associated PCR products
  my $pcr_result = "N/A";
  my @PCRs = $gene->Corresponding_PCR_product;
  foreach my $item (@PCRs){
    $pcr_result = "Does not amplify";
    my $pcr = $db->fetch(PCR_product => "$item");
    my $amplify_status = $pcr->Amplified;
    if ((defined($amplify_status)) && ($amplify_status == 1)){
      $pcr_result = "Amplifies";
      $pcr->DESTROY();
      last;
    }
    $pcr->DESTROY();
  }


  # flag output status depending on how bad we think the prediction is...
  # 1 - no cDNA evidence AND two types of evidence AGAINST gene prediction (RNAi AND PCR)
  # 2 - no cDNA evidence AND one type of evidence AGAINST gene prediction (RNAi OR PCR) and
  #     other type of evidence not available
  # 3 - no cDNA evidence AND no available RNAi AND PCR information
  # 4 - no cDNA evidence AND contradictory RNAi AND PCR information (wild type AND amplifies, OR 
  #     non-wild type AND doesn't amplify)
  # 5 - no cDNA evidence BUT both RNAi and PCR information confirm it is a real gene
  # 6 - cDNA evidence


  my $evidence;

  if($status eq "Predicted"){
    if(($rnai_result eq "WT") && ($pcr_result eq "Does not amplify")){
      $evidence = "1";
    }
    elsif((($rnai_result eq "WT")  && ($pcr_result eq "N/A")) ||
	  (($rnai_result eq "N/A") && ($pcr_result eq "Does not amplify"))){
      $evidence = "2";
    }
    elsif(($rnai_result eq "N/A") && ($pcr_result eq "N/A")){
      $evidence = "3";
    }
    elsif((($rnai_result eq "WT")     && ($pcr_result eq "Amplifies")) || 
	  (($rnai_result eq "non-WT") && ($pcr_result eq "Does not amplify"))){
      $evidence = "4";
    }
    else{ 
      $evidence = "5";
    }
  }
  else{
    $evidence = "6";
  }
  
  print OUT "EVIDENCE $evidence,$lab,$gene,$length,$status,$rnai_result,$pcr_result\n";
  print     "EVIDENCE $evidence,$lab,$gene,$length,$status,$rnai_result,$pcr_result\n" if ($verbose);



  # kill AcePerl objects
  $gene->DESTROY();
  $protein->DESTROY();
}


################################
# Tidy up and exit             #
################################

close(OUT);
$db->close;
$log->write_to(&runtime,": finishing script\n");
# only mail if running as part of build
$log->mail("$maintainers") if ($build);
exit;





