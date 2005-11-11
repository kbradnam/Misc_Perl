#!/usr/local/bin/perl5.6.1 -w
# 
# A simple script to find the best human, fly, and yeast matches to each wormpep protein
# and also extract the blast score
#
# by Keith Bradnam, March 2003


use Ace;
use strict;
use Getopt::Long;

my ($help, $database);

GetOptions ("help"        => \$help,
            "database=s"  => \$database);


#$|=1;

die "\nNeed to specify valid database path using --database option\n\n" if !($database);

# open a local database connection
# Need to have a database which has a Wormpep class/subclass
# use fetch_many rather than fetch as it is big class to grab

my $db = Ace->connect(-path  =>  "$database")|| do { 
  print "Connection failure: $database not a valid path",Ace->error; die();
};

my $wormpep = $db->fetch_many(-class => "Wormpep");

# print header line for output
print "#WormPep protein, best human match, best human score, best yeast match, best yeast score, best fly match, best fly score, best other match, best other score\n";


# Loop through each protein
while (my $protein = $wormpep->next){

  # flag variable for checking if there is any ENSEMBL match
  my $human_flag = 0;
  my $yeast_flag = 0;
  my $fly_flag = 0;
  my $worm_flag = 0;
  my $other_flag = 0;

  # Will get homologies from a column then row
  my @homol_column = ();
 
  # keep track of best score
  my $best_human_score = 0;
  my $best_yeast_score = 0;
  my $best_fly_score = 0;
  my $best_worm_score = 0;
  my $best_other_score = 0;
 
  # hash for matching score to protein name
  my %human_match = ();
  my %yeast_match = ();
  my %fly_match = ();
  my %worm_match = ();
  my %other_match = ();

  # grab all peptide hits
  @homol_column = $protein->at('Homol.Pep_homol');

  foreach my $entry (@homol_column){

    # grab homology details
    my @homol_row = (); 
    @homol_row = $protein->at("Homol.Pep_homol.$entry")->row;
    my $match = $homol_row[0];

    # Does it have a human hit?
    if ($entry =~ m/^ENSEMBL/){
      $human_flag = 1;
      my $score = $homol_row[2];
      
      if ($score > $best_human_score){
	($best_human_score = $score);
	$human_match{$score} = "$match";
      }
      elsif($score == $best_human_score){
	# need to cope with multiple matching proteins with same score
	$human_match{$score} .= " $match";
      }  
    }

    # Does it have a yeast hit
    if ($entry =~ m/^SGD/){
      $yeast_flag = 1;
      my $score = $homol_row[2];
      
      if ($score > $best_yeast_score){
	($best_yeast_score = $score);
	$yeast_match{$score} = "$match";
      }
      elsif($score == $best_yeast_score){
	# need to cope with multiple matching proteins with same score
	$yeast_match{$score} .= " $match";
      }  
    }
    # Does it have a fly hit
    if ($entry =~ m/^GADFLY/){
      $fly_flag = 1;
      my $score = $homol_row[2];
      
      if ($score > $best_fly_score){
	($best_fly_score = $score);
	$fly_match{$score} = "$match";
      }
      elsif($score == $best_fly_score){
	# need to cope with multiple matching proteins with same score
	$fly_match{$score} .= " $match";
      }  
    }

    # Does it have a worm hit
    if ($entry =~ m/^WP:/ && defined($homol_row[1])){
      $worm_flag = 1;
      my $score = $homol_row[2];
	
      if ($score > $best_worm_score){
	($best_worm_score = $score);
	$worm_match{$score} = "$match";
      }
      elsif($score == $best_worm_score){
	# need to cope with multiple matching proteins with same score
	$worm_match{$score} .= " $match";
      }   
    }


    # Does it have an other hit (SwissProt/Trembl)?
    if ($entry =~ m/^SW:/ || $entry =~ m/^TR:/){
      $other_flag = 1;
      my $score = $homol_row[2];
	
      if ($score > $best_other_score){
	($best_other_score = $score);
	$other_match{$score} = "$match";
      }
      elsif($score == $best_other_score){
	# need to cope with multiple matching proteins with same score
	$other_match{$score} .= " $match";
      }   
    }

  }
  # Print Output
  print "$protein,";

  # print results if there were human matches  
  if($human_flag){
    print "$human_match{$best_human_score},$best_human_score,";
  }
  else{
    print ",,";
  }

  # print results if there were yeast matches
  if($yeast_flag){
    print "$yeast_match{$best_yeast_score},$best_yeast_score,";
  }
  else{
    print ",,";
  }

  # print results if there were fly matches
  if($fly_flag){
    print "$fly_match{$best_fly_score},$best_fly_score,";
  }
  else{
    print ",,";
  }

  # print results if there were worm matches
  if($worm_flag){
    print "$worm_match{$best_worm_score},$best_worm_score,";
  }
  else{
    print ",,";
  }

  # print results if there were worm matches
  if($other_flag){
    print "$other_match{$best_other_score},$best_other_score,";
  }
  else{
    print ",,";
  }

  print "\n";
}

print "\n";
$db->close;
exit(0);





