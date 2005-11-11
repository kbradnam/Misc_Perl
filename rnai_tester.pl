#!/usr/local/bin/perl5.8.0 -w
#
# find_short_genes.pl  
# 
# by Keith Bradnam, aged 12 and half
#
# A script to find (and classify) potentially short, spurious genes (default = <50 aa)
#
# Last updated by: $Author$     
# Last updated on: $Date$     


use strict;
use Ace;
use Getopt::Long;


##############################
# misc variables             #
##############################

my $file         = "rnai_data.csv"; # initial output file


my $verbose;             # turn on extra output
my $database;            # which database to use?


GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose
            );



# open initial output file
open(OUT,">$file") || die "Can't open output file\n";


# open a database connection and grab list of genes
my $db = Ace->connect(-path  =>  "$database")  || die "Cannot open $database",Ace->error;
my @genes = $db->fetch(-class => "elegans_CDS");


# used for testing purposes
my $counter = 0;


################################################
#
# M A I N   L O O P  -  loop through each gene
#
################################################

foreach my $gene (@genes){

  # when debugging can reduce this if you want less output to speed up script
  last if ($counter > 50000);  $counter++;


  # get protein information, and make sure that there is a protein object!
  my ($protein) = $gene->at("Visible.Corresponding_protein");
  next if !($protein);  
  $protein = $db->fetch(Protein => "$protein") || die "Cannot fetch protein\n";
  my $length = $protein->at("Peptide")->right(2);

  my $dna_length = $length * 3;
  my $dna = $gene->asDNA;
  $dna =~ s/\n//g;
  my $length_gene_name = length($gene)+1;
  $dna = substr($dna, $length_gene_name);

  my $C = $dna =~ tr/c/c/;
  my $G = $dna =~ tr/g/g/;

  my $p_GC = sprintf "%.2f",((($C+$G)/$dna_length)*100);

  # Get confirmed/partially_confirmed etc. status
  my $status;
  if   ($gene->at("Properties.Coding.Confirmed_by")){$status = "Confirmed";}	
  elsif($gene->at("Visible.Matching_cDNA"))         {$status = "Partially_confirmed";}
  else                                              {$status = "Predicted";}			
  

  # get RNAi info, status defaults to 'N/A' when there are no RNAi results
  # otherwise status is set to WT or non-WT
  my $rnai_result = "N/A";
  my @RNAi = $gene->RNAi_result;
  my $rnai_counter = 0;
  ($rnai_counter) = scalar(@RNAi) if ($gene->RNAi_result);

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

  
  print "$gene,$dna_length,$p_GC,$protein,$length,$status,$rnai_result,$rnai_counter,$pcr_result\n";



  # kill AcePerl objects
  $gene->DESTROY();
  $protein->DESTROY();
}

close(OUT);
$db->close;


exit;





__END__

=pod

=head2   NAME - find_short_genes.pl

=head1 USAGE

=over 4

=item find_short_genes.pl -[options]

=back

=head1 DESCRIPTION

This script will check the 'elegans_CDS' subclass in any target database
and find all genes less than a certain size (defaults to 50 amino acids).
It will then check each of those genes for the presence/absence of cDNA
information and also whether there are associated RNAi experiments that give
a wild-type phenotype.  Finally, it also looks to see if there are associated
PCR products that don't amplify.  The various combinations of these different
types of evidence allows the script to categorise and prioritise the genes,
such that curators can inspect these short genes to see if they are real.

Can be run as part of build or not.  If so it will write the output file to
/wormsrv2/autoace/CHECKS/ else it will write to /wormsrv2/tmp

=back

=head1 MANDATORY arguments: -database

=over 4

=item -database <path to valid acedb database>

Database will be expected to have a valid 'elegans_CDS' subclass

=back

=head1 OPTIONAL arguments: -cutoff, -build, -debug, -verbose, -help, -html

=over 4

=item -cutoff <integer>

Specify the length (in amino acids) of genes to consider by this script.  Will
look at every gene equal to or less than the specified length


=item -build

If specified will assume that this script is being run as part of the build,
the only differences being that the output file will be placed in /wormsrv2/autoace/CHECKS
rather than in /wormsrv2/tmp and the log file will be emailed to everyone


=item -html

Additionally produce a html file that can be added to the website as part of the 
data consistency checks page

=item -debug <user>

Send log report to specified user only

=item -help

This help.

=back

=item -verbose

Toggles a little extra output on the command line when running the script

=back


=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk) 

=back

=cut

