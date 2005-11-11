#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Ace;
use Getopt::Long;




  #################################################################################################
  # checks if locus has positive_clone: 
  # if yes: checks to see if derived clone from CDS/Transcript already exists, otherwise adds to it
  # if no: adds derived clone from CDS/Transcript
  # This check is essential for cloned loci to be appeared as yellow highlighted in gmap
  #################################################################################################

  # list of loci linked to CDS/Transcript and have/have no positive clone
  my @loci_pos_clone_seq=`echo "table-maker -p $def_dir/loci_linked-not_linked_2_pos_clone.def" | $tace /wormsrv1/geneace`;

  my (%loci_has_pos_clone, %loci_no_pos_clone, @clone, @parent, %parent);
  
  my ($seq, $clone, $locus);
  my $ace=1;
  foreach (sort @loci_pos_clone_seq){
    chomp;
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      my $locus = $1;
      my $seq = $2;
      my $clone = $3;
      $clone =~ s/\..+//;
      push(@{$loci_has_pos_clone{$locus}}, $seq, $clone);
    }
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+$/){
      my $locus = $1;
      my $seq = $2;
      $loci_no_pos_clone{$locus} = $seq;
    } 
  }

  # also write ace file for added derived positive clone to each locus linked to CDS/Transcript
  
  foreach (keys %loci_no_pos_clone){
    $seq = $loci_no_pos_clone{$_};
    $seq =~ s/\..+//;
    $locus_errors++;
    print LOG "WARNING: $_ has no positive_clone, the derived one $seq can be added to it\n";
    if ($ace){
      print ACE "\n\nLocus\t\"$_\"\n";
      print ACE "Positive_clone\t\"$seq\" Inferred_automatically \"From Genomic_sequence or Transcript info\"\n";
    }							
  }
    
  foreach (keys %loci_has_pos_clone){
    @parent =(); %parent=();
    if (scalar @{$loci_has_pos_clone{$_}} == 2){
      $seq = ${@{$loci_has_pos_clone{$_}}}[0];
      $seq =~ s/\..+//;
      $clone = ${@{$loci_has_pos_clone{$_}}}[1];
      if ($seq ne $clone){
	if ($seq ne $clone){
	  push(@parent, $seq);
	}
      } 
      foreach my $e (@parent){$parent{$e}++};
      foreach my $ea (keys %parent){
	$locus_errors++;
        print LOG "WARNING: $_ already has positive_clone, the derived clone $ea can be added to it\n";
	if ($ace){
	  print ACE "\n\nLocus\t\"$_\"\n";
	  print ACE "Positive_clone\t\"$ea\"\n";
	}
      }	
    }
    if (scalar @{$loci_has_pos_clone{$_}} > 2){
      @clone =();
      for(my $i=0; $i < scalar @{$loci_has_pos_clone{$_}}; $i = $i+2){
	$seq = ${@{$loci_has_pos_clone{$_}}}[$i];
        $seq =~ s/\..+//;
        $clone = ${@{$loci_has_pos_clone{$_}}}[$i+1];
        push(@clone, $clone);
        if ($seq ne $clone){
	  push(@parent, $seq);
	}
      } 
      foreach my $p (@parent){$parent{$p}++};
      @parent = keys %parent; 
      my @comp_result=array_comp(\@parent, \@clone);
      my @same_clone=@{$comp_result[1]}; 
      if (!@same_clone){
	foreach my $ea (@parent){
	  $locus_errors++;
	  print LOG "WARNING: $_ already has positive_clone, the derived clone $ea can be added to it\n"; 
	  if ($ace){
	    print ACE "\n\nLocus\t\"$_\"\n";
	    print ACE "Positive_clone\t\"$ea\"\n";
	  }
	}
      } 
    }
  }
  print LOG "\nThere are $locus_errors errors in $size loci.\n";

sub array_comp{

  my(@union, @isect, @diff, %union, %isect, %count, $e);
  my ($ary1_ref, $ary2_ref)=@_;
  @union=@isect=@diff=();
  %union=%isect=();
  %count=();
  foreach $e(@$ary1_ref, @$ary2_ref){
    $count{$e}++;
  }
  foreach $e (keys %count){
    push (@union, $e);
    if ($count{$e}==2){push @isect, $e;}
    else {push @diff, $e;}
  } 
  return \@diff,
}
