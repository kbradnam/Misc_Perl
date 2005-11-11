#!/usr/local/bin/perl5.6.0 -w


while(<>){
#  print;
  /([a-z]+)\s+\[phenotype:\s+(.*)\]/;
  ($locus,$desc) = ($1,$2);

  print "\nGene_class : \"$locus\"\nCGC_approved\nPhenotype \"$desc\"\n\n";
  
}
  

  

