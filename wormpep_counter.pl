#!/usr/local/bin/perl -w


for ($i=0; $i<29945;$i++){
  $protein = sprintf "%05d" , $i;
  $pepace_protein = "WP:CE$protein";
  print "Protein : \"$pepace_protein\"\nWormpep\n\n";

}



