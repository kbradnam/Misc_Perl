#!/usr/local/bin/perl5.6.0 -w

my @bits;
while(<>){

  if ($_ =~ m/^Subsequence/){
    @bits = split(/\s+/,$_);
    if (($bits[2] > 1283934) || ($bits[3] > 1283934)){
      my $new_coord1 = $bits[2] -5;
      my $new_coord2 = $bits[3] -5;
      print "$bits[0] $bits[1] $new_coord1 $new_coord2\n";
    }
    else{
      print;
    }
      

  }
  else{
    print;

  }

}
  

  

