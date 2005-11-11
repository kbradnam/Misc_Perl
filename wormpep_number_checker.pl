#!/usr/local/bin/perl5.6.1 -w


while(<>){
  chomp;
  next unless $_ =~ m/^>/;

  $protein = $_;
  $protein =~ s/>CE0*//;
  push(@wormpep,$protein);


}


for ($counter=1;$counter<=33551;$counter++){
  $match = 0;
  for ($i=0; $i< @wormpep; $i++){
    if ($counter == $wormpep[$i]){
      $match =1;
#      splice(@wormpep,$i,1);
    }
  }
  print "NOMATCH: $counter\n" if ($match ==  0);
}  

  

