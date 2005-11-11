#!/usr/local/bin/perl5.6.0 -w


$release = $ARGV[0];

#&part1;
#&part1;
&part3;


sub part1{
  open(IN,"</wormsrv2/WORMPEP/wormpep$release/wormpep.history$release") || die "Couldn't open file\n";

  my ($cds,$id,$in,$out) = "";
  while(<IN>){
    chomp;
    
    my @list = split(/\t/,$_);
    if(defined($list[3])){  
      my $temp = <IN>;
    my @line2 = split(/\t/,$temp);
      if ($list[0] eq $line2[0]){
	#      print "It's a match: \n$_\n$temp\n";
	print "Protein : \"WP:$list[1]\"\n-D Corresponding_DNA\nReplaced_by WP:$line2[1]\n\n";
      }
      
    }
    
  }
  close(IN);
}


sub part2{
  open(SMALL,"</wormsrv2/WORMPEP/wormpep$release/small") || die "Couldn't open file\n";
  
  
  while($small = <SMALL>){
    $match = 0;
    chomp($small);
    #  print "$small:\n";
    open(LARGE,"</wormsrv2/WORMPEP/wormpep$release/large") || die "Couldn't open file\n";
    while($large = <LARGE>){
      chomp($large);
      if ($small eq $large){
	#      print "$small matches $large\n";
	$match =1;
      }
    }
    print "No match for $small" if ($match == 0);
    close(LARGE);
  }

  
  close(SMALL);

}

sub part3{  
  open(IN,"</wormsrv2/WORMPEP/wormpep$release/fix4.ace") || die "Couldn't open file\n";

#  open(OUT,">/wormsrv2/WORMPEP/wormpep$release/fix5.ace") || die "Couldn't open file\n";
  
  $/="";
  while(<IN>){
#    print "*$_*\n";
    
    $protein = $1 if ($_ =~ /(\".*\")/);
    $dna = $1 if ($_ =~ m/Corresponding_DNA (\".*\")/);

    $next = <IN>;
    $protein2 = $1 if ($next =~ /(\".*\")/);
    $dna2 = $1 if ($next =~ m/Corresponding_DNA (\".*\")/);

    if ($protein eq $protein2){
#      print "\nMATCH!\nProtein : \"$protein\"\nCorresponding_DNA $dna\nCorresponding_DNA $dna2\n\n"; 
    }
    else{
      print "$_";
      print "$next";
      ;
    }
   } 
  close(IN);
#  close(OUT);

  }
