#!/usr/local/bin/perl5.6.0 -w

$/ = "";

while(<>){
  $name = $sequence = $locus = $text = $longtext = "";
  chomp;
  if($_ =~ m/^Sequence/){
    $sequence = $_;
    $sequence =~ s/Sequence :\s+\"(.*)\"//;
    $name = $1;
  }
  elsif($_ =~ m/^Locus/){
    $locus = $_;
    $locus =~ s/Locus :\s+\"(.*)\"//;
    $name = $1;
  }
  else{
    print "NO Locus or Sequence name!\n\n***$_***\n\n\n";
  }
#  print "\nName is $name\n--------\n";

  if($_ =~ m/Detailed_description\s+\"(.*?)\"\n(PMID_evidence|Paper_evidence)/s){
    $text = $1;
    $longtext = "\nLongText : \"$name\"\n".$text."\n\*\*\*LongTextEnd\*\*\*\n\n";
    $_ =~ s/Detailed_description.*?(PMID_evidence|Paper_evidence)/$1/sg;

    if($_ =~ m/PMID_evidence\s+\"\d+\"/sg){
      $_ =~ s/PMID_evidence\s+\"(\d+)\"/Detailed_description \"$name\" PMID_evidence \"$1\"/sg;    
    }

    if($_ =~ m/Paper_evidence\s+\"\[.*?\]\"/sg){
      $_ =~ s/Paper_evidence\s+\"(.*?)\"/Detailed_description \"$name\" Paper_evidence \"$1\"/sg;    
    }
    print "$_\n";
   print "$longtext\n";

  }
  elsif($_ =~ m/Provisional_description\s+\"(.*?)\"/s){
#	print "Name is $name\n---------\n$_\n\n\n";
    $text = $1;
    $longtext = "\nLongText : \"$name\"\n".$text."\n\*\*\*LongTextEnd\*\*\*\n\n";

    $text =~ s/\n/ /g;
    $_ =~ s/Provisional_description.*?(PMID_evidence|Paper_evidence)/$1/sg;

    if($_ =~ m/PMID_evidence\s+\"\d+\"/sg){
      $_ =~ s/PMID_evidence\s+\"(\d+)\"/Provisional_description \"$name\" PMID_evidence \"$1\"/sg;    
    }

    if($_ =~ m/Paper_evidence\s+\"\[.*?\]\"/sg){
      $_ =~ s/Paper_evidence\s+\"(.*?)\"/Provisional_description \"$name\" Paper_evidence \"$1\"/sg;    
    }
    print "\n$_\n";
    print "$longtext\n";
  }
}
  

  

