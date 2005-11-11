#!/usr/local/bin/perl5.6.0 -w


while(<>){
  chomp;
  
  @list = split(/,/,$_);
  next if ($list[2] eq "");

  $list[0] =~ s/\"//g;
  $list[1] =~ s/\"//g;
  $list[2] =~ s/\"//g;
  $list[3] =~ s/\"//g;


  $seq_map = $list[1];
  $seq_map =~ s/Sequence-//;
  $chr = $list[3];

  if($list[3]  =~ m/CHROMOSOME/){
    $chr  =~ s/CHROMOSOME_//;
    if(($seq_map ne $list[2]) || ($seq_map ne $chr) || ($list[2] ne $chr)){
      print "$list[0]\t$list[1]\tgMap-$list[2]\t$list[3]\n";
    }
  }

  elsif($list[3]  =~ m/SUPERLINK/){
    $chr  =~ s/SUPERLINK_CB_//;
    chop($chr) if ($chr =~ m/L$/);
    chop($chr) if ($chr =~ m/R$/);

    if(($seq_map ne $list[2]) || ($seq_map ne $chr) || ($list[2] ne $chr)){
      print "$list[0]\t$list[1]\tgMap-$list[2]\t$list[3]\n";
    }
  }
  else{
    print "$list[0]\t$list[1]\tgMap-$list[2]\t$list[3]\n";
  }
}
  

  

