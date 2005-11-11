#!/usr/local/bin/perl5.6.0 -w

use Ace;


$|=1;

# open a local database connection
$db = Ace->connect(-path  =>  '/wormsrv2/pepace');

$prot = "WP:CE26491";

open(IN,"</wormsrv2/WORMPEP/wormpep64/fix3.ace") || die "Could not open file\n";

while(<IN>){
  chomp;
  if ($_ =~ /Protein/){
    $prot = $_;
    $prot =~ s/Protein : //;
    $prot =~ s/\"//g;
#    print "$prot\n";
    $sequence = $db->fetch(Protein => "$prot");
    @list = $sequence->Corresponding_DNA;
    $length = scalar(@list);
    if($length > 1){      
      print "Protein : \"$prot\"\n-D Inactive\n\n";

    }
  }

}
exit;





