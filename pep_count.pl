#!/usr/bin/perl -w

use strict;

open(IN,$ARGV[0]) || die "Can't open file\n";

my ($a,$c,$d,$e,$f,$g,$h,$i,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$v,$w,$y) = 0;

my $seq;

while(my $temp=<IN>){
    if($temp !~ />/){
        chomp($temp);
        $temp =~  tr/a-z/A-Z/;
        $seq .= $temp;
    }
}

$a = $seq =~ tr/A//; 
$c = $seq =~ tr/C//;
$d = $seq =~ tr/D//;
$e = $seq =~ tr/E//;
$f = $seq =~ tr/F//; 
$g = $seq =~ tr/G//;
$h = $seq =~ tr/H//;
$i = $seq =~ tr/I//;
$k = $seq =~ tr/K//;
$l = $seq =~ tr/L//; 
$m = $seq =~ tr/M//;
$n = $seq =~ tr/N//;
$p = $seq =~ tr/P//;
$q = $seq =~ tr/Q//;
$r = $seq =~ tr/R//; 
$s = $seq =~ tr/S//;
$t = $seq =~ tr/T//;
$v = $seq =~ tr/V//;
$w = $seq =~ tr/W//;
$y = $seq =~ tr/Y//; 



# now calculate totals and percentages
my $total = $a+$c+$d+$e+$f+$g+$h+$i+$k+$l+$m+$n+$p+$q+$r+$s+$t+$v+$w+$y;
my ($pa,$pc,$pd,$pe,$pf,$pg,$ph,$pi,$pk,$pl,$pm,$pn,$po,$pp,$pq,$pr,$ps,$pt,$pv,$pw,$py) = 0;

$pa = sprintf("%.4f",($a/$total));
$pc = sprintf("%.4f",($c/$total));
$pd = sprintf("%.4f",($d/$total));
$pe = sprintf("%.4f",($e/$total));
$pf = sprintf("%.4f",($f/$total));
$pg = sprintf("%.4f",($g/$total));
$ph = sprintf("%.4f",($h/$total));
$pi = sprintf("%.4f",($i/$total));
$pk = sprintf("%.4f",($k/$total));
$pl = sprintf("%.4f",($l/$total));
$pm = sprintf("%.4f",($m/$total));
$pn = sprintf("%.4f",($n/$total));
$pp = sprintf("%.4f",($p/$total));
$pq = sprintf("%.4f",($q/$total));
$pr = sprintf("%.4f",($r/$total));
$ps = sprintf("%.4f",($s/$total));
$pt = sprintf("%.4f",($t/$total));
$pv = sprintf("%.4f",($v/$total));
$pw = sprintf("%.4f",($w/$total));
$py = sprintf("%.4f",($y/$total));

print "A,$a,$pa\n";
print "C,$c,$pc\n";
print "D,$d,$pd\n";
print "E,$e,$pe\n";
print "F,$f,$pf\n";
print "G,$g,$pg\n";
print "H,$h,$ph\n";
print "I,$i,$pi\n";
print "K,$k,$pk\n";
print "L,$l,$pl\n";
print "M,$m,$pm\n";
print "N,$n,$pn\n";
print "P,$p,$pp\n";
print "Q,$q,$pq\n";
print "R,$r,$pr\n";
print "S,$s,$ps\n";
print "T,$t,$pt\n";
print "V,$v,$pv\n";
print "W,$w,$pw\n";
print "Y,$y,$py\n";
print "total,$total\n";

close(IN);
