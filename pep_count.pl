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

$a = $seq =~ tr/A/A/; 
$c = $seq =~ tr/C/C/;
$d = $seq =~ tr/D/D/;
$e = $seq =~ tr/E/E/;
$f = $seq =~ tr/F/F/; 
$g = $seq =~ tr/G/G/;
$h = $seq =~ tr/H/H/;
$i = $seq =~ tr/I/I/;
$k = $seq =~ tr/K/K/;
$l = $seq =~ tr/L/L/; 
$m = $seq =~ tr/M/M/;
$n = $seq =~ tr/N/N/;
$p = $seq =~ tr/P/P/;
$q = $seq =~ tr/Q/Q/;
$r = $seq =~ tr/R/R/; 
$s = $seq =~ tr/S/S/;
$t = $seq =~ tr/T/T/;
$v = $seq =~ tr/V/V/;
$w = $seq =~ tr/W/W/;
$y = $seq =~ tr/Y/Y/; 

print "a,$a\n";
print "c,$c\n";
print "d,$d\n";
print "e,$e\n";
print "f,$f\n";
print "g,$g\n";
print "h,$h\n";
print "i,$i\n";
print "k,$k\n";
print "l,$l\n";
print "m,$m\n";
print "n,$n\n";
print "p,$p\n";
print "q,$q\n";
print "r,$r\n";
print "s,$s\n";
print "t,$t\n";
print "v,$v\n";
print "w,$w\n";
print "y,$y\n";

close(IN);
