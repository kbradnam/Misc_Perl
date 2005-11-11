#!/usr/local/bin/perl -w

use strict;

die "The end\n" if (@ARGV < 2);

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

open(FILE1,"<$file1") || die "Couldn't open input file 1\n";



my $counter;

while(<FILE1>){
  my $line = $_;

  chomp($line);
 open(FILE2,"<$file2") || die "Couldn't open input file 2\n";
  while(<FILE2>){
    chomp;
    if (/$line/){
#      print "$line\t$_, Match\n";
      print "\n-R Peptide \"$line\" \"$_\"\n\n";
      print "Protein : \"$_\"\nPeptide $_\n\n";
      print "-D Protein : $line\n\n";
      $counter++;
    }    
  }
  close(FILE2);
}
close(FILE1);



print "Number of matches: $counter\n";
