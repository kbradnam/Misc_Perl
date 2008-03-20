#!/usr/bin/perl
#
# Usage: rename.pl perlexpr [files]

($regexp = shift @ARGV) || die "Usage:  rename.pl perlexpr [filenames]\n";

if (!@ARGV) {
   @ARGV = <STDIN>;
   chomp(@ARGV);
}


foreach $_ (@ARGV) {
   $old_name = $_;
   eval $regexp;
   die $@ if $@;
   rename($old_name, $_) unless $old_name eq $_;
}

exit(0);