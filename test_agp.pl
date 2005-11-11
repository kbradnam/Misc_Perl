#!/usr/local/bin/perl5.6.1 -w

use strict;

 ##############################
 # getz query to build hashes #
 ##############################

# ID   AF016448   standard; DNA; INV; 44427 BP.
# SV   AF016448.1
# AL031629

open (SEQ, "getz -f \'id seqversion\' \"([emblnew-org:Caenorhabditis elegans] \& [emblnew-div:RNA]) |");
while (<SEQ>) {
    print "$_";
}
close(SEQ);

