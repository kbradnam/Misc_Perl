#!/usr/bin/perl
# rows to columns script
# by Keith Bradnam 6/17/09

use warnings;
use strict;

# multidimensionl hash to store input
my @input;

# loop through input file (specify on command line)
while(<>){
	chomp;
	# add whole line (after splitting) to @input
	push(@input,[split(/,/)]);
	#push(@input,[split(/\t/)]);

}

my @output;

# first loop through row input array
for (my $y = 0; $y < @input; $y++) {

	# now lop through each column of each row
    for (my $x = 0; $x <  @{$input[$y]}; $x++) {

		# and reverse output into a new array
        $output[$x][$y] = $input[$y][$x];
    }
}

# loop through output array to print
for (my $line =0; $line < @output; $line++) {
	print "@{$output[$line]}\n";
}