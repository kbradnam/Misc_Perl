package Keith;

# a module of frequently used bioinformatics subroutines
# Last updated by: $Author$
# Last updated on: $Date$

use strict;

sub base_composition{
	my $seq = uc(shift);
	my $length = length($seq);
	my ($a,$t,$c,$g,$n) = 0;
	return($a,$c,$g,$t,$n);
}
    
1;


__END__

=head1 NAME

Keith;