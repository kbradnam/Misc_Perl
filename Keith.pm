package Keith;

# a module of frequently used bioinformatics subroutines
# Last updated by: $Author$
# Last updated on: $Date$

use strict;

sub base_composition{
	
	# recieve sequence, make sure it is all uppercase
	my $seq = uc(shift);

	# how many dp to calculate percentage?
	my $precision = shift;
	$precision = 2 if (!defined($precision));
	
	# hash to count frequencies
	my %mono = ('A' => 0,'C' => 0,'G' => 0,'T' => 0,'N' => 0);
	
	my $length = length($seq);

	# count mononucleotide frequencies
    $mono{"A"} += ($seq =~ tr/A/A/); 
    $mono{"C"} += ($seq =~ tr/C/C/); 
    $mono{"G"} += ($seq =~ tr/G/G/); 
    $mono{"T"} += ($seq =~ tr/T/T/); 
    $mono{"N"} += ($seq =~ tr/N/N/);

	
	my $a_percent  = sprintf("%3.${precision}f",($mono{'A'}/$length)*100);
	my $c_percent  = sprintf("%3.${precision}f",($mono{'C'}/$length)*100);
	my $g_percent  = sprintf("%3.${precision}f",($mono{'G'}/$length)*100);
	my $t_percent  = sprintf("%3.${precision}f",($mono{'T'}/$length)*100);
	my $n_percent  = sprintf("%3.${precision}f",($mono{'N'}/$length)*100);

	return($a_percent,$c_percent,$g_percent,$t_percent,$n_percent,);
}

sub gc_composition {
	my ($seq) = @_;
	my $gc = 0;
	my $tot = 0;
	for (my $i = 0; $i < length($seq); $i++) {
		my $nt = substr($seq, $i, 1);
		if    ($nt =~ /[GgCc]/) {$gc++; $tot++}
		elsif ($nt =~ /[AaTt]/) {$tot++}
		# ignore non-[ACGT]
	}
	return $gc/$tot;
}
sub gc_composition2 {
	my ($seq) = lc(shift);
	my $gc = 0;
	my ($a,$t,$c,$g) = (0,0,0,0);
	$c = $seq =~ tr/c/c/;
	$g = $seq =~ tr/g/g/;
	$a = $seq =~ tr/a/a/; 
	$t = $seq =~ tr/t/t/;
	return (($g+$c)/($c+$g+$a+$t));
}
    
1;


__END__

=head1 NAME

Keith;
