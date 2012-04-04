package Tmstat; 

#-------------------------------------------------------------------------------
# Tmstat.pm
#-------------------------------------------------------------------------------
# 
# Exported functions: 
# tmstat: calculates mean and sd of hydrophobicity of all sequences of length n
# 
# created 28 February 2011
# modified 4 March 2011			Lenna Xiao Ping Peterson

require Exporter;
@ISA = "Exporter";
@EXPORT = "tmstat";
use strict;

# standard hydrophobicity scale (octanol-water)
my %hy = (	A =>  0.50,
			R =>  1.81,
			N =>  0.85,
			D =>  3.64,
			C => -0.02,
			Q =>  0.77,
			E =>  3.63,
			G =>  1.15,
			H =>  2.33,
			I => -1.12,
			L => -1.25,
			K =>  2.80,
			M => -0.67,
			F => -1.71,
			P =>  0.14,
			S =>  0.46,
			T =>  0.25,
			W => -2.09,
			Y => -0.71,
			V => -0.46,
);

# tmstat ----------------------------------------------------------------------
# Calculates mean and sd of hydrophobicity of all sequences of specified length
# in the given sequence (rough estimate of baseline hydrophobicity)
# @sequence is array of hashes with sequence key 'seq'
# USAGE: 
# 		my ( $mean, $sd ) = tmstat( @sequence, $tmlen );
sub tmstat {
	my ( @sequence ) = @_;
	my $tmlen = pop @sequence; 	# removes $tmlen from end of array
	my @all_scores;
	foreach my $seq_ref ( @sequence ) { 	# sequence is array of hashes
		my $seq = $$seq_ref{'seq'}; 
		my @seq_arr = split "", $seq; 
		foreach my $i ( 0 .. @seq_arr - $tmlen ) {
			my $end = $i + $tmlen - 1;
			my @slice = @seq_arr[$i .. $end]; 	# slice length $tmlen
			my $score;
			foreach my $aa ( @slice ) {
				$score += $hy{$aa}; 
			}
			$score /= $tmlen; 	# calculate average
			push @all_scores, $score; 
		}
	}
	my $total; 
	foreach my $scores ( @all_scores ) {
		$total += $scores; 
	}
	my $num_scores = @all_scores;
	my $mean = $total / $num_scores;
	my $square_diff;
	foreach my $scores ( @all_scores ) {
		$square_diff += ($scores - $mean)**2;
	}
	my $sd = sqrt( $square_diff / ( $num_scores - 1 ) );
	return ( $mean, $sd ); 
}
# end function tmstat ---------------------------------------------------------

1; # end package Tmstat.pm ----------------------------------------------------