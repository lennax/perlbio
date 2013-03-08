package Hydrophobic;

#-------------------------------------------------------------------------------
# Hydrophobic.pm
#-------------------------------------------------------------------------------
# 
# Exported functions:
# getFasta: reads FASTA from stdin, returns array of hashes of seq data
# tmBest: finds n most hydrophobic sequences in given protein
#
# created 26 February 2011
# modified 4 March 2011			Copyright Lenna Xiao Ping Peterson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


require Exporter;
@ISA = "Exporter";
@EXPORT = qw ( getFasta tmBest );
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

# getFasta --------------------------------------------------------------------
# reads FASTA file from STDIN, separates FASTA sequences
# returns array of hashes containing 'name', 'doc'umentation, and 'seq'uence
# USAGE: 
# 		my @sequence = getFasta();
sub getFasta {
	my ( @sequences, $seq ); 
	my $seq_count = 0; 
	while ( my $line = <> ) {
		$line =~ s/\r\n$//; 	# remove trailing CR-LF newline
		chomp $line; 
		if ( $line =~ s/^>// ) { 	# remove leading '>'
			if ( $seq && $seq_count > 0 ) {
				$sequences[$seq_count-1]{'seq'} = $seq; 	# in prev element
				$seq = "";
			}
			my ( $seq_name, $seq_doc ) = split " ", $line, 2; 
			$sequences[$seq_count]{'name'} = $seq_name;
			$sequences[$seq_count]{'doc'} = $seq_doc;
			$seq_count++;
		}
		else { $seq .= $line; } 	# concatenate sequence
	}
	$sequences[$seq_count-1]{'seq'} = $seq if $seq; 
	return @sequences;
}
# end function getFasta -------------------------------------------------------


# tmBest ----------------------------------------------------------------------
# finds n most hydrophobic regions of given length in given protein
# returns array of hashes: 'seq'uence, 'begin', 'end', hydrophobic 'score'
# $seq is hash reference with sequence key 'seq'
# USAGE:
# 		my @tmseq = tmBest ( $seq, $tmlen, $nbest );
{ 	# naked block to share array for sorting
	my @tm_seq; 
	sub tmBest {
		my ( $seq_hash_ref, $tmlen, $nbest ) = @_;
		my $seq = $$seq_hash_ref{'seq'};
		@tm_seq = (); 
		foreach my $i ( 0 .. length ($seq) - $tmlen ) {
			my $score;
			my $testregion = substr($seq, $i, $tmlen); 
			my @region = split "", $testregion;
			foreach my $aa ( @region ) {
				$score += $hy{$aa}; 	# sum hydrophobicity
			}
			my $fudge = 100000; 	# fudge factor to avoid float pt errors
			# calculate average, truncate, correct
			$score = (int($score * $fudge / $tmlen)) / $fudge; 	
			my %stored = (seq => $testregion, begin => $i+1, end => $i+$tmlen, 
						score => $score); 
			if ( $i < $nbest ) {
				push @tm_seq, \%stored; 	# store hash reference 
			} else {
				@tm_seq = sort byScore @tm_seq; 
				# check if current region has lower score than highest in array
				if ( $score < $tm_seq[-1]{score} ) {
					( $tm_seq[-1] ) = \%stored; 	# store hash reference
				}
			}
		}
		@tm_seq = sort byScore @tm_seq;
		return @tm_seq;
	}
	# sort function for tmBest - sort by 'score' then by 'begin'
	sub byScore { $$a{score} <=> $$b{score} || $$a{begin} <=> $$b{begin} }
}
# end function tmBest ---------------------------------------------------------

1; # end package Hydrophobic.pm -----------------------------------------------
