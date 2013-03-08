#!/usr/bin/perl -w
# seq_func.pl
# a variety of functions for use in manipulating DNA and protein sequences
# Created: 16 February 2011
# Updated: 18 February 2011			Copyright Lenna Xiao Ping Peterson
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

use strict;

my %aa_list = (
	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
	TGT => "C", TGC => "C", TGA => "*", TGG => "W",
	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	CCT => "P", CCC => "P", CCA => "P", CCG => "P",
	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	CGT => "R", CGC => "R", CGA => "R", CGG => "R",
	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	AGT => "S", AGC => "S", AGA => "R", AGG => "R",
	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	GCT => "A", GCC => "A", GCA => "A", GCG => "A", 
	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	GGT => "G", GGC => "G", GGA => "G", GGG => "G"
);

### Amino acids
# A Ala Alanine
# R Arg Arginine
# N Asn Asparagine
# D Asp "Aspartic acid"
# C Cys Cysteine
# Q Gln Glutamine
# E Glu "Glutamic acid"
# G Gly Glycine
# H His Histidine
# I Ile Isoleucine
# L Leu Leucine
# K Lys Lysine
# M Met MEthionine
# F Phe Phenylalanine
# P Pro Proline
# S Ser Serine
# T Thr Threonine
# W Trp Tryptophan
# Y Tyr Tyrosine
# V Val Valine

### Extended genetic alphabet
# R		A G		puRine
# Y		C T		pYrimidine
# M		A C		aMino
# K		G T		Keto
# B		C G T	not-a (B)
# H		A C T	not-g (H)
# D		A G T	not-c (D)
# V		A C G	not-u (V)
# S		C G		Strong (3 H-bond)
# W		A T		Weak (2 H-bond)
# N 	ACGT	aNy

# extended complements 
# A T
# T A 
# G C 
# C G
# R	Y
# Y	R
# M	K
# K	M
# H	D
# D	H
# B	V
# V	B

# S	S
# W	W
# N N*


# countSeq --------------------------------------------------------------------
# returns hash of counts of each letter in sequence 
# USAGE: 
# 		%count = countSeq($DNA);
sub countSeq {
	my ( $seq ) = @_;
	$seq = uc $seq;
	my %letter_count;
	foreach my $letter ( "A" .. "Z" ) {
		my $count = $seq =~ s/$letter/$letter/g;
		$letter_count{$letter} = $count;
	}
	return %letter_count;
}
# end function countSeq--------------------------------------------------------

# getRange --------------------------------------------------------------------
# returns specified range of sequence
# USAGE: 
# 		$sub_seq = getRange ($seq, $start_pos, $end_pos );
sub getRange {
	my ( $seq, $begin, $end ) = @_;
	my $start = $begin - 1;
	my $range = $end - $begin + 1;
	my $subseq = substr($seq, $start, $range);
	return $subseq;
}
# end function getRange--------------------------------------------------------

# isDna -----------------------------------------------------------------------
# DEPENDENCY: countSeq
# guesses whether a sequence is DNA based on % of ATGC
# optional ATGC % threshold (default .9)
# USAGE: 
# 		if ( isDna($seq, [$threshold]) ) { ... }
sub isDna {
	my ( $seq, $thresh ) = @_;
	$thresh = .9 unless $thresh;
	$seq = uc $seq;
	$seq =~ tr/U/T/; # switch RNA to DNA
	my %count = countSeq($seq);
	my $atgc = 0;
	while ( my ($base, $count ) = each %count ) {
		my $dna_bases = "ATGC";
		$atgc += $count if ($base =~ /[$dna_bases]/);
	}
	my $total = length $seq;
	my $percent = $atgc/$total;
	if ( $percent >= $thresh ) { return 1; }
	else { return 0; }
}
# end function isDna-----------------------------------------------------------

# reverseComplement -----------------------------------------------------------
# finds reverse complement of a DNA sequence
# note for sequence searching: less memory to find revComp of pattern, not seq!
# USAGE:
# 		$reverse = reverseComplement( $dna );
sub reverseComplement {
	my ($seq) = @_;
	my $rev_seq = uc reverse $seq;
	# $rev_seq =~ tr/ATGC/TACG/; 	# find complement
	if ( $rev_seq =~ /U/ && $rev_seq !~ /T/ ) { 	# yes U and no T
		$rev_seq =~ tr/AU/UA/;
	} else { $rev_seq =~ tr/AT/TA/; } 
	$rev_seq =~ tr/GCRYMKHDBV/CGYRKMDHVB/; # including extended alphabet
	return $rev_seq;
}
# end function reverseComplement-----------------------------------------------

# translate -------------------------------------------------------------------
# DEPENDENCY: %aa_list defined above
# translates a DNA sequence into 1-letter amino acid 
# returns array of the 3 reading frames
# USAGE:
# 		@proteins = translate( $dna );
sub translate {
	my ($seq) = @_;
	my @proteins;
	my $length = length $seq;
	for my $frame ( 0 .. 2 ) {
		my $protein = "";
		for ( my $k = 0; $k < $length; $k += 3 ) {
			my $codon = substr $seq, $k + $frame, 3;
			if ( length $codon == 3 ) {
				$protein .= $aa_list{$codon};
			}
		}
		push @proteins, $protein;
	}
	return @proteins;
}
# end function translate-------------------------------------------------------

