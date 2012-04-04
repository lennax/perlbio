#!/usr/bin/perl -w

# translate_fasta.pl
# 
# Specify a FASTA file to produce all 6 translations of the sequences
# 
# main program begins line 92
# 
# Revision 4 		seq in hash, hash->hash->array, simplified loop length
# 11 February 2011			Lenna Xiao Ping Peterson
#

use strict;

# reverseComplement -----------------------------------------------------------
# finds reverse complement of a DNA sequence
# USAGE:
# 		$reverse = reverseComplement( $dna );
sub reverseComplement {
	my ($seq) = @_;
	my $rev_seq = reverse $seq;
	$rev_seq =~ tr/atgcATGC/TACGTACG/; 	# force UC, find complement
	return $rev_seq;
}

# translate -------------------------------------------------------------------
# translates a DNA sequence into 1-letter amino acid 
# returns array of the 3 reading frames
# USAGE:
# 		@proteins = translate( $dna );
sub translate {
	my ($seq) = @_;
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

# formatSeq -------------------------------------------------------------------
# formats a sequence to be printed with a specified block size and line length
# if line length is not an even multiple of block size, rounds down
# if line length is less than block size, treated as # of blocks
# USAGE:
# 		$output = formatSeq( $proteins[0], 10, 50 );
sub formatSeq {
	my ($seq, $block_len, $line_len) = @_;
	# force block_len to be logical
	if ($line_len < $block_len) { $line_len *= $block_len; } 
	elsif ($line_len % $block_len) { $line_len -= $line_len % $block_len; }
	my @sequence = split "", $seq; 	# break sequence into characters
	my $out;
	foreach my $char (0 .. $#sequence) {
		unless ($char % $block_len) {
			unless ($char % $line_len) {
				$out .= "\n" if $char; # \n if even block and line, not first
			} else {
				$out .= " "; 	# space if even block but not line
			}
		}
		$out .= $sequence[$char]; 	# print a letter every loop
	}
	return $out;
}

### grab sequence name and full sequence from input
my ($seq, $key, %sequences, %proteins);
while ( my $line = <> ) {
	$line =~ s/\r\n$//; 	# remove trailing CR-LF newline
	chomp $line; 
	if ( $line =~ s/^>// ) { 	# remove leading >
		$sequences{$key} = $seq, $seq = "" if $seq; 	# store and clear seq
		( my $seq_name ) = split " ", $line, 2; 	# part before space
		( $key ) = split/\|/, $seq_name, 2; 	# part before pipe
	}
	else { $seq .= $line; } 	# concatenate the sequence
}
$sequences{$key} = $seq if $seq; 	# store last sequence

### generate 3D hash of proteins, struct: {$key}{$direction}[$frame]
for my $k ( sort keys %sequences ) {
	$proteins{$k} = { forward => [ translate($sequences{$k}) ],
		 reverse => [ translate(reverseComplement $sequences{$k}) ] };
}

#### format and print sequences in FASTA format
my $block_len = 10;
my $line_len = $block_len * 5;
for my $key ( sort keys %proteins ) {  # loop through sequences
	for my $direction ( sort keys %{$proteins{$key}}) {  # forward and reverse
		my $dir = substr $direction, 0, 1;
		for my $frame ( 0 .. 2 ) { 	# 2 = $#{$proteins{$key}{$direction}}
			my $number = $frame + 1;
			print ">$key", "_$dir$number ", 
			"$key $direction reading frame $number\n",
			formatSeq($proteins{$key}{$direction}[$frame], 
			$block_len, $line_len), "\n";
		}	
	}	
}
### end of translate_fasta.pl ----------------------------------------------------