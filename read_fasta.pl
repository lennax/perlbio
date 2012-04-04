#!/usr/bin/perl -w

# read_fasta.pl
# 
# Specify a FASTA file to output a list of sequence names and lengths
# and show full name and residues 50-99 of the following sequences: 
# 11667.m00003, 11667.m00009, 11667.m00015, 11667.m00031 
#
# Usage example:
# > perl read_fasta.pl data/read_fasta_in.txt
# 
# Revision 4 		works with MS-DOS and *nix newlines
# 					handles lack of seq terminator, streamlined
# 27 January 2011			Lenna Xiao Ping Peterson
# 

# process file line by line
# build hash of array of each sequence and its name
# each sequence begins with >
while  ( $line = <> ) {
	$line =~ s/\r\n$//; 	# remove CR-LF newline
	chomp $line; 		# remove *nix newline
	
	# first line: extract the sequence name for the key, 
	# separate part 1 name and details
	# start building the hash
	if ( $line =~ m/^>/ ) { 	# line starts with >
		if ( $seq_full ) { 	# add stored sequence to hash and clear var
			$sequences{$seq_key}[2] = $seq_full;
			$seq_full = "";		
		}
		
		( $seq_name, $seq_details ) = split " ", $line, 2; 	# pt 1 name
		( $seq_key ) = split/\|/, $seq_name, 2; 	# pt 2 name & hash key
		$seq_key =~ s/^>//; 		# removes leading >
		$sequences{$seq_key} = [ $seq_name, $seq_details ]; 	# build hash
	}
	else { $seq_full .= $line; } # concatenate the sequence
}
# adds last sequence to hash
if ( $seq_full ) { $sequences{$seq_key}[2] = $seq_full; }

### part 1 ###
# sort, hash == pseudorandom
# pull info to be printed from hash (sequence name and length)
# sequence name is [0] and full sequence is [2]
$number = 1;
foreach $seq1 ( sort { $a cmp $b } keys %sequences ) {
	$tmp_name = $sequences{$seq1}[0];
	$tmp_length = length($sequences{$seq1}[2]);
	print "$number\t$tmp_name\t$tmp_length\n";
	$number++;
}

### part 2 ###
# list of desired sequences, sort to be sure
@desired = sort { $a cmp $b } qw(11667.m00015 11667.m00003 11667.m00031 11667.m00009);

# pull stuff from %sequences
# desired sequences are conveniently equal to keys (almost like I planned it!)
# full name: sequence name [0], sequence details [1], note of range
# full sequence [2]; 50th base is index 49, 50-99 is 50 bases 
foreach $k ( 0 .. $#desired ) {
	$key = $desired[$k];
	$name = $sequences{$key}[0];
	$detail = $sequences{$key}[1]; 
	$seq2 = $sequences{$key}[2];
	$seq_slice = substr($seq2, 49, 50); 
	print "$name $detail [residues 50-99]\n";
	print "$seq_slice\n";
}
