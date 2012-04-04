#!/usr/bin/perl -w

# pdb_distance.pl
# 
# Specify a .PDB file of a molecule complexed with ADP to list all every atom
# of any protein residue of which any atom is within 6 angstroms of the ADP
# 
# Revision 3		store target atoms in hash by (chain and residue number)
# 17 February 2011			Lenna Xiao Ping Peterson

use strict; 
# use Time::HiRes ("gettimeofday", "tv_interval"); 
# my $start = [gettimeofday]; 

my ( @adp, %residues, %close ); 
### loop through input and extract needed fields 
# store ADP atoms in an array, store protein atoms in a hash by residue
# need: chain, residue type, residue number, atomtype, xyz position
# input structure (human numbers): 
#  1    2     *3*  *4*      *5*   *6*    *7*8*9*  10  11
# type number atom res-type chain res-num x y z  occ temp 
print "Target atoms\n"; 
while ( my $line = <> ) {
	$line =~ s/\r\n$//; 
	chomp $line; 
	my ( $type, $nottype ) = split /\s+/, $line, 2; 
	if ( $type eq "HETATM" || $type eq "ATOM" ) {
		my ( $id, $atom_type, $residue_type, $tmp ) = split /\s+/, $nottype, 4; 
		# fields 5 and 6 may not have a space, add one if second char is digit
		$tmp =~ s/(^\w)(\d)/$1 $2/; # match: first two char are word & digit
		my ( $chain, $residue_n, $x, $y, $z ) = split /\s+/, $tmp, 6; 
		my $res_id = $chain."_".$residue_n;  # create unique ID of residue
		my @store = ( $id, $residue_type, $atom_type, $x, $y, $z ); 
		
		if ( $residue_type eq "ADP") { 	# store ADP atom coords in an array
			push @adp, [ @store[-3..-1] ]; 	# coords are last 3 elements
			print "$chain\t$residue_type\t$residue_n\t"; 
			print join("\t", @store[-4..-1]), "\n"; 
		} else { 	# store non-ADP atoms in a HoAoA
			push @{$residues{$res_id}}, [ @store ]; 
		}
	} 
}

### check whether each protein atom is within 6 angstroms of any ADP atom
# store residue information in a hash (unique)
# structure: %residues{$res_id}[$i] = @array:
#  0   1             2          3  4  5
# $id $residue_type $atom_type $x $y $z
foreach my $res_id ( keys %residues ) {
	foreach my $atom ( @{$residues{$res_id}}) {
		unless ( $close{$res_id} ) {	# checks if residue is already flagged
			foreach my $adp_atom ( @adp ) {
				my ( $ax, $ay, $az ) = @$adp_atom;
				my ( $x, $y, $z ) = @$atom[-3..-1]; 	# coords are last 3
				my $distance = sqrt( ($x-$ax)**2 + ($y-$ay)**2 + ($z-$az)**2 ); 
				if ( $distance <= 6 ) {
					$close{$res_id}++; 
				}
			}
		}
	}
}

### print all atoms of any residue within 6 angstroms of ADP
# look up close residues in atom hash, sort by atom ID, print
print "Neighbor residues\n";
# byId: sort by atom ID of first atom in residue
sub byId { $residues{$a}[0][0] <=> $residues{$b}[0][0] }
foreach my $res_id ( sort byId keys %close) {
	foreach my $atom ( @{$residues{$res_id}} ) {
		my ( $chain, $residue_n ) = split "_", $res_id; 
		my $residue_type = $atom->[1]; 
		print "$chain\t$residue_type\t$residue_n\t"; 
		print join("\t", @$atom[-4..-1]), "\n"; 
	}
}

# my $clock = tv_interval($start); 
# warn "time: $clock\n"; 
# end pdb_distance.pl ---------------------------------------------------------