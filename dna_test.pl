#!/usr/bin/perl -w

# Copyright 2011 Lenna X Peterson
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

require 'seq_func.pl';

my $dwf5_rev2 = "NNNNNNNNTNNNTCCGGATGANCCTGTATTTNACTTTCTCACAATAGGGCTTCCAATATTTCCCANAN";

my $thresh = .77;
if ( isDna($dwf5_rev2, $thresh) ) {
	print "dwf5_rev2 appears to be DNA\n";
} else {
	my $percent = $thresh * 100;
	print "dwf5_rev2 is less than $percent% ATGC\n";
}

my $test_rna = "aucauguaucuugua";

if ( isDna ($test_rna) ) {
	# print "seq appears to be DNA\n";
	print reverseComplement($test_rna), "\n";
} else {
	print "seq is less than 90% ATGC\n";
}
