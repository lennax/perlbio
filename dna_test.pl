#!/usr/bin/perl -w

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