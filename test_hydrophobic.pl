#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# test_hydrophobic.pl
#
# test the Hydrophobic package by doing a trial calculation of the most
# hydrophobic regions of all proteins in a FASTA file
#
# Copyright 2011 Lenna X. Peterson
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
#-------------------------------------------------------------------------------

use strict;
use Hydrophobic;
use Tmstat;

my $tmlen = 20;     # length of hydrophobic region
my $nbest = 10;     # number of best regions to report

my @sequence = getFasta();      # @sequence is an array of hashes
my ( $mean, $sd ) = tmstat( @sequence, $tmlen );
warn "$mean\n$sd\n";

# calculate most hydrophobic regions for each sequence

foreach my $seq ( @sequence ) {
    print "$$seq{name}: $$seq{doc}\n";
    
    print "Most hydrophobic regions ($tmlen residues)\n";
    my @tmseq = tmBest( $seq, $tmlen, $nbest );
    foreach my $tm ( @tmseq ) {
        my $zscore = ( $$tm{score} - $mean ) / $sd;
        printf "%20s  %4d  %4d  %8.2f  %6.3f\n",
               $$tm{seq}, $$tm{begin}, $$tm{end}, $$tm{score}, $zscore;
    }
    print "\n";
}
