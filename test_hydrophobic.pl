#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# test_hydrophobic.pl
#
# test the Hydrophobic package by doing a trial calculation of the most
# hydrophobic regions of all proteins in a FASTA file
#
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
