## Perl bioinformatics exercises ##

Written in 2011 by Lenna X. Peterson

`dna_test.pl` attempts to guess whether a sequence is DNA or RNA. It depends on `seq_func.pl`.

The following files have sample input and output located in `data/`:

`read_fasta.pl` reads a multi-record FASTA file and outputs a list of sequence names and lengths.

`translate_fasta.pl` reads a FASTA file and outputs the 6 possible translations (3 reading frames, forward and reverse complement).

`pdb_distance.pl` reads a PDB file of a molecule complexed with ADP and outputs a list of every atom within 6 angstroms of the ADP.

`test_hydrophobic.pl` calls a package that finds the most hydrophobic sequences in a protein. It depends on `Hydrophobic.pm` and `Tmstat.pm`.

*Another project around the same time used CGI to perform batch BLAST searches of a given sequence against multiple genomes. However, it was a group project so I lack authorization to upload it to github.*
