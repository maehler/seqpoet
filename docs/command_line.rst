Command line arguments
======================

::

    seqpoet [options] genomedir probe

Mandatory arguments
-------------------

=============   =======================================================
genomedir       directory containing the genome files to use (FASTA or
                GenBank format) or a single GenBank or FASTA file
probe           file containing either a single sequence (probe) or a
                pair of sequences (primer pair; one sequence per line)
=============   =======================================================

Optional arguments
------------------

-h, --help            show this help message and exit
--pcr                 only perform in silico PCR. Requires that the probe
                      file contains a primer pair (default: perform operon
                      extraction)
-m int, --mismatches int
                      the maximum number of mismatches allowed when aligning
                      probe/primer to the genome (default: 2)
-d int, --max-distance int
                      the maximum intergenic distance allowed when
                      assembling operons (default: 500)
--min-product int     minimum PCR product length to consider (default: 0)
--max-product int     maximum PCR product length to consider (default: 3000)
--no-revcomp          don't reverse complement results on the minus strand
                      (default: do reverse complementation)
--downstream int      extend probe/primer match int bases downstream for
                      operon finding (default: 0)
--upstream int        extend probe/primer match int bases upstream for
                      operon finding (default: 0)
-o file, --out file   file for output (default: stdout)
--version             print version and exit
