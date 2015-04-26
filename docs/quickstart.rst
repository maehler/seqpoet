Quick start
===========

Operon extraction
-----------------

.. code-block:: bash

	seqpoet --out output.fa input_directory probe.txt
	seqpoet --out output.fa input_directory primers.txt
	seqpoet --out output.fa input.gb probe.txt
	seqpoet --out output.fa input.gb primers.txt

The file ``input.gb`` should be a valid GenBank file (possibly with
multiple loci) and ``probe.txt`` should contain either a single nucleotide
sequence (*i.e.* a probe) or two nucleotide sequences (*i.e.* a primer pair).
Instead of supplying a single file, a directory of sequence files can be used
as the first argument.

Annotations are needed for the operon extraction, and currently GenBank
is the only supported format for this. The FASTA file ``output.fa`` will
contain the extracted sequences. If the ``--out`` argument is not supplied,
the results are written to stdout.

*In silico* PCR
---------------

.. code-block:: bash

	seqpoet --pcr --out output.fa input.gb primers.txt
	seqpoet --pcr --out output.fa input.fa primers.txt


For *in silico* PCR, only primer pairs are supported, but the sequence input
can be either FASTA or GenBank. The FASTA file ``output.fa`` will contain the
predicted PCR products. If the ``--out`` argument is not supplied,
the results are written to stdout.

Output
------

The output from both operon extraction and *in silico* PCR will be a FASTA
file. The header line for each result sequence is a colon separated string
and will look something like this:

::

	>input.gb:locus:3451:3812:28:+

- ``input.gb``: the original file where the sequence originates from
- ``locus``: the name of the sequence in the original file, either from a FASTA
  header or a GenBank locus name
- ``3451``: the position in the original sequence of the first nucleotide in
  the result sequence
- ``3812``: the position in the original sequence of the last nucleotide in the
  result sequence
- ``28``: the length of the original match
- ``+``: the sequence was found on the plus strand (otherwise ``-``)
