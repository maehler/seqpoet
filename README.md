# seqpoet

The main purpose of seqpoet is to provide a simple interface for in silico PCR
and operon extraction in prokaryotes. The secondary purpose of seqpoet is to be
a Python package that can be used for handling sequence data in the form of
FASTA and GenBank files.

## Requirements

Currently, the only requirement is Python >= 2.7, but Python 3 is not supported
at the moment.

## Installation

Currently, the easiest way to install seqpoet is by downloading the source
files (or clone the repository) and (in the console) run

    python setup.py install

This should make the seqpoet script available system wide, and you can test
this by running

    seqpoet --version
