#-*- encoding: utf-8 -*-
"""Classes and functions for representing DNA sequences.

.. module:: sequence
.. moduleauthor:: Niklas Mähler <niklas.mahler@gmail.com>
"""

import re
import string

class DNA(object):
    """DNA alphabet.
    """

    bases = 'ACGT'
    complement = 'TGCA'
    transtable = string.maketrans(bases, complement)

    @classmethod
    def revcomp(self, s):
        return s.upper().translate(self.transtable)[::-1]

class IUPACDNA(DNA):
    """DNA alpahbet with IUPAC ambiguity bases.
    """

    bases = DNA.bases + 'MRWSYKVHDBN'
    complement = DNA.complement + 'KYWSRMBDHVN'
    transtable = string.maketrans(bases, complement)

class Sequence(object):
    """Represent a DNA sequence.

    :param seq: a string representing a DNA sequence.
    :param alphabet: alphabet to use, see :py:class:`.DNA`
                     and :py:class:`.IUPACDNA`.
    :raises: ValueError if the sequence contains illegal characters.
    """

    #: Reverse complement translation table.
    _revcomp_trans = string.maketrans('acgt', 'tgca')

    def __init__(self, seq, alphabet=IUPACDNA()):
        """Sequence constructor.

        To prevent confusion in sequence comparisons, sequences
        are represented as all lower case.

        Args:
            seq: a string representing a DNA sequence.
            alphabet: alphabet to use.
        Raises:
            ValuError: if the sequence contains illegal characters.
        """
        self.alphabet = alphabet
        self.seq = seq.lower()
        if not re.match(r'^[{0}]*$'.format(alphabet.bases), self.seq, re.I):
            raise ValueError('illegal characters in sequence, '
                'not part of {0} class'.format(alphabet.__class__.__name__))

    def revcomp(self):
        """Get the reverse complement of the sequence.

        :returns:
            a sequence object representing the reverse complement
            of the sequence.
        """
        return Sequence(self.alphabet.revcomp(self.seq))

    def __getitem__(self, key):
        return Sequence(self.seq[key])

    def __eq__(self, seq2):
        if isinstance(seq2, basestring):
            return self.seq == seq2
        return self.seq == seq2.seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __repr__(self):
        return '<Sequence: {0}...>'.format(self.seq[:5])
