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

    base_dict = {
        'A': { 'eq': set('A'), 'comp': 'T' },
        'C': { 'eq': set('C'), 'comp': 'G' },
        'G': { 'eq': set('G'), 'comp': 'C' },
        'T': { 'eq': set('T'), 'comp': 'A' }
    }
    bases      = ''.join(base_dict.keys())
    complement = ''.join(x['comp'] for x in base_dict.values())
    transtable = string.maketrans(bases + bases.lower(),
        complement + complement.lower())

    @classmethod
    def revcomp(self, s):
        return s.translate(self.transtable)[::-1]

class IUPACDNA(DNA):
    """DNA alpahbet with IUPAC ambiguity bases.
    """

    base_dict = {
        'M': { 'eq': set('MAC'), 'comp': 'K' },
        'R': { 'eq': set('RAG'), 'comp': 'Y' },
        'W': { 'eq': set('WAT'), 'comp': 'W' },
        'S': { 'eq': set('SCG'), 'comp': 'S' },
        'Y': { 'eq': set('YCT'), 'comp': 'R' },
        'K': { 'eq': set('KGT'), 'comp': 'M' },
        'V': { 'eq': set('VMRSACG'), 'comp': 'B' },
        'H': { 'eq': set('HMWYACT'), 'comp': 'D' },
        'D': { 'eq': set('DRWKAGT'), 'comp': 'H' },
        'B': { 'eq': set('BSYKCGT'), 'comp': 'V' },
        'N': { 'eq': set('NMRWSYKVHDBACGT'), 'comp': 'N' }
    }
    base_dict.update(DNA.base_dict)
    bases      = ''.join(base_dict.keys())
    complement = ''.join(x['comp'] for x in base_dict.values())
    transtable = string.maketrans(bases + bases.lower(),
        complement + complement.lower())

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
        self.seq = seq
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
            return self.seq.lower() == seq2.lower()
        return self.seq.lower() == seq2.seq.lower()

    def __ne__(self, seq2):
        return not self.__eq__(seq2)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __repr__(self):
        return '<Sequence: {0}...>'.format(self.seq[:5])
