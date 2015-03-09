"""Classes and functions for representing DNA sequences.
"""

import re
import string

class Sequence(object):
    """Represent a DNA sequence.
    """

    revcomp_trans = string.maketrans('acgt', 'tgca')

    def __init__(self, seq):
        """Sequence constructor.

        To prevent confusion in sequence comparisons, sequences
        are represented as all lower case.

        Args:
            seq: a string representing a DNA sequence. Bases A, C,
                 G, T and N are allowed.
        Raises:
            ValuError: if the sequence contains illegal characters.
        """
        self.seq = seq.lower()
        if not re.match(r'^[acgtn]*$', self.seq):
            raise ValueError('illegal characters in sequence, currently ',
                             'only supports DNA sequences')

    def revcomp(self):
        """Get the reverse complement of the sequence.

        Returns:
            a sequence object representing the reverse complement
            of the sequence.
        """
        return Sequence(self.seq.translate(Sequence.revcomp_trans)[::-1])

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
