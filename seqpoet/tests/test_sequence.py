import os
import re

from nose.tools import raises

import seqpoet
from seqpoet.sequence import Sequence, DNA, IUPACDNA

class TestAlphabet:

    def test_dna_bases(self):
        assert DNA.bases == 'ACGT'

    def test_iupac_dna_bases(self):
        assert IUPACDNA.bases == 'ACGTMRWSYKVHDBN'

class TestSequenceAlphabet:

    @raises(ValueError)
    def test_iupac_seq_wrong_alphabet(self):
        assert Sequence('CAGWRSYKHVH', DNA())

    def test_iupac_seq_right_alphabet(self):
        assert isinstance(Sequence('CAGWRSYKHVH', IUPACDNA()), Sequence)

class TestSequence:

    def setup(self):
        self.seq1 = 'ACATacacagaATAgagaCacata'
        self.illegal = 'agagcatgcacthisisnotcorrect'

    def test_sequence_length(self):
        s = seqpoet.Sequence(self.seq1)
        assert len(s) == len(self.seq1)

    def test_casing(self):
        s = seqpoet.Sequence(self.seq1)
        assert re.match('^[acgt]+$', str(s))

    def test_reverse_complement(self):
        s = seqpoet.Sequence(self.seq1)
        s2 = seqpoet.Sequence('acct')
        assert s.revcomp() == 'tatgtgtctctattctgtgtatgt', \
            '"{0}" is not "tatgtgtctctattctgtgtatgt"'.format(s.revcomp().seq)
        assert s2.revcomp() == 'aggt', \
            '"{0}" is not "aggt"'.format(s2.revcomp().seq)

    def test_str(self):
        s = seqpoet.Sequence(self.seq1)
        assert str(s) == self.seq1.lower()

    def test_repr(self):
        s = seqpoet.Sequence(self.seq1)
        assert repr(s) == '<Sequence: acata...>'
        assert repr(s.revcomp()) == '<Sequence: tatgt...>'

    def test_indexing(self):
        s = seqpoet.Sequence(self.seq1)
        assert s[4] == 'a'
        assert s[:5] == 'acata'
        assert s[-6:] == 'cacata'
        assert s[4:8] == 'acac'

    def test_equality(self):
        s = seqpoet.Sequence(self.seq1)
        assert s == self.seq1.lower()
        assert s[:3] == seqpoet.Sequence(self.seq1[:3])

    @raises(ValueError)
    def test_illegal_characters(self):
        s = seqpoet.Sequence(self.illegal)
