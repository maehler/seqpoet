from __future__ import print_function

import os
import re

from nose.tools import raises

import seqpoet
from seqpoet.sequence import Sequence, DNA, IUPACDNA

class TestAlphabet:

    def test_dna_bases(self):
        assert re.match(r'^[acgtn]+$', DNA.bases, re.I)

    def test_iupac_dna_bases(self):
        assert re.match(r'^[ACGTMRWSYKVHDBN]+$', IUPACDNA.bases, re.I)

    def test_dna_revcomp(self):
        assert DNA.revcomp('ACGT') == 'ACGT'

    def test_iupacdna_revcomp(self):
        assert IUPACDNA.revcomp('NACRYSW') == 'WSRYGTN', \
            'expected WSRYGTN, got {0}'.format(IUPACDNA.revcomp('NACRYSW'))

    def test_equals(self):
        assert IUPACDNA.equals('nnnn', 'acgt')
        assert not IUPACDNA.equals('acgt', 'nnnn')
        assert IUPACDNA.equals('n', 'w')
        assert not IUPACDNA.equals('w', 'n')
        assert IUPACDNA.equals('CWGC', 'CAGC')
        assert DNA.equals('GTACACTATGTTTATAACAATCATCCAG',
                          'gtacactatgtttataacaatcatccag')

class TestSequenceAlphabet:

    @raises(ValueError)
    def test_iupac_seq_wrong_alphabet(self):
        assert Sequence('CAGWRSYKHVH', DNA)

    def test_iupac_seq_right_alphabet(self):
        assert isinstance(Sequence('CAGWRSYKHVH', IUPACDNA), Sequence)

class TestSequence:

    def setup(self):
        self.seq1 = 'ACATacacagaATAgagaCacaTA'
        self.illegal = 'agagcatgcacthisisnotcorrect'

    def test_sequence_length(self):
        s = seqpoet.Sequence(self.seq1)
        assert len(s) == len(self.seq1)

    def test_casing(self):
        s = seqpoet.Sequence(self.seq1)
        assert s == self.seq1

    def test_reverse_complement(self):
        s = seqpoet.Sequence(self.seq1)
        s2 = seqpoet.Sequence('acct')
        assert s.revcomp() == 'tatgtgtctctattctgtgtatgt', \
            '"{0}" is not "tatgtgtctctattctgtgtatgt"'.format(s.revcomp().seq)
        assert s2.revcomp() == 'aggt', \
            '"{0}" is not "aggt"'.format(s2.revcomp().seq)

    def test_str(self):
        s = seqpoet.Sequence(self.seq1)
        assert str(s) == self.seq1

    def test_repr(self):
        s = seqpoet.Sequence(self.seq1)
        assert repr(s) == '<Sequence: ACATa...>'
        assert repr(s.revcomp()) == '<Sequence: TAtgt...>'

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

    def test_example_primers(self):
        assert seqpoet.Sequence('GTACACTATGTTTATAACAATCATCCAG', DNA) == \
               seqpoet.Sequence('gtacactatgtttataacaatcatccag', DNA)
        assert seqpoet.Sequence('GCAAACCAGATTCAAAGTCAGTATG', DNA).revcomp() == \
               seqpoet.Sequence('catactgactttgaatctggtttgc', DNA)

    def test_inequality(self):
        s = seqpoet.Sequence(self.seq1)
        assert not s != self.seq1.lower()
        assert seqpoet.Sequence('acgtt') != seqpoet.Sequence('acgta')
        assert seqpoet.Sequence('a') != seqpoet.Sequence('t')

    def test_zip(self):
        for a, b in zip(self.seq1, self.seq1):
            assert a == b

    def test_iter(self):
        seq_iter = iter(self.seq1)
        seq_list = list(seq_iter)
        assert len(seq_list) == len(self.seq1)

    @raises(ValueError)
    def test_illegal_characters(self):
        s = seqpoet.Sequence(self.illegal)
