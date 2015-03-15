from nose.tools import raises
import os

from propex.search import search, hamming_distance
from propex import Sequence

class TestHammingDistance:

    @raises(ValueError)
    def test_sequences_of_different_length(self):
        hamming_distance('gattaca', 'gatt')

    def test_wikipedia_examples(self):
        assert hamming_distance('karolin', 'kathrin') == 3
        assert hamming_distance('karolin', 'kerstin') == 3
        assert hamming_distance('1011101', '1001001') == 2
        assert hamming_distance('2173896', '2233796') == 3

    def test_exact_matches(self):
        assert hamming_distance('karolin', 'karolin') == 0
        assert hamming_distance('gattaca', 'gattaca') == 0

    def test_one_mismatch(self):
        assert hamming_distance('niklas', 'niclas') == 1

class TestSearch:

    def setUp(self):
        self.haystack = 'accgtgacgggcacgaggcatcattatctagcagcacatg'
        self.needle = 'gaggcat'

    def test_exact_match(self):
        res = search(self.needle, self.haystack)
        assert res == [14], 'expected one match in pos 14, found {0}' \
            .format(str(res))

    def test_one_mismatch(self):
        res = search(self.needle, self.haystack, mismatches=1)
        assert res == [14], 'expected one match in pos 14, found {0}' \
            .format(str(res))

        res = search('ggg', self.haystack, mismatches=1)
        assert res == [3, 7, 8, 9, 14, 15, 16], 'found {0}'.format(str(res))
